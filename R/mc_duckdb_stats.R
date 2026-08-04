# R/mc_duckdb_stats.R
#
# ===========================================================================
# WHY THIS IS A HYBRID (DuckDB + R), NOT PURE SQL
#
# Earlier versions of this file tried to do the ENTIRE risk computation --
# including the final reduction across 10,000 vessel iterations -- as SQL
# joins + GROUP BY. That doesn't work, and it's not a tuning problem: a
# relational join fundamentally can't replicate matrix multiplication
# efficiently for this workload. To sum a trip's contribution across every
# vessel iteration that selected it, a join has to materialize one physical
# row per (trip, whale_iter, vessel_iteration-that-selected-it) combination
# before it can aggregate them away. Matrix multiplication (D %*% S)
# computes the exact same result by accumulating directly into a bounded
# output matrix, with no such intermediate ever materialized. That's
# exactly why the original R implementation used %*% in the first place --
# and why the SQL-only approach kept blowing out RAM/disk regardless of how
# much it was tuned (approximate quantiles, batching, memory limits, etc.
# all treated symptoms of this, not the cause).
#
# So this file goes back to genuine sparse matrix multiplication for the
# reduction step, and only uses DuckDB for what it's actually good at:
# collapsing the (up to 800GB) event-level data down to small, per-group
# chunks on disk, via ONE full scan.
#
#   Stage 1 (materialize_per_trip_grp, DuckDB): collapse event-level D to
#     (group, trip, whale_iter) totals via a 1:1 join with the key table --
#     no fan-out here, this was never what crashed. Written to disk,
#     partitioned by group (Hive-style), so Stage 2 can read one group's
#     data at a time without rescanning everything.
#
#   Stage 2 (compute_group_stats_matmul, R): loop over groups ONE AT A TIME.
#     For each group: read its small collapsed table into R, build a sparse
#     matrix D_group (whale_iter x trip), multiply against the
#     vessel-iteration selection matrix S (trip x vessel_iter) -- this is
#     the same ~800MB-per-group dense result the original R implementation
#     produced -- reduce immediately via matrixStats, discard, move to the
#     next group. Structurally identical to the original per-grid-cell R
#     loop; the only difference is D is streamed from parquet instead of
#     pre-built as one giant in-memory matrix via cbind().
#
# Because each group's reduction is a real bounded in-memory matrix again,
# there's no more need for approximate quantiles -- matrixStats gives exact
# quantiles cheaply per group.
# ===========================================================================

#' Stage 1: collapse event-level D to (grp, trip, whale_iter) totals via a
#' full DuckDB scan (the expensive part, but it never fans out and happens
#' exactly once), written to disk partitioned by grp.
#'
#' IMPORTANT: d_paths and key_paths are STORE-RELATIVE paths (same
#' convention as the rest of the pipeline) -- resolved against store_root
#' here, at point of use, never stored resolved.
#'
#' @param d_paths Character vector of store-relative paths to long-format D
#'   parquet files (columns: event_key, whale_iter, value).
#' @param key_paths Character vector of store-relative paths to key parquet
#'   files (columns: event_key, trip_id, zone1, GRID_ID).
#' @param group_col "zone1" or "GRID_ID".
#' @param store_root Resolved get_store() root for THIS machine.
#' @param local_tmp_dir LOCAL (never NAS/store-relative) scratch directory
#'   for DuckDB's spill storage and the Stage 1 output itself. See
#'   _targets.R's duckdb_local_tmp_dir comments for why this must be local.
#' @param memory_limit,threads,max_temp_directory_size DuckDB tuning
#'   pragmas, forwarded as-is. See _targets.R for guidance on setting these.
#'
#' @return A list: stage1_dir (the partitioned parquet directory),
#'   instance_tmp_dir (DuckDB's spill dir, for cleanup), groups (character
#'   vector of distinct group values found), max_whale_iter (largest
#'   whale_iter value observed, needed to size Stage 2's sparse matrices
#'   consistently). NULL if there's no input data.
materialize_per_trip_grp <- function(d_paths,
                                     key_paths,
                                     group_col,
                                     store_root,
                                     local_tmp_dir = tempdir(),
                                     memory_limit = NULL,
                                     threads = NULL,
                                     max_temp_directory_size = NULL) {
  d_paths   <- d_paths[!is.na(d_paths)]
  key_paths <- key_paths[!is.na(key_paths)]

  if (length(d_paths) == 0 || length(key_paths) == 0) return(NULL)

  stopifnot(group_col %in% c("zone1", "GRID_ID"))

  d_full   <- file.path(store_root, d_paths)
  key_full <- file.path(store_root, key_paths)
  dir.create(local_tmp_dir, recursive = TRUE, showWarnings = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Isolated per-call spill directory -- see prior fix history: a shared
  # fixed temp_directory across concurrent crew workers caused spill-file
  # collisions and corruption.
  instance_tmp_dir <- tempfile(pattern = "duckdb_instance_", tmpdir = local_tmp_dir)
  dir.create(instance_tmp_dir, recursive = TRUE, showWarnings = FALSE)

  DBI::dbExecute(con, sprintf("PRAGMA temp_directory='%s'", instance_tmp_dir))
  if (!is.null(memory_limit)) {
    DBI::dbExecute(con, sprintf("PRAGMA memory_limit='%s'", memory_limit))
  }
  if (!is.null(threads)) {
    DBI::dbExecute(con, sprintf("PRAGMA threads=%d", as.integer(threads)))
  }
  if (!is.null(max_temp_directory_size)) {
    DBI::dbExecute(con, sprintf("PRAGMA max_temp_directory_size='%s'", max_temp_directory_size))
  }

  d_glob   <- paste0("['", paste(d_full, collapse = "','"), "']")
  key_glob <- paste0("['", paste(key_full, collapse = "','"), "']")

  stage1_dir <- tempfile(pattern = paste0("per_trip_grp_", group_col, "_"), tmpdir = local_tmp_dir)
  dir.create(stage1_dir, recursive = TRUE, showWarnings = FALSE)

  DBI::dbExecute(con, sprintf("
    COPY (
      SELECT k.%s AS grp, k.trip_id, d.whale_iter, SUM(d.value) AS value
      FROM read_parquet(%s) d
      JOIN read_parquet(%s) k USING (event_key)
      GROUP BY k.%s, k.trip_id, d.whale_iter
    ) TO '%s' (FORMAT PARQUET, PARTITION_BY (grp))
  ", group_col, d_glob, key_glob, group_col, stage1_dir))

  groups <- list.dirs(stage1_dir, recursive = FALSE, full.names = FALSE)
  groups <- sub("^grp=", "", groups)

  if (length(groups) == 0) {
    unlink(stage1_dir, recursive = TRUE, force = TRUE)
    unlink(instance_tmp_dir, recursive = TRUE, force = TRUE)
    return(NULL)
  }

  max_wi <- DBI::dbGetQuery(con, sprintf(
    "SELECT MAX(whale_iter) AS max_wi FROM read_parquet('%s')",
    file.path(stage1_dir, "*", "*.parquet")
  ))$max_wi

  list(
    stage1_dir       = stage1_dir,
    instance_tmp_dir = instance_tmp_dir,
    groups           = groups,
    max_whale_iter   = max_wi
  )
}

#' Stage 2: loop over groups one at a time, doing genuine sparse matrix
#' multiplication for each -- the direct replacement for the original R
#' per-grid-cell loop, sourcing D from disk instead of a pre-built dense
#' matrix.
#'
#' @param stage1_info Return value of materialize_per_trip_grp(). Its
#'   on-disk contents are deleted at the end of this function.
#' @param trip_lookup data.frame/data.table with trip_id and trip_index
#'   columns -- the SAME trip_index values used to build S.
#' @param S Sparse Matrix (trip_index x vessel_iter), e.g. the
#'   mc_sampled_trips target already in R memory. Kept in memory rather
#'   than round-tripped through parquet -- it's small (trip-level, not
#'   event-level), so there's no benefit to persisting it for this purpose.
#'
#' @return A data.frame with columns grp, whale_iter, delta_mean,
#'   delta_median, delta_variance, delta_q025, delta_q975, n_delta_gt0 --
#'   or NA if there was nothing to compute.
compute_group_stats_matmul <- function(stage1_info, trip_lookup, S) {
  if (is.null(stage1_info) || length(stage1_info$groups) == 0) return(NA)
  on.exit({
    unlink(stage1_info$stage1_dir, recursive = TRUE, force = TRUE)
    unlink(stage1_info$instance_tmp_dir, recursive = TRUE, force = TRUE)
  })

  n_trip_ids <- nrow(trip_lookup)
  trip_index_by_id <- setNames(trip_lookup$trip_index, as.character(trip_lookup$trip_id))
  max_wi <- stage1_info$max_whale_iter

  if (is.na(max_wi) || max_wi <= 0) return(NA)

  results <- vector("list", length(stage1_info$groups))

  for (i in seq_along(stage1_info$groups)) {
    grp_val  <- stage1_info$groups[i]
    part_dir <- file.path(stage1_info$stage1_dir, paste0("grp=", grp_val))

    dt <- tryCatch(
      arrow::open_dataset(part_dir, format = "parquet") |>
        dplyr::collect() |>
        data.table::as.data.table(),
      error = function(e) NULL
    )

    if (is.null(dt) || nrow(dt) == 0) next

    trip_idx <- trip_index_by_id[as.character(dt$trip_id)]
    keep <- !is.na(trip_idx)
    if (!any(keep)) next

    # D_group: whale_iter x trip_index, sparse. Same dimensional
    # conformity as the original D_group %*% S -- trip_index here lines up
    # exactly with S's row dimension (both built from the same
    # trip_lookup$trip_index).
    Dg <- Matrix::sparseMatrix(
      i    = dt$whale_iter[keep],
      j    = trip_idx[keep],
      x    = dt$value[keep],
      dims = c(max_wi, n_trip_ids)
    )
    rm(dt)

    # The actual reduction: real sparse matrix multiplication, not a SQL
    # join -- this is what bounds memory to one group's ~800MB result
    # instead of the combinatorial blowup a join produces.
    Rg <- as.matrix(Dg %*% S)  # whale_iter x vessel_iter, dense
    rm(Dg); gc()

    quants <- matrixStats::rowQuantiles(Rg, probs = c(0.025, 0.975))

    results[[i]] <- data.frame(
      grp            = grp_val,
      whale_iter     = seq_len(nrow(Rg)),
      delta_mean     = matrixStats::rowMeans2(Rg),
      delta_median   = matrixStats::rowMedians(Rg),
      delta_variance = matrixStats::rowVars(Rg),
      delta_q025     = quants[, 1],
      delta_q975     = quants[, 2],
      n_delta_gt0    = matrixStats::rowSums2(Rg > 0)
    )

    rm(Rg); gc()
  }

  dplyr::bind_rows(results)
}

#' Convenience wrapper running both stages. Prefer calling this from
#' _targets.R rather than the two stages separately, unless you need to
#' inspect the Stage 1 intermediate for debugging.
#'
#' @param d_paths,key_paths,group_col,store_root,local_tmp_dir,memory_limit,
#'   threads,max_temp_directory_size See materialize_per_trip_grp().
#' @param trip_lookup,S See compute_group_stats_matmul().
compute_stats_duckdb <- function(d_paths,
                                 key_paths,
                                 group_col,
                                 store_root,
                                 trip_lookup,
                                 S,
                                 local_tmp_dir = tempdir(),
                                 memory_limit = NULL,
                                 threads = NULL,
                                 max_temp_directory_size = NULL) {
  stage1_info <- materialize_per_trip_grp(
    d_paths, key_paths, group_col, store_root,
    local_tmp_dir = local_tmp_dir,
    memory_limit = memory_limit,
    threads = threads,
    max_temp_directory_size = max_temp_directory_size
  )
  compute_group_stats_matmul(stage1_info, trip_lookup, S)
}

#' Melt a list of whale-iteration matrices (whale_iter x event) to a long
#' data.table and write to parquet, filtering out zero/NA cells for free
#' space savings. Called once per source .RData file within a zone/chunk
#' target, so peak memory is one file's worth of matrices, not the whole
#' zone's.
#'
#' IMPORTANT: rel_path is a STORE-RELATIVE path (e.g.
#' "whale_long/D_2024-08_NEGSL_1_1.parquet"), not a full filesystem path.
#' The function resolves it against store_root to know where to actually
#' write, but always returns the relative path -- that's what should end up
#' in a target's return value, so the pipeline looks unchanged when re-run
#' from a machine where get_store() resolves to a different root.
#'
#' @param mats List of whale_iter x event matrices to melt.
#' @param store_root Resolved root of the data store for THIS machine --
#'   pass get_store() from the calling target.
#' @param rel_path Store-relative output path.
melt_matrices_to_parquet <- function(mats, store_root, rel_path) {
  if (length(mats) == 0) return(NA_character_)

  long_dt <- data.table::rbindlist(lapply(mats, function(m) {
    v    <- as.vector(m)
    keep <- !is.na(v) & v != 0
    data.table::data.table(
      whale_iter = rep.int(seq_len(nrow(m)), ncol(m))[keep],
      event_key  = rep(colnames(m), each = nrow(m))[keep],
      value      = v[keep]
    )
  }))

  if (nrow(long_dt) == 0) return(NA_character_)

  out_path <- file.path(store_root, rel_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(out_path, ".tmp")
  arrow::write_parquet(long_dt, tmp_path)
  file.rename(tmp_path, out_path)
  rel_path
}

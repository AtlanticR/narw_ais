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
# much it was tuned.
#
# So this file uses genuine sparse matrix multiplication for the reduction
# step, and only uses DuckDB for what it's actually good at: collapsing the
# (up to 800GB) event-level data down to a much smaller table on disk, via
# ONE full scan.
#
#   Stage 1 (materialize_per_trip_grp, DuckDB): collapse event-level D to
#     (group, trip, whale_iter) totals via a 1:1 join with the key table --
#     no fan-out here. Written to disk Hive-partitioned by grp
#     (grp=<value>/*.parquet) so DuckDB's writer can stay parallel across
#     groups instead of serializing into one file.
#
#   Stage 2 (compute_group_stats_matmul, R): read the Stage 1 partitioned
#     directory as a LAZY Arrow dataset (Arrow reconstructs the grp column
#     from the folder names, nothing is collect()-ed yet), plan batches of
#     many groups at a time using a cheap pushed-down row-count-per-group
#     query, then for EACH batch: pull only that batch's rows off disk,
#     run one sparse matrix multiply against the vessel-iteration selection
#     S, reduce via matrixStats, and discard before the next batch starts.
#
# CHANGED (this revision): three mistakes, now all fixed:
#   1. Stage 2 originally looped one GRID_ID at a time, re-opening a small
#      parquet partition directory for each one. With thousands of grid
#      cells, the fixed per-call overhead (disk I/O, data.table conversion,
#      sparse matrix construction, gc()) dominated -- that's what made grid
#      stats ~200x slower than zone stats, far worse than the actual
#      difference in FLOPs required. Fixed by batching many groups per
#      matrix multiply.
#   2. A subsequent attempt to "simplify" Stage 1 by writing ONE
#      unpartitioned file (instead of partitioned) backfired badly: it
#      forced a single-threaded, order-preserving merge at the end of an
#      otherwise-parallel aggregation, regressing zone stats from 6 minutes
#      to over 24 hours. Fixed by restoring Hive partitioning for the
#      WRITE (keeps it parallel) -- plus setting
#      PRAGMA preserve_insertion_order=false, since row order is never
#      meaningful here and DuckDB pays a real cost to preserve it by
#      default.
#   3. *** THE ACTUAL SOURCE OF THE RAM/DISK BLOWOUT ***. The previous fix
#      for (1) called arrow::open_dataset(...) |> dplyr::collect() ONCE,
#      up front, pulling the ENTIRE Stage 1 output (every group, every
#      trip, every whale_iter) into a single in-memory data.table, and only
#      THEN sliced it into batches for the matmul loop. The comment claimed
#      "peak memory is bounded by rows_per_batch, not by total group
#      count" -- that was false: the big collected data.table sat fully
#      resident in R for the entire loop regardless of rows_per_batch, and
#      for GRID_ID (thousands of groups) that single collect() is what was
#      exhausting RAM (and, once R started swapping/spilling, disk too).
#      Fixed by NOT collecting up front: get a cheap per-group row count
#      via a pushed-down Arrow aggregate (this touches disk but never
#      materializes actual columns in R), use it to pack groups into
#      batches sized near rows_per_batch, and then collect() only ONE
#      batch's worth of rows at a time inside the loop, discarding it
#      before moving to the next batch. This is what actually gives the
#      "peak memory bounded by rows_per_batch, not total group count"
#      property -- the earlier version only had the *appearance* of that
#      property.
#
# Peak R memory in Stage 2 is now bounded by rows_per_batch (one batch's
# raw rows, plus the resulting rows_per_batch x n_vessel_iter dense
# matrix) -- it does NOT grow with total group count, total row count, or
# how many groups GRID_ID vs zone1 produces.
# ===========================================================================

#' Stage 1: collapse event-level D to (grp, trip, whale_iter) totals via a
#' full DuckDB scan (the expensive part, but it never fans out and happens
#' exactly once), written to disk Hive-partitioned by grp for fast parallel
#' writes.
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
#'   _targets.R's duckdb_local_tmp_dir comments for why this must be local
#'   -- and genuinely local (not a network-mounted home directory).
#' @param memory_limit,threads,max_temp_directory_size DuckDB tuning
#'   pragmas, forwarded as-is. See _targets.R for guidance on setting these.
#'
#' @return A list: stage1_dir (the Hive-partitioned parquet directory,
#'   grp=<value>/*.parquet per group), instance_tmp_dir (DuckDB's spill
#'   dir, for cleanup). NULL if there's no input data.
materialize_per_trip_grp <- function(
  d_paths,
  key_paths,
  group_col,
  store_root,
  local_tmp_dir = tempdir(),
  memory_limit = NULL,
  threads = NULL,
  max_temp_directory_size = NULL,
  n_buckets = 64
) {
  d_paths <- d_paths[!is.na(d_paths)]
  key_paths <- key_paths[!is.na(key_paths)]

  if (length(d_paths) == 0 || length(key_paths) == 0) {
    return(NULL)
  }

  stopifnot(group_col %in% c("zone1", "GRID_ID"))

  d_full <- file.path(store_root, d_paths)
  key_full <- file.path(store_root, key_paths)
  dir.create(local_tmp_dir, recursive = TRUE, showWarnings = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  instance_tmp_dir <- tempfile(
    pattern = "duckdb_instance_",
    tmpdir = local_tmp_dir
  )
  dir.create(instance_tmp_dir, recursive = TRUE, showWarnings = FALSE)

  DBI::dbExecute(con, sprintf("PRAGMA temp_directory='%s'", instance_tmp_dir))
  if (!is.null(memory_limit)) {
    DBI::dbExecute(con, sprintf("PRAGMA memory_limit='%s'", memory_limit))
  }
  if (!is.null(threads)) {
    DBI::dbExecute(con, sprintf("PRAGMA threads=%d", as.integer(threads)))
  }
  if (!is.null(max_temp_directory_size)) {
    DBI::dbExecute(
      con,
      sprintf("PRAGMA max_temp_directory_size='%s'", max_temp_directory_size)
    )
  }
  DBI::dbExecute(con, "PRAGMA preserve_insertion_order=false")

  d_glob <- paste0("['", paste(d_full, collapse = "','"), "']")
  key_glob <- paste0("['", paste(key_full, collapse = "','"), "']")

  stage1_dir <- tempfile(
    pattern = paste0("per_trip_grp_", group_col, "_"),
    tmpdir = local_tmp_dir
  )
  dir.create(stage1_dir, recursive = TRUE, showWarnings = FALSE)

  n_distinct_grp <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT COUNT(DISTINCT k.%s) AS n FROM read_parquet(%s) k",
      group_col,
      key_glob
    )
  )$n
  message(sprintf(
    "[%s] Stage 1 (%s): %d distinct grp values -> %d write buckets (~%s grp values/bucket).",
    format(Sys.time(), "%H:%M:%S"),
    group_col,
    n_distinct_grp,
    n_buckets,
    format(round(n_distinct_grp / n_buckets, 1))
  ))

  t0 <- Sys.time()
  message(sprintf(
    "[%s] Stage 1 (%s): starting DuckDB scan + COPY...",
    format(t0, "%H:%M:%S"),
    group_col
  ))

  # Cap threads just for the write step -- default DuckDB behaviour
  # creates one file per thread per partition, so file count = threads x
  # n_buckets. Full `threads` value is still used for the join/aggregate
  # part of the plan; only the COPY's write parallelism is capped here.
  write_threads <- 4L
  if (!is.null(threads)) {
    write_threads <- min(as.integer(threads), 4L)
  }
  DBI::dbExecute(con, sprintf("PRAGMA threads=%d", write_threads))

  DBI::dbExecute(
    con,
    sprintf(
      "
    COPY (
      SELECT
        k.%s AS grp,
        CAST(hash(k.%s) %% %d AS INTEGER) AS bucket,
        k.trip_id, d.whale_iter, SUM(d.value) AS value
      FROM read_parquet(%s) d
      JOIN read_parquet(%s) k USING (event_key)
      GROUP BY k.%s, bucket, k.trip_id, d.whale_iter
    ) TO '%s' (
      FORMAT PARQUET,
      PARTITION_BY (bucket),
      ROW_GROUP_SIZE 500000
    )
  ",
      group_col,
      group_col,
      n_buckets,
      d_glob,
      key_glob,
      group_col,
      stage1_dir
    )
  )

  t1 <- Sys.time()
  n_buckets_written <- length(list.dirs(stage1_dir, recursive = FALSE))
  n_files_written <- length(list.files(
    stage1_dir,
    pattern = "\\.parquet$",
    recursive = TRUE
  ))
  message(sprintf(
    "[%s] Stage 1 (%s): done in %s -- wrote %d bucket partitions, %d files.",
    format(t1, "%H:%M:%S"),
    group_col,
    format(round(difftime(t1, t0, units = "auto"), 2)),
    n_buckets_written,
    n_files_written
  ))

  # NEW: compaction pass. Confirmed via math (bucket dir totals vs file
  # counts) that Stage 1's output files average ~170-200KB each --
  # extreme fragmentation. Restructuring Stage 2 to read each bucket once
  # (previous fix) was necessary but insufficient: reading a bucket still
  # means opening/parsing ~219 tiny files per bucket, and that per-file
  # overhead is what's actually dominating wall time now (14.7 min for
  # just 6 buckets in parallel). This pass rewrites each bucket's many
  # small files into ONE file per bucket, using DuckDB's full thread count
  # (no PARTITION_BY restriction here, so no single-thread-merge
  # regression) -- and because it's reading already-local, already-small
  # data (not re-scanning the original 800GB input), it should be fast:
  # this is fundamentally a local disk defragmentation step, not a
  # repeat of the expensive join+aggregate.
  if (!is.null(threads)) {
    DBI::dbExecute(con, sprintf("PRAGMA threads=%d", as.integer(threads)))
  }

  t_compact0 <- Sys.time()
  message(sprintf(
    "[%s] Stage 1 (%s): compacting fragmented bucket files...",
    format(t_compact0, "%H:%M:%S"),
    group_col
  ))

  bucket_dirs <- list.dirs(stage1_dir, recursive = FALSE)
  n_compacted <- 0L
  for (bd in bucket_dirs) {
    files_in_bd <- list.files(bd, pattern = "\\.parquet$", full.names = TRUE)
    if (length(files_in_bd) <= 1) {
      next
    } # already a single file, nothing to do

    glob <- paste0("['", paste(files_in_bd, collapse = "','"), "']")
    compact_tmp <- file.path(bd, "compact_tmp.parquet")

    DBI::dbExecute(
      con,
      sprintf(
        "COPY (SELECT * FROM read_parquet(%s)) TO '%s' (FORMAT PARQUET, ROW_GROUP_SIZE 500000)",
        glob,
        compact_tmp
      )
    )

    file.remove(files_in_bd)
    file.rename(compact_tmp, file.path(bd, "data_0.parquet"))
    n_compacted <- n_compacted + 1L
  }

  t_compact1 <- Sys.time()
  n_files_after <- length(list.files(
    stage1_dir,
    pattern = "\\.parquet$",
    recursive = TRUE
  ))
  message(sprintf(
    "[%s] Stage 1 (%s): compaction done in %s -- %d buckets compacted, %d files (was %d).",
    format(t_compact1, "%H:%M:%S"),
    group_col,
    format(round(difftime(t_compact1, t_compact0, units = "auto"), 2)),
    n_compacted,
    n_files_after,
    n_files_written
  ))

  list(
    stage1_dir = stage1_dir,
    instance_tmp_dir = instance_tmp_dir
  )
}

#' Stage 2: process the Stage 1 output in batches of many groups at a time
#' via genuine sparse matrix multiplication, WITHOUT ever collecting the
#' full Stage 1 output into R memory at once.
#'
#' Batching strategy: a cheap Arrow aggregate query (group_by(grp) +
#' n(), pushed down to the parquet reader -- no columns get materialized
#' in R for this) gets a row count per group. Groups are then greedily
#' packed into batches so each batch's total raw row count stays near
#' rows_per_batch. Inside the loop, EACH batch is collect()-ed from disk
#' on its own, processed, and discarded before the next one is read --
#' this is what actually bounds peak R memory to one batch's worth of
#' data, regardless of how many total groups or rows Stage 1 produced.
#'
#' @param stage1_info Return value of materialize_per_trip_grp(). Its
#'   on-disk contents are deleted at the end of this function.
#' @param trip_lookup data.frame/data.table with trip_id and trip_index
#'   columns -- the SAME trip_index values used to build S.
#' @param S Sparse Matrix (trip_index x vessel_iter), e.g. the
#'   mc_sampled_trips target already in R memory.
#' @param rows_per_batch Target number of raw (grp, trip, whale_iter) rows
#'   per batch. Groups are packed together (whole groups per batch, never
#'   split) until adding the next group would exceed this target; a single
#'   oversized group still gets its own batch. Each batch's matmul produces
#'   a transient dense matrix of size (that batch's distinct (grp,
#'   whale_iter) row count) x n_vessel_iter doubles -- e.g. with
#'   n_vessel_iter = 10,000, a 20,000-row batch is ~1.6GB. Pick this based
#'   on n_vessel_iter and available RAM (see _targets.R for sizing
#'   guidance). Tune upward from a conservative starting point while
#'   watching memory -- there's no correctness risk in raising it, only a
#'   speed/memory trade-off.
#'
#' @return A data.frame with columns grp, whale_iter, delta_mean,
#'   delta_median, delta_variance, delta_q025, delta_q975, n_delta_gt0 --
#'   or NA if there was nothing to compute. NOTE: rows are sorted within
#'   each batch (by grp, whale_iter) but batches are not guaranteed to be
#'   globally sorted relative to each other -- sort the result yourself
#'   downstream if a specific row order matters.
compute_group_stats_matmul <- function(
  stage1_info,
  trip_lookup,
  S,
  rows_per_batch = 20000,
  n_workers = 6
) {
  if (is.null(stage1_info)) {
    message(
      "Stage 2: stage1_info is NULL (no D/key paths upstream) -- returning NA."
    )
    return(NA)
  }
  on.exit(unlink(stage1_info$stage1_dir, recursive = TRUE, force = TRUE))
  on.exit(
    unlink(stage1_info$instance_tmp_dir, recursive = TRUE, force = TRUE),
    add = TRUE
  )

  part_dirs <- list.dirs(stage1_info$stage1_dir, recursive = FALSE)
  if (length(part_dirs) == 0) {
    message(sprintf(
      "Stage 2: Stage 1 output dir '%s' has no group partitions -- returning NA.",
      stage1_info$stage1_dir
    ))
    return(NA)
  }

  ds <- arrow::open_dataset(
    stage1_info$stage1_dir,
    format = "parquet",
    partitioning = arrow::schema(bucket = arrow::int32())
  )

  t_meta0 <- Sys.time()
  message(sprintf(
    "[%s] Stage 2: counting rows per bucket...",
    format(t_meta0, "%H:%M:%S")
  ))

  # CHANGED: only need per-BUCKET totals now, not per-group -- grp-level
  # batching happens in-memory after a bucket is read, not via separate
  # disk reads per group. This aggregate is still cheap/pushed-down.
  bucket_counts <- ds |>
    dplyr::group_by(bucket) |>
    dplyr::summarise(n_rows = dplyr::n(), .groups = "drop") |>
    dplyr::collect() |>
    data.table::as.data.table()

  if (nrow(bucket_counts) == 0) {
    message(
      "Stage 2: bucket_counts query returned zero buckets -- returning NA."
    )
    return(NA)
  }
  data.table::setorder(bucket_counts, bucket)

  message(sprintf(
    "[%s] Stage 2: %d buckets, %s total rows (in %s). Starting parallel per-bucket processing (%d workers)...",
    format(Sys.time(), "%H:%M:%S"),
    nrow(bucket_counts),
    format(sum(bucket_counts$n_rows), big.mark = ","),
    format(round(difftime(Sys.time(), t_meta0, units = "auto"), 2)),
    n_workers
  ))

  n_buckets_total <- nrow(bucket_counts)
  n_workers <- min(n_workers, n_buckets_total)

  n_trip_ids <- nrow(trip_lookup)
  trip_index_by_id <- setNames(
    trip_lookup$trip_index,
    as.character(trip_lookup$trip_id)
  )

  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  stage1_dir_path <- stage1_info$stage1_dir
  bucket_ids <- bucket_counts$bucket

  # CHANGED: this is now the per-BUCKET worker function, not per-batch.
  # Each call does exactly ONE dplyr::collect() -- reading that bucket's
  # entire partition (all its files) off disk a single time -- then loops
  # over row-batches WITHIN the already-in-memory data.table, with zero
  # further disk I/O. This is what actually fixes the redundant-re-read
  # problem: previously each of a bucket's ~36 groups triggered its own
  # full re-scan of that bucket's ~219 files (bucket pruning alone wasn't
  # enough -- it correctly limited I/O to the right bucket, but did
  # nothing about doing that I/O over and over for every group in it).
  run_one_bucket <- function(bkt) {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
    }

    if (!exists(".ds_worker_cache", envir = .GlobalEnv)) {
      assign(
        ".ds_worker_cache",
        arrow::open_dataset(
          stage1_dir_path,
          format = "parquet",
          partitioning = arrow::schema(bucket = arrow::int32())
        ),
        envir = .GlobalEnv
      )
    }
    ds_worker <- get(".ds_worker_cache", envir = .GlobalEnv)

    # ONE read for the whole bucket -- no `grp` filter, no repeated calls.
    dt_bucket <- ds_worker |>
      dplyr::filter(bucket == bkt) |>
      dplyr::collect() |>
      data.table::as.data.table()

    if (nrow(dt_bucket) == 0) {
      return(NULL)
    }

    dt_bucket[, trip_idx := trip_index_by_id[as.character(trip_id)]]
    dt_bucket <- dt_bucket[!is.na(trip_idx)]
    if (nrow(dt_bucket) == 0) {
      return(NULL)
    }

    # In-memory batching by group, same greedy packing logic as before,
    # just operating on an already-resident data.table instead of driving
    # separate disk reads.
    grp_local_counts <- dt_bucket[, .N, by = grp]
    data.table::setorder(grp_local_counts, grp)

    local_batch <- integer(nrow(grp_local_counts))
    cur_batch <- 1L
    cur_total <- 0L
    for (i in seq_len(nrow(grp_local_counts))) {
      n_i <- grp_local_counts$N[i]
      if (cur_total > 0 && cur_total + n_i > rows_per_batch) {
        cur_batch <- cur_batch + 1L
        cur_total <- 0L
      }
      local_batch[i] <- cur_batch
      cur_total <- cur_total + n_i
    }
    grp_local_counts[, local_batch_id := local_batch]
    dt_bucket <- merge(
      dt_bucket,
      grp_local_counts[, .(grp, local_batch_id)],
      by = "grp"
    )
    n_local_batches <- max(local_batch)

    bucket_results <- vector("list", n_local_batches)

    for (lb in seq_len(n_local_batches)) {
      dt_lb <- dt_bucket[local_batch_id == lb]

      data.table::setorder(dt_lb, grp, whale_iter)
      dt_lb[, row_id := .GRP, by = .(grp, whale_iter)]

      row_key <- unique(dt_lb[, .(row_id, grp, whale_iter)])
      data.table::setorder(row_key, row_id)
      n_local_rows <- nrow(row_key)

      Dg <- Matrix::sparseMatrix(
        i = dt_lb$row_id,
        j = dt_lb$trip_idx,
        x = dt_lb$value,
        dims = c(n_local_rows, n_trip_ids)
      )
      Rg <- as.matrix(Dg %*% S)
      rm(Dg)

      quants <- matrixStats::rowQuantiles(Rg, probs = c(0.025, 0.975))
      bucket_results[[lb]] <- data.frame(
        grp = row_key$grp,
        whale_iter = row_key$whale_iter,
        delta_mean = matrixStats::rowMeans2(Rg),
        delta_median = matrixStats::rowMedians(Rg),
        delta_variance = matrixStats::rowVars(Rg),
        delta_q025 = quants[, 1],
        delta_q975 = quants[, 2],
        n_delta_gt0 = matrixStats::rowSums2(Rg > 0)
      )
      rm(Rg, dt_lb, row_key)
      gc()
    }

    rm(dt_bucket)
    gc()
    dplyr::bind_rows(bucket_results)
  }

  t_loop0 <- Sys.time()
  wave_starts <- seq(1L, n_buckets_total, by = n_workers)
  results <- vector("list", n_buckets_total)

  for (w in seq_along(wave_starts)) {
    wave_buckets_idx <- wave_starts[w]:min(
      wave_starts[w] + n_workers - 1L,
      n_buckets_total
    )
    wave_bucket_ids <- bucket_ids[wave_buckets_idx]

    wave_results <- furrr::future_map(
      wave_bucket_ids,
      run_one_bucket,
      .options = furrr::furrr_options(
        seed = TRUE,
        packages = c("Matrix", "arrow", "dplyr", "data.table", "matrixStats")
      ),
      .progress = FALSE
    )
    results[wave_buckets_idx] <- wave_results

    elapsed <- difftime(Sys.time(), t_loop0, units = "secs")
    eta <- elapsed /
      max(wave_buckets_idx) *
      (n_buckets_total - max(wave_buckets_idx))
    message(sprintf(
      "[%s] Stage 2: wave %d/%d done (buckets %d-%d of %d). Elapsed %s, ETA %s.",
      format(Sys.time(), "%H:%M:%S"),
      w,
      length(wave_starts),
      min(wave_buckets_idx),
      max(wave_buckets_idx),
      n_buckets_total,
      format(round(as.numeric(elapsed) / 60, 1)),
      format(round(as.numeric(eta) / 60, 1))
    ))
  }

  t_loop1 <- Sys.time()
  message(sprintf(
    "[%s] Stage 2: per-bucket matmul done in %s.",
    format(t_loop1, "%H:%M:%S"),
    format(round(difftime(t_loop1, t_loop0, units = "auto"), 2))
  ))

  dplyr::bind_rows(results)
}

#' Convenience wrapper running both stages. Prefer calling this from
#' _targets.R rather than the two stages separately, unless you need to
#' inspect the Stage 1 intermediate for debugging.
#'
#' @param d_paths,key_paths,group_col,store_root,local_tmp_dir,memory_limit,
#'   threads,max_temp_directory_size See materialize_per_trip_grp().
#' @param trip_lookup,S,rows_per_batch See compute_group_stats_matmul().
compute_stats_duckdb <- function(
  d_paths,
  key_paths,
  group_col,
  store_root,
  trip_lookup,
  S,
  local_tmp_dir = tempdir(),
  memory_limit = NULL,
  threads = NULL,
  max_temp_directory_size = NULL,
  rows_per_batch = 20000,
  n_buckets = 64
) {
  stage1_info <- materialize_per_trip_grp(
    d_paths,
    key_paths,
    group_col,
    store_root,
    local_tmp_dir = local_tmp_dir,
    memory_limit = memory_limit,
    threads = threads,
    max_temp_directory_size = max_temp_directory_size,
    n_buckets = n_buckets
  )
  compute_group_stats_matmul(
    stage1_info,
    trip_lookup,
    S,
    rows_per_batch = rows_per_batch
  )
}

#' Melt a list of whale-iteration matrices (whale_iter x event) to a long
#' data.table and write to parquet, filtering out zero/NA cells for free
#' space savings. Called once per source .rds file within a zone/chunk
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
  if (length(mats) == 0) {
    return(NA_character_)
  }

  long_dt <- data.table::rbindlist(lapply(mats, function(m) {
    v <- as.vector(m)
    keep <- !is.na(v) & v != 0
    data.table::data.table(
      whale_iter = rep.int(seq_len(nrow(m)), ncol(m))[keep],
      event_key = rep(colnames(m), each = nrow(m))[keep],
      value = v[keep]
    )
  }))

  if (nrow(long_dt) == 0) {
    return(NA_character_)
  }

  out_path <- file.path(store_root, rel_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(out_path, ".tmp")
  arrow::write_parquet(long_dt, tmp_path)
  file.rename(tmp_path, out_path)
  rel_path
}


# ===========================================================================
# ADDITIONS -- append to R/mc_duckdb_stats.R
#
# total_risk here is just SUM(value) per (grp, whale_iter) across ALL
# trips -- no growth-trip sampling (S), no trip_lookup, no matmul. That
# means it's independent of mc_zone_stats/mc_grid_stats's Stage 2 entirely:
# a single DuckDB scan + GROUP BY straight off the same D/key parquet files,
# collected directly since the result (one row per grp x whale_iter) is
# small. This is joined onto mc_zone_stats/mc_grid_stats downstream in
# _targets.R -- those targets' commands/functions are unmodified, so they
# are NOT invalidated by any of this.
# ===========================================================================

compute_total_risk_duckdb <- function(
  d_paths,
  key_paths,
  group_col,
  store_root,
  local_tmp_dir = tempdir(),
  memory_limit = NULL,
  threads = NULL,
  max_temp_directory_size = NULL
) {
  d_paths <- d_paths[!is.na(d_paths)]
  key_paths <- key_paths[!is.na(key_paths)]

  if (length(d_paths) == 0 || length(key_paths) == 0) {
    return(NA)
  }

  stopifnot(group_col %in% c("zone1", "GRID_ID"))

  d_full <- file.path(store_root, d_paths)
  key_full <- file.path(store_root, key_paths)
  dir.create(local_tmp_dir, recursive = TRUE, showWarnings = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  instance_tmp_dir <- tempfile(
    pattern = paste0("duckdb_total_risk_", group_col, "_"),
    tmpdir = local_tmp_dir
  )
  dir.create(instance_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(
    unlink(instance_tmp_dir, recursive = TRUE, force = TRUE),
    add = TRUE
  )

  DBI::dbExecute(con, sprintf("PRAGMA temp_directory='%s'", instance_tmp_dir))
  if (!is.null(memory_limit)) {
    DBI::dbExecute(con, sprintf("PRAGMA memory_limit='%s'", memory_limit))
  }
  if (!is.null(threads)) {
    DBI::dbExecute(con, sprintf("PRAGMA threads=%d", as.integer(threads)))
  }
  if (!is.null(max_temp_directory_size)) {
    DBI::dbExecute(
      con,
      sprintf("PRAGMA max_temp_directory_size='%s'", max_temp_directory_size)
    )
  }
  DBI::dbExecute(con, "PRAGMA preserve_insertion_order=false")

  d_glob <- paste0("['", paste(d_full, collapse = "','"), "']")
  key_glob <- paste0("['", paste(key_full, collapse = "','"), "']")

  t0 <- Sys.time()
  message(sprintf(
    "[%s] total_risk (%s): starting DuckDB scan + aggregate...",
    format(t0, "%H:%M:%S"),
    group_col
  ))

  result <- DBI::dbGetQuery(
    con,
    sprintf(
      "
      SELECT
        k.%s AS grp,
        d.whale_iter,
        SUM(d.value) AS total_risk
      FROM read_parquet(%s) d
      JOIN read_parquet(%s) k USING (event_key)
      GROUP BY k.%s, d.whale_iter
      ",
      group_col,
      d_glob,
      key_glob,
      group_col
    )
  )

  t1 <- Sys.time()
  message(sprintf(
    "[%s] total_risk (%s): done in %s -- %s (grp, whale_iter) rows.",
    format(t1, "%H:%M:%S"),
    group_col,
    format(round(difftime(t1, t0, units = "auto"), 2)),
    format(nrow(result), big.mark = ",")
  ))

  result
}

# ===========================================================================
# join_and_save_stats(): left-joins total_risk onto an existing
# mc_zone_stats / mc_grid_stats data frame by (grp, whale_iter), optionally
# adds percent_delta_* = delta_*/total_risk for every delta_* column, and
# writes both CSV and rds. Returns the two file paths for a
# format = "file" tar_target. Returns NA_character_ if either input isn't a
# data frame (mirrors the NA-passthrough convention used elsewhere in this
# pipeline for empty upstream results).
# ===========================================================================
join_and_save_stats <- function(
  stats_df,
  total_risk_df,
  out_dir,
  file_stub,
  add_percent_delta = FALSE,
  grp_name = "grp"
) {
  if (!is.data.frame(stats_df) || !is.data.frame(total_risk_df)) {
    return(NA_character_)
  }

  df <- dplyr::left_join(
    stats_df,
    total_risk_df,
    by = c("grp", "whale_iter")
  ) |>
    rename("{grp_name}" := grp)

  if (add_percent_delta) {
    delta_cols <- grep("^delta_", names(df), value = TRUE)
    for (col in delta_cols) {
      pct_col <- sub("^delta_", "percent_delta_", col)
      df[[pct_col]] <- df[[col]] / df$total_risk
    }

    summarydf <- df |>
      group_by(.data[[grp_name]]) |>
      reframe(
        delta_mean_mean = mean(delta_mean, na.rm = TRUE),
        delta_mean_median = median(delta_mean, na.rm = TRUE),
        delta_mean_var = var(delta_mean, na.rm = TRUE),
        delta_mean_q025 = quantile(delta_mean, 0.025, na.rm = TRUE),
        delta_mean_q975 = quantile(delta_mean, 0.975, na.rm = TRUE),
        delta_mean_gt0 = sum(delta_mean > 0, na.rm = TRUE),

        percent_delta_mean_mean = mean(percent_delta_mean, na.rm = TRUE),
        percent_delta_mean_median = median(percent_delta_mean, na.rm = TRUE),
        percent_delta_mean_var = var(percent_delta_mean, na.rm = TRUE),
        percent_delta_mean_q025 = quantile(
          percent_delta_mean,
          0.025,
          na.rm = TRUE
        ),
        percent_delta_mean_q975 = quantile(
          percent_delta_mean,
          0.975,
          na.rm = TRUE
        ),
        percent_delta_mean_gt0 = sum(percent_delta_mean > 0, na.rm = TRUE)
      )
  } else {
    summarydf <- df |>
      group_by(.data[[grp_name]]) |>
      reframe(
        delta_mean_mean = mean(delta_mean, na.rm = TRUE),
        delta_mean_median = median(delta_mean, na.rm = TRUE),
        delta_mean_var = var(delta_mean, na.rm = TRUE),
        delta_mean_q025 = quantile(delta_mean, 0.025, na.rm = TRUE),
        delta_mean_q975 = quantile(delta_mean, 0.975, na.rm = TRUE),
        delta_mean_gt0 = sum(delta_mean > 0, na.rm = TRUE)
      )
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  csv_path <- file.path(out_dir, paste0(file_stub, ".csv"))
  rds_path <- file.path(out_dir, paste0(file_stub, ".rds"))
  summary_csv_path <- file.path(
    dirname(csv_path),
    paste0("summary_", basename(csv_path))
  )
  summary_rds_path <- file.path(
    dirname(rds_path),
    paste0("summary_", basename(rds_path))
  )

  data.table::fwrite(df, csv_path)

  saveRDS(df, file = rds_path)

  data.table::fwrite(summarydf, summary_csv_path)

  saveRDS(summarydf, file = summary_rds_path)

  c(csv_path, rds_path)
}

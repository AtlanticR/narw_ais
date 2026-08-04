if (!require("librarian")) {
  install.packages("librarian")
}
pkgs <- c(
  "data.table",
  "purrr",
  "sf",
  "matrixStats",
  "targets",
  "tarchetypes",
  "tibble",
  "qs2",
  "crew",
  "tidyr",
  "dplyr",
  "arrow",
  "duckdb",
  "DBI",
  "Matrix"
)


shelf(pkgs)

tar_option_set(
  controller = crew_controller_group(
    local = crew_controller_local(
      workers = 1
    )
  ),
  packages = basename(pkgs),
  format = "qs"
)

# get the path to the data store
get_store <- local({
  f <- source("R/get_store.R", local = TRUE)$value
  function(...) f(...)
})

tar_config_set(store = file.path(get_store(), "targets"))

tar_source("R/assign_trip_ids.R")
# CHANGED: replaces R/mc_group_stats.R. No more matrix-based
# compute_group_stats()/compute_group_raw_long() -- everything downstream of
# the .RData load is melted to long parquet and aggregated in DuckDB.
tar_source("R/mc_duckdb_stats.R")

# ===========================================================================
# DuckDB tuning for the stats stage. Adjust to fit the machine this runs on;
# these are passed straight to compute_stats_duckdb().
#
# IMPORTANT: duckdb_local_tmp_dir is a LOCAL filesystem path, NOT
# store-relative and NOT on the NAS/Samba share. DuckDB's spill-to-disk
# storage writes/reads in ways network filesystems often don't support --
# pointing it at get_store() produced "Bad address" IO errors. This path
# never needs to be portable across machines (it's transient scratch space,
# cleaned up automatically after each target finishes), so it's fine -- and
# necessary -- for it to differ machine to machine.
#
# DO NOT trust tempdir()'s default blindly: on shared servers and LXD
# containers, /tmp is frequently on its own small-quota partition,
# separate from the machine's real local disk capacity -- this is what
# produced the "Disk quota exceeded" error. Before running at scale, on
# EACH machine that will run this pipeline:
#   df -h /tmp                 # check /tmp's own capacity/quota
#   df -h /some/other/local/mount   # find a roomier LOCAL (non-NAS) path
# and set duckdb_local_tmp_dir below to that path explicitly, rather than
# relying on tempdir()'s default. A few x the size of your largest zone's
# melted long-parquet output is a reasonable amount of headroom to look for.
# ===========================================================================
duckdb_local_tmp_dir <- tempdir()  # OVERRIDE this to an explicit, roomy,
# confirmed-local path before running at
# scale -- see guidance above.
duckdb_memory_limit  <- "170GB"  # Set explicitly so DuckDB spills to disk
# (slower) instead of growing unbounded
# (crashes) -- tune to your machine's
# available RAM, leaving headroom for other
# crew workers running concurrently.
duckdb_max_temp_directory_size <- "250GB"  # OVERRIDE to match confirmed free
# space on duckdb_local_tmp_dir's disk (with
# margin). Left unset, DuckDB auto-sizes this
# to whatever's free on that disk at
# connection time -- if that's small, you
# get an opaque "Out of Memory" error citing
# a limit you never chose (this is what
# happened at "1.6 GiB/1.6 GiB used").
# Check with `df -h` on duckdb_local_tmp_dir
# first and set this explicitly below that.
duckdb_threads         <- NULL
# NOTE: no grid batch-size setting needed anymore -- Stage 2 now loops over
# groups one at a time in R (compute_group_stats_matmul in
# R/mc_duckdb_stats.R), which is what actually bounds memory; there's no
# SQL-side batch parameter to tune.


data_ais_layers <- tibble(
  year_month = c("2023-11", paste0("2024-", sprintf("%02d", 4:10)))
  # year_month = c(paste0("2024-", sprintf("%02d", 8:9))) #for small scale testing
)

large_zones    <- c("NEGSL", "SLE", "SGSL")
n_chunks_large <- 5  # tune to fit worker memory; with the matrix multiply gone
# this now only bounds how many source-file iterations
# get melted into one long_dt in RAM at once per target,
# not a downstream dense-matrix size. Worth re-tuning
# down (or to 1) once you've confirmed melt-and-flush
# memory use on one large zone.

data_ais_zones <- expand.grid(
  year_month = data_ais_layers$year_month,
  zones = sub(".*SEGS_([^_]+)_Unaffected.*", "\\1",
              list.files(file.path(get_store(), "Increased_Traffic"),
                         recursive = TRUE, full.names = TRUE,
                         pattern = paste0("resultsIterations_2024-08.*\\.RData"))),
  # zones = c("SM1","SM2","SGSL"), #for small scale testing
  stringsAsFactors = FALSE
) |>
  mutate(
    n_chunks = ifelse(zones %in% large_zones, n_chunks_large, 1L),
    data_ais = rlang::syms(paste0("data_ais_", gsub("-", ".", year_month)))
  ) |>
  uncount(n_chunks, .id = "chunk", .remove = FALSE)

# CHANGED: growth is now applied to ALL trips together -- no per-Type
# stratification. target_growth still crosses per year_month so sampling
# can be mapped by (year_month, tg_label).
data_mc_sampling_values <- tidyr::crossing(
  data_ais_layers,
  target_growth = c(0.05, 0.10)
) |>
  mutate(
    data_ais = rlang::syms(
      paste0("data_ais_", gsub("-", ".", year_month))
    ),
    trip_nums = rlang::syms(
      paste0("trip_nums_", gsub("-", ".", year_month))
    ),
    tg_label = sprintf("%02d", target_growth * 100)
  )

mapped_data_ais_layers <- tar_map(
  values = data_ais_layers,
  tar_target(
    name = data_ais,
    command = {
      ais <- st_read(
        file.path(get_store(), "data", "SEGS_INT_Mortality.gpkg"),
        year_month
      ) |>
        assign_trip_ids(grid, data_ports) |>
        mutate(Type = ifelse(is.na(Type), "OTHER", Type)) |>
        filter(!(GRID_ID %in% remove_grid_ids))

      gaps <- ais |>
        group_by(trip_id) |>
        reframe(geom = st_combine(geom)) |>
        mutate(
          merged_geom  = st_line_merge(geom),
          has_gap      = st_geometry_type(merged_geom) == "MULTILINESTRING",
          internal_gap = if_else(
            has_gap,
            vapply(merged_geom, \(geom) {
              coords <- st_coordinates(geom)
              l1     <- coords[, "L1"]
              ends   <- coords[!duplicated(l1, fromLast = TRUE), c(1, 2), drop = FALSE]
              starts <- coords[!duplicated(l1),                  c(1, 2), drop = FALSE]
              n <- nrow(ends)
              if (n <= 1L) return(0)
              gaps <- vapply(seq_len(n - 1L), \(i) {
                p <- rbind(ends[i, ],      starts[i, ])
                q <- rbind(ends[i + 1L, ], starts[i + 1L, ])
                min(
                  sqrt((p[, 1] - q[1, 1])^2 + (p[, 2] - q[1, 2])^2),
                  sqrt((p[, 1] - q[2, 1])^2 + (p[, 2] - q[2, 2])^2),
                  sqrt((p[1, 1] - q[, 1])^2 + (p[1, 2] - q[, 2])^2),
                  sqrt((p[2, 1] - q[, 1])^2 + (p[2, 2] - q[, 2])^2)
                )
              }, numeric(1))
              max(gaps)
            }, numeric(1)),
            0
          ),
          complete = internal_gap < 10000
        ) |>
        select(-merged_geom, -geom)

      ais |> left_join(gaps, by = "trip_id")
    }
  )
)

# ===========================================================================
# CHANGED: this tar_map() no longer builds a dense data_whaleresults_matrix.
# Each source .RData file is still load()'ed in full (unavoidable -- load()
# always deserializes the whole file), but everything after that is melted
# to long format and flushed to parquet immediately, one source file at a
# time, so peak memory is one file's matrices, not n_chunks x that.
#
# data_whaleresults_long_path and data_whaleskey_path are the two artifacts
# every large zone now writes to disk instead of holding in R:
#   - data_whaleresults_long_path: (whale_iter, event_key, value) rows
#   - data_whaleskey_path:         (event_key, trip_id, Type, zone1, GRID_ID)
# ===========================================================================
mapped_data_ais_zones <- tar_map(
  values = data_ais_zones,
  names  = c(year_month, zones, chunk),

  tar_target(
    name = data_whaleresults_files,
    command = {
      files <- gsub(get_store(), "",
                    list.files(file.path(get_store(), "Increased_Traffic"),
                               recursive = TRUE, full.names = TRUE,
                               pattern = paste0("resultsIterations_", year_month, "INT_SEGS_", zones, ".*\\.RData")))
      tibble::tibble(file = files)
    }
  ),

  # NOTE: returns STORE-RELATIVE paths (e.g. "whale_long/D_2024-08_..._1.parquet"),
  # never anything containing the resolved get_store() root -- so this
  # target's hash/return value is identical whether it was built on the
  # Ubuntu server, a laptop, or wherever the NAS happens to be mounted.
  tar_target(
    name = data_whaleresults_long_path,
    command = {
      if (length(data_whaleresults_files$file) > 0) {
        out_paths <- character(0)

        for (fi in seq_along(data_whaleresults_files$file)) {
          load(file.path(get_store(), data_whaleresults_files$file[fi]))  # 20-45GB transient, unavoidable

          idx <- seq_along(results)
          if (n_chunks > 1) {
            idx <- idx[((idx - 1L) %% n_chunks) == (chunk - 1L)]
          }

          mats <- map(results[idx], ~ map(.x, "MatrixOri")) |>
            list_flatten() |>
            compact()
          rm(results); gc()

          rel_path <- file.path("whale_long",
                                paste0("D_", year_month, "_", zones, "_", chunk, "_", fi, ".parquet"))
          written <- melt_matrices_to_parquet(mats, get_store(), rel_path)
          rm(mats); gc()

          if (!is.na(written)) out_paths <- c(out_paths, written)
        }
        if (length(out_paths) > 0) out_paths else NA_character_
      } else {
        NA_character_
      }
    }
  ),

  # Small (one row per event) -- kept in memory here, same as data_whaleskey
  # was originally, and used by data_ais / trip-count targets downstream.
  # Derives event_key from the parquet just written rather than from
  # colnames() of an in-memory matrix.
  tar_target(
    name = data_whaleskey,
    command = {
      if (length(data_whaleresults_long_path) > 0 && !all(is.na(data_whaleresults_long_path))) {
        con <- DBI::dbConnect(duckdb::duckdb())
        on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
        # resolve store-relative -> full path here, at point of use, on
        # whichever machine is actually running this target
        full_paths <- file.path(get_store(), data_whaleresults_long_path)
        glob <- paste0("['", paste(full_paths, collapse = "','"), "']")
        keys <- DBI::dbGetQuery(con, sprintf(
          "SELECT DISTINCT event_key FROM read_parquet(%s)", glob
        ))

        keys |>
          separate(event_key, into = c("mmsi","UNIX_start","GRID_ID","zone1","zone2","ym","extra"),
                   sep = "_", remove = FALSE) |>
          mutate(uniqID = paste(mmsi, UNIX_start, GRID_ID, sep = "_")) |>
          left_join(data_ais |> select(uniqID, trip_id, complete, Type) |> st_drop_geometry(), by = "uniqID")
      } else {
        NA
      }
    }
  ),

  # NEW: persist the key table to parquet so the final DuckDB stats query
  # doesn't need R at all -- pure disk-to-disk join + aggregation.
  tar_target(
    name = data_whaleskey_path,
    command = {
      if (is.data.frame(data_whaleskey)) {
        rel_path <- file.path("whale_key",
                              paste0("key_", year_month, "_", zones, "_", chunk, ".parquet"))
        out_path <- file.path(get_store(), rel_path)
        dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
        arrow::write_parquet(
          data_whaleskey |> select(event_key, trip_id, Type, zone1, GRID_ID),
          out_path
        )
        rel_path
      } else {
        NA_character_
      }
    }
  )
)

# Recursively flatten a tar_map() result into a flat list of tar_target objects
flatten_targets <- function(x) {
  if (inherits(x, "tar_target")) {
    return(list(x))
  }
  if (is.list(x)) {
    return(do.call(c, lapply(x, flatten_targets)))
  }
  list()
}

# Given a tar_map() result, its "values" tibble, a grouping column, and a
# target-name prefix, build one tar_combine() per group value. Used where
# the combined thing is a data.frame (row-bound).
build_group_combine_targets <- function(mapped_result, values, group_col, name_prefix, out_prefix, extra_cols = character(0)) {
  targets_flat <- flatten_targets(mapped_result)
  target_lookup <- setNames(
    targets_flat,
    vapply(targets_flat, function(x) x$settings$name, character(1))
  )

  lapply(
    unique(values[[group_col]]),
    function(grp) {

      sub <- values[values[[group_col]] == grp, , drop = FALSE]

      target_names <- apply(
        sub[, c(group_col, extra_cols), drop = FALSE],
        1,
        function(row) {
          paste0(
            name_prefix, "_",
            gsub("-", ".", row[[group_col]]),
            if (length(extra_cols) > 0) paste0("_", paste(row[extra_cols], collapse = "_")) else ""
          )
        }
      )

      missing <- setdiff(target_names, names(target_lookup))
      if (length(missing) > 0) {
        stop(
          "build_group_combine_targets: these target names aren't in the mapped result: ",
          paste(missing, collapse = ", ")
        )
      }

      do.call(
        tar_combine,
        c(
          list(name = paste0(out_prefix, "_", gsub("-", ".", grp))),
          target_lookup[target_names],
          list(
            command = quote(
              dplyr::bind_rows(purrr::keep(list(!!!.x), is.data.frame))
            )
          )
        )
      )
    }
  )
}

# NEW: same pattern, but for combining file-PATH targets into one character
# vector per group, instead of row-binding data.frames.
build_group_path_combine_targets <- function(mapped_result, values, group_col, name_prefix, out_prefix, extra_cols = character(0)) {
  targets_flat  <- flatten_targets(mapped_result)
  target_lookup <- setNames(targets_flat, vapply(targets_flat, function(x) x$settings$name, character(1)))

  lapply(unique(values[[group_col]]), function(grp) {
    sub <- values[values[[group_col]] == grp, , drop = FALSE]
    target_names <- apply(sub[, c(group_col, extra_cols), drop = FALSE], 1, function(row) {
      paste0(name_prefix, "_", gsub("-", ".", row[[group_col]]),
             if (length(extra_cols) > 0) paste0("_", paste(row[extra_cols], collapse = "_")) else "")
    })
    missing <- setdiff(target_names, names(target_lookup))
    if (length(missing) > 0) {
      stop("build_group_path_combine_targets: missing targets: ", paste(missing, collapse = ", "))
    }

    do.call(tar_combine, c(
      list(name = paste0(out_prefix, "_", gsub("-", ".", grp))),
      target_lookup[target_names],
      list(command = quote(unlist(list(!!!.x))))
    ))
  })
}

# ===========================================================================
# NOTE: allwhales_targets / data_allwhaleskey, mapped_shared_grid, and
# shared_grid_ids are RETIRED. They existed to identify grid cells shared
# across zones so their stats could be combined separately from "unique"
# grid cells. Grouping by GRID_ID in the DuckDB stats query below produces
# correct cross-zone sums automatically -- a grid cell touched by two zones
# just accumulates rows from both zones' parquet files during Stage 1
# aggregation, with no special-casing needed.
# ===========================================================================

# CHANGED: no more Type stratification -- just total unique trips per
# zone/chunk.
data_trip_counts_values <- data_ais_zones |>
  mutate(
    data_whaleskey = rlang::syms(paste0("data_whaleskey_", gsub("-", ".", year_month), "_", zones, "_", chunk))
  )

mapped_trip_counts <- tar_map(
  values = data_trip_counts_values,
  names  = c(year_month, zones, chunk),
  tar_target(
    name = unique_trip_counts,
    command = {
      if (is.data.frame(data_whaleskey)) {
        tibble::tibble(
          zone         = zones,
          unique_trips = n_distinct(data_whaleskey$trip_id, na.rm = TRUE)
        )
      } else {
        NA
      }
    }
  )
)

combined_trip_counts_targets <- build_group_combine_targets(
  mapped_result = mapped_trip_counts,
  values        = data_trip_counts_values,
  group_col     = "year_month",
  name_prefix   = "unique_trip_counts",
  out_prefix    = "unique_trip_counts_combined",
  extra_cols    = c("zones", "chunk")
)

data_trip_nums_values <- data_ais_layers |>
  mutate(
    unique_trip_counts_combined = rlang::syms(paste0("unique_trip_counts_combined_", gsub("-", ".", year_month)))
  )

# CHANGED: trip_nums is now just total unique trips x {5%, 10%}, no Type
# breakdown. NOTE: this preserves the pre-existing behaviour of summing
# unique_trips across zones (a trip appearing in multiple zones' data gets
# counted once per zone it appears in) -- that was already the case in the
# original pipeline and isn't something this refactor changes.
mapped_trip_nums <- tar_map(
  values = data_trip_nums_values,
  names  = year_month,
  tar_target(
    name = trip_nums,
    command = {
      if (is.data.frame(unique_trip_counts_combined)) {
        total_trips <- sum(unique_trip_counts_combined$unique_trips, na.rm = TRUE)
        tidyr::crossing(
          total_unique_trips = total_trips,
          growth_rate = c(0.05, 0.10)
        ) |>
          mutate(trip_growth_nums = round(total_unique_trips * growth_rate))
      } else {
        NA
      }
    }
  )
)

# CHANGED: mc_sampled_trips now samples uniformly from ALL trips regardless
# of Type -- no more per-Type sample_sizes / split(). trip_lookup keeps
# Type as descriptive metadata (still useful for downstream reporting) but
# it no longer drives sampling.
mapped_mc_sampling <- tar_map(
  values = data_mc_sampling_values,
  names = c(year_month, tg_label),

  tar_target(
    mc_sampled_trips,
    {
      n_iter      <- 10000L
      n_trip_ids  <- nrow(trip_lookup)
      growth_n    <- trip_nums$trip_growth_nums[trip_nums$growth_rate == target_growth]
      sample_size <- min(growth_n, n_trip_ids)

      sampled_indices <- lapply(
        seq_len(n_iter),
        function(i) sample.int(n_trip_ids, size = sample_size)
      )

      Matrix::sparseMatrix(
        i = unlist(sampled_indices),
        j = rep.int(seq_len(n_iter), lengths(sampled_indices)),
        x = 1,
        dims = c(n_trip_ids, n_iter)
      )
    }
  ),

  tar_target(
    trip_lookup,
    {
      as.data.table(data_ais)[
        complete &
          !is.na(trip_id) &
          !is.na(Type),
        .(trip_id, Type)
      ] |>
        unique() |>
        arrange(trip_id) |>
        mutate(trip_index = row_number())
    }
  ),

  # NEW: persist the vessel-iteration trip selection to parquet -- small
  # (trip-level, not event-level), but this is what lets the DuckDB stats
  # query run entirely off disk.
  tar_target(
    mc_sampled_trips_long_path,
    {
      sel <- Matrix::summary(mc_sampled_trips)  # i = trip_index, j = vessel_iter
      dt <- data.table(
        trip_id     = trip_lookup$trip_id[sel$i],
        vessel_iter = sel$j
      )
      rel_path <- file.path("mc_selection",
                            paste0("sel_", year_month, "_", tg_label, ".parquet"))
      out_path <- file.path(get_store(), rel_path)
      dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
      arrow::write_parquet(dt, out_path)
      rel_path
    }
  )
)

# ===========================================================================
# Gather every zone/chunk's D-path and key-path outputs into one vector per
# year_month, so the stats stage can glob across all zones for that month in
# a single DuckDB query.
# ===========================================================================
combined_D_paths_targets <- build_group_path_combine_targets(
  mapped_result = mapped_data_ais_zones,
  values        = data_ais_zones,
  group_col     = "year_month",
  name_prefix   = "data_whaleresults_long_path",
  out_prefix    = "data_whaleresults_long_paths",
  extra_cols    = c("zones", "chunk")
)

combined_key_paths_targets <- build_group_path_combine_targets(
  mapped_result = mapped_data_ais_zones,
  values        = data_ais_zones,
  group_col     = "year_month",
  name_prefix   = "data_whaleskey_path",
  out_prefix    = "data_whaleskey_paths",
  extra_cols    = c("zones", "chunk")
)

# ===========================================================================
# CHANGED: this is what replaces mapped_outputs / mc_results /
# mc_single_zone_stats / mc_unique_grid_stats / mc_grid_raw_path /
# mc_shared_grid_stats_targets / mc_grid_stats_targets / mc_zone_stats_targets
# in their entirety. compute_stats_duckdb() (see R/mc_duckdb_stats.R) is a
# hybrid: DuckDB collapses event-level D to (group, trip, whale_iter)
# totals via one full scan (Stage 1), then R does genuine sparse matrix
# multiplication against the vessel-iteration selection, one group at a
# time (Stage 2) -- SQL joins can't efficiently replicate D %*% S (see the
# header comment in R/mc_duckdb_stats.R for why), so Stage 2 goes back to
# the same linear-algebra approach and memory bound the original R
# implementation used, just sourcing D from disk instead of a pre-built
# dense matrix. Zone stats and grid stats are two calls with a different
# group_col -- both read the same D/key parquet files, just grouped
# differently.
#
# NOTE: trip_lookup and mc_sampled_trips are passed as the actual in-memory
# R objects (not paths) -- S (mc_sampled_trips) is small (trip-level, not
# event-level), so there's no benefit to round-tripping it through parquet
# for this. mc_sampled_trips_long_path (built earlier) is kept in the
# pipeline as a portable/inspectable artifact, but mapped_stats doesn't
# depend on it.
# ===========================================================================
mc_stats_values <- data_mc_sampling_values |>
  select(year_month, target_growth, tg_label) |>
  mutate(
    d_paths_target       = rlang::syms(paste0("data_whaleresults_long_paths_", gsub("-", ".", year_month))),
    key_paths_target     = rlang::syms(paste0("data_whaleskey_paths_", gsub("-", ".", year_month))),
    trip_lookup_target   = rlang::syms(paste0("trip_lookup_", gsub("-", ".", year_month), "_", tg_label)),
    mc_sampled_trips_target = rlang::syms(paste0("mc_sampled_trips_", gsub("-", ".", year_month), "_", tg_label))
  )

mapped_stats <- tar_map(
  values = mc_stats_values,
  names  = c(year_month, tg_label),

  tar_target(
    mc_zone_stats,
    compute_stats_duckdb(
      d_paths_target, key_paths_target,
      group_col = "zone1",
      store_root = get_store(),  # resolved fresh here, on whichever machine runs this
      trip_lookup = trip_lookup_target,
      S = mc_sampled_trips_target,
      local_tmp_dir = duckdb_local_tmp_dir,
      memory_limit = duckdb_memory_limit,
      max_temp_directory_size = duckdb_max_temp_directory_size,
      threads = duckdb_threads
    )
  ),

  tar_target(
    mc_grid_stats,
    compute_stats_duckdb(
      d_paths_target, key_paths_target,
      group_col = "GRID_ID",
      store_root = get_store(),  # resolved fresh here, on whichever machine runs this
      trip_lookup = trip_lookup_target,
      S = mc_sampled_trips_target,
      local_tmp_dir = duckdb_local_tmp_dir,
      memory_limit = duckdb_memory_limit,
      max_temp_directory_size = duckdb_max_temp_directory_size,
      threads = duckdb_threads
    )
  )
)

list(
  tar_target(
    name = grid,
    command = st_read(file.path(get_store(), "data", "Grid.gpkg"), "GridGulf")
  ),

  tar_target(
    name = remove_grid_ids,
    command = c(
      'ATW-1190','ATX-1190','ATV-1189','ATW-1189','ATX-1189','ATV-1188','ATW-1188',
      'ATX-1188','ATY-1188','ATV-1187','ATW-1187','ATX-1187','ATY-1187','ATU-1186',
      'ATV-1186','ATW-1186','ATX-1186','ATY-1186','ATT-1185','ATU-1185','ATV-1185',
      'ATW-1185','ATX-1185','ATY-1185','ATZ-1185','ATT-1184','ATU-1184','ATV-1184',
      'ATW-1184','ATX-1184','ATY-1184','ATZ-1184','ATU-1183','ATV-1183','ATW-1183',
      'ATX-1183','ATY-1183','ATZ-1183','AUA-1183','ATU-1182','ATV-1182','ATW-1182',
      'ATX-1182','ATY-1182','ATZ-1182','AUA-1182','ATU-1181','ATV-1181','ATW-1181',
      'ATX-1181','ATY-1181','ATZ-1181','AUA-1181','AUB-1181','ATU-1180','ATV-1180',
      'ATW-1180','ATX-1180','ATY-1180','ATZ-1180','AUA-1180','AUB-1180','ATU-1179',
      'ATV-1179','ATW-1179','ATX-1179','ATY-1179','ATZ-1179','AUA-1179','AUB-1179',
      'ATU-1178','ATV-1178','ATW-1178','ATX-1178','ATY-1178','ATZ-1178','AUA-1178',
      'AUB-1178','ATU-1177','ATV-1177','ATW-1177','ATX-1177','ATY-1177','ATZ-1177',
      'AUA-1177','AUB-1177','AUC-1177','ATV-1176','ATW-1176','ATX-1176','ATY-1176',
      'ATZ-1176','AUA-1176','AUB-1176','AUC-1176','ATV-1175','ATW-1175','ATX-1175',
      'ATY-1175','ATZ-1175','AUA-1175','AUB-1175','AUC-1175','ATV-1174','ATW-1174',
      'ATX-1174','ATY-1174','ATZ-1174','AUA-1174','AUB-1174','ATU-1173','ATV-1173',
      'ATW-1173','ATX-1173','ATY-1173','ATZ-1173','AUA-1173','AUB-1173','AUC-1173',
      'ATU-1172','ATV-1172','ATW-1172','ATX-1172','ATY-1172','ATZ-1172','AUA-1172',
      'AUB-1172','AUC-1172','ATU-1171','ATV-1171','ATW-1171','ATX-1171','ATY-1171',
      'ATZ-1171','AUA-1171','AUB-1171','AUC-1171','ATV-1170','ATW-1170','ATX-1170',
      'ATY-1170','ATZ-1170','AUA-1170','AUB-1170','AUC-1170','ATU-1169','ATV-1169',
      'ATW-1169','ATX-1169','ATY-1169','ATZ-1169','AUA-1169','AUB-1169','AUC-1169'
    )
  ),

  tar_target(
    name = zone_ids,
    command = {
      path <- file.path(get_store(), "data", "Grid.gpkg")

      layers <- sf::st_layers(path)$name
      read_ids <- function(lyr, path) {
        sf::st_read(path, lyr, quiet = TRUE)$GRID_ID
      }

      list(
        North_Static_Zone = read_ids("Grid_StaticNorth", path),
        South_Static_Zone = c(
          read_ids("Grid_StaticSouthxR", path),
          read_ids("Grid_RestrictedArea", path)
        ),
        DSZA = read_ids("Grid_Dynamic_A", path),
        DSZB = read_ids("Grid_Dynamic_B", path),
        DSZC = read_ids("Grid_Dynamic_C", path),
        DSZD = read_ids("Grid_Dynamic_D", path),
        DSZE = read_ids("Grid_Dynamic_E", path),
        SM1 = read_ids("Grid_SM1", path),
        SM2 = read_ids("Grid_SM2", path),
        RA = read_ids("Grid_RestrictedArea", path),
        CS = read_ids("Grid_VoluntarySlowdownZone", path)
      )
    }
  ),

  tar_target(
    name = data_ports,
    command = st_read(file.path(get_store(), "data", "Grid.gpkg"), "Ports_100mBuffer")
  ),

  mapped_data_ais_layers,
  mapped_data_ais_zones,
  mapped_trip_counts,
  combined_trip_counts_targets,
  mapped_trip_nums,
  mapped_mc_sampling,
  combined_D_paths_targets,
  combined_key_paths_targets,
  mapped_stats
)

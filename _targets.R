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
  "DBI"
)


shelf(pkgs)

tar_option_set(
  # controller = crew_controller_group(
  #   local = crew_controller_local(
  #     workers = 2
  #   )
  # ),
  packages = basename(pkgs),
  format = "qs"
)

# get the path to the data store
source("R/get_store.R")

tar_config_set(store = file.path(get_store(), "targets"))

tar_source("R/assign_trip_ids.R")

data_ais_layers <- tibble(
  year_month = c("2023-11", paste0("2024-", sprintf("%02d", 4:10)))
  # year_month = c(paste0("2024-", sprintf("%02d", 8:9)))
)

data_ais_zones <- expand.grid(
  year_month    = data_ais_layers$year_month,
  zones = sub(".*SEGS_([^_]+)_Unaffected.*", "\\1",list.files(file.path(get_store(),"Increased_Traffic"), recursive = TRUE, full.names = TRUE, pattern = paste0("resultsIterations_2024-08.*\\.RData"))),
  stringsAsFactors = FALSE
) |>
  mutate(
    data_ais = rlang::syms(paste0("data_ais_", gsub("-",".",year_month)))
  )

data_ais_growth <- expand.grid(
  year_month    = data_ais_layers$year_month,
  zones = unique(data_ais_zones$zones),
  target_growth = c(0.05, 0.1),
  stringsAsFactors = FALSE
) |>
  mutate(
    tg_label     = sprintf("%02d", target_growth * 100),  # "05", "10"
    data_allwhaleskey = rlang::syms(paste0("data_allwhaleskey_", gsub("-", ".", year_month))),
    data_whaleskey = rlang::syms(paste0("data_whaleskey_", gsub("-", ".", year_month), "_",zones)),
    data_whaleresults_matrix = rlang::syms(paste0("data_whaleresults_matrix_", gsub("-", ".", year_month),"_",zones)),
    trip_nums = rlang::syms(paste0("trip_nums_", gsub("-", ".", year_month))),
    trip_lookup = rlang::syms(
      paste0(
        "trip_lookup_",
        gsub("-", ".", year_month),
        "_",
        tg_label
      )
    ),
    mc_sampled_trips = rlang::syms(
      paste0(
        "mc_sampled_trips_",
        gsub("-", ".", year_month),
        "_",
        tg_label
      )
    ),
    # NEW: each zone branch gets a reference to its month's shared-grid-id vector
    shared_grid_ids = rlang::syms(paste0("shared_grid_ids_", gsub("-", ".", year_month)))
  )

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

mapped_data_ais_zones <- tar_map(
  values = data_ais_zones,
  names  = c(year_month, zones),
  tar_target(
    name = data_whaleresults_files,
    command = {
      files <- gsub(get_store(),"",list.files(file.path(get_store(),"Increased_Traffic"), recursive = TRUE, full.names = TRUE, pattern = paste0("resultsIterations_",year_month,"INT_SEGS_",zones,".*\\.RData")))

      tibble::tibble(
        file = files
      )
    }
  ),

  tar_target(
    name = data_whaleresults_matrix,
    command = {
      if(length(data_whaleresults_files$file)>0){
        all_mats <- map(data_whaleresults_files$file, function(file) {
          load(file.path(get_store(),file))
          map(results, ~ map(.x, "MatrixOri")) |>
            list_flatten() |>
            compact()
        }) |>
          list_flatten()

        combined <- do.call(cbind, all_mats)
      } else {
        NA
      }
    }
  ),

  tar_target(
    name = data_whaleskey,
    command = {
      if(is.matrix(data_whaleresults_matrix)){
        data.frame(
          col_name = colnames(data_whaleresults_matrix)
        ) |>
          separate(col_name, into = c("mmsi","UNIX_start","GRID_ID","zone1","zone2","ym","extra"),
                   sep = "_", remove = FALSE) |>
          mutate(uniqID = paste(mmsi, UNIX_start, GRID_ID, sep = "_")) |>
          left_join(data_ais |> select(uniqID,trip_id,complete,Type) |> st_drop_geometry(), by = "uniqID")
      } else {
        NA
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
# target-name prefix, build one tar_combine() per group value.
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

# NEW: combine shared-grid-cell raw parquet files across zones via duckdb
combine_grid_raw_duckdb <- function(paths) {
  paths <- paths[!is.na(paths)]
  if (length(paths) == 0) return(NA)

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  glob <- paste0("['", paste(paths, collapse = "','"), "']")

  DBI::dbGetQuery(con, sprintf(
    "SELECT GRID_ID, vessel_iter, iteration, SUM(value) AS value
     FROM read_parquet(%s)
     GROUP BY GRID_ID, whale_iter",
    glob
  ))
}

allwhales_targets <- build_group_combine_targets(
  mapped_result = mapped_data_ais_zones,
  values        = data_ais_zones,
  group_col     = "year_month",
  name_prefix   = "data_whaleskey",
  out_prefix    = "data_allwhaleskey",
  extra_cols    = "zones"
)

# NEW: identify grid cells that appear in more than one zone, per year_month
data_shared_grid_values <- data_ais_layers |>
  mutate(
    data_allwhaleskey = rlang::syms(paste0("data_allwhaleskey_", gsub("-", ".", year_month)))
  )

mapped_shared_grid <- tar_map(
  values = data_shared_grid_values,
  names  = year_month,
  tar_target(
    shared_grid_ids,
    command = {
      as.data.table(data_allwhaleskey)[
        , .(n_zones = uniqueN(zone1)), by = GRID_ID
      ][n_zones > 1, GRID_ID]
    }
  )
)

data_trip_counts_values <- data_ais_zones |>
  mutate(
    data_whaleskey = rlang::syms(paste0("data_whaleskey_", gsub("-", ".", year_month), "_", zones))
  )

mapped_trip_counts <- tar_map(
  values = data_trip_counts_values,
  names  = c(year_month, zones),
  tar_target(
    name = unique_trip_counts,
    command = {
      if (is.data.frame(data_whaleskey)) {
        data_whaleskey |>
          group_by(Type) |>
          summarise(unique_trips = n_distinct(trip_id, na.rm = TRUE)) |>
          mutate(zone = zones) |>
          as.data.frame()
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
  extra_cols    = "zones"
)

data_trip_nums_values <- data_ais_layers |>
  mutate(
    unique_trip_counts_combined = rlang::syms(paste0("unique_trip_counts_combined_", gsub("-", ".", year_month)))
  )

mapped_trip_nums <- tar_map(
  values = data_trip_nums_values,
  names  = year_month,
  tar_target(
    name = trip_nums,
    command = {
      if (is.data.frame(unique_trip_counts_combined)) {
        tidyr::crossing(
          unique_trip_counts_combined,
          growth_rate = target_growth
        ) |>
          group_by(Type,growth_rate) |>
          reframe(unique_trips =sum(unique_trips ,na.rm = TRUE)) |>
          mutate(
            trip_growth_nums = round(unique_trips * growth_rate)
          )
      } else {
        NA
      }
    }
  )
)

mapped_mc_sampling <- tar_map(
  values = data_mc_sampling_values,
  names = c(year_month, tg_label),

  tar_target(
    mc_sampled_trips,
    {
      n_iter <- 10000L

      trip_ids_by_type <- split(
        trip_lookup$trip_index,
        trip_lookup$Type
      )

      trip_nums_filtered <- trip_nums |>
        filter(growth_rate == target_growth,
               !is.na(Type))

      sample_sizes <- setNames(
        trip_nums_filtered$trip_growth_nums,
        trip_nums_filtered$Type
      )

      trip_ids_by_type <- trip_ids_by_type[trip_nums_filtered$Type]

      sampled_indices <- lapply(
        seq_len(n_iter),
        function(i) {
          unlist(
            mapply(
              function(ids, n) {
                sample(
                  ids,
                  size = min(n, length(ids))
                )
              },
              trip_ids_by_type,
              sample_sizes[names(trip_ids_by_type)],
              SIMPLIFY = FALSE
            ),
            use.names = FALSE
          )
        }
      )

      Matrix::sparseMatrix(
        i = unlist(sampled_indices),
        j = rep.int(
          seq_len(n_iter),
          lengths(sampled_indices)
        ),
        x = TRUE,
        dims = c(
          nrow(trip_lookup),
          n_iter
        )
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
  )
)

# computes zone stats, finalizes unique (non-shared) grid cell stats
# in-memory, and writes only shared grid cells' raw values to parquet.
mapped_outputs <- tar_map(
  values = data_ais_growth,
  names  = c(year_month, tg_label, zones),

  tar_target(
    mc_sampled_mask,
    {
      if (is.data.frame(data_whaleskey)) {

        whales_dt <- as.data.table(data_whaleskey)

        stopifnot(
          nrow(whales_dt) ==
            ncol(data_whaleresults_matrix)
        )

        trip_to_rows <- split(
          seq_len(nrow(whales_dt)),
          whales_dt$trip_id
        )

        selected <- Matrix::summary(mc_sampled_trips)

        trip_ids_selected <- trip_lookup$trip_id[selected$i]

        rows_selected <- trip_to_rows[trip_ids_selected]

        row_lengths <- lengths(rows_selected)

        keep <- row_lengths > 0

        i_idx <- unlist(
          rows_selected[keep],
          use.names = FALSE
        )

        j_idx <- rep(
          selected$j[keep],
          row_lengths[keep]
        )

        stopifnot(length(i_idx) == length(j_idx))

        Matrix::sparseMatrix(
          i = i_idx,
          j = j_idx,
          x = 1,
          dims = c(
            nrow(whales_dt),
            ncol(mc_sampled_trips)
          )
        )
      } else {
        NA
      }
    }
  ),

  tar_target(
    mc_results,
    command = {
      if (is.data.frame(data_whaleskey)) {
        # mask <- as.matrix(mc_sampled_mask)
        idx <- Matrix::summary(mc_sampled_mask)

        zone_groups <- as.factor(data_whaleskey$zone1)
        grid_groups <- as.factor(data_whaleskey$GRID_ID)
        grid_id_levels <- levels(grid_groups)
        is_shared <- grid_id_levels %in% shared_grid_ids

        n_iter <- nrow(data_whaleresults_matrix)

        zone_results_list    <- vector("list", n_iter)
        unique_grid_raw_list <- vector("list", n_iter)
        shared_grid_raw_list <- vector("list", n_iter)

        for (i in seq_len(n_iter)) {
          if (i %% floor(n_iter / 10L) == 0L) {
            message("iteration ", i, " of ", n_iter)
          }

          # maskedvessels <- mask * data_whaleresults_matrix[i, ]
          maskedvessels <- matrix(NA_real_,
                                  nrow = nrow(mc_sampled_mask),
                                  ncol = ncol(mc_sampled_mask))

          maskedvessels[cbind(idx$i, idx$j)] <- data_whaleresults_matrix[i, ][idx$i]

          # zone stats -- unchanged logic, computed immediately per iteration
          vesselsbygroup <- rowsum(maskedvessels, group = zone_groups, na.rm = TRUE)
          zone_results_list[[i]] <- tibble(
            whale_iter = i,
            group = rownames(vesselsbygroup),
            vessel_iter_delta_mean = matrixStats::rowMeans2(vesselsbygroup, na.rm = TRUE),
            vessel_iter_delta_median = matrixStats::rowMedians(vesselsbygroup, na.rm = TRUE),
            vessel_iter_delta_variance = matrixStats::rowVars(vesselsbygroup, na.rm = TRUE),
            vessel_iter_delta_q025 = matrixStats::rowQuantiles(vesselsbygroup, probs = 0.025, na.rm = TRUE),
            vessel_iter_delta_q975 = matrixStats::rowQuantiles(vesselsbygroup, probs = 0.975, na.rm = TRUE),
            vessel_iter_delta_n_delta_gt0 = matrixStats::rowCounts(vesselsbygroup > 0, value = TRUE, na.rm = TRUE)
          )

          # grid -- one rowsum, then split by shared vs. unique
          gridsum <- rowsum(maskedvessels, group = grid_groups, na.rm = TRUE)

          unique_grid_raw_list[[i]] <- gridsum[!is_shared, , drop = FALSE]

          # if (any(is_shared)) {
          #   shared_grid_raw_list[[i]] <- tibble(
          #     GRID_ID   = sub("^grid_", "", rownames(gridsum)[is_shared]),
          #     iteration = i,
          #     value     = gridsum[is_shared, 1]
          #   )
          # }

          if (any(is_shared)) {
            shared_grid_raw_list[[i]] <- data.table::as.data.table(
              gridsum[is_shared, , drop = FALSE],
              keep.rownames = "GRID_ID"
            ) |>
              tidyr::pivot_longer(
                cols = -GRID_ID,
                names_to = "vessel_iter",
                values_to = "value"
              ) |>
              mutate(
                vessel_iter = as.numeric(sub("V", "", vessel_iter)),
                iteration = i # whale iterations
              )
          }
        }

        zone_stats <- data.table::rbindlist(zone_results_list) |>
          mutate(type = sub("^zone_", "", group)) |>
          select(-group)

        # unique grid cells: build matrix across iterations, stats now, no disk round-trip
        unique_mat <- do.call(cbind, unique_grid_raw_list)
        unique_grid_stats <- if (!is.null(unique_mat) && nrow(unique_mat) > 0) {
          tibble(
            GRID_ID = sub("^grid_", "", rownames(unique_mat)),
            vessel_iter_delta_mean = matrixStats::rowMeans2(unique_mat, na.rm = TRUE),
            vessel_iter_delta_median = matrixStats::rowMedians(unique_mat, na.rm = TRUE),
            vessel_iter_delta_variance = matrixStats::rowVars(unique_mat, na.rm = TRUE),
            vessel_iter_delta_q025 = matrixStats::rowQuantiles(unique_mat, probs = 0.025, na.rm = TRUE),
            vessel_iter_delta_q975 = matrixStats::rowQuantiles(unique_mat, probs = 0.975, na.rm = TRUE),
            vessel_iter_delta_n_delta_gt0 = matrixStats::rowCounts(unique_mat > 0, value = TRUE, na.rm = TRUE)
          )
        } else {
          NA
        }

        # shared grid cells only: write raw values to parquet for later duckdb combine
        shared_raw_long <- data.table::rbindlist(shared_grid_raw_list)
        grid_raw_path <- NA_character_
        if (nrow(shared_raw_long) > 0) {
          out_path <- file.path(
            get_store(), "grid_raw",
            paste0("mc_grid_raw_", year_month, "_", tg_label, "_", zones, ".parquet")
          )
          tmp_path <- paste0(out_path, ".tmp")
          dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
          arrow::write_parquet(shared_raw_long, tmp_path)
          file.rename(tmp_path, out_path)
          grid_raw_path <- out_path
        }

        list(
          zone_stats = zone_stats,
          unique_grid_stats = unique_grid_stats,
          grid_raw_path = grid_raw_path
        )

      } else {
        list(zone_stats = NA, unique_grid_stats = NA, grid_raw_path = NA_character_)
      }
    }
  ),

  # thin extraction targets -- no recomputation, just indexing the cached list
  tar_target(mc_zone_stats, mc_results$zone_stats),
  tar_target(mc_unique_grid_stats, mc_results$unique_grid_stats),
  tar_target(mc_grid_raw_path, mc_results$grid_raw_path)
)

# NEW: combine unique (non-shared) grid stats across zones, per year_month/tg_label
# -- these are already final, just need bind_rows across zones.
data_grid_combine_values <- data_ais_growth |>
  mutate(ym_tg = paste0(gsub("-", ".", year_month), "_", tg_label))

mc_unique_grid_stats_targets <- flatten_targets(mapped_outputs)
mc_unique_grid_stats_lookup <- setNames(
  mc_unique_grid_stats_targets,
  vapply(mc_unique_grid_stats_targets, function(x) x$settings$name, character(1))
)

mc_unique_grid_combined_targets <- lapply(
  unique(data_grid_combine_values$ym_tg),
  function(grp) {
    sub <- data_grid_combine_values[data_grid_combine_values$ym_tg == grp, ]
    target_names <- paste0("mc_unique_grid_stats_", grp, "_", sub$zones)

    do.call(
      tar_combine,
      c(
        list(name = paste0("mc_unique_grid_combined_", grp)),
        mc_unique_grid_stats_lookup[target_names],
        list(command = quote(dplyr::bind_rows(purrr::keep(list(!!!.x), is.data.frame))))
      )
    )
  }
)

# NEW: combine shared grid cells' raw parquet files across zones via duckdb,
# per year_month/tg_label
mc_grid_raw_path_lookup <- setNames(
  mc_unique_grid_stats_targets,
  vapply(mc_unique_grid_stats_targets, function(x) x$settings$name, character(1))
)

mc_grid_combined_targets <- lapply(
  unique(data_grid_combine_values$ym_tg),
  function(grp) {
    sub <- data_grid_combine_values[data_grid_combine_values$ym_tg == grp, ]
    target_names <- paste0("mc_grid_raw_path_", grp, "_", sub$zones)

    do.call(
      tar_combine,
      c(
        list(name = paste0("mc_grid_combined_", grp)),
        mc_grid_raw_path_lookup[target_names],
        list(command = quote(combine_grid_raw_duckdb(c(!!!.x))))
      )
    )
  }
)

# NEW: final shared-grid stats -- pivot the combined long table to a matrix,
# run matrixStats once
mc_shared_grid_stats_targets <- lapply(
  unique(data_grid_combine_values$ym_tg),
  function(grp) {
    tar_target_raw(
      name = paste0("mc_shared_grid_stats_", grp),
      command = bquote({
        long <- .(as.name(paste0("mc_grid_combined_", grp)))
        if (is.data.frame(long) && nrow(long) > 0) {
          wide <- tidyr::pivot_wider(long, names_from = iteration, values_from = value, values_fill = 0)
          grid_mat <- as.matrix(wide[, -1])
          rownames(grid_mat) <- wide$GRID_ID

          tibble(
            GRID_ID = rownames(grid_mat),
            vessel_iter_delta_mean        = matrixStats::rowMeans2(grid_mat, na.rm = TRUE),
            vessel_iter_delta_median      = matrixStats::rowMedians(grid_mat, na.rm = TRUE),
            vessel_iter_delta_variance    = matrixStats::rowVars(grid_mat, na.rm = TRUE),
            vessel_iter_delta_q025        = matrixStats::rowQuantiles(grid_mat, probs = 0.025, na.rm = TRUE),
            vessel_iter_delta_q975        = matrixStats::rowQuantiles(grid_mat, probs = 0.975, na.rm = TRUE),
            vessel_iter_delta_n_delta_gt0 = matrixStats::rowCounts(grid_mat > 0, value = TRUE, na.rm = TRUE)
          )
        } else {
          NA
        }
      })
    )
  }
)

# NEW: final grid stats = union of unique-cell stats (already final) and
# shared-cell stats (combined then finalized)
mc_grid_stats_targets <- lapply(
  unique(data_grid_combine_values$ym_tg),
  function(grp) {
    tar_target_raw(
      name = paste0("mc_grid_stats_", grp),
      command = bquote({
        unique_part <- .(as.name(paste0("mc_unique_grid_combined_", grp)))
        shared_part <- .(as.name(paste0("mc_shared_grid_stats_", grp)))
        dplyr::bind_rows(
          purrr::keep(list(unique_part, shared_part), is.data.frame)
        )
      })
    )
  }
)

mc_zone_stats_targets <- lapply(
  unique(data_grid_combine_values$ym_tg),
  function(grp) {
    sub <- data_grid_combine_values[data_grid_combine_values$ym_tg == grp, ]
    target_names <- paste0("mc_results_", grp, "_", sub$zones)

    do.call(
      tar_combine,
      c(
        list(name = paste0("mc_zone_stats_", grp)),
        mc_unique_grid_stats_lookup[target_names],
        list(command = quote(dplyr::bind_rows(purrr::keep(list(!!!.x), is.data.frame))))
      )
    )
  }
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

  tar_target(
    name = target_growth,
    command = c(0.05, 0.10) # this should only be 2 values or will break downstream code
  ),

  mapped_data_ais_layers,
  mapped_data_ais_zones,
  allwhales_targets,
  mapped_shared_grid,               # NEW: identifies multi-zone grid cells per month
  mapped_trip_counts,
  combined_trip_counts_targets,
  mapped_trip_nums,
  mapped_mc_sampling,
  mapped_outputs,                   # now produces mc_results -> mc_zone_stats, mc_unique_grid_stats, mc_grid_raw_path
  mc_unique_grid_combined_targets,  # NEW: bind_rows of already-final unique-cell stats
  mc_grid_combined_targets,         # NEW: duckdb-combined shared-cell raw values
  mc_shared_grid_stats_targets,     # NEW: matrixStats on combined shared cells
  mc_grid_stats_targets,            # NEW: final union of unique + shared grid stats
  mc_zone_stats_targets
)

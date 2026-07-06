if (!require(librarian)) {
  install.packages("librarian")
}
pkgs <- c(
  "data.table",
  "purrr",
  "sf",
  "targets",
  "tarchetypes",
  "tibble",
  "qs2",
  "tidyr",
  "dplyr"
)


shelf(pkgs)


tar_option_set(
  packages = basename(pkgs),
  format = "qs"
)

# get the path to the data store
if(dir.exists("//ci-WPNSBIO9039519-smb-1.mar.dfo-mpo.ca/ocean_data/")) {
  # on the NAS
  store = "//ci-WPNSBIO9039519-smb-1.mar.dfo-mpo.ca/ocean_data/SPA/NARW_AIS"
  } else if (dir.exists("/srv/sambashare/NARW")) {
  # on a linux machine on the sambashare
  store = "/srv/sambashare/NARW"
} else if (
  dir.exists(
    "//wpnsbio9039519.mar.dfo-mpo.ca/sambashare/NARW"
  )
) {
  # on a windows machine on the sambashare
  store <- "//wpnsbio9039519.mar.dfo-mpo.ca/sambashare/NARW"
} else {
  store <-  getwd()
}

tar_config_set(store = file.path(store, "targets"))

tar_source("R/assign_trip_ids.R")

data_ais_layers <- tibble(
  # year_month = c("2023-11", paste0("2024-", sprintf("%02d", 4:10)))
  year_month = c(paste0("2024-", sprintf("%02d", 8:9)))
)

data_ais_growth <- expand.grid(
  year_month    = data_ais_layers$year_month,
  target_growth = c(0.05, 0.1),
  stringsAsFactors = FALSE
) |>
  mutate(
    tg_label     = sprintf("%02d", target_growth * 100),  # "05", "10"
    # data_ais = rlang::syms(paste0("data_ais_", gsub("-",".",year_month))),
    data_whaleskey = rlang::syms(paste0("data_whaleskey_", gsub("-", ".", year_month))),
    data_whaleresults_matrix = rlang::syms(paste0("data_whaleresults_matrix_", gsub("-", ".", year_month)))
  )

mapped_data_ais_layers <- tar_map(
  values = data_ais_layers,
  tar_target(
    name = data_ais,
    command = {
      ais <- st_read(
        file.path(store, "data", "SEGS_INT_Mortality.gpkg"),
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
  ),


  tar_target(
    name = data_whaleresults_files,
    command = {
      files <- list.files(file.path(store,"Increased_Traffic"), recursive = TRUE, full.names = TRUE, pattern = paste0("resultsIterations_",year_month,"INT_SEGS_.*\\.RData"))

      tibble::tibble(
        file = files,
        size = file.info(files)$size
      ) |>
        filter(size <1000000000) #TODO remove 1gb filter

      }
  ),

  tar_target(
    name = data_whaleresults_matrix,
    command = {
      if(length(data_whaleresults_files$file)>0){
        all_mats <- map(data_whaleresults_files$file, function(file) {
          load(file)
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

mapped_outputs <- tar_map(
  values = data_ais_growth,
  names  = c(year_month, tg_label),

  tar_target(
    name = trip_nums,
    command = {
      if(is.data.frame(data_whaleskey)){
        data_whaleskey |>
          group_by(Type) |>
          summarise(unique_trips = n_distinct(trip_id, na.rm = TRUE)) |>
          as.data.frame() |>
          mutate(
            growth_rate        = target_growth,
            trip_growth_nums   = round(unique_trips * target_growth)
          )
      } else {
        NA
      }
    }
  ),


  tar_target(
    name = mc_sampled_mask,
    command = {
      if (is.data.frame(data_whaleskey)) {

        n_iter <- 10000L
        data_whaleskey_dt <- as.data.table(data_whaleskey)
        stopifnot(nrow(data_whaleskey_dt) == ncol(data_whaleresults_matrix))

        complete_dt <- data_whaleskey_dt[complete & !is.na(trip_id) & !is.na(Type)]
        trip_lookup <- unique(complete_dt, by = c("Type", "trip_id"))

        # unique trip_id per Type -- sampling must be at trip level, not
        # AIS-record level, or long trips get oversampled
        trip_ids_by_type <- split(trip_lookup$trip_id, trip_lookup$Type)

        sample_sizes <- setNames(trip_nums$trip_growth_nums, trip_nums$Type)

        # precompute AIS-row indices per trip (columns of data_whaleresults_matrix)
        ais_rows_by_trip <- split(seq_len(nrow(data_whaleskey_dt)), data_whaleskey_dt$trip_id)

        sampled_ais_rows <- lapply(seq_len(n_iter), function(i) {
          sampled_trips <- unlist(mapply(
            function(ids, n) sample(ids, size = min(n, length(ids))),
            trip_ids_by_type,
            sample_sizes[names(trip_ids_by_type)],
            SIMPLIFY = FALSE
          ), use.names = FALSE)
          unlist(ais_rows_by_trip[sampled_trips], use.names = FALSE)
        })

        Matrix::sparseMatrix(
          i = unlist(sampled_ais_rows, use.names = FALSE),
          j = rep.int(seq_len(n_iter), lengths(sampled_ais_rows)),
          x = 1,
          dims = c(nrow(data_whaleskey_dt), n_iter)
        )

      } else {
        NA
      }
    }
  ),


  tar_target(
    mc_results,
    command = {
      mask <- as.matrix(mc_sampled_mask)
      zones <- as.factor(paste0("zone_",data_whaleskey$zone1))
      gridids <- as.factor(paste0("grid_",data_whaleskey$GRID_ID))

      n_iter <- nrow(data_whaleresults_matrix)

      results_list <- vector("list", n_iter)

      for(i in seq_len(n_iter)) {
                if (i %% floor(n_iter / 10L) == 0L) {
                  message(
                    "iteration ", i, " of ", n_iter
                  )
                }

        maskedvessels <- mask * data_whaleresults_matrix[i]


        vesselsbygroup <- rbind(rowsum(maskedvessels, group = zones, na.rm = TRUE),
                                rowsum(maskedvessels, group = gridids, na.rm = TRUE))

        results_list[[i]] <- tibble(
          whale_iter = i,
          group = rownames(vesselsbygroup),

          vessel_iter_delta_mean = matrixStats::rowMeans2(
            vesselsbygroup,
            na.rm = TRUE
          ),

          vessel_iter_delta_median = matrixStats::rowMedians(
            vesselsbygroup,
            na.rm = TRUE
          ),

          vessel_iter_delta_variance = matrixStats::rowVars(
            vesselsbygroup,
            na.rm = TRUE
          ),

          vessel_iter_delta_q025 = matrixStats::rowQuantiles(
            vesselsbygroup,
            probs = 0.025,
            na.rm = TRUE
          ),

          vessel_iter_delta_q975 = matrixStats::rowQuantiles(
            vesselsbygroup,
            probs = 0.975,
            na.rm = TRUE
          ),

          vessel_iter_delta_n_delta_gt0 = matrixStats::rowCounts(
            vesselsbygroup > 0,
            value = TRUE,
            na.rm = TRUE
          )
        ) |>
          separate(group, into = c("type", "id"), sep = "_", remove = TRUE)
      }

      data.table::rbindlist(results_list)

    }
  )
)


list(
  tar_target(
    name = grid,
    command = st_read(file.path(store, "data", "Grid.gpkg"), "GridGulf")
  ),

  tar_target(
    name = remove_grid_ids,
    command = c(
      'ATW-1190',
      'ATX-1190',
      'ATV-1189',
      'ATW-1189',
      'ATX-1189',
      'ATV-1188',
      'ATW-1188',
      'ATX-1188',
      'ATY-1188',
      'ATV-1187',
      'ATW-1187',
      'ATX-1187',
      'ATY-1187',
      'ATU-1186',
      'ATV-1186',
      'ATW-1186',
      'ATX-1186',
      'ATY-1186',
      'ATT-1185',
      'ATU-1185',
      'ATV-1185',
      'ATW-1185',
      'ATX-1185',
      'ATY-1185',
      'ATZ-1185',
      'ATT-1184',
      'ATU-1184',
      'ATV-1184',
      'ATW-1184',
      'ATX-1184',
      'ATY-1184',
      'ATZ-1184',
      'ATU-1183',
      'ATV-1183',
      'ATW-1183',
      'ATX-1183',
      'ATY-1183',
      'ATZ-1183',
      'AUA-1183',
      'ATU-1182',
      'ATV-1182',
      'ATW-1182',
      'ATX-1182',
      'ATY-1182',
      'ATZ-1182',
      'AUA-1182',
      'ATU-1181',
      'ATV-1181',
      'ATW-1181',
      'ATX-1181',
      'ATY-1181',
      'ATZ-1181',
      'AUA-1181',
      'AUB-1181',
      'ATU-1180',
      'ATV-1180',
      'ATW-1180',
      'ATX-1180',
      'ATY-1180',
      'ATZ-1180',
      'AUA-1180',
      'AUB-1180',
      'ATU-1179',
      'ATV-1179',
      'ATW-1179',
      'ATX-1179',
      'ATY-1179',
      'ATZ-1179',
      'AUA-1179',
      'AUB-1179',
      'ATU-1178',
      'ATV-1178',
      'ATW-1178',
      'ATX-1178',
      'ATY-1178',
      'ATZ-1178',
      'AUA-1178',
      'AUB-1178',
      'ATU-1177',
      'ATV-1177',
      'ATW-1177',
      'ATX-1177',
      'ATY-1177',
      'ATZ-1177',
      'AUA-1177',
      'AUB-1177',
      'AUC-1177',
      'ATV-1176',
      'ATW-1176',
      'ATX-1176',
      'ATY-1176',
      'ATZ-1176',
      'AUA-1176',
      'AUB-1176',
      'AUC-1176',
      'ATV-1175',
      'ATW-1175',
      'ATX-1175',
      'ATY-1175',
      'ATZ-1175',
      'AUA-1175',
      'AUB-1175',
      'AUC-1175',
      'ATV-1174',
      'ATW-1174',
      'ATX-1174',
      'ATY-1174',
      'ATZ-1174',
      'AUA-1174',
      'AUB-1174',
      'ATU-1173',
      'ATV-1173',
      'ATW-1173',
      'ATX-1173',
      'ATY-1173',
      'ATZ-1173',
      'AUA-1173',
      'AUB-1173',
      'AUC-1173',
      'ATU-1172',
      'ATV-1172',
      'ATW-1172',
      'ATX-1172',
      'ATY-1172',
      'ATZ-1172',
      'AUA-1172',
      'AUB-1172',
      'AUC-1172',
      'ATU-1171',
      'ATV-1171',
      'ATW-1171',
      'ATX-1171',
      'ATY-1171',
      'ATZ-1171',
      'AUA-1171',
      'AUB-1171',
      'AUC-1171',
      'ATV-1170',
      'ATW-1170',
      'ATX-1170',
      'ATY-1170',
      'ATZ-1170',
      'AUA-1170',
      'AUB-1170',
      'AUC-1170',
      'ATU-1169',
      'ATV-1169',
      'ATW-1169',
      'ATX-1169',
      'ATY-1169',
      'ATZ-1169',
      'AUA-1169',
      'AUB-1169',
      'AUC-1169'
    )
  ),

  tar_target(
    name = zone_ids,
    command = {
      path <- file.path(store, "data", "Grid.gpkg")

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
    command = st_read(file.path(store, "data", "Grid.gpkg"), "Ports_100mBuffer")
  ),

  tar_target(
    name = target_growth,
    command = c(0.05, 0.10) # this should only be 2 values or will break downstream code
  ),

  mapped_data_ais_layers,
  mapped_outputs
)

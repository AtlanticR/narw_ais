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
if (dir.exists("/srv/sambashare/NARW")) {
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
  year_month = c("2023-09","2023-11", paste0("2024-", sprintf("%02d", 4:10)))
)

data_ais_growth <- expand.grid(
  year_month    = data_ais_layers$year_month,
  target_growth = c(0.05, 0.1),
  stringsAsFactors = FALSE
) |>
  mutate(
    tg_label     = sprintf("%02d", target_growth * 100),  # "05", "10"
    # data_ais = rlang::syms(paste0("data_ais_", gsub("-",".",year_month))),
    data_whalestats = rlang::syms(paste0("data_whalestats_", gsub("-", ".", year_month)))
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
    command = list.files("data", full.names = TRUE, pattern = paste0("resultsIterations_",year_month,".*ffected\\.RData")),
    format = "file"
  ),

  tar_target(
    name = data_whaleresults_matrix,
    command = {
      if(length(data_whaleresults_files)>0){
        all_mats <- map(data_whaleresults_files, function(file) {
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
    name = data_whalestats,
    command = {
      if(is.matrix(data_whaleresults_matrix)){
        data.frame(
          col_name = colnames(data_whaleresults_matrix),
          MortMedian = as.numeric(apply(data_whaleresults_matrix, 2, median)),
          MortVar = as.numeric(apply(data_whaleresults_matrix, 2, var))
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
      if(is.data.frame(data_whalestats)){
        data_whalestats |>
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
    name = mc_sampled_uniqID,
    command = {
      if(is.data.frame(data_whalestats)){
        n_iter    <- 100L  # TODO: change to 10000
        Grid_IDs  <- unique(data_whalestats$GRID_ID)
        data_whalestats_dt <- as.data.table(data_whalestats)

        trip_ids_by_type <- split(
          data_whalestats_dt$trip_id[data_whalestats_dt$complete],
          data_whalestats_dt$Type[data_whalestats_dt$complete]
        )

        sample_sizes <- setNames(trip_nums$trip_growth_nums, trip_nums$Type)


        sampled_trips <- lapply(seq_len(n_iter), function(i) {
          unlist(mapply(
            sample,
            trip_ids_by_type,
            sample_sizes[names(trip_ids_by_type)],
            SIMPLIFY = FALSE
          ))

        })

      } else {
        NA
      }

    }
  ),

  tar_target(
    name = mc_result_matrix,
    command = {
      if(is.data.frame(data_whalestats)){
        Grid_IDs    <- unique(data_whalestats$GRID_ID)
        data_whalestats_dt <- as.data.table(st_drop_geometry(data_whalestats))
        n_iter      <- length(mc_sampled_uniqID)

        result_matrix <- matrix(
          NA_real_,
          nrow     = length(Grid_IDs),
          ncol     = n_iter,
          dimnames = list(Grid_IDs = Grid_IDs, n_iter = seq_len(n_iter))
        )
        for (i in seq_len(n_iter)) {
          if (i %% (n_iter / 10) == 0)
            message("MC iteration: ", i, " of ", n_iter,
                    " | growth: ", target_growth,
                    " | ", year_month)

          temp <- data_whalestats_dt[
            trip_id %in% mc_sampled_uniqID[[i]],
            .(MortMedian = sum(MortMedian, na.rm = TRUE)),
            by = GRID_ID
          ]
          result_matrix[temp$GRID_ID, i] <- temp$MortMedian
        }
        result_matrix
      } else {
        NA
      }


    }
  ),

  tar_target(
    name = mc_stats,
    command = {
      if(is.data.frame(data_whalestats)){
        # na.rm = TRUE handles grid cells with no vessel traffic in some iterations
        list(
          median = apply(mc_result_matrix, 1, median, na.rm = TRUE),
          mean   = apply(mc_result_matrix, 1, mean,   na.rm = TRUE),
          p10    = apply(mc_result_matrix, 1, quantile, probs = 0.10, na.rm = TRUE),
          p90    = apply(mc_result_matrix, 1, quantile, probs = 0.90, na.rm = TRUE),
          cv     = apply(mc_result_matrix, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
        )
      } else {
        NA
      }

    }
  ),


  tar_target(
    name = Grid_sums,
    command = {
      grid_base <- as.data.table(st_drop_geometry(data_ais))[,
                                                             .(
                                                               MortMedian = sum(MortMedian, na.rm = TRUE),
                                                               MortVar    = sum(MortVar,    na.rm = TRUE)
                                                             ),
                                                             by = GRID_ID
      ]

      data.frame(
        GRID_ID          = rownames(mc_stats),
        MortMedian_MC    = apply(mc_stats, 1, median, na.rm = TRUE),
        MortVar_MC       = apply(mc_stats, 1, var,    na.rm = TRUE)
      ) |>
        left_join(grid_base, by = "GRID_ID") |>
        mutate(
          MortMedian_sd_MC = sqrt(MortVar_MC),
          FMortMedian      = MortMedian + MortMedian_MC,
          Residual_risk_perc = (MortMedian_MC / MortMedian) * 100
        )
    }
  ),

  tar_target(
    name = delta_risk,
    command = {
      gs <- Grid_sums
      for (zone in names(zone_ids)) {
        gs[[zone]] <- as.integer(gs$GRID_ID %in% zone_ids[[zone]])
      }
      gs[["Study_Area"]] <- 1L

      lapply(c("Study_Area", names(zone_ids)), function(zone) {
        mask  <- gs[[zone]]
        old   <- sum(gs$MortMedian   * mask, na.rm = TRUE)
        delta <- sum(gs$MortMedian_MC * mask, na.rm = TRUE)
        total <- sum(gs$FMortMedian  * mask, na.rm = TRUE)
        v_old <- sum(gs$MortVar      * mask, na.rm = TRUE)
        v_mc  <- sum(gs$MortVar_MC   * mask, na.rm = TRUE)

        data.frame(
          Area              = zone,
          target_growth     = target_growth,
          delta_risk        = delta,
          delta_risk_perc   = delta / old * 100,
          error             = 100 * sqrt(
            ((total / old^2)^2) * v_old +
              ((1     / v_old  )^2) * v_mc
          )
        )
      }) |>
        bind_rows()
    }
  ),

  tar_target(
    name = outputs,
    command = {
      outdir <- file.path(store, "final")
      dir.create(outdir, showWarnings = FALSE)

      growth_tag <- gsub("\\.", "p", as.character(target_growth))  # "0.05" -> "0p05"

      paths <- list(
        gpkg = file.path(outdir, paste0("Increased_Traffic_FResidual_Risk_RGrid_",  year_month, "_", growth_tag, ".gpkg")),
        grid = file.path(outdir, paste0("Increased_Traffic_FGrid_sums_RGrid_",      year_month, "_", growth_tag, ".csv")),
        risk = file.path(outdir, paste0("Increased_Traffic_FDelta_Risk_RGrid_",     year_month, "_", growth_tag, ".csv"))
      )

      dplyr::right_join(Grid_sums, grid, by = "GRID_ID") |>
        sf::st_as_sf() |>
        sf::st_write(paths$gpkg, layer = "polygons", delete_dsn = TRUE)

      write.csv(Grid_sums,   paths$grid, row.names = FALSE)
      write.csv(delta_risk,  paths$risk, row.names = FALSE)

      unlist(paths)
    },
    format = "file"
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

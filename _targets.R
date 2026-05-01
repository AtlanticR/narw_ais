if (!require(librarian)) {
  install.packages("librarian")
}
pkgs <- c(
  "data.table",
  "sf",
  "targets",
  "tarchetypes",
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
}

tar_config_set(store = file.path(store, "targets"))

tar_source("R/assign_trip_ids.R")

data_ais_layers <- list(
  year_month = c("2023-11", paste0("2024-", sprintf("%02d", 4:10)))
)

mapped_data_ais_layers <- tar_map(
  values = data_ais_layers,
  tar_target(
    name = data_ais,
    command = st_read(
      file.path(store, "data", "SEGS_INT_Mortality.gpkg"),
      year_month
    ) |>
      assign_trip_ids(grid, data_ports) |>
      mutate(Type = ifelse(is.na(Type), "OTHER", Type)) |>
      filter(!(GRID_ID %in% remove_grid_ids))
  ),

  tar_target(
    name = trip_nums,
    command = {
      trip_nums <- data_ais |>
        st_drop_geometry() |>
        group_by(Type) |>
        summarise(unique_trips = n_distinct(trip_id, na.rm = TRUE)) |>
        as.data.frame() |>
        left_join(conversion_rates, by = "Type") |>
        mutate(
          growth_rate = Annual_increase / 100,
          scaled_growth_rate1 = growth_rate *
            (target_growth[1] /
              sum(unique_trips / sum(unique_trips) * growth_rate)),
          scaled_growth_rate2 = growth_rate *
            (target_growth[2] /
              sum(unique_trips / sum(unique_trips) * growth_rate)),
          trip_growth_Nums1 = round(unique_trips * scaled_growth_rate1),
          trip_growth_Nums2 = round(unique_trips * scaled_growth_rate2)
        )
    }
  ),

  tar_target(
    name = mc_results,
    command = {
      n_iter = 10000 #TODO change to 10000
      Grid_IDs <- unique(data_ais$GRID_ID)

      # data.table is much faster for this type of operation than data.frame, so convert to data.table
      data_ais_dt <- as.data.table(st_drop_geometry(data_ais))

      # split the trip_ids by Type for sampling
      trip_ids_by_type <- split(data_ais_dt$trip_id, data_ais_dt$Type)

      # prepare array
      resultarray <- array(
        NA_real_,
        dim = c(length(target_growth), length(Grid_IDs), ncol = n_iter),
        dimnames = list(
          target_growth = target_growth,
          Grid_IDs = Grid_IDs,
          n_iter = seq_len(n_iter)
        )
      )

      for (p in c(0.05, 0.10)) {
        sample_sizes <- setNames(
          if (p == target_growth[1]) {
            trip_nums$trip_growth_Nums1
          } else {
            trip_nums$trip_growth_Nums2
          },
          trip_nums$Type
        )
        for (i in seq_len(n_iter)) {
          if (i %% n_iter / 10 == 0) {
            message(
              "MC iteration: ",
              i,
              " of ",
              n_iter,
              " for target growth: ",
              p
            )
          }

          sampled_ids <- unlist(
            mapply(
              sample,
              trip_ids_by_type,
              sample_sizes[names(trip_ids_by_type)],
              SIMPLIFY = FALSE
            )
          )

          temp <- data_ais_dt[
            trip_id %in% sampled_ids,
            .(MortMedian = sum(MortMedian, na.rm = TRUE)),
            by = GRID_ID
          ]

          resultarray[target_growth == p, temp$GRID_ID, i] <- temp$MortMedian
        }
      }

      resultarray
    }
  ),

  tar_target(
    name = Grid_sums,
    command = {
      grid_base <- as.data.table(st_drop_geometry(data_ais))[,
        .(
          MortMedian = sum(MortMedian, na.rm = TRUE),
          MortVar = sum(MortVar, na.rm = TRUE)
        ),
        by = GRID_ID
      ]

      Grid_sums <- data.frame(
        GRID_ID = attr(mc_results, "dimnames")$Grid_IDs,
        MortMedian_5p = apply(mc_results[1, , ], 1, median, na.rm = TRUE),
        MortMedian_10p = apply(mc_results[2, , ], 1, median, na.rm = TRUE),
        MortMedian_5p_var_MC = apply(mc_results[1, , ], 1, var, na.rm = TRUE),
        MortMedian_10p_var_MC = apply(mc_results[2, , ], 1, var, na.rm = TRUE)
      ) |>
        left_join(grid_base, by = "GRID_ID") |>
        mutate(
          MortMedian_5p_sd_MC = sqrt(MortMedian_5p_var_MC),
          MortMedian_10p_sd_MC = sqrt(MortMedian_10p_var_MC),
          FMortMedian_5p = rowSums(
            cbind(MortMedian, MortMedian_5p),
            na.rm = TRUE
          ),
          FMortMedian_10p = rowSums(
            cbind(MortMedian, MortMedian_10p),
            na.rm = TRUE
          ),
          Residual_risk_5p_perc2 = (MortMedian_5p / MortMedian) * 100, # GRID-cell percentages (new-old)/old
          Residual_risk_10p_perc2 = (MortMedian_10p / MortMedian) * 100 # GRID-cell percentages (new-old)/old
        )
    }
  ),

  tar_target(
    name = delta_risk,
    command = {
      gs <- Grid_sums
      # Tag each grid cell with zone membership
      for (zone in names(zone_ids)) {
        gs[[zone]] <- as.integer(gs$GRID_ID %in% zone_ids[[zone]])
      }
      gs[["Study_Area"]] <- 1L # all cells

      zones <- c("Study_Area", names(zone_ids))

      lapply(zones, function(zone) {
        mask <- gs[[zone]]

        old <- sum(gs$MortMedian * mask, na.rm = TRUE)
        d5 <- sum(gs$MortMedian_5p * mask, na.rm = TRUE)
        d10 <- sum(gs$MortMedian_10p * mask, na.rm = TRUE)
        tot5 <- sum(gs$FMortMedian_5p * mask, na.rm = TRUE)
        tot10 <- sum(gs$FMortMedian_10p * mask, na.rm = TRUE)

        v_old <- sum(gs$MortVar * mask, na.rm = TRUE)
        v5 <- sum(gs$MortMedian_5p_var_MC * mask, na.rm = TRUE)
        v10 <- sum(gs$MortMedian_10p_var_MC * mask, na.rm = TRUE)

        data.frame(
          Area = zone,
          delta_risk_5p = d5,
          delta_risk_5p_perc2 = d5 / old * 100,
          error_5p = 100 *
            sqrt(((tot5 / old^2)^2) * v_old + ((1 / v_old)^2) * v5),
          delta_risk_10p = d10,
          delta_risk_10p_perc2 = d10 / old * 100,
          error_10p = 100 *
            sqrt(((tot10 / old^2)^2) * v_old + ((1 / v_old)^2) * v10)
        )
      }) |>
        bind_rows()
    }
  ),

  tar_target(
    name = outputs,
    command = {
      ym <- gsub("-", "_", year_month)
      outdir <- file.path(store, "final")
      dir.create(outdir, showWarnings = FALSE)

      paths <- list(
        gpkg = file.path(
          outdir,
          paste0("Increased_Traffic_FResidual_Risk_RGrid_", ym, ".gpkg")
        ),
        grid = file.path(
          outdir,
          paste0("Increased_Traffic_FGrid_sums_RGrid_", ym, ".csv")
        ),
        risk = file.path(
          outdir,
          paste0("Increased_Traffic_FDelta_Risk_RGrid_", ym, ".csv")
        )
      )

      geo <- dplyr::right_join(
        Grid_sums,
        grid,
        by = "GRID_ID"
      ) |>
        sf::st_as_sf()

      sf::st_write(geo, paths$gpkg, layer = "polygons", delete_dsn = TRUE)
      write.csv(Grid_sums, paths$grid, row.names = FALSE)
      write.csv(delta_risk, paths$risk, row.names = FALSE)

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
    name = conversion_rates,
    command = data.frame(
      Type = c(
        "CARGO",
        "FERRY",
        "FISHING",
        "GOV_RES",
        "OTHER",
        "PASSENGER",
        "PLEASURE",
        "TANKER",
        "TUG"
      ),
      Annual_increase = c(
        2.094334579,
        2.978065767,
        2.978065767,
        2.978065767,
        2.978065767,
        2.978065767,
        2.978065767,
        1.409169349,
        2.978065767
      )
    )
  ),

  tar_target(
    name = target_growth,
    command = c(0.05, 0.10) # this should only be 2 values or will break downstream code
  ),

  mapped_data_ais_layers
)

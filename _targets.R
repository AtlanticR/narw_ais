if (!require(librarian)) {
  install.packages("librarian")
}
pkgs <- c(
  "sf",
  "targets",
  "qs2",
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

list(
  tar_target(
    name = grid,
    command = st_read(file.path(store, "data", "Grid.gpkg"), "GridGulf")
  ),

  tar_target(
    name = data_ports,
    command = st_read(file.path(store, "data", "Grid.gpkg"), "Ports_100mBuffer")
  ),

  tar_target(
    name = data_ais,
    command = st_read(
      file.path(store, "data", "SEGS_INT_Mortality.gpkg"),
      "2023-11"
    )
  ),

  tar_target(
    name = processed_ais,
    command = assign_trip_ids(data_ais, grid, data_ports)
  )
)

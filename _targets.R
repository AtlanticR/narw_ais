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

# function to assign trip ids
assign_trip_ids <- function(ais_data, grid, ports, min_gap_records = 3) {
  # common CRS
  if (st_crs(ais_data) != st_crs(ports)) {
    ports <- st_transform(ports, st_crs(ais_data))
  }
  if (st_crs(grid) != st_crs(ports)) {
    grid <- st_transform(grid, st_crs(ports))
  }

  # Build port grid cells, tagged with which port they belong to
  # Each grid cell gets the id of the port it intersects (or is nearest to)
  # so that adjacent ports map to distinct grid cells

  ports_intersecting_idx <- lengths(st_intersects(ports, grid)) > 0

  # grid cells that intersect a port — join port id
  grid_intersects <- grid[lengths(st_intersects(grid, ports)) > 0, ]
  grid_intersects$port_id <- ports$id[
    st_nearest_feature(grid_intersects, ports) # for cells touching multiple ports, assign nearest
  ]

  # ports with no intersecting grid cell — find their nearest grid cell
  ports_no_grid <- ports[!ports_intersecting_idx, ]
  grid_nearest <- grid[st_nearest_feature(ports_no_grid, grid), ] |>
    mutate(port_id = ports_no_grid$id)
  portgrid <- rbind(grid_intersects, grid_nearest) |>
    group_by(GRID_ID) |>
    slice(1) |> # a grid cell does not need 2 IDs
    ungroup()

  port_grid_ids <- unique(portgrid$GRID_ID)

  # lookup: GRID_ID -> port_id
  grid_port_lookup <- portgrid |>
    st_drop_geometry() |>
    select(GRID_ID, port_id)

  # ── 2. Sort + flag at-port records, tagged with which port ──────────────────
  result <- ais_data |>
    arrange(MMSI, UNIX_start) |>
    mutate(in_port_grid = GRID_ID %in% port_grid_ids) |>
    left_join(grid_port_lookup, by = "GRID_ID") |>
    # Smooth yo-yo + label discrete port visits
    # A new visit only starts when BOTH: entering port zone AND the port_id changes
    # (or re-enters after a genuine departure) — prevents neighbouring ports
    # from merging into one visit
    group_by(MMSI) |>
    mutate(
      .consec_nonport = sequence(rle(in_port_grid)$lengths) * !in_port_grid,
      .departed = .consec_nonport >= min_gap_records,
      .in_port_smooth = in_port_grid | (!in_port_grid & !.departed),
      # new visit when: entering port smooth zone AND (different port_id OR re-entry)
      .port_block_change = .in_port_smooth &
        (!lag(.in_port_smooth, default = FALSE) |
          (port_id != lag(port_id, default = NA_real_) & !is.na(port_id))),
      ,
      .visit_id = {
        r <- rle(.in_port_smooth)
        rep(cumsum(r$values) * r$values, r$lengths)
      }
    ) |>
    ungroup() |>
    group_by(.in_port_smooth) |>
    mutate(
      .dist_to_port = ifelse(
        .in_port_smooth,
        apply(st_distance(geom, ports), 1, min),
        Inf
      )
    ) |>
    ungroup() |>

    # Split point = closest record to its specific port per visit
    group_by(MMSI, .visit_id) |>
    mutate(
      .is_split = .visit_id > 0 & row_number() == which.min(.dist_to_port)
    ) |>
    ungroup() |>
    group_by(MMSI) |>
    mutate(
      .has_departure_after = rev(cumsum(rev(
        !.in_port_smooth & !in_port_grid
      ))) >
        0,
      .is_split = .is_split & .has_departure_after
    ) |>
    ungroup() |>

    # Assign trip numbers
    group_by(MMSI) |>
    mutate(
      .trip_num = as.integer(cumsum(.is_split))
    ) |>
    ungroup() |>
    mutate(trip_id = paste0(MMSI, "_T", sprintf("%05d", .trip_num))) |>
    select(-starts_with("."), -port_id, -in_port_grid)

  return(result)
}

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

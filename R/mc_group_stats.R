# Replaces the old whale-iteration loop in mc_results. D is the zone's
# whale_iter x n_events matrix; S is the event x vessel_iter selection mask
# (same S regardless of which group -- zone or grid -- is being summarized).
#
# For a single group g, R_g <- D[, cols in g] %*% S[cols in g, ] gives a
# whale_iter x vessel_iter matrix in one BLAS call, equivalent to what the
# old loop built one whale-iteration-slice at a time via NA-masking +
# rowsum(). Row-wise matrixStats across the vessel_iter axis then reproduce
# the same per-whale-iteration statistics the old loop computed.

#' @param D whale_iter x n_events dense matrix (data_whaleresults_matrix)
#' @param S n_events x n_vessel_iter sparse 0/1 Matrix (mc_sampled_mask)
#' @param group_vec length n_events vector assigning each column of D / row
#'   of S to a group (zone1 or GRID_ID) -- must stay full-length and aligned
#'   to D's columns; do not pre-filter it, filter via `groups` instead
#' @param groups which group ids to compute stats for (subset of unique(group_vec))
#' @param group_col_name name of the output column holding the group id
compute_group_stats <- function(D, S, group_vec, groups, group_col_name = "group") {
  purrr::map_dfr(groups, function(g) {
    cols <- which(group_vec == g)
    if (length(cols) == 0) return(NULL)

    Dg <- D[, cols, drop = FALSE]
    Sg <- S[cols, , drop = FALSE]

    R <- as.matrix(Dg %*% Sg)  # whale_iter x vessel_iter

    tibble::tibble(
      !!group_col_name := g,
      whale_iter = seq_len(nrow(R)),
      vessel_iter_delta_mean        = matrixStats::rowMeans2(R, na.rm = TRUE),
      vessel_iter_delta_median      = matrixStats::rowMedians(R, na.rm = TRUE),
      vessel_iter_delta_variance    = matrixStats::rowVars(R, na.rm = TRUE),
      vessel_iter_delta_q025        = matrixStats::rowQuantiles(R, probs = 0.025, na.rm = TRUE, drop = TRUE),
      vessel_iter_delta_q975        = matrixStats::rowQuantiles(R, probs = 0.975, na.rm = TRUE, drop = TRUE),
      vessel_iter_delta_n_delta_gt0 = matrixStats::rowCounts(R > 0, na.rm = TRUE)
    )
  })
}

#' Raw (unreduced) whale_iter x vessel_iter values for one group, melted to
#' long format. Used only for shared grid cells, where stats aren't valid
#' until this zone's contribution is summed with every other zone's
#' contribution to the same grid cell (see combine_and_stat_shared_grids_duckdb
#' in _targets.R).
#'
#' @param group_id single group id (one shared GRID_ID) to extract
compute_group_raw_long <- function(D, S, group_vec, group_id) {
  cols <- which(group_vec == group_id)
  if (length(cols) == 0) return(NULL)

  Dg <- D[, cols, drop = FALSE]
  Sg <- S[cols, , drop = FALSE]
  R <- as.matrix(Dg %*% Sg)  # whale_iter x vessel_iter

  dt <- data.table::as.data.table(R)
  data.table::setnames(dt, as.character(seq_len(ncol(R))))
  dt[, whale_iter := seq_len(.N)]

  data.table::melt(
    dt,
    id.vars = "whale_iter",
    variable.name = "vessel_iter",
    value.name = "value"
  )[, `:=`(
    vessel_iter = as.integer(as.character(vessel_iter)),
    GRID_ID = group_id
  )]
}

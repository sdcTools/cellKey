#' Utility measures for perturbed counts
#'
#' This function computes utility/information loss measures
#' based on two numeric vectors (original and perturbed)
#'
#' @param orig a numeric vector holding original values
#' @param pert a numeric vector holding perturbed values
#' @param exclude_zeros a scalar logical value; if `TRUE` (the default), all only cells
#' with counts `> 0` are used when computing distances `d1`, `d2` and `d3`. If this
#' argument is `FALSE`, the complete vector is used.
#' @return a `list` containing the following elements:
#' - `overview`: a `data.table` with the following three columns:
#'    * `noise`: amount of noise computed as `orig` - `pert`
#'    * `cnt`: number of cells perturbed with the value given in column `noise`
#'    * `pct`: percentage of cells perturbed with the value given in column `noise`
#' - `measures`: a `data.table` containing measures of the distribution
#' of three different distances between original and perturbed values
#' of the unweighted counts. Column `what` specifies the computed measure.
#' The three distances considered are:
#'    * `d1`: absolute distance between original and masked values
#'    * `d2`: relative absolute distance between original and masked values
#'    * `d3`: absolute distance between square-roots of original and perturbed
#' values
#'
#' - `cumdistr_d1`, `cumdistr_d2` and `cumdistr_d3`: for each distance `d1`, `d2`
#' and `d3`, a `data.table` with the following three columns:
#'    * `cat`: a specific value (for `d1`) or interval (for distances `d2` and `d3`)
#'    * `cnt`: number of records smaller or equal the value in column `cat` for the
#'    given distance
#'    * `pct`: proportion of records smaller or equal the value
#'    in column `cat` for the selected distance
#' - `false_zero`: number of cells that were perturbed to zero
#' - `false_nonzero`: number of cells that were initially zero but
#' have been perturbed to a number different from zero
#' - `exclude_zeros`: were empty cells exluded from computation or not
#' @export
#' @md
#' @examples
#' orig <- c(1:10, 0, 0)
#' pert <- orig; pert[c(1, 5, 7)] <- c(0, 6, 9)
#'
#' # ignore empty cells when computing measures `d1`, `d2`, `d3`
#' ck_cnt_measures(orig = orig, pert = pert, exclude_zeros = TRUE)
#'
#' # use all cells
#' ck_cnt_measures(orig = orig, pert = pert, exclude_zeros = FALSE)
#'
#' # for an application on a perturbed object, see ?cellkey_pkg
ck_cnt_measures <- function(orig, pert, exclude_zeros = TRUE) {
  pct <- cnt <- NULL

  # for measures, see
  # https://ec.europa.eu/eurostat/cros/system/files/methods_for_protecting_census_data.pdf

  if (!is.numeric(orig) | !is.numeric(pert)) {
    stop("Arguments `orig` and `pert` must be numeric vectors.", call. = FALSE)
  }
  if (length(orig) != length(pert)) {
    stop("Number of elements in `orig` and `pert` differs.", call. = FALSE)
  }
  if (!rlang::is_scalar_logical(exclude_zeros)) {
    stop("Argument `exclude_zeros` needs to be a scalar logical value.", call. = FALSE)
  }

  dt_overview <- as.data.table(as.data.frame.table(table(orig - pert)))
  dt_overview$pct <- dt_overview$Freq / length(orig)
  setnames(dt_overview, c("noise", "cnt", "pct"))

  quantvals <- c(0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, Inf)

  # distances per cell
  # absolute distances

  if (exclude_zeros) {
    ii <- which(pert != 0)
    if (length(ii) == 0) {
      dist_d1 <- dist_d2 <- dist_d3 <- rep(NA, length(orig))
    }
    ovec <- orig[ii]
    pvec <- pert[ii]
  } else {
    ovec <- orig
    pvec <- pert
  }

  dist_d1 <- abs(pvec - ovec)

  # relative absolute distance
  dist_d2 <- (abs(pvec - ovec)) / ovec
  dist_d2[is.nan(dist_d2)] <- 0

  # absolute distances of square roots
  dist_d3 <- abs(sqrt(pvec) - sqrt(ovec))

  ## cumdistr
  vv <- get_distr_vals(dist_d1)
  dt_measures <- data.table(what = names(vv), d1 = vv)
  dt_measures$d2 <- get_distr_vals(dist_d2)
  dt_measures$d3 <- get_distr_vals(dist_d3)

  # absolute distances
  cumsum_d1 <- cumsum(table(dist_d1))
  cumdistr_d1 <- data.table(cat = names(cumsum_d1), cnt = cumsum_d1)
  cumdistr_d1[, pct := cnt / length(dist_d1)]

  # relative absolute distances
  cumsum_d2 <- cumsum(table(cut(dist_d2, quantvals, include.lowest = TRUE)))
  cumdistr_d2 <- data.table(cat = names(cumsum_d2), cnt = cumsum_d2)
  cumdistr_d2[, pct := cnt / length(dist_d2)]

  # absolute distance between square roots
  cumsum_d3 <- cumsum(table(cut(dist_d3, quantvals, include.lowest = TRUE)))
  cumdistr_d3 <- data.table(cat = names(cumsum_d3), cnt = cumsum_d3)
  cumdistr_d3[, pct := cnt / length(dist_d3)]

  # false zero cells: cells that were perturbed to zero
  false_zero <- sum(pert == 0 & orig != 0)
  # false positive: zeros that were perturbed to a value !=0
  false_nonzero <- sum(orig == 0 & pert != 0)

  return(
    list(
      overview = dt_overview[],
      measures = dt_measures[],
      cumdistr_d1 = cumdistr_d1,
      cumdistr_d2 = cumdistr_d2,
      cumdistr_d3 = cumdistr_d3,
      false_zero = false_zero,
      false_nonzero = false_nonzero,
      exclude_zeros = exclude_zeros
    )
  )
}

#' ck_cnt_measures
#'
#' compute utility/information loss measures for count variables
#' of a perturbed table.
#' @param x input object of class \code{\linkS4class{pert_table}}
#' @param vname a specific tabulated count variable
#' @return a \code{list} containing the following elements:
#' \itemize{
#' \item \strong{measures}: a \code{data.table} containing measures
#' of the distribution of three different distances between original
#' and perturbed values of the unweighted counts. The three distances are:
#' \itemize{
#' \item absolute distance (column \code{vals_abs}), \code{D1}
#' \item relative distance (column \code{vals_rel}), \code{D2}
#' \item absolute distance between square-roots of original and perturbed
#' values (column \code{vals_r}), \code{D3}
#' }
#' \item \strong{cumDistrA}: a \code{data.table} with 3 columns
#' showing the cummulative sum of absolute distances
#' \itemize{
#' \item \code{kat}: possible values for distance \code{D1}
#' \item \code{val_abs}: number of records smaller or equal the
#' value in column \code{kat}
#' \item \code{prop_abs}: proportion or records smaller or equal
#' the value in column \code{kat}
#' }
#' \item \strong{cumDistrB}: a \code{data.table} with 5 columns
#' which show the cummulative distributions for distances
#' \code{D2} (columns \code{val_rel} and \code{prop_rel})
#' and \code{D3} (columns \code{val_r} and \code{prop_r}) for a set
#' of intervals shown in column \code{kat}.
#' \itemize{
#' \item \code{kat}: a specific interval
#' \item \code{val_rel}: number of records smaller or equal the
#' value in column \code{kat} for distance \code{D2}
#' \item \code{prop_rel} proportion of records smaller or equal the
#' value in column \code{kat} for distance \code{D2}
#' \item \code{val_r}: number of records smaller or equal the value
#' in column \code{kat} for distance \code{D3}
#' \item \code{prop_r}: proportion of records smaller or equal the
#' value in column \code{kat} for distance \code{D3}
#' }
#' \item \strong{false_zero}: number of cells that were perturbed to zero
#' \item \strong{false_positives}: number of cells that were initially
#' zero but have been perturbed to a number different from zero
#' }
#' @export
#' @examples
#' ## see example in perturbTable
ck_cnt_measures <- function(x, vname="Total") {
  stopifnot(isS4(x), "pert_table" %in% class(x))
  stopifnot(is_scalar_character(vname))
  stopifnot(vname %in% slot(x, "countVars"))


  # unweighted
  tab <- slot(x, "tab")

  # unweighted
  v_orig <- paste0("UWC_", vname)
  v_pert <- paste0("pUWC_", vname)
  return(ck_cnt_measures_basic(tab[[v_orig]], tab[[v_pert]]))
}

#' ck_cnt_measures_basic
#'
#' compute utility/information loss measures two numeric vectors
#' (original and perturbed)
#' @param orig a numeric vector holding original values
#' @param pert a numeric vector holding perturbed values
#' @return a \code{list} containing the following elements:
#' \itemize{
#' \item \strong{measures}: a \code{data.table} containing measures
#' of the distribution of three different distances between
#' original and perturbed values of the unweighted counts.
#' The three distances are:
#' \itemize{
#' \item absolute distance (column \code{vals_abs}), \code{D1}
#' \item relative distance (column \code{vals_rel}), \code{D2}
#' \item absolute distance between square-roots of original and perturbed
#' values (column \code{vals_r}), \code{D3}
#' }
#' \item \strong{cumDistrA}: a \code{data.table} with 3 columns
#' showing the cummulative sum of absolute distances
#' \itemize{
#' \item \code{kat}: possible values for distance \code{D1}
#' \item \code{val_abs}: number of records smaller or equal the
#' value in column \code{kat}
#' \item \code{prop_abs}: proportion or records smaller or equal
#' the value in column \code{kat}
#' }
#' \item \strong{cumDistrB}: a \code{data.table} with 5 columns
#' which show the cummulative distributions for distances
#' \code{D2} (columns \code{val_rel} and \code{prop_rel}) and
#' \code{D3} (columns \code{val_r} and \code{prop_r}) for a set
#' of intervals shown in column \code{kat}.
#' \itemize{
#' \item \code{kat}: a specific interval
#' \item \code{val_rel}: number of records smaller or equal
#' the value in column \code{kat} for distance \code{D2}
#' \item \code{prop_rel} proportion of records smaller or
#' equal the value in column \code{kat} for distance \code{D2}
#' \item \code{val_r}: number of records smaller or equal
#' the value in column \code{kat} for distance \code{D3}
#' \item \code{prop_r}: proportion of records smaller or equal the
#' value in column \code{kat} for distance \code{D3}
#' }
#' \item \strong{false_zero}: number of cells that were perturbed to zero
#' \item \strong{false_positives}: number of cells that were initially zero but
#' have been perturbed to a number different from zero
#' }
#' @export
#' @examples
#' orig <- 1:10
#' pert <- orig; pert[c(1,5,7)] <- c(0,6,9)
#' ck_cnt_measures_basic(orig=orig, pert=pert)
#'
ck_cnt_measures_basic <- function(orig, pert) {
  stopifnot(is.numeric(orig), is.numeric(pert), length(orig) == length(pert))

  val_abs <- prop_abs <- val_rel <- vals_rel <- prop_rel <- NULL
  val_r <- vals_r <- prop_r <- NULL
  quantvals <- c(
    seq(from = 0, to = 0.1, by = 0.01),
    seq(from = 0.15, to = 1, by = 0.05),
    Inf
  )

  dabs <- abs(pert - orig)
  drel <- (abs(pert - orig)) / orig
  drel[is.nan(drel)] <- 0

  dr <- abs(sqrt(pert) - sqrt(orig))

  # absolute
  cs_abs <- cumsum(table(dabs))
  dt_a <- data.table(kat = names(cs_abs), val_abs = cs_abs)
  dt_a[, prop_abs := val_abs / length(dabs)]

  # relative and squared distances
  drelf <- cut(drel, quantvals, include.lowest = TRUE)
  drf <- cut(dr, quantvals, include.lowest = TRUE)

  cs_rel <- cumsum(table(drelf))
  cs_r <- cumsum(table(drf))

  dt_b <- data.table(kat = names(cs_rel), val_rel = cs_rel)
  dt_b[, prop_rel := val_rel / length(drel)]

  dt_b[, val_r := cs_r]
  dt_b[, prop_r := val_rel / length(dr)]

  ## cumdistr
  vv <- get_distr_vals(dabs)
  dt2 <- data.table(what = names(vv), vals_abs = vv)
  dt2[, vals_rel := get_distr_vals(drel)]
  dt2[, vals_r := get_distr_vals(dr)]

  # false zero cells: cells that were perturbed to zero
  false_zero <- sum(pert == 0 & orig != 0)
  # false positive: zeros that were perturbed to a value !=0
  false_positives <- sum(orig == 0 & pert != 0)

  return(
    list(
      measures = dt2,
      cumdistrA = dt_a,
      cumdistrB = dt_b,
      false_zero = false_zero,
      false_positives = false_positives
    )
  )
}

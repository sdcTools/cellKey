#' Generate random record keys
#'
#' This function allows to create random record keys from a uniform distribution. If no seed is 
#' specified, a seed value is computed from the input data set to allow for reproducability depending
#' on the input data set.
#'
#' @param dat microdata used to generated hash for random seed
#' @param nr_digits maximum number of digits in the record keys.
#' @param seed if not `NULL`, a number specifying the initial seed value 
#' for the random number generator. If `NULL`, a seed is computed from `dat` itself.
#' @return a numeric vector with `nrow(dat)` record keys
#' @export
#' @md
#' @examples
#' dat <- ck_create_testdata()
#' dat$rkeys <- ck_generate_rkeys(dat = ck_create_testdata(), nr_digits = 6)
ck_generate_rkeys <- function(dat, nr_digits=10, seed=NULL) {
  if (!is.null(seed)) {
    stopifnot(is_scalar_integerish(seed))
  } else {
    seed <- ck_create_seed_from_hash(dat)
  }

  stopifnot(is_scalar_integerish(nr_digits))
  stopifnot(nr_digits >= 5 & nr_digits <= 20)
  set.seed(seed)
  round(runif(n = nrow(dat), min = 0, max = 1), digits = nr_digits)
}

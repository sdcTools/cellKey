#' ck_generate_rkeys
#'
#' creates random record keys by sampling between 1:max_val
#'
#' @param dat microdata used to generated hash for random seed
#' @param max_val maximum number for record-key
#' @return a numeric vector with \code{N} record keys
#' @export
#'
#' @examples
#' dat <- ck_create_testdata()
#' dat$rkeys <- ck_generate_rkeys(dat=dat, max_val=(2^6)-1)
ck_generate_rkeys <- function(dat, max_val=(2^6)-1) {
  N <- nrow(dat)
  seed <- ck_create_seed_from_hash(dat)
  set.seed(seed)
  rkeys <- sample(1:max_val, size=N, replace=TRUE)
  rkeys
}

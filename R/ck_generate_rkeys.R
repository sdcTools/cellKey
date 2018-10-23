#' ck_generate_rkeys
#'
#' creates random record keys by sampling between 1:max_val
#'
#' @param dat microdata used to generated hash for random seed
#' @param max_val maximum number for record-key
#' @param max_digits maximum number of digits in the record keys.
#' This parameter is only used if \code{type} equals \code{abs}.
#' This parameter is only used if \code{type} equals \code{destatis}.
#' @param type (character) choice how to compute record keys. Valid choices
#' are \code{abs} and \code{destatis}.
#' @param verbose (logical) if \code{TRUE}, print additional information
#' @return a numeric vector with \code{N} record keys given a seed depending on \code{dat}
#' @export
#' @examples
#' dat <- ck_create_testdata()
#' dat$rkeys <- ck_generate_rkeys(dat=dat, max_val=2*nrow(dat), type="abs")
#' ## rkeys for destatis method
#' # ck_generate_rkeys(dat=dat, max_digits=6, type="destatis")
ck_generate_rkeys <- function(dat, max_val=nrow(dat), max_digits=10, type="abs", verbose=TRUE) {
  gen_key_abs <- function(N, max_val, seed) {
    set.seed(seed)
    rkeys <- sample(1:max_val, size=N, replace=TRUE)
    rkeys
  }
  gen_key_destatis <- function(N, max_digits, seed) {
    set.seed(seed)
    rkeys <- round(runif(N, min=0, max=1), digits=max_digits)
    rkeys
  }

  stopifnot(is_scalar_logical(verbose))
  stopifnot(is_scalar_character(type))
  stopifnot(type %in% c("abs","destatis"))

  seed <- ck_create_seed_from_hash(dat)
  if (type=="abs") {
    stopifnot(is_scalar_integerish(max_val))
    stopifnot(max_val >= nrow(dat))
    if (verbose) {
      mf <- match.call(expand.dots=FALSE)
      if (!is.null(mf$max_digits)) {
        message("Note: Argument 'max_digits' is ignored!")
      }
    }
    rkeys <- gen_key_abs(N=nrow(dat), max_val=max_val, seed=seed)
  }
  if (type=="destatis") {
    stopifnot(is_scalar_integerish(max_digits))
    stopifnot(max_digits>=5 & max_digits<=20)
    if (verbose) {
      mf <- match.call(expand.dots=FALSE)
      if (!is.null(mf$max_val)) {
        message("Note: Argument 'max_val' is ignored!")
      }
    }
    rkeys <- gen_key_destatis(N=nrow(dat), max_digits=max_digits, seed=seed)
  }
  return(rkeys)
}
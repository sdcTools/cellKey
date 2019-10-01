#' @import data.table
#' @importFrom rlang is_double is_integerish is_na is_scalar_integerish is_scalar_integer
#' @importFrom rlang is_scalar_character is_bare_integerish is_scalar_atomic
#' @importFrom rlang is_scalar_logical is_character is_scalar_double
#' @importFrom methods new validObject is slot slot<-
#' @importFrom stats addmargins formula weights xtabs runif median na.omit quantile
#' @importFrom utils data RShowDoc
#' @importFrom digest sha1
#' @importFrom sdcTable makeProblem sdcProb2df
#' @importFrom ptable pt_create_pParams pt_create_pTable
#' @import sdcHierarchies
#'

Sys.setenv(".ck_nr_digits" = 8)
.ck_digits <- function(digits = 8) {
  d <- as.numeric(Sys.getenv(".ck_nr_digits"))
  if (is.na(d)) {
    return(digits)
  }
  return(d)
}

# number of cores for parallel processing
# number of cores to be used when doing parallel computing
.ck_cores <- function() {
  if (Sys.info()["sysname"] == "Windows") {
    return(1)
  }
  available_cores <- parallel::detectCores() - 2
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    # see https://stackoverflow.com/a/50571533
    num_workers <- 1L
  } else {
    num_workers <- max(1, available_cores)
  }
  max(1, min(num_workers, available_cores))
}
NULL

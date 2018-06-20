.onLoad <- function(libname, pkgname) {
  data("ptable_destatis", "ptable_destatis", package=pkgname, envir=parent.env(environment()))
}

#' @import data.table
#' @importFrom rlang is_double is_integerish is_na is_scalar_integerish is_scalar_integer
#' @importFrom rlang is_scalar_character is_bare_integerish is_scalar_atomic
#' @importFrom rlang is_scalar_logical is_character
#' @importFrom primes is_prime
#' @importFrom methods new validObject is slot slot<-
#' @importFrom stats addmargins formula weights xtabs runif
#' @importFrom utils data
#' @importFrom digest sha1
#' @importFrom sdcTable makeProblem sdcProb2df
#' @importFrom data.tree Node Prune Traverse FindNode
NULL

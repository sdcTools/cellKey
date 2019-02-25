#' An S4 class to represent perturbation parameters for continous tables
#' @slot big_n (integer) large prime number used to derive
#' cell keys from record keys
#' @slot small_n (integer) parameter for parameter `small_n` from as
#' introduced in the ABS-specification of the algorithm.
#' @slot ptab (data.table) perturbation table with 256 rows
#' @slot ptab_size (integer) number of columns of \code{pTable}
#' @slot mtab numeric vector specifying parameter mTable
#' for continous perturbation
#' @slot small_c (integer) specifying parameter smallC for
#' continous perturbation
#' the row in the lookup-table that was used to get the
#' perturbation value
#' @slot stab numeric vector specifying parameter sTable
#' for continous perturbation
#' @slot top_k (integer) specifiying the number of units in
#' each cell whose values
#' will be perturbed differently
#' @slot type (character) specifying the type of pTable
#' (either "destatis", "abs" or "free")
#' @name pert_params-class
#' @rdname pert_params-class
#' @export
setClass("pert_params",
representation = list(
  big_n = "integer",
  small_n = "integer",
  ptab = "data.table",
  ptab_size = "integer",
  stab = "data.table",
  mtab = "numeric",
  small_c = "integer",
  top_k = "integer",
  type = "character"
),
prototype = list(
  big_n = integer(),
  small_n = integer(),
  ptab = data.table(),
  ptab_size = integer(),
  mtab = numeric(),
  small_c = integer(),
  stab = data.table(),
  top_k = integer(),
  type = character()
),
validity = function(object) {
  if (length(slot(object, "mtab")) > 0) {
    stopifnot(all(slot(object, "mtab") > 0))
    stopifnot(is_scalar_integer(slot(object, "top_k")))
  }
  if (nrow(slot(object, "stab")) > 0) {
    stopifnot(is_scalar_integer(slot(object, "small_c")))
  }
  stopifnot(is_scalar_character(slot(object, "type")))
  if (!slot(object, "type") %in% c("destatis", "abs", "free")) {
    stop("`type` must be either `destatis`, `abs` or `free`", call. = FALSE)
  }
  stopifnot(is_scalar_integer(slot(object, "big_n")))
  stopifnot(is_scalar_integer(slot(object, "small_n")))

  if (slot(object, "type") == "abs") {
    if (!is_prime(slot(object, "big_n"))) {
      stop("`big_n` must be a prime number!", call. = FALSE)
    }
  } else if (slot(object, "type") == "abs") {
    if (slot(object, "big_n") != 1) {
      stop("`big_n` must equal 1 for method `destatis`.", call. = FALSE)
    }
  }
  return(TRUE)
})
NULL

#' An S4-class to represent input data for applying the cell-key method
#'
#' @slot microdat data.table containing microdata
#' @slot rkeys (integer) vector specifying record keys
#' @slot pert_params pert_params information about perturbation parameters
#' @name pert_inputdat-class
#' @rdname pert_inputdat-class
#' @export
setClass("pert_inputdat",
representation = list(
  microdat = "data.table",
  rkeys = "numeric",
  pert_params = "pert_params"
),
prototype = list(
  microdat = data.table(),
  rkeys = numeric(),
  pert_params = NULL
),
validity = function(object) {
  if (!is.null(slot(object, "rkeys"))) {
    stopifnot(length(slot(object, "rkeys")) == nrow(slot(object, "microdat")))
    stopifnot(all(slot(object, "rkeys") >= 0))
  }
  return(TRUE)
})
NULL

#' An S4 class to represent a perturbed table
#'
#' @slot tab a \code{data.table} containing original and
#' pertubed values. The following variables are always present:
#' \itemize{
#' \item \strong{UWC: } unweighted counts
#' \item \strong{pUWC: } perturbed unweighted counts
#' \item \strong{WC: } weighted counts
#' \item \strong{pWC: } perturbed weighted counts
#' \item \strong{WCavg: } average weight for each cell
#' }
#' Additionally, for each numerical variable the perturbed variable
#' named {vName.pert} is included in this table.
#' @slot count_modifications a \code{data.table} with 3 columns
#' (\code{row_indices}, \code{col_indices}
#' and \code{applied_perturbation}) that contains information
#' for each cell, where the applied perturbation
#' was extracted from the perturbation table (useful for debugging)
#' @slot numvars_modifications (list) containing for each
#' numerical variable that was tabulated a list element with
#' a \code{data.table} showing which values have been modified
#' prior to tabulation.
#' @slot dimvars (character) vector of variable names used to define
#' the table structure
#' @slot numvars (character) variable names of numeric variables
#' that have been tabulated
#' @slot countvars (character) vector containing names of variables
#' that have been tabulated using the frequency approach
#' @slot by (character) variable name (that must also be listed
#' in \code{countvars} and that is used to create perturbed tables
#' for \code{numvars} given the groups defined in \code{by}
#' @slot is_weighted (logical) TRUE if sampling weights have been used
#' @slot type (character) either \code{"abs"} or \code{"destatis"}
#' depending on the format
#' of the perturbation table that was used.
#' @name pert_table-class
#' @rdname pert_table-class
#' @export
setClass("pert_table",
representation = list(
  tab = "data.table",
  count_modifications = "data.table",
  numvars_modifications = "data.table",
  dimvars = "character",
  countvars = "character",
  numvars = "character",
  by = "character",
  is_weighted = "logical",
  type = "character"
),
prototype = list(
  tab = data.table(),
  count_modifications = data.table(),
  numvars_modifications = data.table(),
  dimvars = character(),
  countvars = character(),
  numvars = character(),
  by = character(),
  is_weighted = c(),
  type = character()
),
validity = function(object) {
  stopifnot(object@type %in% c("abs", "destatis", "free"))
  by <- slot(object, "by")
  if (by != "") {
    stopifnot(length(by) == 1)
    stopifnot(by %in% slot(object, "countvars"))
  }
})
NULL

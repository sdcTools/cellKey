#' An S4 class to represent perturbation parameters for continous tables
#' @slot bigN (integer) large prime number used to derive cell keys from record keys
#' @slot smallN (integer) parameter for smallN
#' @slot pTable (data.table) perturbation table with 256 rows
#' @slot pTableSize (integer) number of columns of \code{pTable}
#' @slot mTable numeric vector specifying parameter mTable for continous perturbation
#' @slot smallC (integer) specifying parameter smallC for continous perturbation
#' the row in the lookup-table that was used to get the perturbation value
#' @slot sTable numeric vector specifying parameter sTable for continous perturbation
#' @slot topK (integer) specifiying the number of units in each cell whose values
#' will be perturbed differently
#' @slot type (character) specifying the type of pTable (either 'abs' or 'destatis')
#' @name pert_params-class
#' @rdname pert_params-class
#' @export
setClass("pert_params",
representation = list(
  bigN = "integer",
  smallN = "integer",
  pTable = "data.table",
  pTableSize = "integer",
  sTable = "data.table",
  mTable = "numeric",
  smallC = "integer",
  topK = "integer",
  type = "character"
),
prototype = list(
  bigN = integer(),
  smallN = integer(),
  pTable = data.table(),
  pTableSize = integer(),
  mTable = numeric(),
  smallC = integer(),
  sTable = data.table(),
  topK = integer(),
  type = character()
),
validity = function(object) {
  if (length(slot(object, "mTable")) > 0) {
    stopifnot(all(slot(object, "mTable") > 0))
    stopifnot(is_scalar_integer(slot(object, "topK")))
  }
  if (nrow(slot(object, "sTable")) > 0) {
    stopifnot(is_scalar_integer(slot(object, "smallC")))
  }
  stopifnot(is_scalar_character(slot(object, "type")))
  if (!slot(object, "type") %in% c("abs", "destatis", "custom")) {
    stop("type must be either 'abs', 'destatis' or 'custom'\n")
  }
  stopifnot(is_scalar_integer(slot(object, "bigN")))
  stopifnot(is_scalar_integer(slot(object, "smallN")))

  if (slot(object, "type") == "abs") {
    if (!is_prime(slot(object, "bigN"))) {
      stop("bigN must be a prime number!\n")
    }
  } else if (slot(object, "type") == "abs") {
    if (slot(object, "bigN") != 1) {
      stop("bigN must equal 1 for the destatis-method!\n")
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
#' @slot tab a \code{data.table} containing original and pertubed values. The following variables are always present:
#' \itemize{
#' \item \strong{UWC: } unweighted counts
#' \item \strong{pUWC: } perturbed unweighted counts
#' \item \strong{WC: } weighted counts
#' \item \strong{pWC: } perturbed weighted counts
#' \item \strong{WCavg: } average weight for each cell
#' }
#' Additionally, for each numerical variable the perturbed variable named {vName.pert} is included in this table.
#' @slot count_modifications a \code{data.table} with 3 columns (\code{row_indices}, \code{col_indices}
#' and \code{applied_perturbation}) that contains information for each cell, where the applied perturbation
#' was extracted from the perturbation table (useful for debugging)
#' @slot numvars_modifications (list) containing for each numerical variable
#' that was tabulated a list element with a \code{data.table} showing which
#' values have been modified prior to tabulation.
#' @slot dimVars (character) vector of variable names used to define the table structure
#' @slot numVars (character) variable names of numeric variables that have been tabulated
#' @slot countVars (character) vector containing names of variables that have been tabulated
#' using the frequency approach
#' @slot by (character) variable name (that must also be listed in \code{countVars}
#' and that is used to create perturbed tables for \code{numVars} given the groups defined
#' in \code{by}
#' @slot is_weighted (logical) TRUE if sampling weights have been used
#' @slot type (character) either \code{"abs"} or \code{"destatis"} depending on the format
#' of the perturbation table that was used.
#' @name pert_table-class
#' @rdname pert_table-class
#' @export
setClass("pert_table",
representation = list(
  tab = "data.table",
  count_modifications = "data.table",
  numvars_modifications = "data.table",
  dimVars = "character",
  countVars = "character",
  numVars = "character",
  by = "character",
  is_weighted = "logical",
  type = "character"
),
prototype=list(
  tab=data.table(),
  count_modifications=data.table(),
  numvars_modifications=data.table(),
  dimVars=character(),
  countVars=character(),
  numVars=character(),
  by=character(),
  is_weighted=c(),
  type=character()),
validity = function(object) {
  stopifnot(object@type %in% c("abs", "destatis", "custom"))
  by <- slot(object, "by")
  if (by != "") {
    stopifnot(length(by) == 1)
    stopifnot(by %in% slot(object, "countVars"))
  }
})
NULL

#' Method to extract results
#' @exportMethod print
#' @rdname print-methods
#' @param x input object
setMethod(f="print", signature="pert_table",
definition=function(x) {
  cat("still to write!\n\n")
})

#' Method to extract information about modifications for count tables
#' @name mod_counts
#' @rdname mod_counts-methods
#' @param x input object
#' @param ... additional parameters passed to specific methods
#' @exportMethod mod_counts
#' @examples
#' ## see example in ?perturbTable
setGeneric("mod_counts", function(x, ...) {
  standardGeneric("mod_counts")
})

#' @rdname mod_counts-methods
#' @param verbose (logical) if \code{TRUE}, additional information is printed
setMethod(f="mod_counts", signature="pert_table",
definition=function(x, verbose=FALSE) {
  type <- slot(x, "type")
  stopifnot(is_scalar_logical(verbose))
  if (verbose) {
    message("The perturbation was based on the ", shQuote(type), " method.")
  }
  cc <- slot(x, "count_modifications")
  if (type=="destatis") {
    cc$col_indices <- NULL
  }
  cc
})

#' Method to extract information about modifications for counts
#' @name mod_numvars
#' @rdname mod_numvars-methods
#' @param x input object
#' @param ... additional parameters passed to methods
#' @exportMethod mod_numvars
setGeneric("mod_numvars", function(x, ...) {
  standardGeneric("mod_numvars")
})
#' @rdname mod_numvars-methods
setMethod(f="mod_numvars", signature="pert_table",
definition=function(x) {
  slot(x, "numvars_modifications")
})
NULL

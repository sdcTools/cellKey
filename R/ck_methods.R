#' Method to extract results
#' @name results
#' @rdname results-methods
#' @param x input object
#' @param ... additional parameters passed to methods
#' @exportMethod results
setGeneric("results", function(x, ...) {
  standardGeneric("results")
})

#' @rdname results-methods
#' @param meanBeforeSum (logical) defines if for (continuous) variables the weighted means should be computed
#' before weighted sums or the other way round.
#' @seealso \url{https://www.unece.org/fileadmin/DAM/stats/documents/ece/ces/ge.46/2013/Topic_1_ABS.pdf}
setMethod(f="results", signature="pert_table",
definition=function(x, meanBeforeSum=TRUE) {
  numVars <- slot(x, "numVars")
  tab <- slot(x, "tab")

  if (length(numVars)>0) {
    for (i in 1:length(numVars)) {
      nV <- numVars[i]
      nVPert <- paste0(nV,".pert")
      if (meanBeforeSum==TRUE) {
        tmp <- mean_before_sum(tab[,get(nVPert)], pWC=tab[, get("pWC")])
        setnames(tmp, paste0(nV,".", names(tmp)))
        tab <- cbind(tab, tmp)
      } else {
        tmp <- sum_before_mean(tab[,get(nVPert)], pWC=tab[, get("pWC")])
        setnames(tmp, paste0(nV,".", names(tmp)))
        tab <- cbind(tab, tmp)
      }
    }
    tab[,c(paste0(numVars,".pert")):=NULL]
  }
  tab[]
})

#' Method to extract results
#' @exportMethod print
#' @rdname print-methods
#' @param x input object
setMethod(f="print", signature="pert_table",
definition=function(x) {
  cat("please use method results()\n")
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
definition=function(x, verbose=TRUE) {
  type <- slot(x, "type")
  stopifnot(is_scalar_logical(verbose))
  if (verbose) {
    message("The perturbation was based on the ", shQuote(type), " method.")
  }
  cc <- slot(x, "count_modifications")
  if (type=="destatis") {
    cc$col_indices <- NULL
  }
  print(cc)
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

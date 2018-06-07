#' Method to extract results
#' @name results
#' @param x input object
#' @param ... additional parameters passed to methods
#' @rdname results-methods
#' @exportMethod results
setGeneric("results", function(x, ...) {
  standardGeneric("results")
})

#' @rdname pert_table-class
#' @aliases results,pert_table-method
#' @param x an object of class \code{\link{pert_table-class}}
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

#' @rdname pert_table-class
#' @aliases print,pert_table-method
setMethod(f="print", signature="pert_table",
definition=function(x) {
  cat("please use method results()\n")
})

#' Method to extract information about modifications for counts
#' @name mod_counts
#' @param x input object
#' @param ... additional parameters passed to methods
#' @rdname mod_counts-methods
#' @exportMethod mod_counts
setGeneric("mod_counts", function(x, ...) {
  standardGeneric("mod_counts")
})
#' @rdname pert_table-class
#' @aliases mod_counts,pert_table-method
setMethod(f="mod_counts", signature="pert_table",
definition=function(x) {
  type <- slot(x, "type")
  message("The perturbation was based on the ", shQuote(type), " method.")

  cc <- slot(x, "count_modifications")
  if (type=="destatis") {
    cc$col_indices <- NULL
  }
  print(cc)
})

#' Method to extract information about modifications for counts
#' @name mod_numvars
#' @param x input object
#' @param ... additional parameters passed to methods
#' @rdname mod_numvars-methods
#' @exportMethod mod_numvars
setGeneric("mod_numvars", function(x, ...) {
  standardGeneric("mod_numvars")
})
#' @rdname pert_table-class
#' @aliases mod_numvars,pert_table-method
setMethod(f="mod_numvars", signature="pert_table",
definition=function(x) {
  slot(x, "numvars_modifications")
})
NULL

#' Method to extract results
#' @exportMethod print
#' @rdname ck_methods
#' @param x input object of class \code{\linkS4class{pert_table}}
#' @param vname if not \code{NULL}, a character specifying a tabulated count
#' or numeric variable for which for all combinations of dimensional variables the 
#' (perturbed) unweighted and (perturbed) weighted counts are printed.
setMethod(f="print", signature="pert_table",
definition=function(x, vname=NULL) {
  if (!is.null(vname)) {
    stopifnot(is_scalar_character(vname))
    validVars <- c(slot(x, "countVars"), slot(x,"numVars"))
    if (!vname %in% validVars) {
      txt <- paste("specified variable name", shQuote(vname),"is invalid.\n")
      txt <- paste0(txt, "It should be one of ", paste(shQuote(validVars), collapse=", "),".")
      stop(txt)
    }
  }
  
  txt <- "The"
  if (slot(x, "is_weighted")) {
    txt <- paste(txt, "weighted")
  } else {
    txt <- paste(txt, "unweighted")
  }
  txt <- paste0(txt, " ",length(slot(x, "dimVars")),"-dimensional table")
  txt <- paste0(txt, " consists of ", nrow(slot(x, "tab")), " cells.")

  if (length(slot(x, "by"))==0) {
    txt <- paste(txt, "The results are based on all units in the input data.")
  } else {
    txt <- paste(txt, "The results are restricted to units for which by-variable")
    txt <- paste(txt, shQuote(slot(x, "by")),"equals 1.")
  }
  cat(txt, "\n")
  
  txt <- "The dimensions are given by the following variables\no"
  txt <- paste(txt, paste(slot(x, "dimVars"), collapse="\no "))
  cat(txt,"\n\n")
  
  txt <- paste("Type of pTable-used:", shQuote(slot(x, "type")))
  cat(txt,"\n")

  txt <- "The following count-variables have been tabulated/perturbed:\no"
  txt <- paste(txt, paste(slot(x, "countVars"), collapse="\no "))
  cat(txt,"\n")

  if (length(slot(x, "numVars"))==0) {
    txt <- "No numeric variables have been tabulated/perturbed in this table"
  } else {
    txt <- "The following numeric variables have been tabulated/perturbed:\no"
    txt <- paste(txt, paste(slot(x, "numVars"), collapse="\no "))
  } 
  cat(txt,"\n")
  
  if (!is.null(vname)) { 
    dt <- slot(x, "tab")
    
    if (vname %in% slot(x, "countVars")) {
      vnames <- c(slot(x, "dimVars"), paste0(c("UWC","WC","pUWC","pWC"),"_", vname))
    } else {
      vnames <- c(slot(x, "dimVars"), paste0(c("UW","WS","pUW","pWS"),"_", vname))
    }
    dt <- dt[,vnames, with=F]
    cat("Total-Counts\n")
    print(dt)
  }
  return(invisible(NULL))
})

#' Method to summarize perturbation results
#' @exportMethod summary
#' @rdname ck_methods
#' @param object input object of class \code{\linkS4class{pert_table}}
setMethod(f="summary", signature="pert_table",
definition=function(object) {
  vals.pert <- pert <- NULL
  cat("Perturbation statistics on count variables:\n")
  info_cnts <- mod_counts(object)
  out_cnts <- info_cnts[, list(
    min=min(pert), 
    max=max(pert), 
    mean=mean(pert), 
    median=median(pert)), by="countVar"]
  print(out_cnts)

  # overall ratios (weighted_sum / perturbed weighted sum)
  cnt_ratios <- rbindlist(lapply(slot(object, "countVars"), function(x) {
    tmp <- suppressMessages(ck_freq_table(object, x))
    v1 <- paste0("WC_",x)
    v2 <- paste0("pWC_",x)
    res <- na.omit((tmp[[v1]]/tmp[[v2]]))
    data.table(countVar=x, 
      min=round(min(res), digits=2),
      max=round(max(res), digits=2),
      mean=round(mean(res), digits=2), 
      median=round(median(res), digits=2))
  }))
  cat("\nRatios weighted counts / perturbed weighted counts:\n")
  print(cnt_ratios)
  
  info_nums <- mod_numvars(object)
  if (nrow(info_nums)>0) {
    out_nums <- info_nums[,list(
      min=min(vals.pert), 
      max=max(vals.pert),
      mean=mean(vals.pert),
      median=median(vals.pert)), by="numVar"]
    cat("\nPerturbation statistics on numerical variables:\n")
    print(out_nums)
    
    # overall ratios (weighted_sum / perturbed weighted sum)
    out_ratios <- rbindlist(lapply(slot(object, "numVars"), function(x) {
      tmp <- suppressMessages(ck_cont_table(object, x))
      v1 <- paste0("WS_",x)
      v2 <- paste0("pWS_",x)
      res <- na.omit((tmp[[v1]]/tmp[[v2]]))
      data.table(numVar=x, 
        min=round(min(res), digits=2),
        max=round(max(res), digits=2),
        mean=round(mean(res), digits=2), 
        median=round(median(res), digits=2))
    }))
    cat("\nRatios weighted sum / perturbed weighted sum for numerical variables:\n")
    print(out_ratios)
  }
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

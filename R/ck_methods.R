#' Method to extract results
#' @exportMethod print
#' @rdname ck_methods
#' @param x input object of class \code{\linkS4class{pert_table}}
#' @param vname if not \code{NULL}, a character specifying a tabulated count
#' or numeric variable for which for all combinations of dimensional variables
#' the (perturbed) unweighted and (perturbed) weighted counts are printed.
setMethod(f = "print", signature = "pert_table",
definition = function(x, vname=NULL) {
  if (!is.null(vname)) {
    stopifnot(is_scalar_character(vname))
    validVars <- c(slot(x, "countVars"), slot(x, "numVars"))
    if (!vname %in% validVars) {
      vv <- shQuote(vname)
      e <- c(
        "The specified variable name",
        shQuote(vname),
        "is invalid.\n",
        "It shoud be one of:",
        paste(shQuote(validVars), collapse = ", ")
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }
  }

  txt <- "The"
  if (slot(x, "is_weighted")) {
    txt <- paste(txt, "weighted")
  } else {
    txt <- paste(txt, "unweighted")
  }
  txt <- paste0(txt, " ", length(slot(x, "dimVars")), "-dimensional table")
  txt <- paste0(txt, " consists of ", nrow(slot(x, "tab")), " cells.")

  if (slot(x, "by") == "Total") {
    txt <- paste(
      txt,
      "The results are based on all units in the input data."
    )
  } else {
    txt <- paste(
      txt,
      "The results are restricted to units for which by-variable"
    )
    txt <- paste(txt, shQuote(slot(x, "by")), "equals 1.")
  }
  cat(txt, "\n")

  txt <- "The dimensions are given by the following variables\no"
  txt <- paste(txt, paste(slot(x, "dimVars"), collapse = "\no "))
  cat(txt, "\n\n")

  txt <- paste("Type of pTable-used:", shQuote(slot(x, "type")))
  cat(txt, "\n")

  txt <- "The following count-variables have been tabulated/perturbed:\no"
  txt <- paste(txt, paste(slot(x, "countVars"), collapse = "\no "))
  cat(txt, "\n")

  if (length(slot(x, "numVars")) == 0) {
    txt <- "No numeric variables have been tabulated/perturbed in this table"
  } else {
    txt <- "The following numeric variables have been tabulated/perturbed:\no"
    txt <- paste(txt, paste(slot(x, "numVars"), collapse = "\no "))
  }
  cat(txt, "\n")

  if (!is.null(vname)) {
    dt <- slot(x, "tab")

    dv <- slot(x, "dimVars")
    if (vname %in% slot(x, "countVars")) {
      vnames <- c(dv, paste0(c("UWC", "WC", "pUWC", "pWC"), "_", vname))
    } else {
      vnames <- c(dv, paste0(c("UW", "WS", "pUW", "pWS"), "_", vname))
    }
    dt <- dt[, vnames, with = FALSE]
    cat("Total-Counts\n")
    print(dt)
  }
  return(invisible(NULL))
})

#' Method to summarize perturbation results
#' @exportMethod summary
#' @rdname ck_methods
#' @param object input object of class \code{\linkS4class{pert_table}}
setMethod(f = "summary", signature = "pert_table",
definition = function(object) {
  vals.pert <- pert <- NULL
  cat("Perturbation statistics count variables:\n")
  info_cnts <- mod_counts(object)

  cnt_info <- info_cnts[, as.list(get_distr_vals(pert)), by = "countVar"]
  print(cnt_info)
  cat("\n")

  # unweighted only!
  cnt_measures <- lapply(slot(object, "countVars"), function(x) {
    ck_cnt_measures(object, vname = x)
  })
  names(cnt_measures) <- slot(object, "countVars")

  for (vv in names(cnt_measures)) {
    cat("Distance-Based measures for count-variable", shQuote(vv), "\n")
    print(cnt_measures[[vv]]$measures)
    cat("\n")
  }

  info_nums <- mod_numvars(object)
  num_info <- num_ratios <- NULL
  if (nrow(info_nums) > 0) {
    num_info <- info_nums[, as.list(get_distr_vals(vals.pert)), by = "numVar"]
    cat("\nPerturbation statistics on numerical variables:\n")
    print(num_info)
  }

  return(invisible(
    list(cnt_info = cnt_info,
    cnt_measures = cnt_measures,
    num_info = num_info)
  ))
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
setMethod(f = "mod_counts", signature = "pert_table",
definition = function(x, verbose=FALSE) {
  type <- slot(x, "type")
  stopifnot(is_scalar_logical(verbose))
  if (verbose) {
    message("The perturbation was based on the ", shQuote(type), " method.")
  }
  cc <- slot(x, "count_modifications")
  if (type == "destatis") {
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
setMethod(f = "mod_numvars", signature = "pert_table",
definition = function(x) {
  slot(x, "numvars_modifications")
})
NULL

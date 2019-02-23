#' ck_cont_table
#'
#' this function allows to return perturbed magnitude-tables from oututs of
#' \code{\linkS4class{pert_table}} objects typically generated
#' using \code{\link{perturbTable}}.
#'
#' @param inp an object of class \code{\linkS4class{pert_table}} generated with \code{\link{perturbTable}}
#' @param vname either \code{NULL} or the name of a variable for which a perturbed
#' frequeny table has been generated. If \code{NULL}, the function prints the variables for
#' which a perturbation has been applied.
#' @param meanBeforeSum (logical) defines if the weighted means should be computed
#' before weighted sums or the other way round.
#' @return a \code{data.table} containing all combinations of the dimensional variables in
#' the first n columns. Afterwards, the following columns are shown:
#' \itemize{
#' \item \code{UW_{vname}}: unweighted sum for variable specified in \code{vname}
#' \item \code{pUW_{vname}}: perturbed unweighted sum for given variable
#' \item \code{WS_{vname}}: weighted sum for given variable
#' \item \code{pWS_{vname}}: perturbed weighted sum for variable specified in \code{vname}
#' \item \code{pWM_{vname}}: perturbed weighted mean for variable specified in \code{vname}
#' }
#' The return argument has a attribute called "modifications" which contains a \code{data.table} holding all
#' the specific perturbation on unit-level for reproducibility. This can be extracted using \code{attributes(ck_cont_table(...), "modifications")}.
#' .
#' @details The computation of the perturbed weighted mean and the perturbed weighted sum depends
#' on argument \code{meanBeforeSum}. If \code{meanBeforeSum} is \code{TRUE}, the perturbed mean is computed first
#' and the perturbed sum is derived from the perturbed mean an the perturbed weighted counts. If
#' \code{meanBeforeSum} is \code{FALSE}, the computation is the other way round. Details are available in the ABS
#' paper \emph{Methodology for the Automatic Confidentialisation of Statistical Outputs from Remote Servers
#' at the Australian Bureau of Statistics} (Thompson, Broadfoot, Elazar).
#' @author Bernhard Meindl
#' @seealso ck_freq_table perturbTable
#' @export
#'
#' @examples
#' # see examples in perturbTable()
ck_cont_table <- function(inp, vname=NULL, meanBeforeSum=TRUE) {
  stopifnot(isS4(inp))
  stopifnot(class(inp) == "pert_table")

  id <- magnitude <- noise <- numVar <- NULL
  vals.mod <- vals.orig <- vals.pert <- pWMean <- pWSum <- NULL

  avail <- slot(inp, "numVars")
  if (is.null(vname)) {
    if (length(avail) == 0) {
      cat("No continous variables have been perturbed\n")
    } else {
      cat("The following continous variables have been perturbed:\n")
      cat(paste("  -->", avail), sep = "\n")
    }
    return(invisible(NULL))
  }
  stopifnot(is_scalar_character(vname))
  stopifnot(vname %in% avail)

  byVar <- slot(inp, "by")
  if (byVar != "Total") {
    message(sprintf(
      "Note: this table is restricted to groups defined by %s!\n",
      shQuote(byVar))
    )
  }

  data <- slot(inp, "tab")
  dt <- cbind(
    data[, slot(inp, "dimVars"), with = FALSE],
    data[, grep(paste0("_", vname), names(data)), with = FALSE])

  # perturbed weighted counts
  pwc <- data[, get(paste0("pWC_", byVar))]
  if (meanBeforeSum == TRUE) {
    out <- mean_before_sum(dt[, get(paste0("pWS_", vname))], pWC = pwc)
  } else {
    out <- sum_before_mean(dt[, get(paste0("pWS_", vname))], pWC = pwc)
  }
  out[is.nan(pWSum), pWSum := 0]
  out[is.nan(pWMean), pWMean := 0]

  set(dt, j = paste0("pWS_", vname), value = out$pWSum)
  set(dt, j = paste0("pWM_", vname), value = out$pWMean)


  mods <- mod_numvars(inp)
  mods <- mods[numVar == vname,
    .(id, magnitude, dir, noise, vals.orig, vals.pert, vals.mod)]
  attr(dt, "modifications") <- mods
  dt
}

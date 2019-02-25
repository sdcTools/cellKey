#' ck_cont_table
#'
#' this function allows to return perturbed magnitude-tables from oututs of
#' \code{\linkS4class{pert_table}} objects typically generated
#' using \code{\link{perturb_table}}.
#'
#' @param inp an object of class \code{\linkS4class{pert_table}} generated
#' with \code{\link{perturb_table}}
#' @param vname either \code{NULL} or the name of a variable for which a
#' perturbed frequeny table has been generated. If \code{NULL}, the function
#' prints the variables for which a perturbation has been applied.
#' @param mean_before_sum (logical) defines if the weighted means should
#' be computed before weighted sums or the other way round.
#' @return a \code{data.table} containing all combinations of the dimensional
#' variables in the first n columns. Afterwards, the following columns
#' are shown:
#' \itemize{
#' \item \code{UW_{vname}}: unweighted sum for variable specified in
#' \code{vname}
#' \item \code{pUW_{vname}}: perturbed unweighted sum for given variable
#' \item \code{WS_{vname}}: weighted sum for given variable
#' \item \code{pWS_{vname}}: perturbed weighted sum for variable specified
#' in \code{vname}
#' \item \code{pWM_{vname}}: perturbed weighted mean for variable specified
#' in \code{vname}
#' }
#' The return argument has a attribute called "modifications" which contains a
#' \code{data.table} holding all the specific perturbation on unit-level for
#' reproducibility. This can be extracted using
#' \code{attributes(ck_cont_table(...), "modifications")}.
#' .
#' @details The computation of the perturbed weighted mean and the perturbed
#' weighted sum depends on argument \code{mean_before_sum}. If
#' \code{mean_before_sum} is \code{TRUE}, the perturbed mean is computed first
#' and the perturbed sum is derived from the perturbed mean an the perturbed
#' weighted counts. If \code{mean_before_sum} is \code{FALSE}, the computation
#' is the other way round. Details are available in the ABS paper
#' \emph{Methodology for the Automatic Confidentialisation of Statistical
#' Outputs from Remote Servers at the Australian Bureau of Statistics}
#' (Thompson, Broadfoot, Elazar).
#' @author Bernhard Meindl
#' @seealso ck_freq_table perturb_table
#' @export
#'
#' @examples
#' # see examples in perturb_table()
ck_cont_table <- function(inp, vname=NULL, mean_before_sum=TRUE) {
  stopifnot(isS4(inp))
  stopifnot(class(inp) == "pert_table")

  . <- id <- magnitude <- noise <- numvar <- NULL
  vals_mod <- vals_orig <- vals_pert <- pw_mean <- pw_sum <- NULL

  avail <- slot(inp, "numvars")
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

  by_var <- slot(inp, "by")
  if (by_var != "total") {
    message(sprintf(
      "Note: this table is restricted to groups defined by %s!\n",
      shQuote(by_var))
    )
  }

  data <- slot(inp, "tab")
  dt <- cbind(
    data[, slot(inp, "dimvars"), with = FALSE],
    data[, grep(paste0("_", vname), names(data)), with = FALSE]
  )

  # perturbed weighted counts
  pwc <- data[, get(paste0("pwc_", by_var))]
  if (isTRUE(mean_before_sum)) {
    out <- .mean_before_sum(dt[, get(paste0("pws_", vname))], pwc = pwc)
  } else {
    out <- .sum_before_mean(dt[, get(paste0("pws_", vname))], pwc = pwc)
  }
  out[is.nan(pw_sum), pw_sum := 0]
  out[is.nan(pw_mean), pw_mean := 0]

  set(dt, j = paste0("pws_", vname), value = out$pw_sum)
  set(dt, j = paste0("pwm_", vname), value = out$pw_mean)

  mods <- mod_numvars(inp)
  mods <- mods[numvar == vname,
    .(id, magnitude, dir, noise, vals_orig, vals_pert, vals_mod)]
  attr(dt, "modifications") <- mods
  dt
}

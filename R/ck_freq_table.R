#' ck_freq_table
#'
#' this function allows to return perturbed count-tables from oututs of
#' \code{\linkS4class{pert_table}} objects typically generated
#' using \code{\link{perturbTable}}.
#'
#' @param inp an object of class \code{\linkS4class{pert_table}} generated with \code{\link{perturbTable}}
#' @param vname either \code{NULL} or the name of a variable for which a perturbed
#' frequeny table has been generated. If \code{NULL}, the function prints the variables for
#' which a perturbation has been applied.
#'
#' @return a \code{data.table} containing all combinations of the dimensional variables in
#' the first n columns. Afterwards, the following columns are shown:
#' \itemize{
#' \item \code{UWC_{vname}}: unweighted counts for variable specified in \code{vname}
#' \item \code{WC_{vname}}: weighted counts for given variable
#' \item \code{pUWC_{vname}}: perturbed unweighted counts for given variable
#' \item \code{pWC_{vname}}: perturbed weighted counts for given variable
#' \item \code{WCavg_{vname}}: average cell weight for variable specified in \code{vname}
#' \item \code{row_indices}: row index in respect to the perturbation table used
#' \item \code{col_indices}: column index in respect to the perturbation table used. This
#' value is always \code{NA} in case "destatis"-type lookup tables have been used.
#' \item \code{pert}: actual amount of perturbation \code{vname}
#' \item \code{cellKey}: actual cellKey used in the lookup-procedure \code{vname}
#' }
#' @author Bernhard Meindl
#' @seealso ck_cont_table perturbTable
#' @export
#'
#' @examples
#' ## see examples in perturbTable()
ck_freq_table <- function(inp, vname=NULL) {
  stopifnot(isS4(inp))
  stopifnot(class(inp)=="pert_table")

  . <- col_indices <- countVar <- row_indices <- pert <- cellKey <- NULL

  avail <- slot(inp, "countVars")
  if (is.null(vname)) {
    cat("The following variables have been perturbed using the frequency-approach:\n")
    cat(paste("  -->", avail), sep="\n")
    return(invisible(NULL))
  }
  stopifnot(is_scalar_character(vname))
  stopifnot(vname %in% avail)

  data <- slot(inp, "tab")

  vv <- slot(inp, "dimVars")

  dt <- cbind(
    data[,slot(inp, "dimVars"), with=F],
    data[,grep(vname, names(data)), with=F])

  # mods <-
  type <- slot(inp, "type")
  mods <- mod_counts(inp)
  if (type=="destatis") {
    mods[,col_indices:=NA]
  }
  mods <- mods[countVar==vname, .(row_indices, col_indices, pert, cellKey)]

  dt <- cbind(dt, mods)
  dt[]
}

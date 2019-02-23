#' Method to extract information about modifications for count tables
#' @name ck_export_table
#' @rdname ck_export_table
#' @param x input object
#' @param vname (character) variable name of a perturbed count variable or a continous variable for which the output
#' should be prepared
#' @param type (character) specify the type of output. Allowed values are:
#' \itemize{
#' \item \code{"both"}: unweighted and weighted original and perturbed counts/values are returned
#' \item \code{"unweighted"}: unweighted original and perturbed counts/values are returned
#' \item \code{"weighted"}: weighted and weighted original and perturbed counts/values are returned
#' }
#' @return a \code{data.table} containing all dimensional variables and, depending
#' on the value of argument \code{type} and weather if \code{vname} is a count
#' or a continuously scaled variable, the following columns: For count-variables:
#' \itemize{
#' \item{UWC}: unweighted counts (argument \code{type} is "both" or "unweighted")
#' \item{WC}: perturbed unweighted counts (argument \code{type} is "both" or "unweighted")
#' \item{pUWC}: weighted counts (argument \code{type} is "both" or "weighted")
#' \item{pWC}: perturbed weighted counts (argument \code{type} is "both" or "weighted")
#' }
#' For continuously-scaled variables (magnitude tables):
#' \itemize{
#' \item{UW}: unweighted sum (argument \code{type} is "both" or "unweighted")
#' \item{WS}: perturbed unweighted sum (argument \code{type} is "both" or "unweighted")
#' \item{pUW}: weighted sum (argument \code{type} is "both" or "weighted")
#' \item{pWS}: perturbed weighted sum (argument \code{type} is "both" or "weighted")
#' }
#' @export
#' @examples
#' ## see example in ?perturbTable
ck_export_table <- function(x, vname=NULL, type="both") {
  stopifnot(is_scalar_character(type))
  stopifnot(type %in% c("both", "weighted", "unweighted"))

  allvars <- c(slot(x, "countVars"), slot(x, "numVars"))
  if (!is.null(vname)) {
    stopifnot(is_scalar_character(vname))
    if (!vname %in% allvars) {
      e <- c(
        "The variable was not found in the perturbed dataset.\n",
        "Possible variables are: ",
        paste(shQuote(allvars), collapse = ", ")
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }
  }

  tab <- slot(x, "tab")

  names_count <- function(v, type) {
    if (type == "both") {
      return(c(paste0(c("UWC", "WC", "pUWC", "pWC"), "_", v)))
    }
    if (type == "unweighted") {
      return(c(paste0(c("UWC", "pUWC"), "_", v)))
    }
    if (type == "weighted") {
      return(c(paste0(c("WC", "pWC"), "_", v)))
    }
    stop("non-allowed value in 'type'")

  }
  names_cont <- function(v, type) {
    if (type == "both") {
      return(c(paste0(c("UW", "WS", "pUW", "pWS"), "_", v)))
    }
    if (type == "unweighted") {
      return(c(paste0(c("UW", "pUW"), "_", v)))
    }
    if (type == "weighted") {
      return(c(paste0(c("WS", "pWS"), "_", v)))
    }
    stop("non-allowed value in 'type'")
  }

  if (vname %in% slot(x, "countVars")) {
    get_vars <- names_count(vname, type = type)
  } else {
    get_vars <- names_cont(vname, type = type)
  }

  vdim <- slot(x, "dimVars")
  tt <- tab[, vdim, with = FALSE]
  t2 <- tab[, get_vars, with = FALSE]
  names(t2) <- gsub(paste0("_", vname), "", names(t2))
  tt[, vname := vname]
  tt <- cbind(tt, t2)
  return(tt)
}

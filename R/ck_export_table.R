#' Create tabular output from perturbed table
#'
#' This function allows to extract information from a perturbed table
#' as a \code{data.table}.
#'
#' @name ck_export_table
#' @rdname ck_export_table
#' @param x input object
#' @param vname (character) variable name of a perturbed count variable or a
#' continous variable for which the output should be prepared
#' @param type (character) specify the type of output. Allowed values are:
#' \itemize{
#' \item \code{"both"}: unweighted and weighted original and perturbed
#' counts/values are returned
#' \item \code{"unweighted"}: unweighted original and perturbed counts/values
#' are returned
#' \item \code{"weighted"}: weighted and weighted original and perturbed
#' counts/values are returned
#' }
#' @return a \code{data.table} containing all dimensional variables and,
#' depending on the value of argument \code{type} and weather if \code{vname}
#' is a count or a continuously scaled variable, the following columns: For
#' count-variables:
#' \itemize{
#' \item{uwc}: unweighted counts (argument \code{type} is "both" or
#' "unweighted")
#' \item{wc}: perturbed unweighted counts (argument \code{type} is "both" or
#' "unweighted")
#' \item{puwc}: weighted counts (argument \code{type} is "both" or "weighted")
#' \item{pwc}: perturbed weighted counts (argument \code{type} is "both" or
#' "weighted")
#' }
#' For continuously-scaled variables (magnitude tables):
#' \itemize{
#' \item{uw}: unweighted sum (argument \code{type} is "both" or "unweighted")
#' \item{ws}: perturbed unweighted sum (argument \code{type} is "both" or
#' "unweighted")
#' \item{puw}: weighted sum (argument \code{type} is "both" or "weighted")
#' \item{pws}: perturbed weighted sum (argument \code{type} is "both" or
#' "weighted")
#' }
#' @export
#' @examples
#' ## see example in ?perturb_table
ck_export_table <- function(x, vname=NULL, type="both") {
  stopifnot(is_scalar_character(type))
  stopifnot(type %in% c("both", "weighted", "unweighted"))

  allvars <- c(slot(x, "countvars"), slot(x, "numvars"))
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

  .n_freq <- function(v, type) {
    if (type == "both") {
      return(c(paste0(c("uwc", "wc", "puwc", "pwc"), "_", v)))
    }
    if (type == "unweighted") {
      return(c(paste0(c("wuc", "puwc"), "_", v)))
    }
    if (type == "weighted") {
      return(c(paste0(c("wc", "pwc"), "_", v)))
    }
  }
  .n_cont <- function(v, type) {
    if (type == "both") {
      return(c(paste0(c("uw", "ws", "puw", "pws"), "_", v)))
    }
    if (type == "unweighted") {
      return(c(paste0(c("uw", "puw"), "_", v)))
    }
    if (type == "weighted") {
      return(c(paste0(c("ws", "pws"), "_", v)))
    }
  }

  if (vname %in% slot(x, "countvars")) {
    get_vars <- .n_freq(vname, type = type)
  } else {
    get_vars <- .n_cont(vname, type = type)
  }

  vdim <- slot(x, "dimvars")
  tt <- tab[, vdim, with = FALSE]
  t2 <- tab[, get_vars, with = FALSE]
  names(t2) <- gsub(paste0("_", vname), "", names(t2))
  tt[, vname := vname]
  tt <- cbind(tt, t2)
  return(tt)
}

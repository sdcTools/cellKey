#' Create perturbation parameters for count variables
#'
#' This function allows to generate required perturbation parameters that are used
#' to perturb count variables.
#'
#' @inheritParams ck_params_nums
#' @param ptab an object created with [ptable::create_ptable()],
#' or [ptable::create_cnt_ptable()]
#' @return an object suitable as input to method `$params_cnts_set()` for the perturbation
#' of counts and frequencies.
#' @export
#' @md
#' @seealso This function uses functionality from package
#' `ptable` (https://github.com/sdcTools/ptable), expecially
#' [ptable::create_ptable()] and
#' [ptable::create_cnt_ptable()]. More detailed information on the parameters
#' is available from the respective help-pages of these functions.
#' @inherit cellkey_pkg examples
ck_params_cnts <- function(ptab, path = NULL) {
  if (inherits(ptab, "ptable_params")) {
    if (ptab@table != "cnts") {
      e <- "input from `ptable::create_ptable` or `ptable::create_num_ptable` not suitable for count variables."
      stop(e, call. = FALSE)
    }
    ptab <- ptable::create_ptable(params = ptab)
  }
  out <- list(
    params = list(ptable = .chk_ptab(ptab, type = "cnts")),
    type = "cnts")
  class(out) <- "ck_params"

  if (!is.null(path)) {
    out$version <- paste(utils::packageVersion("cellKey"), collapse = ".")
    out$ptype <- "params_cnts"
    .yaml_write(x = out, path = path)
    message("yaml configuration ", shQuote(path), " successfully written.")
  }
  out$version <- out$ptype <- NULL
  return(out)
}

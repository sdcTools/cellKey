#' Create perturbation parameters for count variables
#'
#' This function allows to generate required perturbation parameters that are used
#' to perturb count variables.
#'
#' @param D perturbation parameter for maximum perturbation (scalar or vector)
#' @param V perturbation parameter for variance (scalar)
#' @param js treshold value for blocking of small frequencies (i.e. there won't occur
#' positive target frequencies below the treshold value)
#' @param pstay optional parameter to set
#' @param optim optimization parameter: `1` standard approach (default)
#' @param mono (logical) vector specifying optimization parameter for monotony condition
#' @param epsilon (double)
#'
#' @return an object suitable as input to [cellkey_pkg()] for the perturbation
#' of counts and frequencies.
#' @export
#' @md
#' @seealso This function uses functionality from package
#' `ptable` (https://github.com/sdcTools/ptable), expecially [ptable::pt_create_pParams()]
#' and [ptable::pt_create_pTable()]. More detailed information on the parameters
#' is available from the respective help-pages of these functions.
#' @inherit cellkey_pkg examples
ck_params_cnts <- function(D, V, js=0, pstay=NULL, optim=1, mono=TRUE, epsilon=1e-07) {
  pt_para <- ptable::pt_create_pParams(
    D = D,
    V = V,
    js = js,
    pstay = pstay,
    optim = optim,
    mono = mono,
    epsilon = epsilon
  )

  ptable = ptable::pt_create_pTable(
    params = pt_para,
    type = "destatis",
    monitoring = FALSE,
    debugging = FALSE
  )

  out <- list(
    params = list(
      ptable = ptable
    ),
    type = "cnts"
  )
  class(out) <- "ck_params"
  return(out)
}

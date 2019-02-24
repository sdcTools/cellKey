#' ck_create_pTable
#'
#' generates perturbation tables using functionality from package ptable
#' @param D perturbation parameter for maximum perturbation (scalar or vector)
#' @param V perturbation parameter for variance (scalar)
#' @param js treshold value for blocking of small frequencies (i.e. there won't occur
#' positive target frequencies below the treshold value)
#' @param pstay optional parameter to set
#' @param optim optimization parameter: \code{1} standard approach (default)
#' @param mono (logical) vector specifying optimization parameter for monotony condition
#' @param epsilon (double)
#' @param pTableSize (number) defining the number of required columns in the "abs"-approach
#' @param type character specifying the format of the input table. Valid choices
#' are \code{"abs"} and \code{"destatis"}.
#'
#' @return an object of class \code{ptable} defined in the ptable-package.
#' @export
#' @seealso This function uses functionality from package ptable (https://github.com/sdcTools/ptable), expecially
#' \code{pt_create_pParams} and \code{pt_create_pTable}. More detailed information on the parameters is available
#' from the respective help-pages of these functions.
#' @examples
#' ck_create_pTable(D=5, V=3, js=2, pTableSize=70, type="abs")
#' ck_create_pTable(D=5, V=3, js=2, pstay=0.5, optim=1, mono=TRUE, typ="destatis")
ck_create_pTable <-
  function(D,
           V,
           type,
           js = 0,
           pstay = NULL,
           optim = 1,
           mono = TRUE,
           epsilon = 1e-07,
           pTableSize = 70) {
    stopifnot(is_scalar_character(type))
    type <- tolower(type)
    stopifnot(type %in% c("abs", "destatis"))

  pt_para <- pt_create_pParams(
    D = D,
    V = V,
    js = js,
    pstay = pstay,
    optim = optim,
    mono = mono,
    epsilon = epsilon,
    pTableSize = pTableSize
  )
  return(pt_create_pTable(
    params = pt_para,
    type = type,
    monitoring = FALSE,
    debugging = FALSE
  ))
}

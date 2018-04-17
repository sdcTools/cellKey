#' ck_create_input
#'
#' creates the required input for \code{\link{perturbTable}}.
#'
#' @param dat a \code{data.table} containing micro data
#' @param def_rkey either a column name within \code{dat} specifying a variable containing record
#' keys or a single number used as upper bound to compute record keys between 1 and def_rkey
#' @param pert_params information on perturbation parameters created with \code{\link{ck_create_pert_params}}
#'
#' @return an object of class \code{\link{pert_inputdat-class}}
#' @export
#'
#' @examples
#' ## define parameters
#' bigN <- 20011
#' smallN <- 12
#' maxV <- 214748365 #(2^31)-1
#' dat <- ck_create_testdata()
#' dat$rkeys <- ck_generate_rkeys(dat=dat, max_val=maxV)
#' mTable <- c(0.6,0.4,0.2)
#' sTable <- ck_generate_sTable(smallC=12)
#' pTable <- ck_create_pTable(pTableSize=70)
#'
#' pert_params <- ck_create_pert_params(
#'   bigN=bigN, smallN=smallN, pTable=pTable, sTable=sTable, mTable=mTable)
#'
#' ## create suitable input data using exising record keys
#' inp <- ck_create_input(dat=dat, def_rkey="rkeys", pert_params=pert_params)
#' print(class(inp))
ck_create_input <- function(dat, def_rkey, pert_params) {
  out <- new("pert_inputdat")
  dat <- as.data.table(dat)

  if (is.character(def_rkey)) {
    stopifnot(is_scalar_character(def_rkey))
    stopifnot(def_rkey %in% names(dat))
    out@rkeys <- dat[[def_rkey]]
  } else if (is_bare_integerish(def_rkey)) {
    stopifnot(is_scalar_atomic(def_rkey))
    stopifnot(def_rkey>=1)
    out@rkeys <- ck_generate_rkeys(dat=dat, max_val=ceiling(def_rkey))
  } else {
    stop("Argument",shQuote("def_rkey"),"must either be character or numeric!\n")
  }

  out@microdat <- dat
  out@pert_params <- pert_params
  validObject(out)
  out
}

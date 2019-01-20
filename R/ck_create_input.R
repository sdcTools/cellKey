#' ck_create_input
#'
#' creates the required input for \code{\link{perturbTable}}.
#'
#' @param dat a \code{data.table} containing micro data
#' @param def_rkey either a column name within \code{dat} specifying a variable
#' containing record keys or a single (integer) number.
#'
#' If it is a number, we have to distinguish two cases:
#'
#' \itemize{
#' \item The perturbation table is provided in "abs"-format:
#'
#' \code{def_rkey} is used as upper bound to compute record keys between 1 and \code{def_rkey}
#' \item The perturbation table is provided in "destatis"-format:
#'
#' \code{def_rkey} is used to specify the maximum number of digits for record keys
#' that are generated from a uniform distribution between 0 and 1.
#' }
#' @param pert_params information on perturbation parameters created
#' using \code{\link{ck_create_pert_params}}
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
#' pTable <- ck_create_pTable(D=5, V=3, pTableSize=70, type="abs")
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

  type <- slot(pert_params, "type")

  if (is.character(def_rkey)) {
    stopifnot(is_scalar_character(def_rkey))
    stopifnot(def_rkey %in% names(dat))

    # check the record keys
    if (check_rkeys(rkeys = dat[[def_rkey]], type = type)) {
      out@rkeys <- dat[[def_rkey]]
    } else {
      stop("invalid records detected in the dataset!\n")
    }
  } else if (is_bare_integerish(def_rkey)) {
    stopifnot(is_scalar_atomic(def_rkey))
    if (type == "abs") {
      slot(out, "rkeys") <- ck_generate_rkeys(
        dat = dat,
        max_val = def_rkey,
        type = type)
    } else {
      slot(out, "rkeys") <- ck_generate_rkeys(
        dat = dat,
        max_digits = def_rkey,
        type = type)
    }
  } else {
    stop(
      sprintf("Argument %s must either be character or numeric!\n",
        shQuote("def_rkey"))
    )
  }

  slot(out, "microdat") <- dat
  slot(out, "pert_params") <- pert_params
  validObject(out)
  out
}

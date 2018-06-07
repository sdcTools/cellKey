#' ck_create_pTable
#'
#' generates a simple perturbation table for developing/testing purposes.
#' Real-world lookup tables will be likely created in another package/software.
#'
#' @param pTableSize (number) defining the number of required columns
#' @param type character specifying the format of the input table. Valid choices are \code{"abs"} and \code{"destatis"}.
#' @param ... additional parameters (currently ignored)
#'
#' @return a \code{data.table} with 256 rows and the specified number of
#' columns given by parameter \code{pTableSize} holding perturbation values in case argument \code{type} was set to \code{"abs"}.
#' \code{data.table} with three columns and 66 rows in case \code{type} was set to \code{"destatis"}.
#' @export
#'
#' @examples
#' ck_create_pTable(pTableSize=10, type="abs")
ck_create_pTable <- function(pTableSize=75, type="abs", ...) {
  gen_test_ptable_abs <- function(pTableSize, ...) {
    V1 <- V2 <- V3 <- NULL

    nrows <- 256
    dt <- as.data.table(matrix(sample(-10:10, 256*pTableSize, replace=TRUE), nrow=256, ncol=pTableSize))

    # Cells <=3 will be perturbed to 0!
    dt[,V1:=-1]
    dt[,V2:=-2]
    dt[,V3:=-3]

    setattr(dt, "type", "abs")
    dt[]
  }
  gen_test_ptable_destatis <- function() {
    dt <- copy(ptable_destatis)
    dt[,c("j","kum_p_u","p"):=NULL]
    setattr(dt, "type", "destatis")
    dt[]
  }

  stopifnot(is.character(type))
  type <- tolower(type)
  stopifnot(length(type)==1)
  stopifnot(type %in% c("abs","destatis"))

  if (type=="abs") {
    return(gen_test_ptable_abs(pTableSize=pTableSize, ...))
  }
  if (type=="destatis") {
    message(paste("Note: argument",shQuote("pTableSize"),"is ignored!"))
    return(gen_test_ptable_destatis())
  }
}

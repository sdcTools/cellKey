#' ck_create_pTable
#'
#' generates a simple perturbation table for developing/testing purposes.
#' Real-world lookup tables will be likely created in another package/software.
#'
#' @param pTableSize (number) defining the number of required columns
#' @param ... additional parameters (currently ignored)
#'
#' @return a \code{data.table} with 256 rows and the specified number of
#' columns given by parameter \code{pTableSize} holding perturbation values.
#' @export
#'
#' @examples
#' ck_create_pTable(pTableSize=10)
ck_create_pTable <- function(pTableSize=75, ...) {
  V1 <- V2 <- V3 <- NULL

  nrows <- 256
  dt <- as.data.table(matrix(sample(-10:10, 256*pTableSize, replace=TRUE), nrow=256, ncol=pTableSize))

  # Cells <=3 will be perturbed to 0!
  dt[,V1:=-1]
  dt[,V2:=-2]
  dt[,V3:=-3]
  dt[]
}

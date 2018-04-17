#' ck_create_testdata
#'
#' this function generates some test-data
#'
#' @return a \code{data.frame}
#' @export
#'
#' @examples
#' dat <- ck_create_testdata(); print(str(dat))
ck_create_testdata <- function() {
  testdata <- sampling_weight <- savings <- expend <- income <- NULL

  data(testdata, envir = environment())

  microdat <- as.data.table(testdata)
  microdat[, sampling_weight:=sample(20:100, nrow(microdat), replace=TRUE)]
  microdat[,savings:=round(savings/10000)]
  microdat[,expend:=round(expend/10000)]
  microdat[,income:=round(income/10000)]
  return(microdat)
}

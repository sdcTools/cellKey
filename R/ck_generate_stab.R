#' Generate a perturbation table for numerical variables
#'
#' This function allows to create a perturbation table (differently for for
#' large and small cell values) for testing purposes.
#'
#' @param small_c number of columns for the perturbation table for large cells
#' @param ...  additional parameters, not yet used
#'
#' @return a \code{data.table} with 256 rows and 32 + \code{small_c} columns
#' @export
#'
#' @examples
#' ck_generate_stab(small_c = 10)
ck_generate_stab <- function(small_c=12, ...) {
  nr <- 256
  nc1 <- 32
  nc2 <- small_c

  # for "large" cells
  s <- seq(from = 0.3, to = 1.7, by = 0.01)
  dt1 <- as.data.table(
    matrix(
      data = sample(s, size = nr * nc1, replace = TRUE),
      nrow = nr,
      ncol = nc1
    )
  )

  # for "small" cells
  s <- seq(from = 0.7, to = 1.3, by = 0.01)
  dt2 <- as.data.table(
    matrix(
      data = sample(s, size = nr * nc2, replace = TRUE),
      nrow = nr,
      ncol = nc2
    )
  )

  dt <- cbind(dt1, dt2)
  setnames(dt, paste0("V", 1:ncol(dt)))
  dt[]
}

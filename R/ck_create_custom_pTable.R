#' Create a custom pTable
#'
#' This function allows to create a custom perturbation table (pTable) that contains 256 rows and
#' \code{pTableSize} columns. Each cell of this perturbation table must be a function that returns
#' a single number. By default, each cell gets initialized with a function that always returns
#' \code{0}, so no perturbation would be applied.
#'
#' To actually modify a custom pTable, function \code{\link{ck_update_custom_pTable}} needs
#' to be used.
#'
#' @param pTableSize an integer specifying the number of columns for the perturbation table
#'
#' @return an object of class \code{ptable} defined in the ptable-package.
#' @export
#' @author Bernhard Meindl
#' @examples
#' # initialize
#' pt_custom <- ck_create_custom_pTable(pTableSize=50)
#'
#' # modify
#' fn1 <- function() round(rnorm(1, mean = 5, sd = 10))
#' fn2 <- function() rpois(1, lambda=5)
#' fn3 <- function() return(1)
#' fn4 <- function() return(-1)
#'
#' # fn1 provides perturbation values from a normal distribution
#' # with mean = 5 and sd = 10
#' # we use this for all cells
#' pt_custom <- ck_update_custom_pTable(pt_custom, fun = fn1)
#'
#' # perturbation values from poisson-distribution with lambda = 5
#' # for some cells (rows)
#' pt_custom <- ck_update_custom_pTable(
#'   pTable = pt_custom,
#'   fun = fn2,
#'   cols = 1:5,
#'   rows = 1:20
#' )
#'
#' # we can of course write functions, that return scalars, such
#' # as fn3() (always returns 1) or fn4() (always returns -1)
#' pt_custom <- ck_update_custom_pTable(
#'   pTable = pt_custom,
#'   fun = fn3,
#'   cols = 10:20
#' )
#' pt_custom <- ck_update_custom_pTable(
#'   pTable = pt_custom,
#'   fun = fn4,
#'   cols = 21:30
#' )
#'
#' ## see also the example in ?perturbTable
ck_create_custom_pTable <- function(pTableSize=70) {
  stopifnot(is_scalar_integerish(pTableSize))
  stopifnot(pTableSize > 0)

  out <- list(); length(out) <- 256
  tmp <- list(); length(tmp) <- pTableSize

  out <- lapply(1:length(out), function(x) {
    out[[x]] <- tmp
  })

  # no perturbation by default
  fn <- function() return(0)

  for (i in 1:length(out)) {
    for (j in 1:length(out[[i]])) {
      out[[i]][[j]] <- fn
    }
  }

  nr <- length(out)
  nc <- length(out[[1]])
  dt <- data.table(matrix(NA_character_, nrow = nr, ncol = nc))
  for (i in 1:ncol(dt)) {
    dt[[i]] <- lapply(out, function(x) x[[i]])
  }

  # prepare output
  ret <- new("ptable")
  params <- new("ptable_params")
  slot(params, "pTableSize") <- as.integer(pTableSize)
  slot(params, "label") <- "custom perturbation"
  slot(ret, "pParams") <- params
  slot(ret, "pTable") <- dt
  slot(ret, "type") <- "custom"
  slot(ret, "tStamp") <- format(Sys.time(), "%Y%m%d%H%M%S")
  return(ret)
}

#' Modify a custom perturbation table
#'
#' @param pTable an object created using \code{\link{ck_create_custom_pTable}}
#' @param fun a function that will be used for given cells to derive the
#' perturbation values. The function must return exactly one number when called
#' without parameters.
#' @param cols numeric index specifying columns for which the given perturbation
#' function should be set. If \code{NULL}, the function will be set on all columns
#' of the perturbation table.
#' @param rows numeric index specifying rows for which the given perturbation
#' function should be set. If \code{NULL}, the function will be set on all
#' rows of the perturbation table.
#'
#' @return an object of class \code{ptable} defined in the ptable-package.
#' @export
#' @author Bernhard Meindl
#' @export
#' @examples
#' ## see examples in ?ck_create_custom_pTable and ?perturbTable
ck_update_custom_pTable <- function(pTable, fun, cols=NULL, rows=NULL) {
  valid_custom_ptable(pTable)
  stopifnot(is.function(fun))

  tab <- slot(pTable, "pTable")
  nr <- nrow(tab)
  nc <- ncol(tab)

  if (!is.null(cols)) {
    stopifnot(is_integerish(cols), all(cols %in% 1:nc))
  } else {
    cols <- 1:ncol(tab)
  }

  if (!is.null(rows)) {
    stopifnot(is_integerish(rows), all(rows %in% 1:nr))
  }

  # check function
  rr <- fun()
  stopifnot(is_scalar_integerish(rr))

  for (column in cols) {
    new <- tab[[column]]

    if (!is.null(rows)) {
      for (row in rows) {
        new[[row]] <- fun
      }
    } else {
      # set for entire column
      new <- lapply(1:nr, function(x) {
        new[[x]] <- fun
      })
    }
    tab[[column]] <- new
  }
  slot(pTable, "pTable") <- tab
  validObject(pTable)
  check_custom_pTable(pTable)
  return(pTable)
}

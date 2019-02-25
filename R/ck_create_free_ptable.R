#' Create a perturbation table in `free` format
#'
#' This function allows to create a perturbation table in "free" format
#' that contains 256 rows and \code{ncol} columns. Each cell of this
#' perturbation table must be a function that returns a single number.
#' By default, each cell gets initialized with a function that always
#' returns \code{0}, so no perturbation is be applied.
#'
#' To actually modify such an perturbation table, function
#' \code{\link{ck_update_free_ptable}} needs to be used.
#'
#' @param nrcol an integer specifying the number of columns for
#' the perturbation table
#'
#' @return an object of class \code{ptable} as defined in the ptable-package.
#' @export
#' @author Bernhard Meindl
#' @examples
#' # initialize
#' pt_free <- ck_create_free_ptable(nrcol=50)
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
#' pt_free <- ck_update_free_ptable(
#'   ptab = pt_free,
#'   fun = fn1
#' )
#'
#' # perturbation values from poisson-distribution with lambda = 5
#' # for some cells (rows)
#' pt_free <- ck_update_free_ptable(
#'   ptab = pt_free,
#'   fun = fn2,
#'   cols = 1:5,
#'   rows = 1:20
#' )
#'
#' # we can of course write functions, that return scalars, such
#' # as fn3() (always returns 1) or fn4() (always returns -1)
#' pt_free <- ck_update_free_ptable(
#'   ptab = pt_free,
#'   fun = fn3,
#'   cols = 10:20
#' )
#' pt_free <- ck_update_free_ptable(
#'   ptab = pt_free,
#'   fun = fn4,
#'   cols = 21:30
#' )
#'
#' ## see also the example in ?perturb_table
ck_create_free_ptable <- function(nrcol=70) {
  stopifnot(is_scalar_integerish(nrcol))
  stopifnot(nrcol > 0)

  out <- list(); length(out) <- 256
  tmp <- list(); length(tmp) <- nrcol

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
  slot(params, "pTableSize") <- as.integer(nrcol)
  slot(params, "label") <- "free_perturbation"
  slot(ret, "pParams") <- params
  slot(ret, "pTable") <- dt
  slot(ret, "type") <- "free"
  slot(ret, "tStamp") <- format(Sys.time(), "%Y%m%d%H%M%S")
  return(ret)
}

#' Modify a perturbation table in "free" format
#'
#' This function allows to change and update a perturbation table that
#' is in "free" format.
#'
#' @param ptab an object created using \code{\link{ck_create_free_ptable}}
#' @param fun a function that will be used for given cells to derive the
#' perturbation values. The function must return exactly one number when called
#' without parameters.
#' @param cols numeric index specifying columns for which the given perturbation
#' function should be set. If \code{NULL}, the function will be set on all
#' columns of the perturbation table.
#' @param rows numeric index specifying rows for which the given perturbation
#' function should be set. If \code{NULL}, the function will be set on all
#' rows of the perturbation table.
#'
#' @return an object of class \code{ptable} defined in the ptable-package.
#' @export
#' @author Bernhard Meindl
#' @export
#' @examples
#' ## see examples in ?ck_create_free_ptable and ?perturb_table
ck_update_free_ptable <- function(ptab, fun, cols=NULL, rows=NULL) {
  .valid_free_ptable(ptab)
  stopifnot(is.function(fun))

  tab <- slot(ptab, "pTable")
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
  slot(ptab, "pTable") <- tab
  validObject(ptab)
  .check_free_ptable(ptab)
  return(ptab)
}

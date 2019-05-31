gen_stab <- function(D = 3, l = 0.5) {
  .gen <- function(i, min_d, max_d, l) {
    pval <- seq(min_d, max_d, by = l)
    x <- data.table(
      i = i,
      j = pval - min(pval),
      p = NA_real_,
      kum_p_u = NA_real_,
      kum_p_o = NA_real_,
      diff = pval
    )

    x$p <- diff(sort(runif(nrow(x) + 1)))
    cs <- cumsum(x$p)
    x$kum_p_u <- c(0, cs)[1:nrow(x)]
    x$kum_p_o <- c(cs[1:(nrow(x) - 1)], 1)
    x
  }

  dt <- data.table(i = 0, j = 0, p = 1, kum_p_u = 0, kum_p_o = 1, diff = 0)
  # i:=1
  dt_a <- .gen(i = 1, min_d = -1, max_d = D, l = l)
  # i == D
  dt_b <- .gen(i = D, min_d = -D, max_d = D, l = l)
  dt <- rbind(dt, dt_a, dt_b)
  dt
}

#' Create perturbation parameters for continuous variables
#'
#' This function allows to define perturbation parameters used to
#' perturb cells in magnitude tables.
#'
#' @param P a character value defining the way to identify the `magnifier`,
#' e.g which contributions/values in a cell should be the used  in the
#' perturbation procedure. Possible choices are:
#' - `top_contr`: the `k` largest contributions are used. In this case,
#' it is also required to specify argument `top_k`
#' - `mean`: the (weighted) cellmean is used as starting point
#' - `range`: the difference between largest and smallest contribution
#' is used.
#' - `sum`: the (weighted) cellvalue itself is used as starting point
#' @param top_k it is ignored if `P` is different from `top_contr`. Otherwise,
#' @param D maximum perturbation value
#' @param l stepwidth parameter for computation of perturbation tables
#' a integerish number >= 1 specifying the number of top contributions whose
#' values should be perturbed.
#' @param same_key (logical) should original cell key (`TRUE`) used for
#' for finding perturbation values of the largest contributor to a
#' cell or should a perturbation of the cellkey itself (`FALSE`) take place.
#' @param use_zero_rkeys (logical) skalar defining if record keys of
#' units not contributing to a specific numeric variables should be
#' used (`TRUE`) or ignored (`FALSE`) when computing cell keys.
#' @return an object suitable as input to [ck_setup()] for the perturbation
#' of magnitude tables
#' @export
#' @md
#' @examples
#' ck_params_nums(D = 3, l = 0.5, P = "top_contr", top_k = 3)
#' ck_params_nums(D = 10, l = .5, P = "mean")
ck_params_nums <-
  function(P = "top_contr",
           top_k = NULL,
           D,
           l,
           same_key = TRUE,
           use_zero_rkeys = FALSE) {


  stopifnot(is_scalar_integerish(D))
  stopifnot(is_scalar_double(l), l > 0, l < 1)

  if (!is_scalar_character(P)) {
    stop("`P` needs to be a scalar character", call. = FALSE)
  }
  if (!P %in% c("top_contr", "mean", "range", "sum")) {
    stop("invalid value in `P` detected.", call. = FALSE)
  }

  if (!is_scalar_logical(same_key)) {
    stop("`same_key` needs to be a scalar logical", call. = FALSE)
  }
  if (!is_scalar_logical(use_zero_rkeys)) {
    stop("`use_zero_rkeys` needs to be a scalar logical", call. = FALSE)
  }

  if (P == "top_contr") {
    if (is.null(top_k)) {
      stop("please provide a value for `top_k`", call. = FALSE)
    }
    if (!is_integerish(top_k) | top_k < 1 | top_k > 6) {
      stop("`top_k` must be an integer(ish) number >= 1 and <= 6")
    }
  } else {
    if (!is.null(top_k)) {
      message("ignoring argument `top_k`")
    }
    top_k <- 1
  }

  stab <- gen_stab(D = D, l = l)
  out <- list(
    params = list(
      P = P,
      top_k = top_k,
      stab = stab,
      same_key = same_key,
      use_zero_rkeys = use_zero_rkeys
    ),
    type = "nums"
  )
  class(out) <- "ck_params"
  out
}

# example stab from destatis paper
stab_paper1 <- function() {
  dt <- gen_stab(D = 3, l = 0.5)

  # paper nachvollziehen
  dt$p <- c(1,0.18078,0.12789,0.50000,0.06400,0.04528,0.03203,0.02266,0.01603,0.01134,0.00850,0.01719,0.03059,0.04788,0.06594,0.07990,0.50000,0.07990,0.06594,0.04788,0.03059,0.01719,0.00850)
  dt$kum_p_u <- c(0,  0,  0.18078,  0.30867,  0.80867,  0.87267,  0.91795,  0.94997,  0.97263,  0.98866,  0,  0.00850,  0.02570,  0.05629,  0.10417,  0.17010,  0.25000,  0.75000,  0.82990,  0.89583,  0.94371,  0.97430,  0.99150)
  dt$kum_p_o <- c(1,  0.18078,  0.30867,  0.80867,  0.87267,  0.91795,  0.94997,  0.97263,  0.98866,  1,  0.00850,  0.02570,  0.05629,  0.10417,  0.17010,  0.25000,  0.75000,  0.82990,  0.89583,  0.94371,0.97430,0.99150,1)
  dt
}


# example stab for even/odd numbers
stab_paper2 <- function() {
  i <- c(0, rep(1, 9), rep(3, 13))
  j <- c(0, seq(0, 4, by = 0.5), seq(0, 6, by = 0.5))
  p_even <- c(1,0.18078,0.12789,0.50000,0.06400,0.04528,0.03203,0.02266,0.01603,0.01134,0.00850,0.01719,0.03059,0.04788,0.06594,0.07990,0.50000,0.07990,0.06594,0.04788,0.03059,0.01719,0.00850)
  kum_p_u_even <- c(0,0,0.18078,0.30867,0.80867,0.87267,0.91795,0.94997,0.97263,0.98866,0,0.00850,0.02570,0.05629,0.10417,0.17010,0.25000,0.75000,0.82990,0.89583,0.94371,0.97430,0.99150)
  kum_p_o_even <- c(1,0.18078,0.30867,0.80867,0.87267,0.91795,0.94997,0.97263,0.98866,1,0.00850,0.02570,0.05629,0.10417,0.17010,0.25000,0.75000,0.82990,0.89583,0.94371,0.97430,0.99150,1)

  p_odd <- c(1,0.28806,0.22003,0.16309,0.11732,0.08189,0.05548,0.03647,0.02326,0.01440,0.00234,0.00908,0.02756,0.06536,0.12111,0.17536,0.19839,0.17536,0.12111,0.06536,0.02756,0.00908,0.00234)
  kum_p_u_odd <- c(0,0,0.28806,0.50808,0.67118,0.78850,0.87039,0.92587,0.96233,0.98560,0,0.00234,0.01141,0.03897,0.10433,0.22544,0.40080,0.59920,0.77456,0.89567,0.96103,0.98859,0.99766)
  kum_p_o_odd <- c(1,0.28806,0.50808,0.67118,0.78850,0.87039,0.92587,0.96233,0.98560,1,0.00234,0.01141,0.03897,0.10433,0.22544,0.40080,0.59920,0.77456,0.89567,0.96103,0.98859,0.99766,1)

  p <- c(0, seq(-1, 3, by = 0.5), seq(-3, 3, by = 0.5))
  data.table(
    i = i, j = j,
    p_even = p_even, kum_p_u_even = kum_p_u_even, kum_p_o_even = kum_p_o_even,
    p_odd = p_odd, kum_p_u_odd = kum_p_u_odd, kum_p_o_odd = kum_p_o_odd,
    diff = p
  )
}

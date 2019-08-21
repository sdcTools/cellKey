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

#' Set perturbation parameters for continuous variables
#'
#' This function allows to define perturbation parameters used to
#' perturb cells in magnitude tables.
#'
#' @param type a character value defining the way to identify the `magnifier`,
#' e.g which contributions/values in a cell should be the used in the
#' perturbation procedure. Possible choices are:
#' - `top_contr`: the `k` largest contributions are used. In this case,
#' it is also required to specify argument `top_k`
#' - `mean`: the (weighted) cellmean is used as starting point
#' - `range`: the difference between largest and smallest contribution
#' is used.
#' - `sum`: the (weighted) cellvalue itself is used as starting point
#' @param top_k it is ignored if `variant` is different from `top_contr`. Otherwise,
#' @param D maximum perturbation value
#' @param l stepwidth parameter for computation of perturbation tables
#' a integerish number >= 1 specifying the number of top contributions whose
#' values should be perturbed.
#' @param mult_params an object derived by either [ck_gridparams()] or [ck_flexparams()]
#' that specified parameters for the computation of the multiplier.
#' @param mu_c fixed extra protection amount (`>= 0)` applied to the absolute of the
#' perturbation value of the first (largest) noise component; defaults to
#' `0` (no additional protection)
#' @param same_key (logical) should original cell key (`TRUE`) used for
#' for finding perturbation values of the largest contributor to a
#' cell or should a perturbation of the cellkey itself (`FALSE`) take place.
#' @param use_zero_rkeys (logical) skalar defining if record keys of
#' units not contributing to a specific numeric variables should be
#' used (`TRUE`) or ignored (`FALSE`) when computing cell keys.
#' @param m_fixed_sq (numeric scalar); fixed noise variance for very small values
#' @param pos_neg_var a number defining the strategy to look up perturbation values in case
#' the observations can take positive and negative values. This setting is ignored if the variable
#' has no negative values. The possible settings for parameter `pos_neg_var` are:
#' - `0`: the variable must be strictly positive
#' - `1`: the perturbation value is always selected from the block of perturbation
#' values referring to the symmetric case (`i` equals `D` in the perturbation table) independent
#' of the actual value of the observation
#' - `2`: the lookup of a perturbation value `v` for a value `x` is done like in the case
#' where variables only take positive values using `abs(x)` to find the relevant block
#' in the perturbation table. The perturbed value is then computed as `sign(x) * (abs(x) + v)`.
#' - `3`:	use variant `2` for all `x != 0` and apply `1` for `x == 0`
#' @return an object suitable as input to method `$params_nums_set()` for the perturbation
#' of continous variables.
#' @export
#' @seealso [ck_gridparams()], [ck_flexparams()]
#' @md
#' @examples
#' ck_params_nums(
#'   D = 3,
#'   l = 0.5,
#'   type = "top_contr",
#'   top_k = 3,
#'   mult_params = ck_flexparams(
#'     flexpoint = 1000,
#'     m_small = 0.20,
#'     m_large = 0.03,
#'     q = 2
#'   ),
#'   use_zero_rkeys = TRUE,
#'   m_fixed_sq = 2,
#'   mu_c = 3,
#'   pos_neg_var = 2
#' )
#' ck_params_nums(
#'   D = 10,
#'   l = .5,
#'   type = "mean",
#'   mult_params = ck_gridparams(
#'     grid = c(0, 10, 100, 10000),
#'     pcts = c(0.25, 0.20, 0.10, 0.05)
#'   ),
#'   use_zero_rkeys = FALSE,
#'   m_fixed_sq = 4,
#'   mu_c = 5,
#'   pos_neg_var = 3
#' )
ck_params_nums <-
  function(type = "top_contr",
           top_k = NULL,
           D,
           l,
           mult_params,
           mu_c = 0,
           same_key = TRUE,
           use_zero_rkeys = FALSE,
           m_fixed_sq = NULL,
           pos_neg_var = 1) {

  stopifnot(is_scalar_integerish(D))
  stopifnot(is_scalar_double(l), l > 0, l < 1)

  if (!is_scalar_character(type)) {
    stop("`type` needs to be a scalar character", call. = FALSE)
  }
  if (!type %in% c("top_contr", "mean", "range", "sum")) {
    stop("invalid value in `type` detected.", call. = FALSE)
  }
  if (!rlang::is_scalar_double(mu_c)) {
    stop("Argument `mu_c` is not a number.", call. = FALSE)
  }
  if (mu_c < 0) {
    stop("Argument `mu_c` is not >= 0.", call. = FALSE)
  }
  if (!inherits(mult_params, "params_m_grid") & !inherits(mult_params, "params_m_flex")) {
    stop("Argument `mult_params` needs to be created via `ck_gridparams()` or `ck_flexparams()`", call. = FALSE)
  }

  if (!is_scalar_logical(same_key)) {
    stop("`same_key` needs to be a scalar logical", call. = FALSE)
  }
  if (!is_scalar_logical(use_zero_rkeys)) {
    stop("`use_zero_rkeys` needs to be a scalar logical", call. = FALSE)
  }

  if (!is_scalar_integerish(pos_neg_var)) {
    stop("`pos_neg_var` needs to be an integer(ish) number", call. = FALSE)
  }
  if (!pos_neg_var %in% 0:3) {
    stop("Argument `pos_neg_var` needs to be `0`, 1`, `2` or `3`.", call. = FALSE)
  }

  if (type == "top_contr") {
    if (is.null(top_k)) {
      stop("please provide a value for `top_k`", call. = FALSE)
    }
    if (!is_integerish(top_k) | top_k < 1 | top_k > 6) {
      stop("`top_k` must be an integer(ish) number >= 1 and <= 6", call. = FALSE)
    }
  } else {
    if (!is.null(top_k)) {
      message("ignoring argument `top_k`")
    }
    top_k <- 1
  }

  if (!is.null(m_fixed_sq)) {
    if (!rlang::is_scalar_double(m_fixed_sq)) {
      stop("Argument `m_fixed_sq` is not a number.", call = FALSE)
    }
    if (m_fixed_sq <= 0) {
      stop("Argument `m_fixed_sq` must be positive.", call. = FALSE)
    }
  }

  stab <- gen_stab(D = D, l = l)
  out <- list(
    params = list(
      type = type,
      top_k = top_k,
      stab = stab,
      mu_c = mu_c,
      m_fixed_sq = m_fixed_sq,
      mult_params = mult_params,
      same_key = same_key,
      use_zero_rkeys = use_zero_rkeys,
      pos_neg_var = pos_neg_var
    ),
    type = "nums"
  )
  class(out) <- "ck_params"
  out
}

#' Define parameters for numeric perturbation magnitudes using a fixed grid
#'
#' [ck_gridparams()] allows to define a grid that is used to lookup perturbation
#' magnitudes (percentages) used when perturbing continuous variables.
#'
#' @details details about the grid function can be found in Deliverable D4.2, Part I in
#' SGA *"Open Source tools for perturbative confidentiality methods"*
#'
#' @inheritParams ck_flexparams
#' @param grid a numeric vector (ascending order) defining the bounds for which
#' a specific percentage (defined in `pcts`) needs to be applied
#' @param pcts a numeric vector defining percentages (between `0` and `1`)
#'
#' @return an object suitable as input for [ck_params_nums()].
#' @export
#' @inherit cellkey_pkg examples
#' @seealso [ck_params_nums()], [ck_flexparams()]
#' @md
ck_gridparams <- function(grid, pcts) {
  if (!is.numeric(grid)) {
    stop("Argument `grid` is not numeric.", call = FALSE)
  }
  if (!is.numeric(pcts)) {
    stop("Argument `pcts` is not numeric.", call = FALSE)
  }
  if (length(pcts) != length(grid)) {
    stop("Arguments `pcts` and `grid` differ in length.", call = FALSE)
  }

  if (any(grid != sort(grid))) {
    stop("Argument `grid` is not in ascending order.", call. = FALSE)
  }
  if (any(pcts != sort(pcts, decreasing = TRUE))) {
    stop("Argument `pcts` is not in decreasing order.", call. = FALSE)
  }

  if (min(pcts) <= 0) {
    stop("Argument `pcts` must contain values > 0", call. = FALSE)
  }

  if (min(pcts) <= 0) {
    stop("Argument `pcts` must contain values > 0", call. = FALSE)
  }
  if (max(pcts) > 1) {
    stop("Argument `pcts` must contain values <= 1", call. = FALSE)
  }
  out <- list(
    grid = grid,
    pcts = pcts
  )
  class(out) <- "params_m_grid"
  out
}


#' Define parameters for numeric perturbation magnitudes using a flex function
#'
#' [ck_gridparams()] allows to define a flex function that is used to lookup perturbation
#' magnitudes (percentages) used when perturbing continuous variables.
#'
#' @details details about the grid function can be found in Deliverable D4.2, Part I in
#' SGA *"Open Source tools for perturbative confidentiality methods"*
#' @param flexpoint (numeric scalar); at which point should the noise coefficient
#' function reaches its desired maximum (defined by `m_small`)
#' @param m_small (numeric scalar); the desired maximum percentage for the function (for small values)
#' @param m_large (numeric scalar); the desired noise percentage for larger values
#' @param q (numeric scalar); Parameter of the function; `q` needs to be `>= 1`
#'
#' @return an object suitable as input for [ck_params_nums()].
#' @export
#' @inherit cellkey_pkg examples
#' @seealso [ck_params_nums()], [ck_gridparams()]
#' @md
ck_flexparams <- function(flexpoint, m_small = 0.25, m_large = 0.05, q = 3) {
  if (!rlang::is_scalar_double(flexpoint)) {
    stop("Argument `flexpoint` is not a number.", call = FALSE)
  }
  if (flexpoint <= 0) {
    stop("Argument `flexpoint` must be positive.", call. = FALSE)
  }
  if (!rlang::is_scalar_double(m_small)) {
    stop("Argument `m_small` is not a number.", call = FALSE)
  }
  if (m_small <= 0 | m_small > 1) {
    stop("Argument `m_small` must be > 0 and <= 1.", call. = FALSE)
  }
  if (!rlang::is_scalar_double(m_large)) {
    stop("Argument `m_large` is not a number.", call = FALSE)
  }
  if (m_large <= 0 | m_large > 1) {
    stop("Argument `m_large` must be > 0 and <= 1.", call. = FALSE)
  }
  if (!rlang::is_scalar_double(q)) {
    stop("Argument `q` is not a number.", call = FALSE)
  }
  if (q < 1) {
    stop("Argument `q` needs to be >= 1", call. = FALSE)
  }
  out <- list(
    flexpoint = flexpoint,
    m_small = m_small,
    m_large = m_large,
    q = q
  )
  class(out) <- "params_m_flex"
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

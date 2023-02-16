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
#' @param top_k it is ignored if `variant` is different from `top_contr`.
#' Otherwise, a scalar number `>0` is expected.
#' @param ptab in this argument, one ore more perturbation tables are given as
#' input. the following choices are possible:
#'
#' - an object derived from [ptable::create_ptable()]
#' or [ptable::create_num_ptable()]: this case is the same as specifying a named
#' list with only a single element `"all"` (as described below)
#' - a named `list` where the allowed names are shown below and each element
#' must be the output of [ptable::create_ptable())] or
#' [ptable::create_num_ptable()]
#'   * `"all"`: this ptable will be used for all cells; if specified, no
#'   elements named `"even"` or `"odd"` must exist.
#'   * `"even"`: will be used to look up perturbation values for cells with an
#'   even number of contributors. if specified, also list-element `"odd"` must
#'   exist.
#'   * `"odd"`: will be used to look up perturbation values for cells with an
#'   odd number of contributors; if specified, also list-element `"even"` must
#'   exist.
#'   * `"small_cells"`: if specified, this ptable will be used to extract
#'   perturbation values for very small cells
#' @param mult_params an object derived with [ck_flexparams()] or
#' [ck_simpleparams()] that contain required parameters for the computation of
#' the perturbation multiplier
#' @param mu_c fixed extra protection amount (`>= 0)` applied to the absolute of
#' the perturbation value of the first (largest) noise component if the cell is
#' sensitive. This value defaults to `0` (no additional protection). Please note
#' that sensitive cells can be defined according using the `supp_freq()`,
#' `supp_val`, `supp_p()`, `supp_nk()` and `supp_pq()` methods. An examples is
#' given in `?cellkey_pkg`.
#' @param same_key (logical) should original cell key (`TRUE`) used for
#' for finding perturbation values of the largest contributor to a
#' cell or should a perturbation of the cellkey itself (`FALSE`) take place.
#' @param use_zero_rkeys (logical) scalar defining if record keys of
#' units not contributing to a specific numeric variables should be
#' used (`TRUE`) or ignored (`FALSE`) when computing cell keys.
#' @param path a scalar character specifying a path to which the parameters
#' created with this functions should be written to (in yaml format)
#' @return an object suitable as input to method `$params_nums_set()` for the
#' perturbation of continous variables.
#' @export
#' @seealso [ck_flexparams()]
#' @md
#' @examples
#' # create a perturbation table using
#' # functionality from ptable-pkg; see help(pa = "ptable")
#' # this returns an extra ptable for very small cells
#' ptab <- ptable::pt_ex_nums(separation = TRUE)
#'
#' # create parameters
#' ck_params_nums(
#'   type = "top_contr",
#'   top_k = 3,
#'   ptab = ptab,
#'   mult_params = ck_flexparams(
#'     fp = 1000,
#'     p = c(0.20, 0.03),
#'     epsilon = c(1, 0.5, 0.2),
#'     q = 2),
#'   use_zero_rkeys = TRUE,
#'   mu_c = 3)
ck_params_nums <-
  function(type = "top_contr",
           top_k = NULL,
           ptab,
           mult_params,
           mu_c = 0,
           same_key = TRUE,
           use_zero_rkeys = FALSE,
           path = NULL) {

  if (!is_scalar_character(type)) {
    stop("`type` needs to be a scalar character", call. = FALSE)
  }
  if (!type %in% c("top_contr", "mean", "range", "sum")) {
    stop("invalid value in `type` detected.", call. = FALSE)
  }

  if (inherits(ptab, "ptable_params")) {
    if (ptab@table != "nums") {
      e <- "input from `ptable::create_ptable` or `ptable::create_cnt_ptable` not suitable for continuous variables."
      stop(e, call. = FALSE)
    }
    ptab <- ptable::create_ptable(params = ptab)
  }

  ptab <- .chk_ptab(ptab, type = "nums")
  parity <- even_odd <- ifelse("even" %in% ptab$type, TRUE, FALSE)
  separation <- ifelse("small_cells" %in% ptab$type, TRUE, FALSE)

  if (!is.numeric(mu_c) | !rlang::is_scalar_atomic(mu_c)) {
    stop("Argument `mu_c` is not a number.", call = FALSE)
  }

  if (mu_c < 0) {
    stop("Argument `mu_c` is not >= 0.", call. = FALSE)
  }

  if (!is_scalar_logical(same_key)) {
    stop("`same_key` needs to be a scalar logical", call. = FALSE)
  }
  if (!is_scalar_logical(use_zero_rkeys)) {
    stop("`use_zero_rkeys` needs to be a scalar logical", call. = FALSE)
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

  if (parity & top_k > 1) {
    e <- "different ptables for even/odd cases are not possible if `top_k` is > 1."
    stop(e, call. = FALSE)
  }

  if (!class(mult_params) %in% c("params_m_flex", "params_m_simple")) {
    e <- "`mult_params` needs to be creating with ck_flexparams() or ck_simpleparams()"
    stop(e, call. = FALSE)
  }

  if (length(mult_params$epsilon) != top_k) {
    stop("Invalid length or argument `epsilon` in `mult_params` detected.", call. = FALSE)
  }

  # separation: treat very small values in magnitude tables
  # should be treated like counts
  # if we have a value for `m_fixed_sq`; we have a suitable ptable for small cells

  # compute parameter m1^2 from ptable
  # used in `ck_params_nums()` in case separation is TRUE
  .compute_m1sq <- function(ptab) {
    type <- i <- NULL
    dt <- ptab[type == "small_cells"]
    dt <- dt[i == max(i)]
    if (nrow(dt) == 0) {
      return(stop("invalid ptab found when computing parameter m1sq.", call. = FALSE))
    }
    return(sum((dt$i - dt$j) ^ 2 * dt$p)) # nolint
  }

  m_fixed_sq <- NA
  if (separation) {
    m_fixed_sq <- .compute_m1sq(ptab)
  } else {
    ptab <- ptab[type != "small_cells"]
  }

  # separation point z_s (previously g1)
  .separation_point <- function(m_fixed_sq, e, p) {
    if (is.na(m_fixed_sq)) {
      return(0)
    }
    sqrt(m_fixed_sq) / (sqrt(e) * p)
  }

  e <- sum(mult_params$epsilon ^ 2)
  if (inherits(mult_params, "params_m_flex")) {
    p <- mult_params$p_small
  } else if (inherits(mult_params, "params_m_simple")) {
    p <- mult_params$p
  } else {
    stop("invalid input.", call. = FALSE)
  }

  zs <- .separation_point(
    m_fixed_sq = m_fixed_sq,
    e = e,
    p = p)

  out <- list(
    params = list(
      type = type,
      top_k = top_k,
      ptab = ptab,
      mu_c = mu_c,
      m_fixed_sq = m_fixed_sq,
      zs = zs,
      E = e, # nolint
      mult_params = mult_params,
      same_key = same_key,
      use_zero_rkeys = use_zero_rkeys,
      even_odd = even_odd,
      separation = separation),
    type = class(mult_params))
  class(out) <- "ck_params"

  if (!is.null(path)) {
    out$version <- paste(utils::packageVersion("cellKey"), collapse = ".")
    out$ptype <- "params_nums"
    .yaml_write(x = out, path = path)
    message("yaml configuration ", shQuote(path), " successfully written.")
  }
  out$version <- out$ptype <- NULL
  out
}

#' Set parameters required to perturb numeric variables using a flex function
#'
#' [ck_flexparams()] allows to define a flex function that is used to lookup perturbation
#' magnitudes (percentages) used when perturbing continuous variables.
#'
#' @details details about the flex function can be found in Deliverable D4.2, Part I in
#' SGA *"Open Source tools for perturbative confidentiality methods"*
#' @param fp (numeric scalar); at which point should the noise coefficient
#' function reaches its desired maximum (defined by the first element of `p`)
#' @param p a numeric vector of length `2` where both elements specify a percentage.
#' The first value refers to the desired maximum perturbation percentage for small
#' cells (depending on `fp`) while the second element refers to the desired maximum
#' perturbation percentage for large cells. Both values must be between `0` and `1` and
#' need to be in descending order.
#' @param q (numeric scalar); Parameter of the function; `q` needs to be `>= 1`
#' @param epsilon a numeric vector in descending order with all values `>= 0` and `<= 1` with the first
#' element forced to equal 1. The length of this vector must correspond with the number `top_k`
#' specified in [ck_params_nums()] when creating parameters for `type == "top_contr"` which is
#' checked at runtime. This setting allows to use different flex-functions for the largest `top_k` contributors.
#' @return an object suitable as input for [ck_params_nums()].
#' @export
#' @inherit cellkey_pkg examples
#' @seealso [ck_simpleparams()], [ck_params_nums()]
#' @md
ck_flexparams <- function(fp, p = c(0.25, 0.05), epsilon = 1, q = 3) {
  if (!is.numeric(fp) | !rlang::is_scalar_atomic(fp)) {
    stop("Argument `fp` is not a number.", call. = FALSE)
  }
  if (fp <= 0) {
    stop("Argument `fp` must be positive.", call. = FALSE)
  }
  if (!is.numeric(p)) {
    stop("`p` is not a numeric vector.", call. = FALSE)
  }
  if (length(p) != 2) {
    stop("`p` does not have 2 elements.", call. = FALSE)
  }
  if (any(diff(p) > 0)) {
    stop("Argument `p` must contain numbers in descending order", call. = FALSE)
  }
  if (max(p) > 1 | min(p) < 0) {
    stop("values `< 0` or `> 1` are not allowed in argument `p`", call. = FALSE)
  }

  if (!(rlang::is_scalar_double(q) | rlang::is_scalar_integer(q))) {
    stop("Argument `q` is not a number", call. = FALSE)
  }

  if (q < 1) {
    stop("Argument `q` needs to be >= 1", call. = FALSE)
  }

  if (!is.numeric(epsilon)) {
    stop("Argument `epsilon` must be a numeric vector", call. = FALSE)
  }
  if (any(epsilon < 0)) {
    stop("Argument `epsilon` must contain only numbers >= 0", call. = FALSE)
  }
  if (any(epsilon > 1)) {
    stop("Argument `epsilon` must contain only numbers <= 1", call. = FALSE)
  }
  if (epsilon[1] != 1) {
    stop("The first element of `epsilon` must equal 1", call. = FALSE)
  }
  if (length(epsilon) > 1) {
    if (any(diff(epsilon) > 0)) {
      stop("Argument `epsilon` must contain numbers in descending order", call. = FALSE)
    }
  }
  out <- list(
    fp = fp,
    p_small = p[1],
    p_large = p[2],
    epsilon = epsilon,
    q = q)
  class(out) <- "params_m_flex"
  out
}

#' Set parameters required to perturb numeric variables using a simple approach
#'
#' [ck_simpleparams()] allows to define parameters for a simple perturbation
#' approach based on a single magnitude parameter (`m`). The values of `epsilon`
#' are used to  `"weight"` parameter `m` in case `type == "top_contr"` is set in
#' [ck_params_nums()].
#'
#' @inherit ck_flexparams details
#' @inherit ck_flexparams return
#' @inherit cellkey_pkg examples
#' @inheritParams ck_flexparams
#' @param p a percentage value used as magnitude for perturbation
#' @export
#' @md
#' @seealso [ck_flexparams()], [ck_params_nums()]
ck_simpleparams <- function(p, epsilon = 1) {
  # basically scaling = FALSE
  if (!rlang::is_scalar_double(p)) {
    stop("Argument `p` is not a number.", call = FALSE)
  }
  if (p <= 0 | p > 1) {
    stop("Argument `p` must be > 0 and <= 1.", call. = FALSE)
  }
  if (!is.numeric(epsilon)) {
    stop("Argument `epsilon` must be a numeric vector", call. = FALSE)
  }
  if (any(epsilon < 0)) {
    stop("Argument `epsilon` must contain only numbers >= 0", call. = FALSE)
  }
  if (any(epsilon > 1)) {
    stop("Argument `epsilon` must contain only numbers <= 1", call. = FALSE)
  }
  if (epsilon[1] != 1) {
    stop("The first element of `epsilon` must equal 1", call. = FALSE)
  }
  if (length(epsilon) > 1) {
    if (any(diff(epsilon) > 0)) {
      stop("Argument `epsilon` must contain numbers in descending order", call. = FALSE)
    }
  }
  out <- list(p = p, epsilon = epsilon)
  class(out) <- "params_m_simple"
  out
}

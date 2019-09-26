gen_stab <- function(D = 3, l = 0.5, add_small_cells = TRUE, even_odd = TRUE) {
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

  if (!even_odd) {
    dt <- data.table(i = 0, j = 0, p = 1, kum_p_u = 0, kum_p_o = 1, diff = 0)
    # i:=1
    dt_a <- .gen(i = 1, min_d = -1, max_d = D, l = l)
    # i == D
    dt_b <- .gen(i = D, min_d = -D, max_d = D, l = l)
    dt <- rbind(dt, dt_a, dt_b)
    dt$type <- "all"
  } else {
    dt <- data.table(i = 0, j = 0, p = 1, kum_p_u = 0, kum_p_o = 1, diff = 0)
    # i:=1
    dt_a <- .gen(i = 1, min_d = -1, max_d = D, l = l)
    # i == D
    dt_b <- .gen(i = D, min_d = -D, max_d = D, l = l)
    dt_even <- rbind(dt, dt_a, dt_b)
    dt_even$type <- "even"

    dt <- data.table(i = 0, j = 0, p = 1, kum_p_u = 0, kum_p_o = 1, diff = 0)
    # i:=1
    dt_a <- .gen(i = 1, min_d = -1, max_d = D, l = l)
    # i == D
    dt_b <- .gen(i = D, min_d = -D, max_d = D, l = l)
    dt_odd <- rbind(dt, dt_a, dt_b)
    dt_odd$type <- "odd"
    dt <- rbind(dt_even, dt_odd)
  }

  if (add_small_cells) {
    dt_small <- data.table(i = 0, j = 0, p = 1, kum_p_u = 0, kum_p_o = 1, diff = 0)
    # i:=1
    dt_a <- .gen(i = 1, min_d = -1, max_d = D, l = l)
    # i == D
    dt_b <- .gen(i = D, min_d = -D, max_d = D, l = l)
    dt_small <- rbind(dt_small, dt_a, dt_b)
    dt_small$type <- "small_cells"
    dt <- rbind(dt, dt_small)
  }
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
#' @param mult_params an object derived with [ck_flexparams()] or [ck_simpleparams()]
#' that contain required parameters for the computation of the perturbation multiplier
#' @param mu_c fixed extra protection amount (`>= 0)` applied to the absolute of the
#' perturbation value of the first (largest) noise component if the cell is sensitive.
#' This value defaults to `0` (no additional protection). Please note that sensitive cells
#' can be defined according using the `supp_freq()`, `supp_val`, `supp_p()`, `supp_nk()`
#' and `supp_pq()` methods. An examples is given in `?cellkey_pkg`.
#' @param same_key (logical) should original cell key (`TRUE`) used for
#' for finding perturbation values of the largest contributor to a
#' cell or should a perturbation of the cellkey itself (`FALSE`) take place.
#' @param use_zero_rkeys (logical) scalar defining if record keys of
#' units not contributing to a specific numeric variables should be
#' used (`TRUE`) or ignored (`FALSE`) when computing cell keys.
#' @param separation (logical) scalar defining if very small cell values should
#' be perturbed like counts (exact lookup). If `TRUE`, the ptable provided must contain
#' values with value `small_cells` in column `type`.
#' @param parity (logical) scalar defining if different perturbation tables should be used
#' for cells with an even or odd number of contributors.
#' @param path a scalar character specifying a path to which the parameters created with this functions
#' should be written to (in yaml format)
#' @return an object suitable as input to method `$params_nums_set()` for the perturbation
#' of continous variables.
#' @export
#' @seealso [ck_flexparams()]
#' @md
#' @examples
#' ck_params_nums(
#'   D = 3,
#'   l = 0.5,
#'   type = "top_contr",
#'   top_k = 3,
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
           D,
           l,
           mult_params,
           mu_c = 0,
           same_key = TRUE,
           use_zero_rkeys = FALSE,
           separation = FALSE,
           parity = FALSE,
           path = NULL) {

  stopifnot(is_scalar_integerish(D))
  stopifnot(is_scalar_double(l), l > 0, l < 1)

  if (!is_scalar_character(type)) {
    stop("`type` needs to be a scalar character", call. = FALSE)
  }
  if (!type %in% c("top_contr", "mean", "range", "sum")) {
    stop("invalid value in `type` detected.", call. = FALSE)
  }

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

  if (!is_scalar_logical(separation)) {
    stop("`separation` needs to be a scalar logical", call. = FALSE)
  }

  if (!is_scalar_logical(parity)) {
    stop("`parity` needs to be a scalar logical", call. = FALSE)
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
    stop("`parity` can only be `TRUE` when `top_k` is > 1.", call. = FALSE)
  }

  if (!class(mult_params) %in% c("params_m_flex", "params_m_simple")) {
    e <- "`mult_params` needs to be creating with ck_flexparams() or ck_simpleparams()"
    stop(e, call. = FALSE)
  }

  if (length(mult_params$epsilon) != top_k) {
    stop("Invalid length or argument `epsilon` in `mult_params` detected.", call. = FALSE)
  }

  # even_odd: different lookups for even/odd unweighted cell values
  # ~ parameter `parity` in the documents and tau-argus
  even_odd <- parity

  # needs to come from ptable-package!
  stab <- gen_stab(
    D = D,
    l = l,
    even_odd = even_odd,
    add_small_cells = separation)

  if (even_odd) {
    if (nrow(stab[type == "even"])== 0 | nrow(stab[type == "odd"]) == 0) {
      e <- "argument `parity` is TRUE but 0 cells with `type` being 'even' or `odd` found in ptable."
      stop(e, call. = FALSE)
    }
  } else {
    if (nrow(stab[type == "all"]) == 0) {
      e <- "argument `parity` is FALSE but 0 cells with `type == 'all'` found in ptable."
      stop(e, call. = FALSE)
    }
  }

  # separation: treat very small values in magnitude tables
  # should be treated like counts
  # if we have a value for `m_fixed_sq`; we have a suitable ptable for small cells

  # compute parameter m1^2 from ptable
  # used in `ck_params_nums()` in case separation is TRUE
  .compute_m1sq <- function(ptab) {
    type <- i <- NULL
    dt <- ptab[type == "small_cells" & i == max(i)]
    if (nrow(dt) == 0) {
      return(stop("invalid ptab found when computing parameter m1sq.", call. = FALSE))
    }
    sum((dt$i - dt$j)^2 * dt$p)
  }

  m_fixed_sq <- NA
  if (separation) {
    m_fixed_sq <- .compute_m1sq(stab)
  } else {
    stab <- stab[type != "small_cells"]
  }


  # separation point z_s (previously g1)
  .separation_point <- function(m_fixed_sq, E, p_large) {
    if (is.na(m_fixed_sq)) {
      return(0)
    }
    sqrt(m_fixed_sq) / (sqrt(E) * p_large)
  }

  E <- sum(mult_params$epsilon^2)
  if (inherits(mult_params, "params_m_flex")) {
    p_large <- mult_params$p_large
  } else if (inherits(mult_params, "params_m_simple")) {
    p_large <- mult_params$p
  #} else if (inherits(mult_params, "params_m_grid")) {
    #E <- top_k
    #p_large <- mult_params[[1]]$pcts[1]
  } else {
    stop("invalid input.", call. = FALSE)
  }

  zs <- .separation_point(
    m_fixed_sq = m_fixed_sq,
    E = E,
    p_large = p_large)

  out <- list(
    params = list(
      type = type,
      top_k = top_k,
      stab = stab,
      mu_c = mu_c,
      m_fixed_sq = m_fixed_sq,
      zs = zs,
      E = E,
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
  out <- list(
    p = p,
    epsilon = epsilon)
  class(out) <- "params_m_simple"
  out
}

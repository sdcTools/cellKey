# perturb cell-keys by removing the first digit and appending it at the back
.scramble_ckeys <- function(ck, top_k, same_key) {
  beg <- substr(ck, 1, 2)
  res <- rep(NA, top_k)

  len <- nchar(ck)

  if (!same_key) {
    res[1] <- paste0(beg, substr(ck, 4, len), substr(ck, 3, 3))
  } else {
    res[1] <- ck
  }

  if (top_k == 1) {
    return(as.numeric(res))
  }

  for (i in 2:top_k) {
    res[i] <- paste0(beg, substr(res[i - 1], 4, len), substr(res[i - 1], 3, 3))
  }
  as.numeric(res)
}

# returns the value of the flex-function defined by
# x: value to be checked
# fp: flexpoint
# p_lg: percentage of perturbation for large cells
# p_sm: percentage of perturbation for small cells
# q: parameter for flex function
.para_m <- function(x, p_sm, p_lg, fp, q) {
  f1 <- ((p_sm * x) - (p_lg * fp)) / (p_lg * fp)
  f2 <- ((2 * fp) / (fp + x)) ^ q
  p_lg * (1 + (f1 * f2))
}

.lookup_v_flex <- function(cellkeys, params) {
  i <- type <- ub <- NULL
  d <- params$max_i # parameter `D`

  v <- rep(NA, length(cellkeys))
  a <- params$a

  # perturbation value should be 0 for all cells below the separation point
  # (these are the ones with lookup == "small_cells") and j > 1
  ind_smallcells <- setdiff(which(params$lookup == "small_cells"), 1)
  if (length(ind_smallcells) > 0) {
    params$lookup[ind_smallcells] <- "_zero_"
  }

  a[a > d] <- d # we cut at the maximum value!

  ptab <- params$ptab
  poss <- sort(unique(ptab$i))

  ind_exact <- which(a %in% poss)
  if (length(ind_exact) > 0) {
    v[ind_exact] <- sapply(ind_exact, function(x) {
      if (params$lookup[x] == "_zero_") {
        return(0)
      }
      ptab[type == params$lookup[x] & i == a[x] & cellkeys[x] < ub, v][1]
    })
  }

  ind_comb <- which(a < d)
  if (length(ind_comb) > 0) {
    v[ind_comb] <- sapply(ind_comb, function(x) {
      if (params$lookup[x] == "_zero_") {
        return(0)
      }
      a0 <- poss[which(a[x] < poss) - 1]
      a1 <- poss[max(which(a[x] > poss)) + 1]
      lambda <- (a[x] - a0) / (a1 - a0)

      # compute two perturbation values
      v_low <- ptab[type == params$lookup[x] & i == a0 & cellkeys[x] < ub, v][1]
      v_up <- ptab[type == params$lookup[x] & i == a1 & cellkeys[x] < ub, v][1]

      # combine to get final perturbation value
      (1 - lambda) * v_low + lambda * v_up
    })
  }
  v
}

# perturb a given cell with the following parameters
# cv: original cell value
# x: top_k largest contributors to cell
# ck: corresponding cell keys
# prot_req: additional protection required? if TRUE, we use mu_c, else 0
# lookup: where to look in the perturbation table
# params: a perturbation parameter object created with `ck_params_nums` and
#         mult_params need to be of class `params_m_flex`
.perturb_cell_flex <- function(cv, x, ck, lookup, prot_req, params) {
  type <- i <- NULL
  dig <- .ck_digits()

  debug <- FALSE
  fp <- params$mult_params$fp
  p_lg <- params$mult_params$p_large
  p_sm <- params$mult_params$p_small
  q <- params$mult_params$q
  mu <- ifelse(prot_req, params$mu_c, 0)
  ptab <- params$ptab
  zs <- params$zs
  deltas <- rep(NA, length(x))

  # lookup
  lookup_params <- list(
    ptab = ptab,
    max_i = max(ptab$i))

  xo <- cv
  sign_xo <- sign(xo)
  x_hats <- rep(NA, length(x))
  for (j in seq_len(length(x))) {
    zero_pert <- FALSE
    lookup_params$lookup <- lookup[j]
    if (debug) {
      message("j:= ", j, " | ck: ", ck[j])
    }

    xj <- x[j]
    abs_xj <- abs(xj)

    if (abs_xj < zs) {
      x_delta <- 1
      lookup_params$lookup <- "small_cells"
    } else {
      if (abs_xj < fp) {
        x_delta <- xj * p_sm
      } else {
        m <- .para_m(x = abs_xj, p_sm = p_sm, p_lg = p_lg, fp = fp, q = q)
        x_delta <-  xj * m * params$mult_params$epsilon[j]
      }
    }
    abs_xdelta <- abs(x_delta)
    if (abs(xo) < abs_xdelta) {
      x_delta <- xo
      xj <- xo / (p_sm * params$mult_params$epsilon[j])
      abs_xj <- abs(xj)
      if (debug) {
        message("abs_x < abs_xdelta --> rescaling")
        message("x_delta (new): ", x_delta)
        message("x (new): ", x)
      }

      if (abs(abs_xj) < zs) {
        x_delta <- 1
        lookup_params$lookup <- "small_cells"
        if (j > 1) {
          zero_pert <- TRUE
        }
      }
    }

    deltas[j] <- x_delta

    # compute perturbation value
    if (!zero_pert) {
      lookup_params$a <- abs(xo / x_delta)
      lookup_params$max_i <- lookup_params$ptab[type == lookup_params$lookup, max(i)]

      v <- .lookup_v_flex(cellkeys = ck[j], params = lookup_params)
      if (debug) {
        message("xj: ", round(xj, digits = dig), appendLF = FALSE)
        message(" | x_delta: ", round(x_delta, digits = dig), appendLF = FALSE)
        message(" | a: ", round(lookup_params$a, digits = dig), appendLF = FALSE)
        message(" | v: ", v)
      }
      sign_v <- sign(v)

      x_hats[j] <- sign_xo * abs(x_delta) * (sign_v * (mu + abs(v)))
      if (sign_v == -1) {
        if (sign_xo == -1) {
          x_hats[j] <- min(x_hats[j], -xo)
        } else if (sign_xo > sign_v) {
          x_hats[j] <- max(x_hats[j], -xo)
        }
      }
    } else {
      x_hats[j] <- 0
    }

    if (debug) {
      message("x_hat: ", x_hats[j])
      message("updating original value: old: ", xo, "", appendLF = FALSE)
    }
    xo <- xo + x_hats[j]
    sign_xo <- sign(xo)
    if (debug) {
      message(" | new: ", xo)
    }
  }
  list(
    x_hats = x_hats,
    cv = cv,
    cv_p = cv + sum(x_hats),
    lookup = lookup_params$lookup,
    x_delta = deltas)
}

# perturb a given cell with the following parameters
# cv: original cell value
# x: top_k largest contributors to cell
# ck: corresponding cell keys
# prot_req: additional protection required? if TRUE, we use mu_c, else 0
# lookup: where to look in the perturbation table
# params: a perturbation parameter object created with `ck_params_nums` and
#         mult_params need to be of class `params_m_simple`
.perturb_cell_simple <- function(cv, x, ck, lookup, prot_req, params) {
  type <- i <- NULL
  debug <- TRUE
  dig <- .ck_digits()


  p <- params$mult_params$p # default percentage
  if (debug) {
    message("default_percentage: ", p)
  }

  # separation_point
  zs <- params$zs
  if (debug) {
    message("separation point: ", zs)
  }

  mu <- ifelse(prot_req, params$mu_c, 0)

  ptab <- params$ptab
  deltas <- rep(NA, length(x))

  # lookup
  lookup_params <- list(ptab = ptab)

  xo <- cv
  sign_xo <- sign(xo)
  x_hats <- rep(NA, length(x))
  cnt_sc <- 0
  for (j in seq_len(length(x))) {
    zero_pert <- FALSE

    lookup_params$lookup <- lookup[j]
    if (debug) {
      message("j:= ", j, " | ck: ", ck[j], " | lookup: ", shQuote(lookup[j]))
    }

    xj <- x[j]
    abs_xj <- abs(xj)

    if (abs_xj < zs) {
      x_delta <- xj # in this case, m equals 1
      lookup_params$lookup <- "small_cells"
      cnt_sc <- cnt_sc + 1
    } else {
      x_delta <- xj * p * params$mult_params$epsilon[j]
    }

    abs_xdelta <- abs(x_delta)
    if (abs(xo) < abs_xdelta) {
      x_delta <- xo
      xj <- xo / (p * params$mult_params$epsilon[j])
      abs_xj <- abs(xj)
      if (debug) {
        message("abs_x < abs_xdelta --> rescaling")
        message("x_delta (new): ", x_delta)
        message("x (new): ", x)
      }

      if (abs(abs_xj) < zs) {
        x_delta <- 1
        lookup_params$lookup <- "small_cells"
        cnt_sc <- cnt_sc + 1
        if (cnt_sc > 1) {
          zero_pert <- TRUE
        }
      }
    }

    deltas[j] <- x_delta

    # compute perturbation value
    if (!zero_pert) {
      lookup_params$a <- abs(xo / x_delta)
      lookup_params$max_i <- lookup_params$ptab[type == lookup_params$lookup, max(i)]
      v <- .lookup_v_flex(cellkeys = ck[j], params = lookup_params)
      if (debug) {
        message("xj: ", round(xj, digits = dig), appendLF = FALSE)
        message(" | x_delta: ", round(x_delta, digits = dig), appendLF = FALSE)
        message(" | a: ", round(lookup_params$a, digits = dig), appendLF = FALSE)
        message(" | v: ", v)
      }
      sign_v <- sign(v)

      x_hats[j] <- sign_xo * abs(x_delta) * (sign_v * (mu + abs(v)))
      if (sign_v == -1) {
        if (sign_xo == -1) {
          x_hats[j] <- min(x_hats[j], -xo)
        } else if (sign_xo > sign_v) {
          x_hats[j] <- max(x_hats[j], -xo)
        }
      }
    } else {
      x_hats[j] <- 0
    }

    if (debug) {
      message("x_hat: ", x_hats[j])
      message("updating original value: old: ", xo, "", appendLF = FALSE)
    }
    xo <- xo + x_hats[j]
    sign_xo <- sign(xo)
    if (debug) {
      message(" | new: ", xo)
    }
  }
  list(
    x_hats = x_hats,
    cv = cv,
    cv_p = cv + sum(x_hats),
    lookup = lookup_params$lookup,
    x_delta = deltas)
}
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
    res[i] <- paste0(beg, substr(res[i - 1], 4, len), substr(res[i- 1], 3, 3))
  }
  as.numeric(res)
}

# compute x_delta  multiplication parameters m_j (x) * x
.get_x_delta <- function(params, x, top_k, m_fixed_sq) {
  UseMethod(".get_x_delta", params)
}
.get_x_delta.default <- function(params, x, top_k, m_fixed_sq) {
  stop("invalid input in `.get_x_delta()` detected!", call. = FALSE)
}
.get_x_delta.params_m_flex <- function(params, x, top_k, m_fixed_sq) {
  # by default: perturbation for small cells -> highest pctg
  m <- rep(params$m_small, length(x))
  fp <- params$flexpoint

  ind_lg <- which(x > fp)
  if (length(ind_lg) > 0) {
    x_lg <- x[ind_lg]

    f1 <- ((params$m_small * x_lg) - (params$m_large * fp)) / (params$m_large * fp)
    f2 <- (2 * fp) / (fp + x_lg)
    m[ind_lg] <- params$m_large * (1 + (f1 * f2 ^ params$q))
  }

  # fixed variance for very small observations
  if (!is.null(m_fixed_sq)) {
    # compute g1 (formula 2.3 on p.9)
    g1 <- sqrt(m_fixed_sq) / (sqrt(top_k) * params$m_large)
    # very small values
    ind_vs <- which(x < g1)
    if (length(ind_vs) > 0) {
      m[ind_vs] <- sqrt(m_fixed_sq)
      x[ind_vs] <- 1
    }
  }
  x * m
}
.get_x_delta.params_m_grid <- function(params, x, top_k, m_fixed_sq) {
  m <- rep(NA, length(x))

  # footnote 3, page 7 -> if x < 0 -> use abs(x)
  rr <- lapply(abs(x), function(y) y <= params$grid)

  # smallest cells, highest perturbation
  ind_sm <- sapply(rr, function(x) all(x == TRUE))
  m[ind_sm] <- params$pcts[1]

  # largest cells, smallest perturbation
  ind_lg <- sapply(rr, function(x) all(x == FALSE))
  m[ind_lg] <- utils::tail(params$pcts, 1)

  ind_mid <- is.na(m)
  m[ind_mid] <- sapply(rr[ind_mid], function(x) {
    params$pcts[min(which(x == TRUE))]
  })

  # fixed variance for very small observations
  if (!is.null(m_fixed_sq)) {
    # compute g1 (formula 2.3 on p.9)
    g1 <- sqrt(m_fixed_sq) / (sqrt(top_k) * params$pcts[1])
    # very small (absolute) values
    ind_vs <- which(abs(x) < g1)
    if (length(ind_vs) > 0) {
      m[ind_vs] <- sqrt(m_fixed_sq)
      x[ind_vs] <- 1
    }
  }
  x * m
}

# .get_x_delta(params = ck_flexparams(flexpoint = 100), x = c(100, 20), top_k = 3, m_fixed_sq = 2)
# .get_x_delta(params = ck_flexparams(flexpoint = 100), x = c(-100, 20), top_k = 3, m_fixed_sq = 2)

# returning perturbation values from `stab` based on cellkeys, x_delta (x * m), the weighted cell-value and
# pos_neg_var (defined in `ck_parmams_nums()`)
.lookup_v <- function(stab, cellkeys, x_delta, cell_sum, pos_neg_var) {
  d <- max(stab$i) # parameter `D`

  v <- rep(NA, length(cellkeys))
  a <- cell_sum / x_delta

  if (pos_neg_var == 1) {
    # we only search in the symetric block
    a[1:length(a)] <- d
  } else if (pos_neg_var == 2) {
    # we use absolute values
    a <- abs(a)
  } else if (pos_neg_var == 3) {
    # we use absolute values for x_delta != 0 and the
    # symetric case for x_delta == 0
    a[a == 0] <- d
    ii <- a != 0
    a[ii] <- abs(a[ii])
  }
  a[a > d] <- d # we cut at the maximum value!

  ind_exact <- which(a %in% unique(stab$i))
  if (length(ind_exact) > 0) {
    v[ind_exact] <- sapply(ind_exact, function(x) {
      stab$diff[min(which(cellkeys[x] < stab$kum_p_o & stab$i == d))]
    })
  }

  ind_comb <- which(a < d)
  if (length(ind_comb) > 0) {
    v[ind_comb] <- sapply(ind_comb, function(x) {
      poss <- sort(unique(stab$i))
      a0 <- poss[which(x < poss) - 1]
      a1 <- poss[max(which(x > poss)) + 1]

      lambda <- (x - a0) / (a1 - a0)

      # perturbation_value (low)
      v0 <- stab$diff[min(which(cellkeys[x] < stab$kum_p_o & stab$i == a0))]

      # perturbation_value (up)
      v1 <- stab$diff[min(which(cellkeys[x] < stab$kum_p_o & stab$i == a1))]

      # combine to get final perturbation value
      (1 - lambda) * v0 + lambda * v1
    })
  }
  v
}

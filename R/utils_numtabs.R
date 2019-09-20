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
# not yet supported
.x_delta_grid <- function(x, inp_params) {
  stop("grids will be available in a future version!", call. = FALSE)

  even_odd <- x$even_odd
  x <- x$x
  top_k <- inp_params$top_k
  m_fixed_sq <- inp_params$m_fixed_sq

  # params is a list of parameters of length top_k
  m <- rep(NA, length(x))

  # footnote 3, page 7 -> if x < 0 -> use abs(x)
  rr <- lapply(1:length(x), function(i) {
    abs(x[i]) <= inp_params[[i]]$grid
  })

  # smallest cells, highest perturbation
  ind_sm <- sapply(rr, function(x) all(x == TRUE))
  m[ind_sm] <- sapply(inp_params[1:length(x)], function(x) x$pcts[1])

  # largest cells, smallest perturbation
  ind_lg <- sapply(rr, function(x) all(x == FALSE))
  m[ind_lg] <- sapply(inp_params[1:length(x)], function(x) utils::tail(x$pcts, 1))

  ind_mid <- which(is.na(m))
  if (length(ind_mid) > 0) {
    m[ind_mid] <- sapply(ind_mid, function(x) {
      inp_params[[x]]$pcts[min(which(rr[[x]] == TRUE))]
    })
  }

  # fixed variance for very small observations
  if (!is.na(m_fixed_sq)) {
    # separation point
    zs <- inp_params$zs

    # very small (absolute) values
    ind_vs <- which(abs(x) < zs)
    if (length(ind_vs) > 0) {
      m[ind_vs] <- sqrt(m_fixed_sq)
      x[ind_vs] <- 1
    }
  }
  x * m
}

# returns a list with the x_delta values (x * m) and
# where to look in the ptable (all_cells or small_cells)
.x_delta_flex <- function(x, inp_params) {
  even_odd <- x$even_odd
  x <- abs(x$x)

  m_fixed_sq <- inp_params$m_fixed_sq

  # flexpoint and separation points
  fp <- inp_params$fp
  zs <- inp_params$zs

  if (is.na(even_odd)) {
    lookup <- rep("all", length(x))
  } else {
    lookup <- ifelse(even_odd, "even", "odd")
  }

  m <- rep(inp_params$p_small, length(x))
  ind_lg <- which(x > fp)
  if (length(ind_lg) > 0) {
    x_lg <- x[ind_lg]

    m_lg <- inp_params$p_large * inp_params$epsilon[ind_lg]
    m_sm <- inp_params$p_small * inp_params$epsilon[ind_lg]

    f1 <- ((m_sm * x_lg) - (m_lg * fp)) / (m_lg * fp)
    f2 <- (2 * fp) / (fp + x_lg)
    m[ind_lg] <- m_lg * (1 + (f1 * f2 ^ inp_params$q))
  }

  # fixed variance for very small observations below separation point
  ind_vs <- which(x < zs)
  if (length(ind_vs) > 0) {
    # formula: correction.docx
    m[ind_vs] <- sqrt(m_fixed_sq) * (inp_params$E * inp_params$epsilon[ind_vs])
    x[ind_vs] <- 1
    lookup[ind_vs] <- "small_cells"
  }
  list(x_delta = x * m, lookup = lookup)
}

# lookup using the simple, non-grid version
.x_delta_simple <- function(x, inp_params) {
  even_odd <- x$even_odd
  x <- x$x

  p <- inp_params$p # default percentage
  if (is.na(even_odd)) {
    lookup <- rep("all", length(x))
  } else {
    lookup <- ifelse(even_odd, "even", "odd")
  }

  # separation_point
  zs <- inp_params$zs

  m <- rep(1, length(x))

  # for details, see scaling.docx (from pp)
  ind_lg <-  which(abs(x) >= zs)
  if (length(ind_lg) > 0) {
    m[ind_lg] <- p * inp_params$epsilon[ind_lg]
  }
  list(x_delta = abs(x * m), lookup = lookup)
}

# returning perturbation values from `stab` based on cellkeys,
# x_delta (x * m), the (absolute) values of the highest contributions
# in future; the lookup may depend on argument `pos_neg_var` that will need
# to be defined in `ck_parmams_nums()`
.lookup_v <- function(cellkeys, params) {
  d <- params$max_i # parameter `D`

  v <- rep(NA, length(cellkeys))
  a <- abs(params$x) / abs(params$x_delta)
  # this parameter is currently not implemented (see delivarable)
  # it may come back in a later version

  #if (pos_neg_var == 1) {
  #  # we only search in the symetric block
  #  a[1:length(a)] <- d
  #} else if (pos_neg_var == 2) {
  #  # we use absolute values
  #  a <- abs(a)
  #} else if (pos_neg_var == 3) {
  #  # we use absolute values for x_delta != 0 and the
  #  # symetric case for x_delta == 0
  #  a[a == 0] <- d
  #  ii <- a != 0
  #  a[ii] <- abs(a[ii])
  #}
  a[a > d] <- d # we cut at the maximum value!

  stab <- params$stab
  poss <- sort(unique(stab$i))

  ind_exact <- which(a %in% poss)
  if (length(ind_exact) > 0) {
    v[ind_exact] <- sapply(ind_exact, function(x) {
      stab[type == params$lookup[x] & i == a[x] & cellkeys[x] < kum_p_o, diff][1]
    })
  }

  ind_comb <- which(a < d)
  if (length(ind_comb) > 0) {
    v[ind_comb] <- sapply(ind_comb, function(x, poss) {
      a0 <- poss[which(a[x] < poss) - 1]
      a1 <- poss[max(which(a[x] > poss)) + 1]
      lambda <- (a[x] - a0) / (a1 - a0)

      # compute two perturbation values
      v_low <- stab[type == params$lookup[x] & i == a0 & cellkeys[x] < kum_p_o, diff][1]
      v_up <- stab[type == params$lookup[x] & i == a1 & cellkeys[x] < kum_p_o, diff][1]

      # combine to get final perturbation value
      (1 - lambda) * v_low + lambda * v_up
    }, poss = poss)
  }
  v
}

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

# g1 value (section 2.3.1 in deliverable 4.2);
# modified according to "correction.docx" from 28.8.19
.g1 <- function(m_fixed_sq, E, p_large) {
  # in case m_fixed_sq is NA (no separation) -> g1 is set to 0
  if (is.na(m_fixed_sq)) {
    return(0)
  }
  sqrt(m_fixed_sq) / (sqrt(E) * p_large)
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
    # compute g1 ("correction.docx" from sarah, 28.8.19)
    # modified formula 2.3 on p.9)
    g1 <- .g1(
      m_fixed_sq = m_fixed_sq,
      E = top_k,
      p_large = inp_params[[1]]$pcts[1])

    # very small (absolute) values
    ind_vs <- which(abs(x) < g1)
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
  x <- x$x

  m_fixed_sq <- inp_params$m_fixed_sq
  fp <- inp_params$fp
  E <- sum(inp_params$epsilon^2)

  if (is.na(even_odd)) {
    lookup <- rep("all", length(x))
  } else {
    lookup <- ifelse(even_odd, "even", "odd")
  }

  # g1 is 0 in case m_fixed_sq was not specified
  g1 <- .g1(
    m_fixed_sq = m_fixed_sq,
    E = E,
    p_large = inp_params$p_large)

  m <- rep(inp_params$p_small, length(x))
  ind_lg <- which(abs(x) > fp)
  if (length(ind_lg) > 0) {
    x_lg <- x[ind_lg]

    m_lg <- inp_params$p_large * inp_params$epsilon[ind_lg]
    m_sm <- inp_params$p_small * inp_params$epsilon[ind_lg]

    f1 <- ((m_sm * x_lg) - (m_lg * fp)) / (m_lg * fp)
    f2 <- (2 * fp) / (fp + x_lg)
    m[ind_lg] <- m_lg * (1 + (f1 * f2 ^ inp_params$q))
  }

  # fixed variance for very small observations
  # very small values
  ind_vs <- which(abs(x) < g1)
  if (length(ind_vs) > 0) {
    m[ind_vs] <- (sqrt(m_fixed_sq) * inp_params$epsilon[ind_vs]) / sqrt(E)
    x[ind_vs] <- 1
    lookup[ind_vs] <- "small_cells"
  }
  list(x_delta = abs(x * m), lookup = lookup)
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

  # g1 is 0 in case we do not have seperation
  # meaning no rows in ptable with type == "small_cells"
  E <- sum(inp_params$epsilon^2)
  g1 <- .g1(
    m_fixed_sq = inp_params$m_fixed_sq,
    E = E,
    p_large = p)

  m <- rep(1, length(x))

  # for details, see scaling.docx (from pp)
  ind_lg <-  which(abs(x) >= g1)
  if (length(ind_lg) > 0) {
    m[ind_lg] <- p * inp_params$epsilon[ind_lg]
  }
  list(x_delta = abs(x * m), lookup = lookup)
}

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

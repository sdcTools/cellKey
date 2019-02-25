# check if input is a valid ptable in "free" format
.valid_free_ptable <- function(ptab) {
  stopifnot(isS4(ptab))
  stopifnot(isS4(ptab), class(ptab) == "ptable")
  stopifnot(slot(ptab, "type") == "free")
  return(TRUE)
}

.check_free_ptable <- function(ptab) {
  .valid_free_ptable(ptab)
  tab <- slot(ptab, "pTable")

  for (j in 1:ncol(tab)) {
    fns <- unique(tab[[j]])
    for (x in seq_along(fns)) {
      r <- fns[[x]]()
      if (!is_scalar_integerish(r)) {
        e <- "Function returns a non-scalar integer value!"
        stop(e, call. = FALSE)
      }
    }
  }
  return(invisible(TRUE))
}

# convert bin to decimal
.bintodec <- function(y) {
  # find the decimal number corresponding to binary sequence 'y'
  if (!(all(y %in% 0:1))) {
    stop("not a binary sequence")
  }
  p <- (length(y):1) - 1
  res <- sum(y * 2 ^ p)
  return(res)
}

# compute cell key from record keys
.calc_cellkey <- function(rec_keys, big_n) {
  sum(rec_keys) %% big_n
}

# +1 if 9th bit of record_key is 1, -1 else
# used for perturbation of numerical variables
.direction <- function(rec_keys, type) {
  .abs <- function(rec_keys) {
    out <- rep(-1, length(rec_keys))
    rr <- sapply(1:length(rec_keys), function(x) {
      as.integer(intToBits(rec_keys[x]))[9] == 1
    })
    out[rr] <- 1
    out
  }
  .destatis <- function(rec_keys) {
    ifelse(rec_keys <= 0.5, -1, 1)
  }
  if (type %in% c("abs", "free")) {
    return(.abs(rec_keys))
  }
  if (type == "destatis") {
    return(.destatis(rec_keys))
  }
}

# in which row should we look in the perturbation table
# used for perturbation of numerical variables
.row_index_cont <- function(rec_keys, type) {
  .abs <- function(rec_keys) {
    sapply(1:length(rec_keys), function(x) {
      .bintodec(as.integer(intToBits(rec_keys[x]))[1:8]) + 1
    })
  }
  .destatis <- function(rec_keys) {
    as.numeric(cut(rec_keys, seq(0, 1, length = 257)))
  }

  if (type %in% c("abs", "free")) {
    return(.abs(rec_keys))
  }
  if (type == "destatis") {
    return(.destatis(rec_keys))
  }
}

# in which row should we look in the perturbation table
# used for perturbation of numerical variables
.col_index_cont <- function(ckey, n, small_c, type) {
  .abs <- function(ckey, n, small_c) {
    if (n <= small_c) {
      return(n + 32)
    }
    .bintodec(as.integer(intToBits(ckey))[1:5]) + 1
  }
  .destatis <- function(ckey, n, small_c) {
    if (n <= small_c) {
      return(n + 32)
    }
    as.numeric(cut(ckey, seq(0, 1, length = small_c + 1)))
  }

  if (type %in% c("abs", "free")) {
    return(.abs(ckey, n, small_c))
  }
  if (type == "destatis") {
    return(.destatis(ckey, n, small_c))
  }
}

# row index for perturbation - see page 11
# used for counts
.row_index_freq <- function(cKey) {
  bytes <- as.integer(intToBits(cKey))
  r <- bitwXor(bytes[1:8], bytes[9:16])
  r <- bitwXor(r, bytes[17:24])
  r <- bitwXor(r, bytes[25:32])

  p <- (length(r):1) - 1
  row_index <- sum(r * 2 ^ p) + 1
  row_index
}

# column index for perturbation - see page 11
# used for counts
.col_index_freq <- function(N, ptab_size, small_n) {
  if (N <= ptab_size - small_n) {
    return(N)
  } else {
    return(ptab_size - small_n + N %% small_n + 1)
  }
}

# perform actual lookup to get perturbation values;
# used in perturb_table()
.lookup <- function(tab, pert_params, ckeyname, freqvarname, type) {
  # lookup perturbation values using the abs-method
  .abs <- function(tab, pert_params, freqs, cellkeys) {
    ck <- row_indices <- col_indices <- pert <- NULL

    row_indices <- sapply(freqs, .row_index_freq)
    col_indices <- sapply(freqs, function(z) {
      .col_index_freq(
        N = z,
        ptab_size = slot(pert_params, "ptab_size"),
        small_n = slot(pert_params, "small_n")
      )
    })

    dt <- data.table(row_indices = row_indices, col_indices = col_indices)

    pert_vals <- lapply(1:nrow(dt), function(z) {
      ptab <- slot(pert_params, "ptab")
      ptab[dt[z, row_indices], dt[z, col_indices], with = FALSE]
    })
    ii <- which(sapply(pert_vals, function(x) nrow(x) != 1))
    if (length(ii) > 0) {
      pert_vals[ii] <- NA
    }
    dt[, pert := unlist(pert_vals)]
    dt[, ck := cellkeys]
    dt
  }

  # ptable in "free" format user defined functions
  .free <- function(tab, pert_params, freqs, cellkeys) {
    ck <- row_indices <- col_indices <- pert <- NULL

    row_indices <- sapply(freqs, .row_index_freq)
    col_indices <- sapply(freqs, function(z) {
      .col_index_freq(
        N = z,
        ptab_size = slot(pert_params, "ptab_size"),
        small_n = slot(pert_params, "small_n")
      )
    })

    if (any(col_indices < 0)) {
      e <- c(
        "Negative column indices were computed.",
        "Please provide a perturbation table with more columns."
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }

    dt <- data.table(row_indices = row_indices, col_indices = col_indices)
    pTab <- slot(pert_params, "ptab")
    pert_vals <- lapply(1:nrow(dt), function(z) {
      set.seed(cellkeys[z]) # reproducibility
      pTab[[dt$col_indices[z]]][[dt$row_indices[z]]]()
    })
    dt[, pert := unlist(pert_vals)]
    dt[, ck := cellkeys]
    dt
  }

  # lookup perturbation values if the perturbation table is in destatis
  # format. This function is based on the fifo()-function provided
  # by Tobias Enderle
  .destatis <- function(tab, pert_params, freqs, cellkeys) {
    i <- kum_p_o <- NULL
    ptab <- slot(pert_params, "ptab")
    symmetry <- max(ptab$i)

    row_indices <- rep(-1, nrow(tab))
    pert_vals <- rep(0L, nrow(tab))

    for (d in 1:(symmetry)) {
      if (d == symmetry) {
        rkind <- freqs >= d
      } else {
        rkind <- freqs == d
      }
      if (sum(rkind) > 0) {
        v <- ptab[i == d, kum_p_o]
        ck <- cellkeys[rkind]
        diffs <- ptab[i == d, diff]

        # row_ind
        rI <- sapply(1:sum(rkind), function(x) {
          which.max(ck[x] <= v)
        })
        pert_vals[rkind] <- as.integer(diffs[rI])
        row_indices[rkind] <- rI
      }
    }
    data.table(
      row_indices = row_indices,
      col_indices = NA,
      pert = pert_vals,
      ck = cellkeys
    )
  }

  stopifnot(is_scalar_character(ckeyname))
  stopifnot(ckeyname %in% names(tab))
  stopifnot(is_scalar_character(freqvarname))
  stopifnot(freqvarname %in% names(tab))

  cellkeys <- tab[, get(ckeyname)]
  freqs <- tab[, get(freqvarname)]

  if (type == "abs") {
    return(.abs(
      tab = tab,
      pert_params = pert_params,
      freqs = freqs,
      cellkeys = cellkeys
    ))
  }
  if (type == "free") {
    return(.free(
      tab = tab,
      pert_params = pert_params,
      freqs = freqs,
      cellkeys = cellkeys
    ))
  }
  if (type == "destatis") {
    return(.destatis(
      tab = tab,
      pert_params = pert_params,
      freqs = freqs,
      cellkeys = cellkeys
    ))
  }
}

# identify top k contributors to each cell and compute
# the required amount of perturbation
.identify_topk_cells <- function(dat, rkeys, dim_list, pert_params, v=v, type) {
  tmprkeyfortabulation <- tmpidforsorting <- NULL
  magnitude <- noise <- NULL

  keys <- names(dim_list)
  tmp <- copy(dat)

  mtab <- slot(pert_params, "mtab")
  stab <- slot(pert_params, "stab")
  big_n <- slot(pert_params, "big_n")
  small_c <- slot(pert_params, "small_c")

  tmp <- tmp[, c(v, keys, "tmpidforsorting"), with = FALSE]
  tmp[, tmprkeyfortabulation := rkeys]

  tmp$magnitude <- NA_real_
  tmp$dir <- NA_real_
  tmp$noise <- NA_real_

  # using absolute values for sorting
  tmp$tmpvarfortabulation <- abs(tmp[[v]])

  v.pert <- paste0(v, "_pert")
  v.mod <- paste0(v, "_mod")

  tmp[[v.pert]] <- 0
  tmp[[v.mod]] <- 0

  setorderv(tmp, c(keys, "tmpvarfortabulation"), order = -1)

  res <- NULL
  spl <- split(tmp, by = keys)

  res <- lapply(1:length(spl), function(i) {
    z <- spl[[i]]

    top_k <- min(nrow(z), slot(pert_params, "top_k"))
    z <- z[1:top_k]

    rkeys <- z$tmprkeyfortabulation

    z[, magnitude := mtab[1:top_k]]
    z[, dir := .direction(rkeys, type = type)]
    rind <- .row_index_cont(
      rec_keys = rkeys,
      type = type
    )
    cind <- .col_index_cont(
      ckey = .calc_cellkey(spl[[i]]$tmprkeyfortabulation, big_n),
      n = nrow(spl[[i]]),
      small_c = small_c,
      type = type
    )
    z[, noise := stab[[cind]][rind]]

    # .pert: the amount of perturbation
    z[, eval(v.pert) := z[[v]] * z$magnitude * z$dir * z$noise]

    # .mod: new values of weighted variable
    z[, eval(v.mod) := z[[v]] + z[[v.pert]]]
    z
  })

  mod <- rbindlist(res)
  setkey(mod, tmpidforsorting)

  vv <- c(
    "tmpidforsorting",
    keys,
    "magnitude",
    "dir",
    "noise",
    v,
    v.pert,
    v.mod
  )
  mod <- mod[, vv, with = FALSE]
  return(mod)
}

# x: values of a perturbed numerical variable
# pWC: perturbed weighted counts
.mean_before_sum <- function(x, pwc) {
  pw_mean <- x / pwc
  pw_sum <- pw_mean * pwc
  pw_mean[is.na(pw_mean)] <- 0
  data.table(pw_mean = pw_mean, pw_sum = pw_sum)
}
# x: values of a perturbed numerical variable
# pwc: perturbed weighted counts
.sum_before_mean <- function(x, pwc) {
  pw_sum <- x
  pw_mean <- x / pwc
  pw_mean[is.na(pw_mean)] <- 0
  data.table(pw_mean = pw_mean, pw_sum = pw_sum)
}

# simple check functions for record keys
.check_rkeys <- function(rkeys, type) {
  .abs <- function(rkeys) {
    stopifnot(is_integerish(rkeys))
    stopifnot(all(rkeys > 0))
    return(TRUE)
  }
  .destatis <- function(rkeys) {
    stopifnot(is_double(rkeys))
    stopifnot(all(rkeys >= 0))
    stopifnot(all(rkeys <= 1))
    return(TRUE)
  }

  if (type == "abs") {
    .abs(rkeys)
  }
  if (type == "destatis") {
    .destatis(rkeys)
  }
  return(TRUE)
}

# statistics of a distribution of a numeric vector
.distr_vals <- function(dd) {
  stopifnot(is.numeric(dd))
  dd <- na.omit(dd)
  vals <- c(
    min(dd),
    quantile(dd, seq(10, 40, by = 10) / 100),
    mean(dd),
    median(dd),
    quantile(dd, seq(60, 90, by = 10) / 100),
    max(dd)
  )
  vals <- round(vals, digits = 3)
  names(vals)[1] <- "Min"
  names(vals)[2:5] <- paste0("Q", seq(10, 40, by = 10))
  names(vals)[6] <- "Mean"
  names(vals)[7] <- "Median"
  names(vals)[8:11] <- paste0("Q", seq(60, 90, by = 10))
  names(vals)[12] <- "Max"
  vals
}

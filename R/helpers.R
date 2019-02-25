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
.calc_cellkey <- function(rec_keys, bigN) {
  sum(rec_keys) %% bigN
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
.col_index_cont <- function(cKey, n, smallC, type) {
  .abs <- function(cKey, n, smallC) {
    if (n <= smallC) {
      return(n + 32)
    }
    .bintodec(as.integer(intToBits(cKey))[1:5]) + 1
  }
  .destatis <- function(cKey, n, smallC) {
    if (n <= smallC) {
      return(n + 32)
    }
    as.numeric(cut(cKey, seq(0, 1, length = smallC + 1)))
  }

  if (type %in% c("abs", "free")) {
    return(.abs(cKey, n, smallC))
  }
  if (type == "destatis") {
    return(.destatis(cKey, n, smallC))
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
.col_index_freq <- function(N, pTableSize, smallN) {
  if (N <= pTableSize - smallN) {
    return(N)
  } else {
    return(pTableSize - smallN + N %% smallN + 1)
  }
}

# perform actual lookup to get perturbation values;
# used in perturb_table()
.lookup <- function(tab, pert_params, ckeyname, freqvarname, type) {
  # lookup perturbation values using the abs-method
  .abs <- function(tab, pert_params, freqs, cellkeys) {
    cK <- row_indices <- col_indices <- pert <- NULL

    row_indices <- sapply(freqs, .row_index_freq)
    col_indices <- sapply(freqs, function(z) {
      .col_index_freq(
        N = z,
        pTableSize = slot(pert_params, "pTableSize"),
        smallN = slot(pert_params, "smallN")
      )
    })

    dt <- data.table(row_indices = row_indices, col_indices = col_indices)

    pert_vals <- lapply(1:nrow(dt), function(z) {
      pert_params@pTable[dt[z, row_indices], dt[z, col_indices], with = FALSE]
    })
    ii <- which(sapply(pert_vals, function(x) nrow(x) != 1))
    if (length(ii) > 0) {
      pert_vals[ii] <- NA
    }
    dt[, pert := unlist(pert_vals)]
    dt[, cK := cellkeys]
    dt
  }

  # ptable in "free" format user defined functions
  .free <- function(tab, pert_params, freqs, cellkeys) {
    cK <- row_indices <- col_indices <- pert <- NULL

    row_indices <- sapply(freqs, .row_index_freq)
    col_indices <- sapply(freqs, function(z) {
      .col_index_freq(
        N = z,
        pTableSize = slot(pert_params, "pTableSize"),
        smallN = slot(pert_params, "smallN")
      )
    })

    if (any(col_indices < 0)) {
      e <- c(
        "Negative column indices were computed.",
        "Please provide a pTable with more columns."
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }

    dt <- data.table(row_indices = row_indices, col_indices = col_indices)
    pTab <- slot(pert_params, "pTable")
    pert_vals <- lapply(1:nrow(dt), function(z) {
      set.seed(cellkeys[z]) # reproducibility
      pTab[[dt$col_indices[z]]][[dt$row_indices[z]]]()
    })
    dt[, pert := unlist(pert_vals)]
    dt[, cK := cellkeys]
    dt
  }

  # lookup perturbation values if the perturbation table is in destatis
  # format. This function is based on the fifo()-function provided
  # by Tobias Enderle
  .destatis <- function(tab, pert_params, freqs, cellkeys) {
    i <- kum_p_o <- NULL
    pTable <- slot(pert_params, "pTable")
    symmetry <- max(pTable$i)

    row_indices <- rep(-1, nrow(tab))
    pert_vals <- rep(0L, nrow(tab))

    for (d in 1:(symmetry)) {
      if (d == symmetry) {
        rkind <- freqs >= d
      } else {
        rkind <- freqs == d
      }
      if (sum(rkind) > 0) {
        v <- pTable[i == d, kum_p_o]
        ck <- cellkeys[rkind]
        diffs <- pTable[i == d, diff]

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
      cK = cellkeys
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
  is_topK <- magnitude <- noise <- NULL

  keys <- names(dim_list)
  tmp <- copy(dat)

  mTable <- slot(pert_params, "mTable")
  sTable <- slot(pert_params, "sTable")
  bigN <- slot(pert_params, "bigN")
  smallC <- slot(pert_params, "smallC")

  tmp <- tmp[, c(v, keys, "tmpidforsorting"), with = FALSE]
  tmp[, tmprkeyfortabulation := rkeys]

  tmp$is_topK <- FALSE
  tmp$magnitude <- NA_real_
  tmp$dir <- NA_real_
  tmp$noise <- NA_real_

  # using absolute values for sorting
  tmp$tmpvarfortabulation <- abs(tmp[[v]])

  v.pert <- paste0(v, ".pert")
  v.mod <- paste0(v, ".mod")

  tmp[[v.pert]] <- 0
  tmp[[v.mod]] <- 0

  setorderv(tmp, c(keys, "tmpvarfortabulation"), order = -1)

  res <- NULL
  spl <- split(tmp, by = keys)

  res <- lapply(1:length(spl), function(i) {
    z <- spl[[i]]

    topK <- min(nrow(z), slot(pert_params, "topK"))
    z <- z[1:topK]

    rkeys <- z$tmprkeyfortabulation

    z[, is_topK := TRUE]
    z[, magnitude := mTable[1:topK]]

    z[, dir := .direction(rkeys, type = type)]
    rind <- .row_index_cont(
      rec_keys = rkeys,
      type = type
    )
    cind <- .col_index_cont(
      cKey = .calc_cellkey(spl[[i]]$tmprkeyfortabulation, bigN),
      n = nrow(spl[[i]]),
      smallC = smallC,
      type = type
    )

    z[, noise := sTable[[cind]][rind]]

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
.mean_before_sum <- function(x, pWC) {
  pWMean <- x / pWC
  pWSum <- pWMean * pWC
  pWMean[is.na(pWMean)] <- 0
  data.table(pWMean = pWMean, pWSum = pWSum)
}
# x: values of a perturbed numerical variable
# pWC: perturbed weighted counts
.sum_before_mean <- function(x, pWC) {
  pWSum <- x
  pWMean <- x / pWC
  pWMean[is.na(pWMean)] <- 0
  data.table(pWMean = pWMean, pWSum = pWSum)
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

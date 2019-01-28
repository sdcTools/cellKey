# convert bin to decimal
bintodec <- function(y) {
  # find the decimal number corresponding to binary sequence 'y'
  if (!(all(y %in% c(0, 1)))) {
    stop("not a binary sequence")
  }
  b <- (length(y):1) - 1
  res <- sum(y * 2 ^ b)
  return(res)
}

# compute cell key from record keys
calc_cKey <- function(rec_keys, bigN) {
  sum(as.numeric(rec_keys)) %% bigN
}

# +1 if 9th bit of record_key is 1, -1 else
# used for perturbation of numerical variables
get_direction <- function(rec_keys, type) {
  get_direction_abs <- function(rec_keys) {
    out <- rep(-1, length(rec_keys))
    rr <- sapply(1:length(rec_keys), function(x) {
      as.integer(intToBits(rec_keys[x]))[9] == 1
    })
    out[rr] <- 1
    out
  }
  get_direction_destatis <- function(rec_keys) {
    ifelse(rec_keys <= 0.5, -1, 1)
  }
  stopifnot(is_scalar_character(type))
  stopifnot(type %in% c("abs", "destatis"))

  if (type == "abs") {
    return(get_direction_abs(rec_keys))
  }
  if (type == "destatis") {
    return(get_direction_destatis(rec_keys))
  }
  stop("error in get_direction()\n")
}



# in which row should we look in the perturbation table
# used for perturbation of numerical variables

get_row_index_cont <- function(rec_keys, type) {
  get_row_index_cont_abs <- function(rec_keys) {
    out <- rep(-1, length(rec_keys))
    sapply(1:length(rec_keys), function(x) {
      bintodec(as.integer(intToBits(rec_keys[x]))[1:8]) + 1
    })
  }
  get_row_index_cont_destatis <- function(rec_keys) {
    out <- rep(-1, length(rec_keys))
    as.numeric(cut(rec_keys, seq(0, 1, length = 257)))
  }

  stopifnot(is_scalar_character(type))
  stopifnot(type %in% c("abs", "destatis"))

  if (type == "abs") {
    return(get_row_index_cont_abs(rec_keys))
  }
  if (type == "destatis") {
    return(get_row_index_cont_destatis(rec_keys))
  }
  stop("error in get_row_index_cont()\n")
}

# in which row should we look in the perturbation table
# used for perturbation of numerical variables
get_col_index_cont <- function(cKey, n, smallC, type) {
  get_col_index_cont_abs <- function(cKey, n, smallC) {
    if (n <= smallC) {
      return(n + 32)
    }
    bintodec(as.integer(intToBits(cKey))[1:5]) + 1
  }
  get_col_index_cont_destatis <- function(cKey, n, smallC) {
    if (n <= smallC) {
      return(n + 32)
    }
    as.numeric(cut(cKey, seq(0, 1, length = smallC + 1)))
  }

  stopifnot(is_scalar_character(type))
  stopifnot(type %in% c("abs", "destatis"))

  if (type == "abs") {
    return(get_col_index_cont_abs(cKey, n, smallC))
  }
  if (type == "destatis") {
    return(get_col_index_cont_destatis(cKey, n, smallC))
  }
  stop("error in get_col_index_cont_destatis()\n")
}

# row index for perturbation - see page 11
# used for counts
get_rowIndex <- function(cKey) {
  bytes <- as.integer(intToBits(cKey))
  r <- bitwXor(bytes[1:8], bytes[9:16])
  r <- bitwXor(r, bytes[17:24])
  r <- bitwXor(r, bytes[25:32])

  b <- (length(r):1) - 1
  row_index <- sum(r * 2 ^ b) + 1
  row_index
}

# column index for perturbation - see page 11
# used for counts
get_colIndex <- function(N, pTableSize, smallN) {
  if (N <= pTableSize - smallN) {
    return(N)
  } else {
    return(pTableSize - smallN + N %% smallN + 1)
  }
}

# perform actual lookup to get perturbation values;
# used in perturbTable()
lookup <- function(tab, pert_params, ckeyname, freqvarname, type) {
  # lookup perturbation values using the abs-method
  # lookup_abs() is used in perturbTable()
  lookup_abs <- function(tab, pert_params, freqs, cellkeys) {
    . <- cK <- N <- row_indices <- col_indices <- pert <- NULL

    row_indices <- sapply(freqs, get_rowIndex)
    col_indices <- sapply(freqs, function(z) {
      get_colIndex(z, slot(pert_params, "pTableSize"), slot(pert_params, "smallN"))
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

  # lookup perturbation values if the perturbation table is in destatis
  # format. This function is based on the fifo()-function provided
  # by Tobias Enderle
  # lookup_destatis() is used in perturbTable()
  lookup_destatis <- function(tab, pert_params, freqs, cellkeys) {
    sumW <- i <- kum_p_o <- tmpcellkey <- col_indices <- pert  <- NULL
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

  stopifnot(is_scalar_character(type))
  stopifnot(is_scalar_character(ckeyname))
  stopifnot(ckeyname %in% names(tab))

  stopifnot(is_scalar_character(freqvarname))
  stopifnot(freqvarname %in% names(tab))

  stopifnot(type %in% c("abs", "destatis"))

  cellkeys <- tab[, get(ckeyname)]
  freqs <- tab[, get(freqvarname)]

  if (type == "abs") {
    return(
      lookup_abs(
        tab = tab,
        pert_params = pert_params,
        freqs = freqs,
        cellkeys = cellkeys
      )
    )
  }
  if (type == "destatis") {
    return(
      lookup_destatis(
        tab = tab,
        pert_params = pert_params,
        freqs = freqs,
        cellkeys = cellkeys
      )
    )
  }
  stop("error in lookup()\n")
}

check_weight <- function(dat, w) {
  tmpweight_for_tabulation <- NULL
  w_new <- "tmpweight_for_tabulation"
  if (!is.null(w)) {
    stopifnot(w %in% names(dat))
    stopifnot(length(w) == 1)
    expr <- paste0("dat[,", w_new, ":=", w, "]")
    eval(parse(text = expr))
    is_weighted <- TRUE
  } else {
    dat[, tmpweight_for_tabulation := 1]
    is_weighted <- FALSE
  }
  list(
    dat = dat,
    is_weighted = is_weighted,
    w_new = w_new
  )
}

identify_topK_cells <- function(dat, rkeys, dimList, pert_params, v=v, type) {
  tmprkeyfortabulation <- is_topK <- magnitude <- noise <- NULL
  tmpidforsorting <- NULL

  keys <- names(dimList)
  tmp <- copy(dat)

  mTable <- slot(pert_params, "mTable")
  sTable <- slot(pert_params, "sTable")
  bigN <- slot(pert_params, "bigN")
  smallC <- slot(pert_params, "smallC")

  tmp <- tmp[, c(v, keys, "tmpidforsorting"), with = FALSE]
  tmp[, tmprkeyfortabulation := rkeys]

  # groups defined by key-variables
  setkeyv(tmp, keys)

  tmp[, is_topK := FALSE]
  tmp[, magnitude := NA_real_]
  tmp[, dir := NA_real_]
  tmp[, noise := NA_real_]

  # using absolute values for sorting
  tmp[, ("tmpvarfortabulation") := abs(get(v))]

  res <- NULL
  spl <- split(tmp, by = keys)
  v.pert <- paste0(v, ".pert")
  v.mod <- paste0(v, ".mod")
  res <- lapply(1:length(spl), function(i) {
    z <- spl[[i]]
    setorderv(z, "tmpvarfortabulation", order = -1)

    topK <- min(nrow(z), slot(pert_params, "topK"))
    mTab <- mTable[1:topK]

    z[1:topK, is_topK := TRUE]
    z[1:topK, magnitude := mTab]

    d <- get_direction(
      z[1:topK, tmprkeyfortabulation],
      type = type
    )
    z[1:topK, dir := d]
    rind <- get_row_index_cont(
      rec_keys = z[1:topK, tmprkeyfortabulation],
      type = type
    )
    cind <- get_col_index_cont(
      cKey = calc_cKey(z[, tmprkeyfortabulation], bigN),
      n = nrow(z),
      smallC = smallC,
      type = type
    )

    z[1:topK, noise := unlist(sTable[rind, cind, with = FALSE])]

    # .pert is the actual noise
    z[, eval(v.pert) := 0]
    z[is_topK == TRUE, eval(v.pert) := get(v) * magnitude * dir * noise]

    # .mod: new weighted values
    z[, eval(v.mod) := get(v) + get(v.pert)]
  })

  res <- rbindlist(res)
  setkey(res, tmpidforsorting)
  # modified
  mod <- res[is_topK == TRUE]
  mod <- mod[,
    c("tmpidforsorting", keys, "magnitude", "dir", "noise", v, v.pert, v.mod), with = FALSE]
  return(mod)
}

# x: values of a perturbed numerical variable
# pWC: perturbed weighted counts
mean_before_sum <- function(x, pWC) {
  pWMean <- x / pWC
  pWSum <- pWMean * pWC
  pWMean[is.na(pWMean)] <- 0
  data.table(pWMean = pWMean, pWSum = pWSum)
}
# x: values of a perturbed numerical variable
# pWC: perturbed weighted counts
sum_before_mean <- function(x, pWC) {
  pWSum <- x
  pWMean <- x / pWC
  pWMean[is.na(pWMean)] <- 0
  data.table(pWMean = pWMean, pWSum = pWSum)
}

# simple check functions for record keys
check_rkeys <- function(rkeys, type) {
  check_rkeys_abs <- function(rkeys) {
    stopifnot(is_integerish(rkeys))
    stopifnot(all(rkeys > 0))
    return(TRUE)
  }
  check_rkeys_destatis <- function(rkeys) {
    stopifnot(is_double(rkeys))
    stopifnot(all(rkeys >= 0))
    stopifnot(all(rkeys <= 1))
    return(TRUE)
  }

  if (type == "abs") {
    check_rkeys_abs(rkeys)
  }
  if (type == "destatis") {
    check_rkeys_destatis(rkeys)
  }
  return(TRUE)
}

# statistics of a distribution of a numeric vector
get_distr_vals <- function(dd) {
  stopifnot(is.numeric(dd))
  dd <- na.omit(dd)
  vals <- c(
    min(dd), quantile(dd, seq(10, 40, by = 10) / 100),
    mean(dd), median(dd),
    quantile(dd, seq(60, 90, by = 10) / 100), max(dd))
  vals <- round(vals, digits = 3)
  names(vals)[1] <- "Min"
  names(vals)[2:5] <- paste0("Q", seq(10, 40, by = 10))
  names(vals)[6] <- "Mean"
  names(vals)[7] <- "Median"
  names(vals)[8:11] <- paste0("Q", seq(60, 90, by = 10))
  names(vals)[12] <- "Max"
  vals
}

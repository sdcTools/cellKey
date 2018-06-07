# convert bin to decimal
bintodec <- function(y) {
  # find the decimal number corresponding to binary sequence 'y'
  if (!(all(y %in% c(0,1)))) {
    stop("not a binary sequence")
  }
  res <- sum(y*2^((length(y):1) - 1))
  return(res)
}

# compute cell key from record keys
calc_cKey <- function(rec_keys, bigN) {
  sum(as.numeric(rec_keys)) %% bigN
}

# +1 if 9th bit of record_key is 1, -1 else
# used for perturbation of numerical variables
get_direction <- function(rec_keys) {
  out <- rep(-1, length(rec_keys))
  rr <- sapply(1:length(rec_keys), function(x) {
    as.integer(intToBits(rec_keys[x]))[9]==1
  })
  out[rr] <- 1
  out
}

# in which row should we look in the perturbation table
# used for perturbation of numerical variables
get_row_index_cont <- function(rec_keys) {
  out <- rep(-1, length(rec_keys))
  sapply(1:length(rec_keys), function(x) {
    bintodec(as.integer(intToBits(rec_keys[x]))[1:8])+1
  })
}

# in which row should we look in the perturbation table
# used for perturbation of numerical variables
get_col_index_cont <- function(cKey, n, smallC) {
  if (n<=smallC) {
    return(n+32)
  }
  bintodec(as.integer(intToBits(cKey))[1:5])+1
}

# row index for perturbation - see page 11
# used for counts
get_rowIndex <- function(cKey) {
  bytes <- as.integer(intToBits(cKey))
  r <- bitwXor(bytes[1:8], bytes[9:16])
  r <- bitwXor(r, bytes[17:24])
  r <- bitwXor(r, bytes[25:32])

  row_index <- sum(r*2^((length(r):1) - 1))+1
  row_index
}

# column index for perturbation - see page 11
# used for counts
get_colIndex <- function(N, pTableSize, smallN) {
  if (N<=pTableSize-smallN) {
    return(N)
  } else {
    return(pTableSize-smallN+N%%smallN + 1)
  }
}

# lookup perturbation values using the abs-method
# lookup_abs() is used in perturbTable()
lookup_abs <- function(tab, pert_params) {
  . <- CKey <- N <- row_indices <- col_indices <- pert <- NULL

  set(tab, j="row_indices", value=sapply(tab[,CKey], get_rowIndex))
  set(tab, j="col_indices", value=sapply(tab[,N], function(z) {
    get_colIndex(z, slot(pert_params, "pTableSize"), slot(pert_params, "smallN"))
  }))

  df <- tab[,.(row_indices, col_indices)]

  pert_vals <- lapply(1:nrow(df), function(z) {
    pert_params@pTable[df[z, row_indices], df[z, col_indices], with=F]
  })
  ii <- which(sapply(pert_vals, function(x) nrow(x)!=1))
  if (length(ii)>0) {
    pert_vals[ii] <- NA
  }
  pert_vals <- unlist(pert_vals)
  tab[,pert:=pert_vals]
  return(tab)
}

# lookup perturbation values if the perturbation table is in destatis
# format. This function is based on the fifo()-function provided
# by Tobias Enderle
# lookup_destatis() is used in perturbTable()
lookup_destatis <- function(tab, pert_params, symmetry=8) {
  sumW <- i <- kum_p_o <- CKey <- col_indices <- pert <- tmpfreqforcalc <- NULL
  pTable <- slot(pert_params, "pTable")
  tab[,tmpfreqforcalc:=sumW]
  row_indices <- rep(-1, nrow(tab))
  pert_vals <- rep(0L, nrow(tab))

  for (d in 1:(symmetry)) {
    if (d==symmetry) {
      rkind <- tab[,tmpfreqforcalc>=d]
    } else {
      rkind <- tab[,tmpfreqforcalc==d]
    }
    if (sum(rkind)>0) {
      v <- pTable[i==d,kum_p_o]
      ck <- tab[rkind, CKey]
      diffs <- pTable[i==d,diff]

      # row_ind
      rI <- sapply(1:sum(rkind), function(x) {
        which.max(ck[x] <= v)
      })
      pert_vals[rkind] <- as.integer(diffs[rI])
      row_indices[rkind] <- rI
    }
  }

  tab[,tmpfreqforcalc:=NULL]
  tab[,row_indices:=row_indices]
  tab[,col_indices:=NA]
  tab[,pert:=pert_vals]
  return(tab)
}

check_weight <- function(dat, w) {
  tmpweight_for_tabulation <- NULL
  w_new <- "tmpweight_for_tabulation"
  if (!is.null(w)) {
    stopifnot(w %in% names(dat))
    stopifnot(length(w)==1)
    expr <- paste0("dat[,",w_new,":=",w,"]")
    eval(parse(text=expr))
    is_weighted <- TRUE
  } else {
    dat[,tmpweight_for_tabulation:=1]
    is_weighted <- FALSE
  }
  return(list(dat=dat, is_weighted=is_weighted, w_new=w_new))
}

identify_topK_cells <- function(dat, rkeys, dimList, pert_params, v=v) {

  tmprkeyfortabulation <- is_topK <- magnitude <- noise <- NULL
  tmpidforsorting <- NULL

  keys <- names(dimList)
  tmp <- copy(dat)

  mTable <- slot(pert_params, "mTable")
  sTable <- slot(pert_params, "sTable")
  bigN <- slot(pert_params, "bigN")
  smallC <- slot(pert_params, "smallC")

  # groups defined by key-variables
  setkeyv(tmp, keys)

  tmp <- tmp[,c(v, keys, "tmpidforsorting"), with=F]
  tmp[,tmprkeyfortabulation:=rkeys]

  tmp[,is_topK:=FALSE]
  tmp[,magnitude:=NA_real_]
  tmp[,dir:=NA_real_]
  tmp[,noise:=NA_real_]

  # using absolute values for sorting
  tmp[,("tmpvarfortabulation"):=abs(get(v))]

  res <- NULL
  spl <- split(tmp, by=keys)
  v.pert <- paste0(v,".pert")
  v.mod <- paste0(v,".mod")
  res <- lapply(1:length(spl), function(i) {
    z <- spl[[i]]
    setorderv(z, "tmpvarfortabulation", order=-1)

    topK <- min(nrow(z), slot(pert_params, "topK"))
    mTab <- mTable[1:topK]

    z[1:topK, is_topK:=TRUE]
    z[1:topK, magnitude:=mTab]
    z[1:topK, dir:=get_direction(z[1:topK,tmprkeyfortabulation])]
    rind <- get_row_index_cont(rec_keys=z[1:topK,tmprkeyfortabulation])
    cind <- get_col_index_cont(cKey=calc_cKey(z[,tmprkeyfortabulation], bigN), n=nrow(z), smallC=smallC)
    z[1:topK, noise:=unlist(sTable[rind, cind, with=F])]

    # .pert die tatsaechliche verschmutzung
    z[,eval(v.pert):=0]
    z[is_topK==TRUE, eval(v.pert):=get(v)*magnitude*dir*noise]

    # .mod: die neuen gewichteten werte der variablen
    z[,eval(v.mod):=get(v)+get(v.pert)]
  })

  res <- rbindlist(res)
  setkey(res, tmpidforsorting)
  # modified
  mod <- res[is_topK==TRUE]
  mod <- mod[,c("tmpidforsorting", keys, "magnitude", "dir","noise",v, v.pert, v.mod), with=F]
  return(mod)
}

# x: values of a perturbed numerical variable
# pWC: perturbed weighted counts
mean_before_sum <- function(x, pWC) {
  pWMean <- x / pWC
  pWSum <- pWMean * pWC
  pWMean[is.na(pWMean)] <- 0
  data.table(pWMean=pWMean, pWSum=pWSum)
}
# x: values of a perturbed numerical variable
# pWC: perturbed weighted counts
sum_before_mean <- function(x, pWC) {
  pWSum <- x
  pWMean <- x / pWC
  pWMean[is.na(pWMean)] <- 0
  data.table(pWMean=pWMean, pWSum=pWSum)
}

# simple check functions for record keys
check_rkeys <- function(rkeys, type) {
  check_rkeys_abs <- function(rkeys) {
    stopifnot(is_integerish(rkeys))
    stopifnot(all(rkeys>0))
    return(TRUE)
  }
  check_rkeys_destatis <- function(rkeys) {
    stopifnot(is_double(rkeys))
    stopifnot(all(rkeys)>=0)
    stopifnot(all(rkeys)<=1)
    return(TRUE)
  }

  if (type=="abs") {
    check_rkeys_abs(rkeys)
  }
  if (type=="destatis") {
    check_rkeys_destatis(rkeys)
  }
  return(TRUE)
}

# fixes tables with negative counts after perturbation
# by computing perturbation value times -1
fix_negative_counts <- function(tab) {
  neg_counts <- pert <- N <- NULL
  tab[,neg_counts:=FALSE]
  tab[pert+N<0, neg_counts:=TRUE]
  if (sum(tab[,neg_counts]) >0) {
    warning("after perturbations we got negative counts!\n")
    warning("for now, we compute pert_val=-1*pert_val for such cases, but this should be fixed in pTables\n")
    tab[neg_counts==TRUE, pert:=-1*pert]
  }
  tab[,neg_counts:=NULL]
  tab
}

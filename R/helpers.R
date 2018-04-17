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

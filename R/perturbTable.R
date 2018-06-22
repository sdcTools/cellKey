#' perturbTable
#'
#' compute perturbed tables containing counts and continous variables.
#'
#' @param inp an object of class \code{\link{pert_inputdat-class}}
#' generated with \code{\link{ck_create_input}}
#' @param dimList a list containing slots for each variable that should be
#' tabulated. Each slot consists of a \code{data.frame} or \code{data.table}
#' with columns \code{levels} and \code{codes} specifiying the
#' hierarchies of the variables.
#' @param countVars (character) vector of numerical variables that should be
#' counted and perturbed; must be 0/1 coded.
#' @param numVars (character) vector of numerical variables that
#' should be tabulated or \code{NULL}
#' @param by \code{NULL} or a scalar character. If specified, it is possible to use a valid variable
#' name being 0/1 coded that will be used to create the perturbation of all variables
#' specified in \code{numVars} only for subgroups defined by the specified count variable. If the
#' variable is not listed in argument \code{countVars} it will be automatically added. If not specified or
#' \code{NULL}, the table will be based on all units in the microdata.
#' @param weightVar (character) vector of variable holding sampling
#' weights or \code{NULL}
#' @return an object of class \code{\link{pert_table-class}}.
#' @seealso \url{https://www.unece.org/fileadmin/DAM/stats/documents/ece/ces/ge.46/2013/Topic_1_ABS.pdf}
#' @export
#'
#' @examples
#' ## loading testdata and adding record keys (for abs-algorithm)
#' dat <- ck_create_testdata()
#'
#' ## create some 0/1 variables that should be perturbed later
#' dat[,cnt_females:=ifelse(sex=="male", 0, 1)]
#' dat[,cnt_males:=ifelse(sex=="male", 1, 0)]
#' dat[,cnt_highincome:=ifelse(income>=9000, 1, 0)]
#'
#' ## create record keys
#' dat$rkeys <- ck_generate_rkeys(dat=dat, max_val=2*nrow(dat), type="abs", verbose=TRUE)
#'
#' ## loading test perturbation in abs- and destatis format
#' pTable_abs <- ck_create_pTable(pTableSize=70, type="abs")
#' pTable_destatis <- ck_create_pTable(type="destatis")
#'
#' ## perturbation tables and parameters for numerical variables
#' sTable <- ck_generate_sTable(smallC=12)
#' mTable <- c(0.6,0.4,0.2)
#'
#' ## creating perturbation objects that are later used
#' ## as inputs in perturbTable()
#' pert_params_abs <- ck_create_pert_params(
#'   bigN=17312941, smallN=12,
#'   pTable=pTable_abs,
#'   sTable=sTable,
#'   mTable=mTable)
#'
#' pert_params_destatis <- ck_create_pert_params(
#'   smallN=12,
#'   pTable=pTable_destatis,
#'   sTable=sTable,
#'   mTable=mTable)
#'
#' ## create suitable input data using existing record keys
#' ## and a perturbation table in "abs-format"
#' inp_abs <- ck_create_input(
#'   dat=dat,
#'   def_rkey="rkeys",
#'   pert_params=pert_params_abs)
#'
#' ## create input using a perturbation table in "destatis-format"
#' ## and no existing record keys. In this case, def_rkey=9 means
#' ## that the generated record keys are rounded to 9 digits
#' inp_destatis <- ck_create_input(
#'   dat=dat,
#'   def_rkey=9,
#'   pert_params=pert_params_destatis)
#'
#' ## definining table hierarchies
#' ## variable 'sex'
#' dim.sex <- ck_create_node(total_lab="Total")
#' dim.sex <- ck_add_nodes(dim.sex, reference_node="Total", node_labs=c("male","female"))
#' print(dim.sex)
#'
#' ## variable 'sex'
#' dim.age <- ck_create_node(total_lab="Total")
#' dim.age <- ck_add_nodes(dim.age,
#'   reference_node="Total",
#'   node_labs=paste0("age_group",1:6))
#' print(dim.age)
#'
#' ## define required inputs
#' dimList <- list(sex=dim.sex, age=dim.age)
#' weightVar <- "sampling_weight"
#' numVars <- c("savings", "income")
#'
#' ## perturbing the table
#' res <- perturbTable(
#'   inp=inp_abs,
#'   dimList=dimList,
#'   weightVar=weightVar, countVars=NULL, numVars=numVars)
#' ck_freq_table(res, vname="Total")
#' ck_cont_table(res, vname="income", meanBeforeSum=TRUE)
#' ck_cont_table(res, vname="savings", meanBeforeSum=TRUE)
#'
#' res <- perturbTable(
#'   inp=inp_destatis,
#'   dimList=dimList,
#'   weightVar=weightVar, countVars=NULL, numVars=numVars)
#' ck_freq_table(res, vname="Total")
#' ck_cont_table(res, vname="income", meanBeforeSum=TRUE)
#' ck_cont_table(res, vname="savings", meanBeforeSum=TRUE)
#'
#' ## information on modifications for count variables
#' mod_counts(res, verbose=TRUE)
#'
#' ## information on modifications for numerical variables
#' mod_numvars(res)
#'
#' ## an example using additional countVars
#' res <- perturbTable(
#'   inp=inp_destatis,
#'   dimList=dimList,
#'   weightVar=weightVar,
#'   countVars=c("cnt_females", "cnt_males","cnt_highincome"),
#'   numVars=numVars)
#'
#' ## show count tables
#' ck_freq_table(res, vname=NULL)
#' ck_freq_table(res, vname="Total")
#' ck_freq_table(res, vname="cnt_females")
#' ck_freq_table(res, vname="cnt_males")
#' ck_freq_table(res, vname="cnt_highincome")
#'
#' ## show magnitude tables alog with its modifications
#' p_inc <- ck_cont_table(res, vname="income"); p_inc
#' attr(p_inc, "modifications")
#'
#' ## create perturbed magnitude table for subgroups of data
#' res <- perturbTable(
#'   inp=inp_destatis,
#'   dimList=dimList,
#'   weightVar=weightVar,
#'   countVars=c("cnt_females", "cnt_males","cnt_highincome"),
#'   numVars=numVars, by="cnt_males")
#'
#' ## the perturbed tables of savings and income are computed only where
#' ## variable cnt_males is 1
#' p_sav <- ck_cont_table(res, vname="savings"); p_sav
#' attr(p_sav, "modifications") # no modifications in cells containing females!
perturbTable <- function(inp, dimList, countVars=NULL, numVars=NULL, by=NULL, weightVar=NULL) {
  # rename variables
  gen_vnames <- function(countVars, prefix="sumRK") {
    if (length(countVars)==0) {
      return(paste0(prefix,"_Total"))
    }
    c(paste0(prefix,"_Total"), paste0(prefix,"_",countVars))
  }

  . <- CKey <- N <- WCavg <- col_indices <- neg_counts <- pUWC <- pWC <-  NULL
  row_indices <- sdcStatus <- strID <- sumRec <- sumW <- tmpfreqvarfortabulation <- NULL
  tmprkeysfortabulation <- tmpweightvarfortabulation <- tmpidforsorting <- NULL
  WCavg_Total <- pert <- countVar <- numVar <- NULL

  stopifnot(isS4(inp))
  stopifnot(class(inp)=="pert_inputdat")
  pert_params <- slot(inp, "pert_params")
  type <- slot(pert_params, "type")

  dat  <- copy(slot(inp, "microdat"))
  dat[,tmprkeysfortabulation:=slot(inp, "rkeys")]
  dat[,tmpidforsorting:=.I]

  for (i in 1:length(dimList)) {
    if ("nodedim" %in% class(dimList[[i]])) {
      dimList[[i]] <- node_to_sdcinput(dimList[[i]], addNumLevels=FALSE)
    }
  }

  dV <- match(names(dimList), names(dat))
  if (any(is.na(dV))) {
    stop("check input 'dimList': --> at least one variable was not found in the input dataset\n")
  }

  freqVar <- "tmpfreqvarfortabulation"
  dat[,c(freqVar):=1]

  dat[,tmpweightvarfortabulation:=1]
  if (!is.null(weightVar)) {
    stopifnot(weightVar %in% names(dat))
    dat[,tmpweightvarfortabulation:=get(weightVar)]
  }

  ## check by-argument
  if (!is.null(by)) {
    stopifnot(is_scalar_character(by))
    if (!by %in% countVars) {
      countVars <- c(countVars, by)
    }
  } else {
    by <- "Total"
  }

  # check and prepare countVars
  if (is.null(countVars)) {
    countVars <- countVars_w <- countVars_rec <- c()
  } else {
    # checking if all count_vars exist!
    cV <- match(countVars, names(dat))
    if (any(is.na(cV))) {
      stop("check input 'countVars': --> at least one variable was not found in the input dataset\n")
    }
    # check if all count-vars are integerish 0/1
    chk <- sapply(countVars, function(x) {
      identical(range(dat[[x]]), c(0,1))
    })
    if (!all(chk)) {
      stop("check input 'countVars': --> at least one variable is not 0/1 coded!\n")
    }

    # compute weighted version!
    countVars_w <- paste0("wc_", countVars)
    dat[,(countVars_w):=lapply(.SD, function(x) { x*tmpweightvarfortabulation}), .SDcols=countVars]

    # record_keys using an indicator
    countVars_rec <- paste0("recKey_", countVars)
    dat[,(countVars_rec):=lapply(.SD, function(x) { x*tmprkeysfortabulation}), .SDcols=countVars]
  }

  countVar_names <- c("Total", countVars)

  fV <- match(freqVar, names(dat))
  nV <- c("tmprkeysfortabulation","tmpweightvarfortabulation", countVars, countVars_w, countVars_rec)
  pert_info_cont <- list()
  if (!is.null(numVars)) {
    nV2 <- match(numVars, names(dat))
    if (any(is.na(nV2))) {
      stop("check numVars!\n")
    }
    nV <- c(nV, numVars)
    length(pert_info_cont) <- length(numVars)
    names(pert_info_cont) <- numVars
    for (i in 1:length(numVars)) {
      v <- numVars[i]
      vnum_orig <- paste0(v, "_tmporigforcalculation")
      nV <- c(nV, vnum_orig)

      # copy numerical variable
      if (by=="Total") {
        vals_orig <- dat[,get(v)]
        set(dat, j=vnum_orig, value=vals_orig)
      } else {
        # we need to set original numerical variable to 0
        # where by variable is also 0!
        dat[,(numVars[i]):=get(numVars[i])*get(by)]
        vals_orig <- dat[,get(v)*get(by)]
        set(dat, j=vnum_orig, value=vals_orig)
      }
      rr <- identify_topK_cells(dat=dat, rkeys=slot(inp, "rkeys"),
        dimList=dimList, pert_params=pert_params, v=v, type=type)
      dat[rr[,tmpidforsorting],c(v):= rr[,get(paste0(v,".mod"))]]
      setnames(rr, "tmpidforsorting", "id")

      setnames(rr, c(v, paste0(v, c(".pert",".mod"))), c(paste0("vals.", c("orig","pert","mod"))))
      rr[,numVar:=v]
      pert_info_cont[[i]] <- rr
    }
    pert_info_cont <- rbindlist(pert_info_cont)
  } else {
    by <- "Total"
    pert_info_cont <- data.table()
  }
  nV <- match(nV, names(dat))
  prob <- makeProblem(
    data=dat, dimList=dimList, dimVarInd=dV, freqVarInd=fV, numVarInd=nV, weightInd=NULL,
    sampWeightInd=NULL)
  tab <- sdcProb2df(prob, addDups=TRUE, addNumVars=TRUE, dimCodes="original")
  tab[,strID:=NULL]
  tab[,sdcStatus:=NULL]

  cols_reckey <- gen_vnames(countVars, prefix="sumRK")
  setnames(tab, c("tmprkeysfortabulation", countVars_rec), cols_reckey)

  cols_n <- gen_vnames(countVars, prefix="UWC")
  setnames(tab, c("freq", countVars), cols_n)

  cols_wc <- gen_vnames(countVars, prefix="WC")
  setnames(tab, c("tmpweightvarfortabulation", countVars_w), cols_wc)

  # compute average weights
  cols_wcavg <- gen_vnames(countVars, prefix="WCavg")
  for (i in seq_along(cols_n)) {
    wcavg <- tab[, get(cols_wc[i]) / get(cols_n[i])]
    wcavg[is.nan(wcavg)] <- 0
    set(tab, j=cols_wcavg[i], value=wcavg)
  }

  # compute cell-keys for each column independently
  cols_ck <- gen_vnames(countVars, prefix="cellKey")
  for (i in seq_along(cols_n)) {
    ck <- tab[, get(cols_reckey[i]) %% slot(pert_params, "bigN")]
    set(tab, j=cols_ck[i], value=ck)
  }

  # list with perturbation results for each countVar
  cols_pert <- gen_vnames(countVars, prefix="pert")
  count_modifications <- lapply(1:length(cols_n), function(x) {
    dt <- lookup(
      tab=tab,
      pert_params=pert_params,
      ckeyname=cols_ck[x],
      freqvarname=cols_n[x], type=type)
    set(tab, j=cols_pert[x], value=dt[,pert])
    dt[,countVar:=countVar_names[x]]
    setnames(dt, "cK", "cellKey")
    dt
  })
  count_modifications <- rbindlist(count_modifications)
  count_modifications <- cbind(tab[,c(names(dimList)), with=F], count_modifications)

  #tab <- fix_negative_counts(tab)

  # compute perturbed unweighted counts
  cols_puwc <- gen_vnames(countVars, prefix="pUWC")
  cols_pwc <- gen_vnames(countVars, prefix="pWC")
  for (i in seq_along(cols_n)) {
    # compute perturbed unweighted counts
    puwc <- tab[,get(cols_pert[i]) + get(cols_n[i])]
    puwc[is.na(puwc)] <- 0
    set(tab, j=cols_puwc[i], value=puwc)

    # compute perturbed weighted counts
    pwc <- tab[, round(get(cols_puwc[i])*get(cols_wcavg[i]))]
    pwc[is.na(pwc) | is.infinite(pwc)] <- 0
    set(tab, j=cols_pwc[i], value=pwc)
  }

  # compute measures for each numerical variable
  if (!is.null(numVars)) {
    # fix names
    # unweighted (orig): UW_{numVar}
    # unweighted (perturbed): pUW_{numVar}
    from <- c(paste0(numVars,"_tmporigforcalculation"), numVars)

    nvOrig <- paste0("UW_",numVars)
    nvPert <- paste0("pUW_",numVars)
    setnames(tab, from, c(nvOrig, nvPert))

    cn1 <- paste0("WS_", numVars)
    cn2 <- paste0("pWS_", numVars)

    by_wgt <- paste0("WCavg_", by)

    res1 <- tab[,lapply(.SD, function(x) {
      x*get(by_wgt)}), .SDcols=nvOrig]
    tab[,c(cn1):=res1]

    res2 <- tab[,lapply(.SD, function(x) {
      x*get(by_wgt)}), .SDcols=nvPert]
    tab[,c(cn2):=res2]

    keepNV <- c(nvOrig, nvPert, cn1, cn2)
  } else {
    keepNV <- numVars <- character()
  }

  # prepare output object
  is_weighted <- ifelse(is.null(weightVar), FALSE, TRUE)

  # keep variables
  vv <- c(names(dimList), cols_n, cols_wc, cols_puwc, cols_pwc, cols_wcavg, keepNV)
  tab <- tab[,vv, with=F]
  res <- new("pert_table",
    tab=tab,
    count_modifications=count_modifications,
    numvars_modifications=pert_info_cont,
    dimVars=names(dimList),
    countVars=countVar_names,
    numVars=numVars,
    by=by,
    is_weighted=is_weighted,
    type=type)
  validObject(res)
  return(res)
}

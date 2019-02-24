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
#' @examples
#' # loading testdata and adding record keys
#' dat <- ck_create_testdata()
#'
#' # create some binary (0,1) variables that should be perturbed later
#' dat[, cnt_females := ifelse(sex=="male", 0, 1)]
#' dat[, cnt_males := ifelse(sex=="male", 1, 0)]
#' dat[, cnt_highincome := ifelse(income>=9000, 1, 0)]
#'
#' ## create record keys
#' dat$rkeys <- ck_generate_rkeys(
#'   dat = dat,
#'   max_val = 2*nrow(dat),
#'   type = "abs",
#'   verbose = TRUE
#' )
#'
#' # creating test perturbation in "abs", "destatis" and "custom" format
#' pTable_abs <- ck_create_pTable(
#'   D = 5,
#'   V = 3,
#'   pTableSize = 70,
#'   type = "abs"
#' )
#' pTable_destatis <- ck_create_pTable(
#'   D = 5,
#'   V = 3,
#'   type = "destatis"
#' )
#'
#' # "free" format is a bit more complicated to set up
#' pTable_custom <- ck_create_custom_pTable(
#'   pTableSize = 70
#' )
#'
#' # create some functions first that return a single number
#' fn1 <- function() { round(rnorm(1, mean=5, sd=10)) }
#' fn2 <- function() { rpois(1, lambda=2) }
#' fn3 <- function() { rpois(1, lambda=5) }
#' fn4 <- function() return(1)
#' fn5 <- function() return(-1)
#'
#' # fn1 provides perturbation values from a normal distribution
#' # with mean = 5 and sd = 10; we use this for all cells
#' pTable_custom <- ck_update_custom_pTable(
#'   pTable = pTable_custom,
#'   fun = fn1
#' )
#'
#' # perturbation values from poisson distributions for some cells
#' pTable_custom <- ck_update_custom_pTable(
#'   pTable = pTable_custom,
#'   fun = fn2,
#'   cols = 1:20,
#'   rows = 1:100
#' )
#' pTable_custom <- ck_update_custom_pTable(
#'   pTable = pTable_custom,
#'   fun = fn3,
#'   cols = 1:20,
#'   rows = 101:256
#' )
#'
#' # we can of course write functions, that return scalars, such
#' # as fn4() and fn5() always return 1 or -1
#' pTable_custom <- ck_update_custom_pTable(
#'   pTable = pTable_custom,
#'   fun=fn4,
#'   cols=10:20
#' )
#' pTable_custom <- ck_update_custom_pTable(
#'   pTable = pTable_custom,
#'   fun=fn5,
#'   cols = 21:30
#' )
#'
#' # set up perturbation tables and parameters for numerical variables
#' sTable <- ck_generate_sTable(smallC = 12)
#' mTable <- c(0.6, 0.4, 0.2)
#'
#' # creating perturbation objects that are later used
#' # as inputs in perturbTable()
#' pert_params_abs <- ck_create_pert_params(
#'   bigN = 17312941,
#'   smallN = 12,
#'   pTable = pTable_abs,
#'   sTable = sTable,
#'   mTable = mTable
#' )
#'
#' pert_params_destatis <- ck_create_pert_params(
#'   smallN = 12,
#'   pTable = pTable_destatis,
#'   sTable = sTable,
#'   mTable = mTable
#' )
#'
#' pert_params_custom <- ck_create_pert_params(
#'   bigN = 17312941,
#'   smallN = 12,
#'   pTable = pTable_custom,
#'   sTable = sTable,
#'   mTable = mTable
#' )
#'
#' # create suitable input data using existing record keys
#' # and a perturbation table in "abs" format"
#' inp_abs <- ck_create_input(
#'   dat = dat,
#'   def_rkey = "rkeys",
#'   pert_params = pert_params_abs
#' )
#'
#' # create input using a perturbation table in "destatis" format"
#' # and no existing record keys. In this case, def_rkey = 9 means
#' # that the generated record keys are rounded to 9 digits
#' inp_destatis <- ck_create_input(
#'   dat = dat,
#'   def_rkey = 9,
#'   pert_params = pert_params_destatis
#' )
#'
#' # input using a perturbation table in "free" format
#' # consisting making use of custom functions
#' inp_custom <- ck_create_input(
#'   dat = dat,
#'   def_rkey = "rkeys",
#'   pert_params = pert_params_custom
#' )
#'
#' # define table hierarchies
#' # variable "sex"
#' dim.sex <- hier_create(root = "Total", nodes = c("male", "female"))
#' hier_display(dim.sex)
#'
#' # variable "sex"
#' dim.age <- hier_create(root = "Total", nodes = paste0("age_group", 1:6))
#' hier_display(dim.age)
#'
#' # define required inputs
#' dimList <- list(sex = dim.sex, age = dim.age)
#' weightVar <- "sampling_weight"
#' numVars <- c("savings", "income")
#'
#' # finally perturb the table
#' res <- perturbTable(
#'   inp = inp_abs,
#'   dimList = dimList,
#'   weightVar = weightVar,
#'   countVars = NULL,
#'   numVars = numVars
#' )
#'
#' # access the results
#' ck_freq_table(res, vname = "Total")
#' ck_cont_table(res, vname = "income", meanBeforeSum = TRUE)
#' ck_cont_table(res, vname = "savings", meanBeforeSum = TRUE)
#'
#' res <- perturbTable(
#'   inp = inp_destatis,
#'   dimList = dimList,
#'   weightVar = weightVar,
#'   countVars = NULL,
#'   numVars = numVars
#' )
#' ck_freq_table(res, vname = "Total")
#' ck_cont_table(res, vname = "income", meanBeforeSum = TRUE)
#' ck_cont_table(res, vname = "savings", meanBeforeSum = TRUE)
#'
#' res <- perturbTable(
#'   inp = inp_custom,
#'   dimList = dimList,
#'   weightVar = weightVar,
#'   countVars = NULL,
#'   numVars = numVars
#' )
#' ck_freq_table(res, vname = "Total")
#' ck_cont_table(res, vname = "income", meanBeforeSum = TRUE)
#' ck_cont_table(res, vname = "savings", meanBeforeSum = TRUE)
#'
#' # information on modifications for count variables
#' mod_counts(res, verbose = TRUE)
#'
#' # information on modifications for numerical variables
#' mod_numvars(res)
#'
#' # measures of information loss and utility statistics
#' # on tabulated count variables
#' ck_cnt_measures(res, vname = "Total")
#'
#' # an example using additional countVars
#' res <- perturbTable(
#'   inp = inp_destatis,
#'   dimList = dimList,
#'   weightVar = weightVar,
#'   countVars = c("cnt_females", "cnt_males", "cnt_highincome"),
#'   numVars = numVars)
#' print(res) # custom print method
#' summary(res) # custom summary method
#'
#' # show count tables
#' ck_freq_table(res, vname = NULL)
#' ck_freq_table(res, vname = "Total")
#' ck_freq_table(res, vname = "cnt_females")
#' ck_freq_table(res, vname = "cnt_males")
#' ck_freq_table(res, vname = "cnt_highincome")
#'
#' # show magnitude tables along with its modifications
#' p_inc <- ck_cont_table(res, vname = "income"); p_inc
#' attr(p_inc, "modifications")
#'
#' # create perturbed magnitude table for subgroups of data
#' res <- perturbTable(
#'   inp = inp_destatis,
#'   dimList = dimList,
#'   weightVar = weightVar,
#'   countVars = c("cnt_females", "cnt_males", "cnt_highincome"),
#'   numVars = numVars,
#'   by = "cnt_males"
#' )
#'
#' # the perturbed tables of savings and income are computed only where
#' # variable cnt_males is 1
#' p_sav <- ck_cont_table(res, vname = "savings"); p_sav
#' attr(p_sav, "modifications") # no modifications in cells containing females!
#'
#' # export table to a simple data.table
#' df_tot <- ck_export_table(res, vname = "Total", type = "both")
#' head(df_tot)
#' df_inc <- ck_export_table(res, vname = "income", type = "weighted")
#' head(df_inc)
perturbTable <-
  function(inp,
           dimList,
           countVars = NULL,
           numVars = NULL,
           by = NULL,
           weightVar = NULL) {
  # rename variables
  gen_vnames <- function(countVars, prefix="sumRK") {
    if (length(countVars) == 0) {
      return(paste0(prefix, "_Total"))
    }
    c(paste0(prefix, "_Total"), paste0(prefix, "_", countVars))
  }

  sdcStatus <- strID <- NULL
  tmprkeysfortabulation <- tmpweightvarfortabulation <- tmpidforsorting <- NULL
  pert <- countVar <- numVar <- NULL

  stopifnot(isS4(inp))
  stopifnot(class(inp) == "pert_inputdat")
  pert_params <- slot(inp, "pert_params")
  type <- slot(pert_params, "type")

  if (!is.null(numVars)) {
    nr <- nrow(slot(pert_params, "sTable"))
    if (nr == 0 || is.null(slot(pert_params, "mTable"))) {
      e <- c(
        "Perturbation of magnitude tables is not posible:\n",
        "`sTable` or `mTable` were not specified"
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }
  }

  dat  <- copy(slot(inp, "microdat"))
  dat[, tmprkeysfortabulation := slot(inp, "rkeys")]
  dat[, tmpidforsorting := .I]

  dV <- match(names(dimList), names(dat))
  if (any(is.na(dV))) {
    e <- c(
      "Check input `dimList`:",
      "--> at least one variable was not found in the input data."
    )
    stop(paste(e, collapse = " "), call. = FALSE)
  }

  freqVar <- "tmpfreqvarfortabulation"
  dat[, c(freqVar) := 1]

  dat[, tmpweightvarfortabulation := 1]
  if (!is.null(weightVar)) {
    stopifnot(weightVar %in% names(dat))
    dat[, tmpweightvarfortabulation := get(weightVar)]
  }

  # check if we have subgroups
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
      e <- c(
        "Check input `countVars`:",
        "at least one variable was not found in the input dataset."
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }
    # check if all count-vars are integerish 0/1
    chk <- sapply(countVars, function(x) {
      all(unique(dat[[x]]) %in% 0:1)
    })
    if (!all(chk)) {
      e <- c(
        "Check input `countVars`:",
        "--> at least one variable is not 0/1 coded."
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }

    # compute weighted version!
    countVars_w <- paste0("wc_", countVars)
    dat[, (countVars_w) := lapply(
      .SD, function(x) {
        x * tmpweightvarfortabulation
      }
    ), .SDcols = countVars]

    # record_keys using an indicator
    countVars_rec <- paste0("recKey_", countVars)
    dat[, (countVars_rec) := lapply(
      .SD, function(x) {
        x * tmprkeysfortabulation
      }
    ), .SDcols = countVars]
  }

  countVar_names <- c("Total", countVars)

  fV <- match(freqVar, names(dat))
  nV <- c(
    "tmprkeysfortabulation",
    "tmpweightvarfortabulation",
    countVars,
    countVars_w,
    countVars_rec
  )
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
      if (by == "Total") {
        vals_orig <- dat[, get(v)]
        set(dat, j = vnum_orig, value = vals_orig)
      } else {
        # we need to set original numerical variable to 0
        # where by variable is also 0!
        dat[, (numVars[i]) := get(numVars[i]) * get(by)]
        vals_orig <- dat[, get(v) * get(by)]
        set(dat, j = vnum_orig, value = vals_orig)
      }
      rr <- identify_topK_cells(
        dat = dat,
        rkeys = slot(inp, "rkeys"),
        dimList = dimList,
        pert_params = pert_params,
        v = v,
        type = type
      )
      dat[rr[, tmpidforsorting], c(v) := rr[, get(paste0(v, ".mod"))]]
      setnames(rr, "tmpidforsorting", "id")

      n_old <- c(v, paste0(v, c(".pert", ".mod")))
      n_new <- c(paste0("vals.", c("orig", "pert", "mod")))
      setnames(rr, n_old, n_new)
      rr[, numVar := v]
      pert_info_cont[[i]] <- rr
    }
    pert_info_cont <- rbindlist(pert_info_cont)
  } else {
    by <- "Total"
    pert_info_cont <- data.table()
  }
  nV <- match(nV, names(dat))
  prob <- makeProblem(
    data = dat,
    dimList = dimList,
    dimVarInd = dV,
    freqVarInd = fV,
    numVarInd = nV,
    weightInd = NULL,
    sampWeightInd = NULL
  )
  tab <- sdcProb2df(
    obj = prob,
    addDups = TRUE,
    addNumVars = TRUE,
    dimCodes = "original"
  )
  tab[, strID := NULL]
  tab[, sdcStatus := NULL]

  cols_reckey <- gen_vnames(countVars, prefix = "sumRK")
  setnames(tab, c("tmprkeysfortabulation", countVars_rec), cols_reckey)

  cols_n <- gen_vnames(countVars, prefix = "UWC")
  setnames(tab, c("freq", countVars), cols_n)

  cols_wc <- gen_vnames(countVars, prefix = "WC")
  setnames(tab, c("tmpweightvarfortabulation", countVars_w), cols_wc)

  # compute average weights
  cols_wcavg <- gen_vnames(countVars, prefix = "WCavg")
  for (i in seq_along(cols_n)) {
    wcavg <- tab[, get(cols_wc[i]) / get(cols_n[i])]
    wcavg[is.nan(wcavg)] <- 0
    set(tab, j = cols_wcavg[i], value = wcavg)
  }

  # compute cell-keys for each column independently
  cols_ck <- gen_vnames(countVars, prefix = "cellKey")
  for (i in seq_along(cols_n)) {
    ck <- tab[, get(cols_reckey[i]) %% slot(pert_params, "bigN")]
    set(tab, j = cols_ck[i], value = ck)
  }

  # list with perturbation results for each countVar
  cols_pert <- gen_vnames(countVars, prefix = "pert")
  count_modifications <- lapply(1:length(cols_n), function(x) {
    dt <- lookup(
      tab = tab,
      pert_params = pert_params,
      ckeyname = cols_ck[x],
      freqvarname = cols_n[x],
      type = type
    )
    set(tab, j = cols_pert[x], value = dt[, pert])
    dt[, countVar := countVar_names[x]]
    setnames(dt, "cK", "cellKey")
    dt
  })
  count_modifications <- rbindlist(count_modifications)
  count_modifications <-
    cbind(tab[, c(names(dimList)), with = F], count_modifications)

  # compute perturbed unweighted counts
  cols_puwc <- gen_vnames(countVars, prefix = "pUWC")
  cols_pwc <- gen_vnames(countVars, prefix = "pWC")
  for (i in seq_along(cols_n)) {
    # compute perturbed unweighted counts
    puwc <- tab[, get(cols_pert[i]) + get(cols_n[i])]
    puwc[is.na(puwc)] <- 0
    set(tab, j = cols_puwc[i], value = puwc)

    # compute perturbed weighted counts
    pwc <- tab[, round(get(cols_puwc[i]) * get(cols_wcavg[i]))]
    pwc[is.na(pwc) | is.infinite(pwc)] <- 0
    set(tab, j = cols_pwc[i], value = pwc)
  }

  # compute measures for each numerical variable
  if (!is.null(numVars)) {
    from <- c(paste0(numVars, "_tmporigforcalculation"), numVars)

    nvOrig <- paste0("UW_", numVars)
    nvPert <- paste0("pUW_", numVars)
    setnames(tab, from, c(nvOrig, nvPert))

    cn1 <- paste0("WS_", numVars)
    cn2 <- paste0("pWS_", numVars)

    by_wgt <- paste0("WCavg_", by)
    res1 <- tab[, lapply(
      .SD, function(x) {
        x * get(by_wgt)
      }
    ), .SDcols = nvOrig]
    tab[, c(cn1) := res1]

    res2 <- tab[, lapply(
      .SD, function(x) {
        x * get(by_wgt)
      }
    ), .SDcols = nvPert]
    tab[, c(cn2) := res2]

    keepNV <- c(nvOrig, nvPert, cn1, cn2)
  } else {
    keepNV <- numVars <- character()
  }

  # prepare output object
  is_weighted <- ifelse(is.null(weightVar), FALSE, TRUE)

  # keep variables
  vv <- c(
    names(dimList),
    cols_n,
    cols_wc,
    cols_puwc,
    cols_pwc,
    cols_wcavg,
    keepNV
  )
  tab <- tab[, vv, with = FALSE]
  res <- new(
    Class = "pert_table",
    tab = tab,
    count_modifications = count_modifications,
    numvars_modifications = pert_info_cont,
    dimVars = names(dimList),
    countVars = countVar_names,
    numVars = numVars,
    by = by,
    is_weighted = is_weighted,
    type = type
  )
  validObject(res)
  return(res)
}

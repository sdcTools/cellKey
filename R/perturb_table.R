#' Perturb a multi-dimensional table using the cell-key method
#'
#' This function allows to compute perturbed tables containing counts as
#' well as continous variables.
#'
#' @param inp an object of class \code{\link{pert_inputdat-class}}
#' generated with \code{\link{ck_create_input}}
#' @param dim_list a list containing slots for each variable that should be
#' tabulated. Each slot consists of a \code{data.frame} or \code{data.table}
#' with columns \code{levels} and \code{codes} specifiying the
#' hierarchies of the variables.
#' @param countvars (character) vector of numerical variables that should be
#' counted and perturbed; must be 0/1 coded.
#' @param numvars (character) vector of numerical variables that
#' should be tabulated or \code{NULL}
#' @param by \code{NULL} or a scalar character. If specified, it is possible
#' to use a valid variable name being 0/1 coded that will be used to create the
#' perturbation of all variables specified in \code{numvars} only for subgroups
#' defined by the specified count variable. If the variable is not listed in
#' argument \code{countvars} it will be automatically added. If not specified or
#' \code{NULL}, the table will be based on all units in the microdata.
#' @param weightvar (character) vector of variable holding sampling
#' weights or \code{NULL}
#' @return an object of class \code{\link{pert_table-class}}.
#' @seealso \url{https://preview.tinyurl.com/yyleyg48}
#' @export
#' @examples
#' # loading testdata and adding record keys
#' dat <- ck_create_testdata()
#'
#' # create some binary (0,1) variables that should be perturbed later
#' dat[, cnt_females := ifelse(sex == "male", 0, 1)]
#' dat[, cnt_males := ifelse(sex == "male", 1, 0)]
#' dat[, cnt_highincome := ifelse(income >= 9000, 1, 0)]
#'
#' ## create record keys
#' dat$rkeys <- ck_generate_rkeys(
#'   dat = dat,
#'   max_val = 2*nrow(dat),
#'   type = "abs",
#'   verbose = TRUE
#' )
#'
#' # creating test perturbation in "destatis", "abs" and "free" format
#' ptab_abs <- ck_create_ptab(
#'   D = 5,
#'   V = 3,
#'   ptab_size = 70,
#'   type = "abs"
#' )
#' ptab_destatis <- ck_create_ptab(
#'   D = 5,
#'   V = 3,
#'   type = "destatis"
#' )
#'
#' # "free" format is a bit more complicated to set up
#' ptab_free <- ck_create_free_ptable(
#'   nrcol = 70
#' )
#'
#' # create some functions first that return a single number
#' fn1 <- function() round(rnorm(1, mean = 5, sd = 10))
#' fn2 <- function() rpois(1, lambda = 2)
#' fn3 <- function() rpois(1, lambda = 5)
#' fn4 <- function() return(1)
#' fn5 <- function() return(-1)
#'
#' # fn1 provides perturbation values from a normal distribution
#' # with mean = 5 and sd = 10; we use this for all cells
#' ptab_free <- ck_update_free_ptable(
#'   ptab = ptab_free,
#'   fun = fn1
#' )
#'
#' # perturbation values from poisson distributions for some cells
#' ptab_free <- ck_update_free_ptable(
#'   ptab = ptab_free,
#'   fun = fn2,
#'   cols = 1:20,
#'   rows = 1:100
#' )
#' ptab_free <- ck_update_free_ptable(
#'   ptab = ptab_free,
#'   fun = fn3,
#'   cols = 1:20,
#'   rows = 101:256
#' )
#'
#' # we can of course write functions, that return scalars, such
#' # as fn4() and fn5() always return 1 or -1
#' ptab_free <- ck_update_free_ptable(
#'   ptab = ptab_free,
#'   fun = fn4,
#'   cols = 10:20
#' )
#' ptab_free <- ck_update_free_ptable(
#'   ptab = ptab_free,
#'   fun = fn5,
#'   cols = 21:30
#' )
#'
#' # set up perturbation tables and parameters for numerical variables
#' stab <- ck_generate_stab(small_c = 12)
#' mtab <- c(0.6, 0.4, 0.2)
#'
#' # creating perturbation objects that are later used
#' # as inputs in perturb_table()
#' pert_params_abs <- ck_create_pert_params(
#'   big_n = 17312941,
#'   small_n = 12,
#'   ptab = ptab_abs,
#'   stab = stab,
#'   mtab = mtab
#' )
#'
#' pert_params_destatis <- ck_create_pert_params(
#'   small_n = 12,
#'   ptab = ptab_destatis,
#'   stab = stab,
#'   mtab = mtab
#' )
#'
#' pert_params_free <- ck_create_pert_params(
#'   big_n = 17312941,
#'   small_n = 12,
#'   ptab = ptab_free,
#'   stab = stab,
#'   mtab = mtab
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
#' inp_free <- ck_create_input(
#'   dat = dat,
#'   def_rkey = "rkeys",
#'   pert_params = pert_params_free
#' )
#'
#' # define table hierarchies
#' # variable "sex"
#' dim_sex <- hier_create(root = "Total", nodes = c("male", "female"))
#' hier_display(dim_sex)
#'
#' # variable "age"
#' dim_age <- hier_create(root = "Total", nodes = paste0("age_group", 1:6))
#' hier_display(dim_age)
#'
#' # define required inputs
#' dim_list <- list(sex = dim_sex, age = dim_age)
#' weightvar <- "sampling_weight"
#' numvars <- c("savings", "income")
#'
#' # finally perturb the table
#' res <- perturb_table(
#'   inp = inp_abs,
#'   dim_list = dim_list,
#'   weightvar = weightvar,
#'   countvars = NULL,
#'   numvars = numvars
#' )
#'
#' # access the results
#' ck_freq_table(res, vname = "total")
#' ck_cont_table(res, vname = "income", mean_before_sum = TRUE)
#' ck_cont_table(res, vname = "savings", mean_before_sum = TRUE)
#'
#' res <- perturb_table(
#'   inp = inp_destatis,
#'   dim_list = dim_list,
#'   weightvar = weightvar,
#'   countvars = NULL,
#'   numvars = numvars
#' )
#' ck_freq_table(res, vname = "total")
#' ck_cont_table(res, vname = "income", mean_before_sum = TRUE)
#' ck_cont_table(res, vname = "savings", mean_before_sum = TRUE)
#'
#' res <- perturb_table(
#'   inp = inp_free,
#'   dim_list = dim_list,
#'   weightvar = weightvar,
#'   countvars = NULL,
#'   numvars = numvars
#' )
#' ck_freq_table(res, vname = "total")
#' ck_cont_table(res, vname = "income", mean_before_sum = TRUE)
#' ck_cont_table(res, vname = "savings", mean_before_sum = TRUE)
#'
#' # information on modifications for count variables
#' mod_counts(res, verbose = TRUE)
#'
#' # information on modifications for numerical variables
#' mod_numvars(res)
#'
#' # measures of information loss and utility statistics
#' # on tabulated count variables
#' ck_cnt_measures(res, vname = "total")
#'
#' # an example using additional count variables
#' res <- perturb_table(
#'   inp = inp_destatis,
#'   dim_list = dim_list,
#'   weightvar = weightvar,
#'   countvars = c("cnt_females", "cnt_males", "cnt_highincome"),
#'   numvars = numvars
#' )
#' print(res) # custom print method
#' summary(res) # custom summary method
#'
#' # show count tables
#' ck_freq_table(res, vname = NULL)
#' ck_freq_table(res, vname = "total")
#' ck_freq_table(res, vname = "cnt_females")
#' ck_freq_table(res, vname = "cnt_males")
#' ck_freq_table(res, vname = "cnt_highincome")
#'
#' # show magnitude tables along with its modifications
#' p_inc <- ck_cont_table(res, vname = "income"); p_inc
#' attr(p_inc, "modifications")
#'
#' # create perturbed magnitude table for subgroups of data
#' res <- perturb_table(
#'   inp = inp_destatis,
#'   dim_list = dim_list,
#'   weightvar = weightvar,
#'   countvars = c("cnt_females", "cnt_males", "cnt_highincome"),
#'   numvars = numvars,
#'   by = "cnt_males"
#' )
#'
#' # the perturbed tables of savings and income are computed only where
#' # variable cnt_males is 1
#' p_sav <- ck_cont_table(res, vname = "savings"); p_sav
#' attr(p_sav, "modifications") # no modifications in cells containing females!
#'
#' # export table to a simple data.table
#' df_tot <- ck_export_table(res, vname = "total", type = "both")
#' head(df_tot)
#' df_inc <- ck_export_table(res, vname = "income", type = "weighted")
#' head(df_inc)
perturb_table <-
  function(inp,
           dim_list,
           countvars = NULL,
           numvars = NULL,
           by = NULL,
           weightvar = NULL) {
  # rename variables
  .vnames <- function(countvars, prefix="sum_rk") {
    if (length(countvars) == 0) {
      return(paste0(prefix, "_total"))
    }
    c(paste0(prefix, "_total"), paste0(prefix, "_", countvars))
  }

  sdcStatus <- strID <- NULL
  tmprkeysfortabulation <- tmpweightvarfortabulation <- tmpidforsorting <- NULL
  pert <- countvar <- numvar <- NULL

  stopifnot(isS4(inp))
  stopifnot(class(inp) == "pert_inputdat")
  pert_params <- slot(inp, "pert_params")
  type <- slot(pert_params, "type")

  if (!is.null(numvars)) {
    nr <- nrow(slot(pert_params, "stab"))
    if (nr == 0 || is.null(slot(pert_params, "mtab"))) {
      e <- c(
        "Perturbation of magnitude tables is not posible:\n",
        "`stab` or `mtab` were not specified"
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }
  }

  dat  <- copy(slot(inp, "microdat"))
  dat[, tmprkeysfortabulation := slot(inp, "rkeys")]
  dat[, tmpidforsorting := .I]

  dV <- match(names(dim_list), names(dat))
  if (any(is.na(dV))) {
    e <- c(
      "Check input `dim_list`:",
      "--> at least one variable was not found in the input data."
    )
    stop(paste(e, collapse = " "), call. = FALSE)
  }

  freq_var <- "tmpfreqvarfortabulation"
  dat[, c(freq_var) := 1]

  dat[, tmpweightvarfortabulation := 1]
  if (!is.null(weightvar)) {
    stopifnot(weightvar %in% names(dat))
    dat[, tmpweightvarfortabulation := get(weightvar)]
  }

  # check if we have subgroups
  if (!is.null(by)) {
    stopifnot(is_scalar_character(by))
    if (!by %in% countvars) {
      countvars <- c(countvars, by)
    }
  } else {
    by <- "total"
  }

  # check and prepare count variables
  if (is.null(countvars)) {
    countvars <- countvars_w <- countvars_rec <- c()
  } else {
    # checking if all count_vars exist!
    cv <- match(countvars, names(dat))
    if (any(is.na(cv))) {
      e <- c(
        "Check input `countvars`:",
        "at least one variable was not found in the input dataset."
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }
    # check if all count-vars are integerish 0/1
    chk <- sapply(countvars, function(x) {
      all(unique(dat[[x]]) %in% 0:1)
    })
    if (!all(chk)) {
      e <- c(
        "Check input `countvars`:",
        "--> at least one variable is not 0/1 coded."
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }

    # compute weighted version!
    countvars_w <- paste0("wc_", countvars)
    dat[, (countvars_w) := lapply(
      .SD, function(x) {
        x * tmpweightvarfortabulation
      }
    ), .SDcols = countvars]

    # record_keys using an indicator
    countvars_rec <- paste0("reckey_", countvars)
    dat[, (countvars_rec) := lapply(
      .SD, function(x) {
        x * tmprkeysfortabulation
      }
    ), .SDcols = countvars]
  }

  countvar_names <- c("total", countvars)

  fV <- match(freq_var, names(dat))
  nV <- c(
    "tmprkeysfortabulation",
    "tmpweightvarfortabulation",
    countvars,
    countvars_w,
    countvars_rec
  )
  pert_info_cont <- list()
  if (!is.null(numvars)) {
    nv2 <- match(numvars, names(dat))
    if (any(is.na(nv2))) {
      stop("check input `numvars`\n", call. = FALSE)
    }
    nV <- c(nV, numvars)
    length(pert_info_cont) <- length(numvars)
    names(pert_info_cont) <- numvars
    for (i in 1:length(numvars)) {
      v <- numvars[i]
      vnum_orig <- paste0(v, "_tmporigforcalculation")
      nV <- c(nV, vnum_orig)

      # copy numerical variable
      if (by == "total") {
        vals_orig <- dat[, get(v)]
        set(dat, j = vnum_orig, value = vals_orig)
      } else {
        # we need to set original numerical variable to 0
        # where by variable is also 0!
        dat[, (numvars[i]) := get(numvars[i]) * get(by)]
        vals_orig <- dat[, get(v) * get(by)]
        set(dat, j = vnum_orig, value = vals_orig)
      }
      rr <- .identify_topk_cells(
        dat = dat,
        rkeys = slot(inp, "rkeys"),
        dim_list = dim_list,
        pert_params = pert_params,
        v = v,
        type = type
      )

      dat[rr[, tmpidforsorting], c(v) := rr[, get(paste0(v, "_mod"))]]
      setnames(rr, "tmpidforsorting", "id")

      n_old <- c(v, paste0(v, c("_pert", "_mod")))
      n_new <- c(paste0("vals_", c("orig", "pert", "mod")))
      setnames(rr, n_old, n_new)
      rr[, numvar := v]
      pert_info_cont[[i]] <- rr
    }
    pert_info_cont <- rbindlist(pert_info_cont)
  } else {
    by <- "total"
    pert_info_cont <- data.table()
  }
  nV <- match(nV, names(dat))
  prob <- makeProblem(
    data = dat,
    dimList = dim_list,
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

  cols_reckey <- .vnames(countvars, prefix = "sumrk")
  setnames(tab, c("tmprkeysfortabulation", countvars_rec), cols_reckey)

  cols_n <- .vnames(countvars, prefix = "uwc")
  setnames(tab, c("freq", countvars), cols_n)

  cols_wc <- .vnames(countvars, prefix = "wc")
  setnames(tab, c("tmpweightvarfortabulation", countvars_w), cols_wc)

  # compute average weights
  cols_wcavg <- .vnames(countvars, prefix = "wc_avg")
  for (i in seq_along(cols_n)) {
    wcavg <- tab[, get(cols_wc[i]) / get(cols_n[i])]
    wcavg[is.nan(wcavg)] <- 0
    set(tab, j = cols_wcavg[i], value = wcavg)
  }

  # compute cell-keys for each column independently
  cols_ck <- .vnames(countvars, prefix = "cellkey")
  for (i in seq_along(cols_n)) {
    ck <- tab[, get(cols_reckey[i]) %% slot(pert_params, "big_n")]
    set(tab, j = cols_ck[i], value = ck)
  }

  # list with perturbation results for each countVar
  cols_pert <- .vnames(countvars, prefix = "pert")
  count_modifications <- lapply(1:length(cols_n), function(x) {
    dt <- .lookup(
      tab = tab,
      pert_params = pert_params,
      ckeyname = cols_ck[x],
      freqvarname = cols_n[x],
      type = type
    )
    set(tab, j = cols_pert[x], value = dt[, pert])
    dt[, countvar := countvar_names[x]]
    setnames(dt, "ck", "cellkey")
    dt
  })
  count_modifications <- rbindlist(count_modifications)
  count_modifications <- cbind(
    tab[, c(names(dim_list)), with = FALSE],
    count_modifications
  )

  # compute perturbed unweighted counts
  cols_puwc <- .vnames(countvars, prefix = "puwc")
  cols_pwc <- .vnames(countvars, prefix = "pwc")
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
  if (!is.null(numvars)) {
    from <- c(paste0(numvars, "_tmporigforcalculation"), numvars)

    nv_orig <- paste0("uw_", numvars)
    nv_pert <- paste0("puw_", numvars)
    setnames(tab, from, c(nv_orig, nv_pert))

    cn1 <- paste0("ws_", numvars)
    cn2 <- paste0("pws_", numvars)

    by_wgt <- paste0("wc_avg_", by)
    res1 <- tab[, lapply(
      .SD, function(x) {
        x * get(by_wgt)
      }
    ), .SDcols = nv_orig]
    tab[, c(cn1) := res1]

    res2 <- tab[, lapply(
      .SD, function(x) {
        x * get(by_wgt)
      }
    ), .SDcols = nv_pert]
    tab[, c(cn2) := res2]

    keep_nv <- c(nv_orig, nv_pert, cn1, cn2)
  } else {
    keep_nv <- numvars <- character()
  }

  # prepare output object
  is_weighted <- ifelse(is.null(weightvar), FALSE, TRUE)

  # keep variables
  vv <- c(
    names(dim_list),
    cols_n,
    cols_wc,
    cols_puwc,
    cols_pwc,
    cols_wcavg,
    keep_nv
  )
  tab <- tab[, vv, with = FALSE]
  res <- new("pert_table",
    tab = tab,
    count_modifications = count_modifications,
    numvars_modifications = pert_info_cont,
    dimvars = names(dim_list),
    countvars = countvar_names,
    numvars = numvars,
    by = by,
    is_weighted = is_weighted,
    type = type
  )
  validObject(res)
  return(res)
}

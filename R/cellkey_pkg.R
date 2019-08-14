cellkey_obj_class <- R6::R6Class("cellkey_obj", cloneable = FALSE,
  public=list(
    #initialize=function(x, rkey, dims, w, countvars = NULL, numvars = NULL, params_cnts, params_nums = NULL) {
    initialize=function(x, rkey, dims, w, countvars = NULL, numvars = NULL) {
      type <- is_perturbed <- NULL

      if (!inherits(x, "data.frame")) {
        stop("an object coercible to a `data.frame` must be provided in `x`.", call. = FALSE)
      }

      x <- as.data.table(x)

      # checking record keys
      rkeyvar <- .tmpvarname("rkey")
      if (is.character(rkey)) {
        stopifnot(is_scalar_character(rkey))
        stopifnot(rkey %in% names(x))

        # check the record keys
        rkeys <- x[[rkey]]
        check_rkeys(rkeys = rkeys)
        x[[rkeyvar]] <- rkeys
      } else if (is_scalar_integerish(rkey)) {
        if (rkey < 5) {
          stop("number of digits for record keys must be at least `5`!", call. = FALSE)
        }
        x[[rkeyvar]] <- ck_generate_rkeys(x, nr_digits = rkey)
      } else {
        e <- "Argument `rkey` must either be character or integer(ish) number > 0!"
        stop(e, call. = FALSE)
      }

      # checking dimensions
      if (!rlang::is_named(dims)) {
        stop("Please provide names to argument `dims`", call. = FALSE)
      }
      if (!is.list(dims)) {
        stop("Argument `dims` must be a list.", call. = FALSE)
      }
      dimvars <- names(dims)
      if (!all(dimvars %in% names(x))) {
        stop("Some variable names provided in `dims` don't exist in `x`.", call. = FALSE)
      }
      chk <- sapply(dims, function(x) {
        !inherits(x, "sdc_hierarchy")
      })
      if (any(chk)) {
        stop("Some elements provided in `dims` are not valid hierarchy objects.", call. = FALSE)
      }

      # checking weights
      wvar <- .tmpvarname("weightvar")
      if (!is.null(w)) {
        if (!w %in% names(x)) {
          stop("weight variable `w` does not exist in `x`.", call. = FALSE)
        }
        if (!is.numeric(x[[w]])) {
          stop("weight variable `w` must be numeric.", call. = FALSE)
        }
        x[[wvar]] <- x[[w]]
        is_weighted <- TRUE
      } else {
        x[[wvar]] <- 1
        is_weighted <- FALSE
      }

      # adding a frequency variable
      freqvar <- .tmpvarname("freq")
      x[[freqvar]] <- 1

      if (!is.null(countvars)) {
        # check countvars
        if (!is.character(countvars)) {
          stop("Argument `countvars` needs to be character specifying variable names.", call. = FALSE)
        }
        if (!all(countvars %in% names(x))) {
          stop("Some variable names provided in `countvars` don't exist in `x`.", call. = FALSE)
        }
        chk <- sapply(numvars, function(nv) {
          !is.numeric(x[[nv]])
        })
        if (any(chk)) {
          stop("Some elements provided in `countvars` are not numeric.", call. = FALSE)
        }
        # check if all count-vars are integerish 0/1
        chk <- sapply(countvars, function(cv) {
          all(unique(x[[cv]]) %in% c(0, 1))
        })
        if (!all(chk)) {
          stop("Some elements provided in `countvars` are not 0/1 coded.", call. = FALSE)
        }

        # compute weighted countvars
        countvars_w <- paste0("wc_", countvars)
        suppressWarnings(
          x[, (countvars_w) := lapply(.SD, function(v) v * get(wvar)), .SDcols = countvars]
        )

        # record_keys using an indicator
        countvars_rec <- paste0("rkey_", countvars)
        x[, (countvars_rec) := lapply(.SD, function(v) v * get(rkeyvar)), .SDcols = countvars]
      } else {
        countvars <- countvars_w <- countvars_rec <- c()
      }

      if (!is.null(numvars)) {
        # check numvars
        if (!is.character(numvars)) {
          stop("Argument `numvars` needs to be character specifying variable names.", call. = FALSE)
        }
        if (!all(numvars %in% names(x))) {
          stop("Some variable names provided in `numvars` don't exist in `x`.", call. = FALSE)
        }
        chk <- sapply(numvars, function(nv) {
          !is.numeric(x[[nv]])
        })
        if (any(chk)) {
          stop("Some elements provided in `numvars` are not numeric.", call. = FALSE)
        }

        # compute weighted variables and index to compute record keys for only
        # units actually contributing
        numvars_w <- paste0("ws_", numvars)
        nozero_reckeys_nums <- paste0("rkey_nozero_", numvars)
        numvars_ind <- paste0("ind_nonzero", numvars)
        for (i in seq_along(numvars)) {
          v <- numvars[i]
          x[[numvars_w[i]]] <- x[[wvar]] * x[[v]]

          # compute index
          index <- rep(1, nrow(x))
          index[is.na(x[[v]])] <- 0
          index[is.nan(x[[v]])] <- 0
          index[x[[v]] == 0] <- 0

          # record keys is 0 for non-contributing units
          x[[nozero_reckeys_nums[i]]] <- x[[rkeyvar]] * index
        }
      } else {
        max_contributions <- numvars_w <- nozero_reckeys_nums <- NULL
      }

      nv <- c(
        rkeyvar, wvar, countvars, countvars_w, countvars_rec,
        numvars, numvars_w, nozero_reckeys_nums
      )

      prob <- makeProblem(
        data = x,
        dimList = dims,
        freqVarInd = match(freqvar, names(x)),
        numVarInd = match(nv, names(x)),
        weightInd = NULL,
        sampWeightInd = NULL
      )

      # with perturbations
      tab <- sdcProb2df(prob, addDups = TRUE, addNumVars = TRUE, dimCodes = "original")

      # duplicated cells (bogus codes) have the same
      # strID as their parent cell!
      strids <- unique(tab$strID)

      # calculate contributing indices
      # we do this in any case!
      contr_indices <- lapply(1:length(strids), function(x) {
        sdcTable:::c_contributing_indices(
          object = prob,
          input = list(strids[x]))
      })
      names(contr_indices) <- strids

      # finding top_k contributors for each cell and numerical variable
      microdat <- prob@dataObj@rawData[, c(names(dims), numvars, wvar), with = FALSE]
      microdat$.tmpid <- 1:nrow(microdat)

      # perhaps c++?
      # for each numerical variable and each cell, get the top_k contributions
      # along with its ids and its values
      .get_max_contributions <- function(indices, microdat, wvar, nv, top_k) {
        res <- vector("list", length = length(indices))
        names(res) <- names(indices)

        for (i in 1:length(res)) {
          out <- vector("list", length = length(nv))
          names(out) <- nv
          xx <- microdat[.tmpid %in% indices[[i]]]
          top_k <- min(top_k, nrow(xx))
          for (v in nv) {
            setorderv(xx, v, order = -1L)

            # unweighted
            out[[v]]$uw_ids <- xx$.tmpid[1:top_k]
            out[[v]]$uw_vals <- xx[[v]][1:top_k]
            out[[v]]$uw_spread <- diff(range(xx[[v]], na.rm = TRUE))
            out[[v]]$uw_sum <- sum(xx[[v]], na.rm = TRUE)
            out[[v]]$uw_mean <- out[[v]]$uw_sum / nrow(xx)

            # weighted
            xx$xw <- xx[[v]] * xx[[wvar]]
            setorderv(xx, "xw", order = -1L)

            out[[v]]$w_ids <- xx$.tmpid[1:top_k]
            out[[v]]$w_vals <- xx$xw[1:top_k]
            out[[v]]$w_cell
            out[[v]]$w_spread <- diff(range(xx[[v]], na.rm = TRUE))
            out[[v]]$w_sum <- sum(xx$xw, na.rm = TRUE)
            out[[v]]$w_mean <- out[[v]]$w_sum / sum(xx[[wvar]], na.rm = TRUE)
          }
          res[[i]] <- out
        }
        res
      }

      # top_k is hardcoded to 6;
      # this is the maximum allowed value for top_k, also in params_nums()
      max_contributions <- .get_max_contributions(
        indices = contr_indices,
        microdat = microdat,
        nv = numvars,
        wvar = wvar,
        top_k = 6
      )

      res <- tab[, c("strID", "freq", dimvars, nv), with = FALSE]
      setnames(res, nv, tolower(nv))

      cols_ck <- gen_vnames(c("total", countvars), prefix = "ck")

      if (length(countvars) > 0) {
        setnames(res, c(rkeyvar, paste0("rkey_", countvars)), cols_ck)
      } else {
        setnames(res, rkeyvar, cols_ck)
      }

      # fix cell keys
      for (i in 1:length(cols_ck)) {
        res[[cols_ck[i]]] <- res[[cols_ck[i]]] %% 1
      }

      cols_n <- gen_vnames(c("total", countvars), prefix = "uwc")
      setnames(res, c("freq", countvars), cols_n)

      cols_wc <- gen_vnames(c("total", countvars), prefix = "wc")
      if (length(countvars) > 0) {
        setnames(res, c(wvar, paste0("wc_", countvars)), cols_wc)
      } else {
        setnames(res, wvar, cols_wc)
      }

      cv <- c("total", countvars)
      rescnts <- vector("list", length = length(cv) + 1)
      names(rescnts) <- c("dims", cv)
      rescnts[[1]] <- res[, c("strID", dimvars), with = FALSE]

      for (i in seq_along(cv)) {
        index <- i + 1
        vv <- c(
          gen_vnames(cv[i], prefix = "ck"),
          gen_vnames(cv[i], prefix = "uwc"),
          gen_vnames(cv[i], prefix = "wc")
        )
        rescnts[[index]] <- res[, vv, with = FALSE]
      }

      if (length(numvars) > 0) {
        resnums <- vector("list", length = length(numvars))
        names(resnums) <- numvars

        cols_uws <- gen_vnames(numvars, prefix = "uws")
        setnames(res, numvars, cols_uws)
        cols_ws <- gen_vnames(numvars, prefix = "ws")
        for (j in 1:length(numvars)) {
          cn1 <- gen_vnames("total", prefix = "ck")
          cn2 <- gen_vnames(numvars[j], prefix = "rkey_nozero")
          vv <- c(cn1, cn2, cols_uws[j], cols_ws[j])
          vv <- res[, vv, with = FALSE]
          vv[[cn2]] <- vv[[cn2]] %% 1 # normalize cellkey

          setnames(vv, cn1, gen_vnames(numvars[j], prefix = "ck"))
          setnames(vv, cn2, gen_vnames(numvars[j], prefix = "ck_nz"))
          resnums[[j]] <- vv
        }
        rescnts <- append(rescnts, resnums)
      }

      varsdt <- data.table(
        vname = c(freqvar, wvar, rkeyvar),
        type = c("freqvar", "wvar", "rkeyvar")
      )
      varsdt <- rbind(varsdt, data.table(
        vname = dimvars, type = "dimvars"
      ))
      varsdt <- rbind(varsdt, data.table(
        vname = c("Total", countvars), type = "countvars"
      ))
      if (length(numvars) > 0) {
        varsdt <- rbind(varsdt, data.table(
          vname = numvars, type = "numvars"
        ))
      }
      varsdt$is_perturbed <- NA
      varsdt[type %in% c("countvars", "numvars"), is_perturbed := FALSE]
      varsdt$vname <- tolower(varsdt$vname)

      # modifications and perturbation parameters
      mods <- pert_params <- list(
        cnts = NULL,
        nums = NULL
      )

      private$.prob <- prob
      private$.max_contributions <- max_contributions
      private$.varsdt <- varsdt
      private$.is_weighted <- is_weighted
      private$.pert_params <- pert_params
      private$.modifications <- mods
      private$.results <- rescnts
      private$.is_initialized <- TRUE
      private$.validate()
      invisible(self)
    },
    # perturb variables
    perturb=function(v) {
      if (!is_character(v)) {
        stop("Argument `v` needs to be a character vector.", call. = FALSE)
      }
      v <- tolower(v)

      # important variables
      countvars <- self$cntvars()
      numvars <- self$numvars()
      avail <- c(countvars, numvars)
      .check_avail(v = v, avail = avail, msg = "Invalid variables specified in `v`:")
      for (i in 1:length(v)) {
        vname <- v[i]
        if (vname %in% countvars) {
          private$.ck_perturb_cnts(v = vname)
        } else {
          private$.ck_perturb_nums(v = vname)
        }
      }
      return(invisible(self))
    },
    # return a table for perturbed count variables
    freqtab=function(v=NULL, path = NULL, type = "both") {
      if (!is_scalar_character(type)) {
        stop("Argument `type` must be a scalar character.", call. = FALSE)
      }
      if (!type %in% c("both", "weighted", "unweighted")) {
        stop("Argument `type` must be either `both`, `weighted` or `unweighted`.", call. = FALSE)
      }
      if (!is.null(path)) {
        if (!is_scalar_character(path)) {
          stop("Argument `path` must be a scalar character.", call. = FALSE)
        }
      }

      avail <- private$.ck_perturbed_vars(what = "countvars")
      if (length(avail) == 0) {
        stop("No perturbed count variables found, please use the `perturb()-method` first.", call. = FALSE)
      }

      if (is.null(v)) {
        v <- avail
      }

      if (!is.character(v)) {
        stop("Argument `v` must be a character vector specifying variable names.", call. = FALSE)
      }

      v <- tolower(v)
      ii <- which(!v %in% avail)
      if (length(ii) > 0) {
        e <- c(
          "Variable(s) ", paste(shQuote(v[ii]), collapse = ", "),
          ", provided in `v`, have not yet been perturbed.")
        stop(e, collapse = " ", call. = FALSE)
      }

      # return a table for a single perturbed count variable
      .onevar_tab_freq <- function(results, vname, type, modifications = NULL) {
        debug <- ifelse(is.null(modifications), FALSE, TRUE)

        dt1 <-  results$dims
        dt1$strID <- NULL
        dt1$vname <- vname
        dt2 <-  results[[vname]]

        cn <- names(dt2)
        setnames(dt2, gsub(paste0("_", vname), "", cn))

        dt2$pert <- NULL

        if (type == "weighted") {
          dt2$uwc <- dt2$puwc <- NULL
        } else if (type == "unweighted") {
          dt2$wc <- dt2$pwc <- NULL
        }

        dt <- cbind(dt1, dt2)
        if (debug) {
          mods <- modifications[countvar == vname, .(row_nr, pert)]
          dt <- cbind(dt, mods)
        } else {
          dt$ck <- dt$wcavg <- NULL
        }
        dt[]
      }

      debug <- FALSE
      if (debug) {
        modifications <- private$.modifications$cnts
      } else {
        modifications <- NULL
      }
      results <- copy(private$.results)
      res <- rbindlist(lapply(v, function(x) {
        .onevar_tab_freq(results = results, vname = x, type = type, modifications = modifications)
      }))

      if (!is.null(path)) {
        f_out <- paste0(path, ".csv")
        res <- tryCatch(fwrite(res, file = f_out, sep = ";"), error = function(e) e)
        if (inherits(res, "error")) {
          warning("File ", shQuote(f_out), " could not be written to disk", call. = FALSE)
        }
        message("File ", shQuote(f_out), " successfully written to disk.")
      }
      return(res)
    },
    # return a table for perturbed numeric variables
    numtab=function(v = NULL, mean_before_sum = FALSE, path = NULL) {
      if (!is.null(path)) {
        if (!is_scalar_character(path)) {
          stop("Argument `path` must be a scalar character.", call. = FALSE)
        }
      }

      if (!rlang::is_scalar_logical(mean_before_sum)) {
        stop("Argument `mean_before_sum` needs to be a logical scalar.", call. = FALSE)
      }

      avail <- private$.ck_perturbed_vars(what = "numvars")
      if (length(avail) == 0) {
        stop("No perturbed numerical variables found, please use the `perturb()-method` first.", call. = FALSE)
      }

      if (is.null(v)) {
        v <- avail
      }
      if (!is.character(v)) {
        stop("Argument `v` must be a character vector specifying variable names.", call. = FALSE)
      }

      v <- tolower(v)
      if (!all(v %in% avail)) {
        e <- "Some provided variable(s) in `v` are not valid, already perturbed numerical variables."
        stop(e, call. = FALSE)
      }

      tab <- rbindlist(lapply(v, function(x) {
        tmp <- private$.results[[x]]
        tmp$vname <- x
        tmp <- tmp[, c(7, 3:4, 6)]
        setnames(tmp, c("vname", "uws", "ws", "pws"))
        cbind(private$.results$dims, tmp)[, -1]
      }))


      if (mean_before_sum) {
        ck_log("applying correction based on argument `mean_before_sum`")
        tab$pws <- tab$pws * (tab$pws / tab$ws)
      }
      tab
    },

    # compute distance-based utility measures
    measures_cnts=function(v, exclude_zeros = TRUE) {
      if (!is_scalar_character(v)) {
        stop("Argument `v` needs to be a scalar character!", call. = FALSE)
      }
      if (!rlang::is_scalar_logical(exclude_zeros)) {
        stop("Argument `exclude_zeros` needs to be a scalar logical value.", call. = FALSE)
      }
      v <- tolower(v)
      if (private$.is_perturbed_countvar(v)) {
        return(private$.ck_utility_cnts(v, exclude_zeros = exclude_zeros))
      } else {
        avail <- private$.ck_perturbed_vars("countvars")
        e <- c(
          "Variable `v` is neither a perturbed count ",
          "nor a perturbed numeric variable. Possible choices are:\n",
          paste("-->", shQuote(avail), collapse = "\n")
        )
        stop(paste(e, collapse = ""), call. = FALSE)
      }
    },

    # compute measures for continuous variables
    measures_nums=function(v) {
      stop("TODO: implement method `measures_nums()`.", call = FALSE)
    },

    # return names of all-, cnt- or numerical variables
    allvars=function() {
      list(
        cntvars = self$cntvars(),
        numvars = self$numvars()
      )
    },
    cntvars=function() {
      private$.ck_vars("countvars")
    },
    numvars=function() {
      private$.ck_vars("numvars")
    },

    # return actual modifications
    mod_cnts=function() {
      private$modifications(type = "cnts")
    },
    mod_nums=function() {
      private$modifications(type = "nums")
    },

    # get/set parameters for cnt-vars
    params_cnts_get=function() {
      return(private$.pert_params$cnts)
    },
    params_cnts_set=function(val, v = NULL) {
      if (!inherits(val, "ck_params")) {
        stop("Please create the input using `ck_params_cnts()`", call. = FALSE)
      }
      if (val$type != "cnts") {
        stop("Please create the input using `ck_params_cnts()`", call. = FALSE)
      }

      cv <- self$cntvars()
      if (is.null(v)) {
        v <- cv
      } else {
        if (!is.character(v)) {
          stop("Argument `v` needs to be a character vector.", call. = FALSE)
        }
        if (!all(v %in% cv)) {
          stop("Please specify only valid countvars in argument `v`.", call. = FALSE)
        }
      }

      ex_params <- private$.pert_params$cnts
      for (curvar in v) {
        if (private$.is_perturbed_countvar(curvar)) {
          message("--> Variable ", shQuote(curvar), " was already perturbed: parameters are not updated.")
        } else {
          if (curvar %in% names(ex_params)) {
            message("--> replacing perturbation parameters for variable ", shQuote(curvar))
          } else {
            message("--> setting perturbation parameters for variable ", shQuote(curvar))
          }
          ex_params[[curvar]] <- val
        }
      }
      private$.pert_params$cnts <- ex_params
      return(invisible(self))
    },

    # reset perturbation and perturbation parameters
    # of the given variables
    reset_cntvars=function(v = NULL) {
      avail <- private$.ck_perturbed_vars(what = "countvars")
      if (length(avail) == 0) {
        message("No perturbed count variables available.")
        return(invisible(self))
      }
      if (is.null(v)) {
        v <- avail
      } else {
        if (!is.character(v)) {
          stop("Argument `v` needs to be a character vector.", call. = FALSE)
        }
        .check_avail(v = v, avail = avail, msg = "Invalid variables specified in `v`:")
      }
      private$.reset_vars(vars = v)
      return(invisible(self))
    },
    reset_numvars=function(v = NULL) {
      avail <- private$.ck_perturbed_vars(what = "numvars")
      if (length(avail) == 0) {
        message("No perturbed numerical variables available.")
        return(invisible(self))
      }
      if (is.null(v)) {
        v <- avail
      } else {
        if (!is.character(v)) {
          stop("Argument `v` needs to be a character vector.", call. = FALSE)
        }
        .check_avail(v = v, avail = avail, msg = "Invalid variables specified in `v`:")
      }
      private$.reset_vars(vars = v)
      return(invisible(self))
    },
    reset_allvars=function() {
      avail <- c(private$.ck_perturbed_vars(what = "countvars"), private$.ck_perturbed_vars(what = "numvars"))
      if (length(avail) == 0) {
        message("No perturbed variables available.")
        return(invisible(self))
      }
      private$.reset_vars(vars = avail)
    },

    # get/set parameters for num-vars
    params_nums_get=function(val) {
      return(private$.pert_params$nums)
    },
    params_nums_set=function(val, v = NULL) {
      nv <- self$numvars()
      if (length(nv) == 0) {
        message("no numeric variables are available!")
        return(invisible(NULL))
      }

      if (!inherits(val, "ck_params")) {
        stop("Please create the input using `ck_params_nums()`", call. = FALSE)
      }
      if (val$type != "nums") {
        stop("Please create the input using `ck_params_nums()`", call. = FALSE)
      }

      if (is.null(v)) {
        v <- nv
      } else {
        if (!is.character(v)) {
          stop("Argument `v` needs to be a character vector.", call. = FALSE)
        }
        if (!all(v %in% nv)) {
          stop("Please specify only valid numvars in argument `v`.", call. = FALSE)
        }
      }
      ex_params <- private$.pert_params$nums
      for (curvar in v) {
        if (private$.is_perturbed_numvar(curvar)) {
          message("--> Variable ", shQuote(curvar), " was already perturbed: parameters are not updated.")
        } else {
          if (curvar %in% names(ex_params)) {
            message("--> replacing perturbation parameters for variable ", shQuote(curvar))
          } else {
            message("--> setting perturbation parameters for variable ", shQuote(curvar))
          }
          ex_params[[curvar]] <- val
        }
      }
      private$.pert_params$nums <- ex_params
      return(invisible(self))
    },

    summary=function() {
      cli::cat_line(cli::boxx("Utility measures for perturbed count variables", padding = 0))

      # perturbed countvars
      cv <- private$.ck_perturbed_vars("countvars")
      if (length(cv) == 0) {
        cli::cat_line(cli::symbol$cross, " no count variables have been perturbed")
      } else {
        cli::cat_rule("Distribution statistics of perturbations")
        info_cnts <- private$modifications("cnts")
        cnt_info <- info_cnts[, as.list(get_distr_vals(pert)), by = "countvar"]
        print(cnt_info)

        # unweighted only!
        cnt_measures <- lapply(cv, function(x, excl) {
          self$measures_cnts(v = x)
        })
        names(cnt_measures) <- cv

        cli::cat_line()
        cli::cat_rule("Distance-based measures")
        for (vv in names(cnt_measures)) {
          cli::cat_line(cli::symbol$tick, " Variable: ", shQuote(vv), "\n")
          print(cnt_measures[[vv]]$measures)
          cli::cat_line()
        }
      }

      numvars_available <- FALSE
      if (numvars_available) {
        cli::cat_line(cli::boxx("Utility measures for perturbed numerical variables", padding = 0))
        nv <- private$.ck_perturbed_vars("numvars")
        if (length(nv) == 0) {
          cli::cat_line(cli::symbol$cross, " no numerical variables have been perturbed")
        } else {
          # info_nums <- mod_numvars(object)
          # num_info <- num_ratios <- NULL
          # if (nrow(info_nums) > 0) {
          #   num_info <- info_nums[, as.list(get_distr_vals(vals.pert)), by = "numVar"]
          #   cat("\nPerturbation statistics on numerical variables:\n")
          #   print(num_info)
          # }
        }
      }
      return(invisible(NULL))
      # return(invisible(
      #   list(cnt_info = cnt_info,
      #        cnt_measures = cnt_measures,
      #        num_info = num_info)
      # ))
    },
    print=function() {
      cli::cat_rule("Table Information")
      nr_dims <- ncol(private$.results$dims) - 1
      nr_cells <- nrow(private$.results$dims)
      dims <- paste(shQuote(private$.ck_vars("dimvars")), collapse = ", ")
      cli::cat_line(cli::symbol$tick, " ", nr_cells, " cells in ", nr_dims, " dimensions (", dims, ")")

      if (private$.is_weighted) {
        cli::cat_line(cli::symbol$tick, " weights: yes")
      } else {
        cli::cat_line(cli::symbol$tick, " weights: no")
      }

      cli::cat_rule("Tabulated / Perturbed countvars")
      countvars <- private$.varsdt[type == "countvars"]
      for (i in 1:nrow(countvars)) {
        v <- countvars$vname[i]
        if (countvars$is_perturbed[i]) {
          cli::cat_line(cli::symbol$checkbox_on, " ", shQuote(v), " (perturbed)")
        } else {
          cli::cat_line(cli::symbol$checkbox_off, " ", shQuote(v))
        }
      }

      numvars <- private$.varsdt[type == "numvars"]
      if (nrow(numvars) > 0) {
        cli::cat_rule("Tabulated / Perturbed numvars")
        for (i in 1:nrow(numvars)) {
          v <- numvars$vname[i]
          if (numvars$is_perturbed[i]) {
            cli::cat_line(cli::symbol$checkbox_on, " ", shQuote(v), " (perturbed)")
          } else {
            cli::cat_line(cli::symbol$checkbox_off, " ", shQuote(v))
          }
        }
      }
      return(invisible(NULL))
    }
    # to delete
    #max_contributions=function() { private$.max_contributions },
    #results=function() { private$.results }
  ),
  private=list(
    .prob = NULL,
    .max_contributions = NULL,
    .varsdt = NULL,
    .is_weighted = NULL,
    .pert_params = list(),
    .modifications = list(),
    .results = list(),
    .is_initialized = FALSE,

    # return important variables
    .ck_vars = function(what) {
      vname <- type <- NULL
      vv <- private$.varsdt
      if (!is.null(what)) {
        stopifnot(length(what) == 1)
        stopifnot(what %in% c("freqvar", "countvars", "numvars", "wvar", "dimvars", "rkeyvar"))
        return(vv[type == what, vname])
      }
      return(vv)
    },
    # names of already perturbed count or continuous variables
    .ck_perturbed_vars=function(what) {
      stopifnot(what %in% c("countvars", "numvars"))
      type <- is_perturbed <- vname <- NULL
      private$.varsdt[type == what & is_perturbed == TRUE, vname]
    },
    # actually perturb count variables
    .ck_perturb_cnts=function(v) {
      lookup_freq <- function(ptab, cellkeys, cellvals) {
        row_nr <- rep(-1, length(cellkeys))
        pert_vals <- rep(0L, length(cellkeys))
        symmetry <- max(ptab$i)
        ptab$ids <- 1:nrow(ptab)

        for (d in 1:(symmetry)) {
          if (d == symmetry) {
            rkind <- cellvals >= d
          } else {
            rkind <- cellvals == d
          }
          if (sum(rkind) > 0) {
            ii <- which(ptab$i == d)
            v <- ptab$p_int_ub[ii]
            ids <- ptab$ids[ii]
            ck <- cellkeys[rkind]
            diffs <- ptab$v[ii]

            # row_ind
            rowind <- sapply(1:sum(rkind), function(x) {
              which.max(ck[x] <= v)
            })
            pert_vals[rkind] <- as.integer(diffs[rowind])
            row_nr[rkind] <- ids[rowind]
          }
        }
        data.table(
          row_nr = row_nr,
          pert = pert_vals,
          ckey = cellkeys
        )
      }

      if (private$.is_perturbed_countvar(v)) {
        message("Variable ", shQuote(v),  " was already perturbed!")
        return(invisible(self))
      }

      countvars <- private$.ck_vars("countvars")
      numvars <- private$.ck_vars("numvars")
      rkeyvar <- private$.ck_vars("rkeyvar")
      freqvar <- private$.ck_vars("freqvar")
      wvar <- private$.ck_vars("wvar")
      dimvars <- private$.ck_vars("dimvars")

      # perturbation parameters
      params <- private$.pert_params$cnts[[v]]
      if (is.null(params)) {
        e <- c(
          "No perturbation parameters found for variable", shQuote(v), "-",
          "Please use method `$params_cnts_set().")
        stop(paste(e, collapse = " "), call. = FALSE)
      }

      newtab <- copy(private$.results[[v]])
      # compute average weights
      col_wcavg <- gen_vnames(v, prefix = "wcavg")
      col_n <- gen_vnames(v, prefix = "uwc")
      col_wc <- gen_vnames(v, prefix = "wc")
      wcavg <- newtab[[col_wc]] / newtab[[col_n]]
      wcavg[is.nan(wcavg)] <- 0
      newtab[[col_wcavg]] <- wcavg

      # compute cell-keys for each column independently
      col_ck <- gen_vnames(v, prefix = "ck")

      # list with perturbation results for each countvar
      col_pert <- gen_vnames(v, prefix = "pert")

      # lookup perturbation values
      mods <- lookup_freq(
        ptab = params$params$ptable@pTable,
        cellkeys = newtab[[col_ck]],
        cellvals = newtab[[col_n]]
      )
      mods$countvar <- v

      mods <- cbind(private$.results$dims[, dimvars, with = FALSE], mods)
      newtab[[col_pert]] <- mods$pert
      modifications <- rbind(private$.modifications$cnts, mods)

      # compute perturbed unweighted counts
      col_puwc <- gen_vnames(v, prefix = "puwc")
      col_pwc <- gen_vnames(v, prefix = "pwc")

      # compute perturbed unweighted counts
      puwc <- newtab[[col_pert]] + newtab[[col_n]]
      puwc[is.na(puwc)] <- 0
      newtab[[col_puwc]] <- puwc

      # compute perturbed weighted counts
      pwc <- newtab[[col_puwc]] + newtab[[col_wcavg]]
      pwc[is.na(pwc) | is.infinite(pwc)] <- 0
      newtab[[col_pwc]] <- pwc

      # restrict table to required columns
      private$.results[[v]] <- newtab
      private$.modifications$cnts <- modifications

      # update vars
      varsdt <- private$.varsdt
      varsdt$is_perturbed[varsdt$vname == v] <- TRUE
      private$.varsdt <- varsdt
      private$.validate()
      message("Count variable ", shQuote(v), " was perturbed.")
      return(invisible(self))
    },
    # actually perturb continuous variables
    .ck_perturb_nums=function(v) {
      if (private$.is_perturbed_numvar(v)) {
        message("Variable ", shQuote(v),  " was already perturbed!")
        return(invisible(self))
      }

      # parameters
      params <- private$.pert_params$nums[[v]]$params
      ck_log("parameters:")
      ck_log("--> type: ", shQuote(params$type))
      ck_log("--> top_k: ", params$top_k)
      ck_log("--> same_key: ", params$same_key)
      ck_log("--> use_zero_rkeys: ", params$use_zero_rkeys)

      mult_params <- params$mult_params
      do_grid <- inherits(mult_params, "params_m_grid")
      if (do_grid) {
        ck_log("--> multiplicator: ", shQuote("grid"))
      } else {
        ck_log("--> multiplicator: ", shQuote("flex"))
      }

      ck_log("--> mu_c: ", params$mu_c)

      # get maximum contributions for each cell and variable!
      # restrict to current variable
      max_contr <- lapply(private$.max_contributions, function(x) {
        x[[v]]
      })

      # compute x_j values that are required input for the computation of
      # the perturbation magnitutes m_j
      w_sums <- lapply(max_contr, function(x) x$w_sum)
      if (params$type == "top_contr") {
        # we may get NAs -> need to ignore them later
        x_vals <- lapply(max_contr, function(x) x$w_vals[1:params$top_k])
      } else if (params$type == "mean") {
        x_vals <- lapply(max_contr, function(x) x$w_mean)
      } else if (params$type == "range") {
        x_vals <- lapply(max_contr, function(x) x$w_spread)
      } else if (params$type == "sum") {
        x_vals <- w_sums
      }

      # the table with cell_keys (non-zero cellkeys and weighted and unweighted sums)
      newtab <- copy(private$.results[[v]])
      dim_dt <- private$.results[["dims"]]

      # compute multiplication parameters m_j (x) * x
      .get_multipliers <- function(params, x, top_k) {
        UseMethod(".get_multipliers", params)
      }
      .get_multipliers.default <- function(params, x, top_k) {
        stop("invalid input in .get_multipliers.default() detected!", call. = FALSE)
      }
      .get_multipliers.params_m_flex <- function(params, x, top_k) {
        # by default: perturbation for small cells -> highest pctg
        m <- rep(params$m_small, length(x))
        fp <- params$flexpoint

        ind_lg <- which(x > fp)
        if (length(ind_lg) > 0) {
          x_lg <- x[ind_lg]

          f1 <- ((params$m_small * x_lg) - (params$m_large * fp)) / (params$m_large * fp)
          f2 <- (2 * fp) / (fp + x_lg)
          m[ind_lg] <- params$m_large * (1 + (f1 * f2 ^ params$q))
        }

        # fixed variance for very small observations
        if (!is.null(params$m_fixed_sq)) {
          # compute g1 (formula 2.3 on p.9)
          g1 <- sqrt(params$m_fixed_sq) / (sqrt(top_k) * params$m_large)
          # very small values
          ind_vs <- which(x < g1)
          if (length(ind_vs) > 0) {
            m[ind_vs] <- sqrt(m_fixed_sq)
            x[ind_vs] <- 1
          }
        }
        #list(x = x, m = m, x_times_m = x * m)
        x * m # x_delta!
      }
      .get_multipliers.params_m_grid <- function(params, x, top_k) {
        m <- rep(NA, length(x))
        rr <- lapply(x, function(y) y <= params$grid)

        # smallest cells, highest perturbation
        ind_sm <- sapply(rr, function(x) all(x == TRUE))
        m[ind_sm] <- params$pcts[1]

        # largest cells, smallest perturbation
        ind_lg <- sapply(rr, function(x) all(x == FALSE))
        m[ind_lg] <- tail(params$pcts, 1)

        ind_mid <- is.na(m)
        m[ind_mid] <- sapply(rr[ind_mid], function(x) {
          params$pcts[min(which(x == TRUE))]
        })
        #list(x = x, m = m, x_times_m = x * m)
        x * m
      }

      # we return x_j * m_j
      # reason: in case of fixed_variance for very small cells, x_j is changed to 1!
      # corresponds to x_delta computation!
      x_delta <- lapply(x_vals, function(x) {
        .get_multipliers(params = params$mult_params, x = x, top_k = params$top_k)
      })

      stab <- params$stab

      ck_log("selecting cell keys based on choice of argument `use_zero_rkey`")
      if (params$use_zero_rkeys) {
        cellkeys <- newtab[[paste0("ck_nz_", v)]]
      } else {
        cellkeys <- newtab[[paste0("ck_", v)]]
      }
      cellvals <- newtab[[paste0("ws_", v)]]

      # we compute top_k-scrambled cell keys
      # required for argument `same_key`
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
      ck_log("scrambling cell keys depending on argument `same_key`")
      cellkeys <- lapply(cellkeys, function(x) {
        .scramble_ckeys(x, top_k = params$top_k, same_key = params$same_key)
      })

      ck_log("computing actual perturbation values")

      # we remove duplicated cells for now and add them later!
      names(cellkeys) <- dim_dt$strID
      index_nondup <- !duplicated(cellkeys)
      cellkeys <- cellkeys[index_nondup]
      cellvals <- cellvals[index_nondup]

      .lookup_v <- function(stab, cellkeys, x_delta, cell_sum) {
        d <- max(stab$i) # parameter `D`

        v <- rep(NA, length(cellkeys))
        a <- cell_sum / x_delta
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

      # additional perturbation
      mu_c <- rep(0, params$top_k)
      mu_c[1] <- params$mu_c
      pvals <- lapply(1:length(cellkeys), function(x) {
        p <- .lookup_v(
          stab = stab,
          cellkeys = cellkeys[[x]],
          x_delta = x_delta[[x]],
          cell_sum = cellvals[x])

        # we add extra perturbation for largest contributor
        # according to formula 2.1, page 5.
        signs <- ifelse(p >= 0, 1, -1)
        (abs(p) + mu_c[1:length(p)]) * signs
      })
      names(pvals) <- names(cellkeys)

      # actual perturbation for cells are sum(pvals * x_delta)
      # these values are added to the weighted cell total to get final perturbed cell values

      # we add some more variables here that we need for the mods-table (for debugging purposes)
      ck_log("compute perturbation values and final perturbed cell values for each cells")
      pert_result <- rbindlist(lapply(names(pvals), function(x) {
        pert <- sum(pvals[[x]] * x_delta[[x]])
        cell_value_pert = w_sums[[x]] + pert
        data.table(
          str_id = x,
          cv = w_sums[[x]],
          pert = pert,
          cv_pert = cell_value_pert,
          ckey  = cellkeys[[x]][1])
      }))

      ck_log("add rows for duplicated cells!")
      pert_result <- merge(dim_dt, pert_result, by.x = "strID", by.y = "str_id", all.x = TRUE)

      ck_log("update private$.results")
      newtab <- private$.results[[v]]
      newtab$pert_total <- pert_result$pert
      newtab[[gen_vnames(v, prefix = "pws")]] <- pert_result$cv_per
      private$.results[[v]] <- newtab

      ck_log("compute and add modifications")
      # note: `ckey`: the first cell key (all others were scrambled)
      mods <- copy(pert_result)
      mods$strID <- mods$cv <- mods$cv_pert <- NULL
      mods$numvar <- v
      private$.modifications$nums <- rbind(private$.modifications$nums, mods)

      # update vars
      varsdt <- private$.varsdt
      varsdt$is_perturbed[varsdt$vname == v] <- TRUE
      private$.varsdt <- varsdt
      private$.validate()
      message("Numeric variable ", shQuote(v), " was perturbed.")
      return(invisible(self))
    },
    # utility measures for count variables
    .ck_utility_cnts=function(v, exclude_zeros) {
      v <- tolower(v)
      dt <- private$.results[[v]]
      return(ck_cnt_measures(
        orig = dt$uwc,
        pert = dt$puwc,
        exclude_zeros = exclude_zeros
      ))
    },
    # utility measures for numeric variables
    .ck_utility_nums=function(v) {
      stop(".ck_utility_nums() needs to be implemented", call. = FALSE)
    },
    # is `v` a perturbed countvar?
    .is_perturbed_countvar=function(v) {
      stopifnot(is_scalar_character(v))
      v %in% private$.ck_perturbed_vars("countvars")
    },
    # is `v` a perturbed numvar?
    .is_perturbed_numvar=function(v) {
      stopifnot(is_scalar_character(v))
      v %in% private$.ck_perturbed_vars("numvars")
    },
    # resets perturbed variables
    .reset_vars=function(vars) {
      for (v in vars) {
        ck_log("removing perturbation results and parameters for ", shQuote(v))
        private$.results[[v]] <- private$.results[[v]][, 1:4]
        if (private$.is_perturbed_countvar(v)) {
          private$.modifications$cnts <- private$.modifications$cnts[countvar != v]
          private$.pert_params$cnts[[v]] <- NULL
        } else {
          private$.modifications$nums <- private$.modifications$nums[numvar != v]
          private$.pert_params$nums[[v]] <- NULL
        }
        private$.varsdt[vname == v, is_perturbed := FALSE]
      }
    },
    # returns modification slot
    modifications=function(type = NULL) {
      if (is.null(type)) {
        return(private$.modifications)
      }
      if (!is_scalar_character(type)) {
        stop("Argument `type` needs to be a scalar characater.", call. = FALSE)
      }
      if (!type %in% c("cnts", "nums")) {
        stop("Argument `type` must be either `cnts` or `nums`.", call. = FALSE)
      }
      return(private$.modifications[[type]])
    },
    .validate = function() {
      if (!private$.is_initialized) {
        return(invisible(TRUE))
      }
      # todo
      return(invisible(TRUE))
    }
  )
)

#' @rdname cellkey_pkg
#' @name cellkey_pkg
#' @md
#' @title cellKey - an post-tabular perturbation method
#' @description This class allows to define statistical tables and perturb both
#' count and numerical variables.
#' @param x an object coercible to a `data.frame`
#' @param rkey either a column name within `x` referring to a variable containing record keys
#' or a single integer(ish) number > `5` that referns to the number of digits for record keys that
#' will be generated internally.
#' @param dims a list containing slots for each variable that should be
#' tabulated. Each slot consists should be created/modified using [sdcHierarchies::hier_create()],
#' [sdcHierarchies::hier_add()] and other functionality from package `sdcHierarchies`.
#' @param w (character) a scalar character referring to a variable in `x` holding sampling
#' weights. If `w` is `NULL` (the default), all weights are assumed to be `1`
#' @param countvars (character) an optional vector containing names of binary (0/1 coded)
#' variables withing `x` that should be included in the problem instance.
#' These variables can later be perturbed.
#' @param numvars (character) an optional vector of numerical variables that can later be tabulated.
#' @format [R6::R6Class] object.
#' @section Usage: For usage details of the implemented methods, see the **Methods** and
#' **Examples** sections.
#' @section Methods: The following methods are available and can be used:
#' - **`print()`**: prints information about the current table
#' - **`summary()`**: some aggregated summary statistics about perturbed variables
#'
#' - **`allvars()`**: a list with elements `cntvars` and `numvars` with each being a character
#' vector containing the pre-defined count- and continuous variables in the current instance.
#'
#' - **`cntvars()`**: returns a character vector of available count variables.
#'
#' - **`numvars()`**: returns a character vector of available numeric variables.
#'
#' - **`perturb(v)`**: Perturb a count- or magnitude variable. The method has the following arguments:
#'    * `v`: name(s) of count or magnitude variables that should be perturbed.
#'
#' - **`freqtab(v)`**: get results from already perturbed count variables as a `data.table`. The required arguments are:
#'    * `v`: a vector of variable names for already perturbed count variables. If `NULL` (the default), the results
#'    are returned for all perturbed count variables.
#'    * `type`: a scalar character depending what variables should be available in the output. Allowed
#'    values are `"both"` (the default), `"weighted"` and `"unweighted"`
#'    * `path`: if not `NULL`, a scalar character defining a (relative or absolute) path to which the result table
#'    should be written. A `csv` file will be generated and `.csv` will be appended to the value provided.
#'
#'    This method returns a `data.table` containing all combinations of the dimensional variables in
#' the first n columns. Additionally, the following columns are shown:
#'    * `vname`: name of the perturbed variable
#'    * `uwc`: unweighted counts (if `type` is `"both"` or `"unweighted"`)
#'    * `wc`: weighted counts (if `type` is `"both"` or `"weighted"`)
#'    * `puwc`: perturbed unweighted counts (if `type` is `"both"` or `"unweighted"`)
#'    * `pwc`: perturbed weighted counts (if `type` is `"both"` or `"weighted"`)
#'
#' - **`numtab(v)`**: get results from already perturbed continuous variables as a `data.table`. The required arguments are:
#'    * `v`: a vector of variable names for already perturbed count variables. If `NULL` (the default), the results
#'    are returned for all perturbed numeric variables.
#'    * `mean_before_sum`: (logical); if `TRUE`, the perturbed values are adjusted by a factor `((n+p))⁄n` with `n`
#'    being the original weighted cellvalue and `p` the perturbed cell value. This makes sense if the the
#'    accuracy of the variable mean is considered more important than accuracy of sums of the variable. The default value
#'    is `FALSE` (no adjustment is done)
#'    * `path`: if not `NULL`, a scalar character defining a (relative or absolute) path to which the result table
#'    should be written. A `csv` file will be generated and `.csv` will be appended to the value provided.
#'
#'    This method returns a `data.table` containing all combinations of the dimensional variables in
#' the first n columns. Additionally, the following columns are shown:
#'    * `vname`: name of the perturbed variable
#'    * `uws`: unweighted sum of the given variable
#'    * `ws`: weighted cellsum
#'    * `pws`: perturbed weighted sum of the given cell
#'
#' - **`measures_cnts(v, exclude_zeros = TRUE)`**: utility measures for perturbed count variables. The required arguments are:
#'    * `v`: name of a count variable for which utility measures should be computed.
#'    * `exclude_zeros`: should empty (zero) cells in the original values be excluded when computing distance measures
#'    This method returns a `list` containing a set of utility measures based on some distance functions.
#'    For a detailed description of the computed measures, see [ck_cnt_measures()]
#'
#' - **`mod_cnts()`**: returns a `data.table` containing modifications applied to count variables
#' - **`mod_nums()`**: returns a `data.table` containing modifications applied to numerical variables
#'
#' - **`params_cnts_get()`**: returns a named list in which each list-element contains the
#' active perturbation parameters for the specific count variable defined by the list-name.
#'
#' - **`params_cnts_set(val, v = NULL)`**: allows to set perturbation parameters for count variables.
#' The following arguments are expected:
#'    * `val`: a perturbation object created with [ck_params_cnts()]
#'    * `v`: a character vector (or `NULL`). If `NULL` (the default), the perturbation parameters
#'    provided in `val` are set for all count variables; otherwise one may specify the names of
#'    the count variables for which the parameters should be set.
#'
#' - **`params_nums_get()`**: returns a named list in which each list-element contains the
#' active perturbation parameters for the specific continuous variable defined by the list-name.
#'
#' - **`params_nums_set(val, v = NULL)`**: allows to set perturbation parameters for continuous variables.
#' The following arguments are expected:
#'    * `val`: a perturbation object created with [ck_params_nums()]
#'    * `v`: a character vector (or `NULL`). If `NULL` (the default), the perturbation parameters
#'    provided in `val` are set for all continuous variables; otherwise one may specify the names of
#'    the numeric variables for which the parameters should be set.
#'
#' - **`reset_cntvars(v = NULL)`**: allows to reset results and parameters for already perturbed count
#' variables specified in `v`. If `v` equals `NULL` (the default), the results are reset for all perturbed
#' count variables.
#' - **`reset_numvars(v = NULL)`**: allows to reset results and parameters for already perturbed numerical
#' variables specified in `v`. If `v` equals `NULL` (the default), the results are reset for all perturbed
#' numerical variables.
#' - **`reset_allvars()`**: allows to reset results and parameters for all already perturbed variables.
#' @export
#' @examples
#' x <- ck_create_testdata()
#'
#' # create some 0/1 variables that should be perturbed later
#' x[, cnt_females := ifelse(sex == "male", 0, 1)]
#' x[, cnt_males := ifelse(sex == "male", 1, 0)]
#' x[, cnt_highincome := ifelse(income >= 9000, 1, 0)]
#'
#' # create record keys
#' x$rkey <- ck_generate_rkeys(dat = x)
#'
#' # define required inputs
#'
#' # hierarchy with some bogus codes
#' d_sex <- hier_create(root = "Total", nodes = c("male", "female"))
#' d_sex <- hier_add(d_sex, root = "female", "f")
#' d_sex <- hier_add(d_sex, root = "male", "m")
#'
#' d_age <- hier_create(root = "Total", nodes = paste0("age_group", 1:6))
#' d_age <- hier_add(d_age, root = "age_group1", "ag1a")
#' d_age <- hier_add(d_age, root = "age_group2", "ag2a")
#'
#' # define the cell key object
#' countvars <- c("cnt_females", "cnt_males", "cnt_highincome")
#' numvars <- c("expend", "income", "savings")
#' tab <- ck_setup(
#'   x = x,
#'   rkey = "rkey",
#'   dims = list(sex = d_sex, age = d_age),
#'   w = "sampling_weight",
#'   countvars = countvars,
#'   numvars = numvars)
#'
#' # show some information about this table instance
#' tab$print() # identical with print(tab)
#'
#' # which variables have been defined?
#' tab$allvars()
#'
#' # count-variables
#' tab$cntvars()
#'
#' # continuous variables
#' tab$numvars()
#'
#' # add perturbation parameters for "total" variable
#' p_cnts1 <- ck_params_cnts(
#'   D = 5,
#'   V = 3,
#'   js = 2,
#'   pstay = 0.5,
#'   optim = 1,
#'   mono = TRUE)
#' tab$params_cnts_set(val = p_cnts1, v = "total")
#'
#' # create alternative perturbation parameters
#' p_cnts2 <- ck_params_cnts(
#'   D = 8,
#'   V = 3,
#'   js = 2,
#'   pstay = 0.5,
#'   optim = 1,
#'   mono = TRUE)
#'
#' # use it for the remaining variables
#' tab$params_cnts_set(val = p_cnts2, v = countvars)
#'
#' # perturb a variable
#' tab$perturb(v = "total")
#'
#' # multiple variables can be perturbed as well
#' tab$perturb(v = c("cnt_males", "cnt_highincome"))
#'
#' # return results (weighted and unweighted)
#' tab$freqtab(v = c("total", "cnt_males"))
#'
#' # only unweighted and perturbed unweighted counts
#' tab$freqtab(v = c("total", "cnt_males"), type = "unweighted")
#'
#' # numerical variables
#' p_nums1 <- ck_params_nums(
#'   D = 10,
#'   l = 0.5,
#'   type = "top_contr",
#'   top_k = 3,
#'   mult_params = ck_flexparams(
#'     flexpoint = 1000,
#'     m_fixed_sq = 2
#'   )
#' )
#'
#' # another set of parameters, using a grid
#' p_nums2 <- ck_params_nums(
#'   D = 10,
#'   l = .5,
#'   type = "mean",
#'   mult_params = ck_gridparams(
#'     grid = c(0, 10, 100, 10000),
#'     pcts = c(0.25, 0.20, 0.10, 0.05)
#'   )
#' )
#'
#' # use `p_nums1` for all variables
#' tab$params_nums_set(p_nums1)
#'
#' # use different parameters for variable `savings`
#' tab$params_nums_set(p_nums2, v = "savings")
#'
#' # perturb variables
#' tab$perturb(v = c("income", "savings"))
#'
#' # extract results
#' tab$numtab("income", mean_before_sum = TRUE)
#' tab$numtab("income", mean_before_sum = FALSE)
#' tab$numtab("savings")
#'
#' # results can be resetted, too
#' tab$reset_cntvars(v = "cnt_males")
#'
#' # we can then set other parameters and perturb again
#' tab$params_cnts_set(val = p_cnts1, v = "cnt_males")
#' tab$perturb(v = "cnt_males")
#'
#' # write to a file "outtab.csv" (.csv is automatically added to the path)
#' \dontrun{
#' tab$freqtab(v = c("total", "cnt_males"), path = "outtab")
#' }
#'
#' # only weighted and perturbed weighted counts
#' tab$freqtab(v = c("total", "cnt_males"), type = "weighted")
#'
#' # utility measures for a count variable
#' tab$measures_cnts(v = "total", exclude_zeros = TRUE)
#'
#' # modifications for perturbed count variables
#' tab$mod_cnts()
#'
#' # display a summary about utility measures
#' tab$summary()
ck_setup <- function(x, rkey, dims, w, countvars = NULL, numvars = NULL) {
  cellkey_obj_class$new(
    x = x,
    rkey = rkey,
    dims = dims,
    w = w,
    countvars = countvars,
    numvars = numvars)
}

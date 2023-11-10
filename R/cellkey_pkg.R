#' R6 Class defining statistical tables that can be perturbed
#'
#' Such objects are typically generated using [ck_setup()].
#' @md
#' @rdname cellkey_pkg
ck_class <- R6::R6Class("cellkey_obj", cloneable = FALSE,
  public = list(
    #' @description Create a new table instance
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
    #' @return A new `cellkey_obj` object. Such objects (internally) contain the fully computed
    #' statistical tables given input microdata (`x`), the hierarchical definitionals (`dims`) as
    #' well as the remaining inputs. Intermediate results are stored internally and can only be
    #' modified / accessed via the exported public methods described below.
    initialize = function(x, rkey, dims, w = NULL, countvars = NULL, numvars = NULL) {
      type <- is_perturbed <- NULL

      if (!inherits(x, "data.frame")) {
        stop("`x` is not coercible to a `data.frame`.", call. = FALSE)
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
          stop("number of digits for record keys must be >= 5!", call. = FALSE)
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

      ck_log("compute dim_info")
      h_info <- lapply(dims, hier_info)
      h_info <- lapply(h_info, function(x) {
        data.table::rbindlist(lapply(x, function(y) {
          data.frame(
            code = y$name,
            level = y$level,
            is_leaf = y$is_leaf,
            parent = y$parent,
            stringsAsFactors = FALSE)
        }))
      })

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
        chk <- sapply(countvars, function(cv) {
          !is.numeric(x[[cv]])
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

        setnames(x, numvars, tolower(numvars))
        numvars <- tolower(numvars)

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

      prob <- sdcTable::makeProblem(
        data = x,
        dimList = dims,
        freqVarInd = match(freqvar, names(x)),
        numVarInd = match(nv, names(x)),
        weightInd = NULL,
        sampWeightInd = NULL)

      # with perturbations
      tab <- sdcTable::sdcProb2df(
        obj = prob,
        addDups = TRUE,
        addNumVars = TRUE,
        dimCodes = "original")

      ck_log("mark duplicated cells")
      dims <- prob@dimInfo@dimInfo
      tab$is_bogus <- FALSE
      for (i in seq_len(length(dims))) {
        dd <- dims[[i]]
        tab$is_bogus[tab[[dd@vName]] %in% dd@dups] <- TRUE
      }

      # duplicated cells (bogus codes) have the same
      # strID as their parent cell!
      strids <- unique(tab[["strID"]])

      # calculate contributing indices
      # we do this in any case!
      ck_log("compute contributing indices")
      # the .contributing_indices_tmp() should be removed
      # once sdcTable 0.32 is on cran as this version has the
      # appropriate fixes;
      .contributing_indices_tmp <- function(prob, ids = NULL) {
        # returns all contributing codes for each dimensions
        # of a sdcProblem-object; dimensions are converted
        # sdcHierarchies-trees and hier_info() is then used.
        .get_all_contributing_codes <- function(x) {
          .sdchier_from_sdc <- function(d) {
            df <- data.frame(
              levels = slot(d, "levels"),
              codes = slot(d, "codesOriginal"),
              stringsAsFactors = FALSE
            )
            df$levels <- sapply(seq_len(nrow(df)), function(x) {
              paste0(rep("@", df$levels[x]), collapse = "")
            })
            hier_import(df, from = "df")
          }

          stopifnot(inherits(x, "sdcProblem"))
          dims_hier <- lapply(x@dimInfo@dimInfo, function(x) {
            .sdchier_from_sdc(x)
          })

          dims_info <- lapply(dims_hier, function(x) {
            hier_info(x)
          })

          all_contr_codes <- lapply(dims_info, function(x) {
            if ("is_rootnode" %in% names(x)) {
              res <- list()
              res[[x$name]] <- list(
                is_root = x$is_rootnode,
                contr_codes = x$contributing_codes
              )
            } else {
              res <- lapply(x, function(y) {
                list(
                  is_root = y$is_rootnode,
                  contr_codes = y$contributing_codes
                )
              })
            }
            return(res)
          })
          all_contr_codes
        }

        . <- NULL
        dt <- sdcProb2df(prob, addDups = FALSE, dimCodes = "original")
        poss_ids <- dt$strID

        if (is.null(ids)) {
          ids <- poss_ids
        } else {
          if (!is.character(ids)) {
            stop("Please provide a character vector in argument `ids`.", call. = FALSE)
          }
          if (!all(ids %in% poss_ids)) {
            e <- c(
              "Some values provided in `ids` are not valid. ",
              "See column `strID` in `sdcProb2df()` for valid ids."
            )
            stop(paste(e, collapse = " "), call. = FALSE)
          }
        }

        dimvars <- slot(prob, "dimInfo")@vNames
        nr_dims <- length(dimvars)

        dt <- dt[, c("strID", "freq", dimvars), with = FALSE]
        data.table::setnames(dt, old = "strID", new = "id")

        # we compute all unique codes once
        unique_codes <- lapply(dt[, dimvars, with = FALSE], function(x) {
          sort(unique(x))
        })

        # get contributing codes
        contr_codes <- .get_all_contributing_codes(prob)

        # positions in strID
        str_info <- prob@dimInfo@strInfo
        names(str_info) <- dimvars

        # merge inner cell-info to data
        dt_inner <- data.table(id = prob@dimInfo@strID, is_inner = TRUE)
        dt_inner$idx <- seq_len(nrow(dt_inner))
        dt_inner <- dt[dt_inner, on = "id"]

        dt_inner$tmp <- apply(dt_inner[, dimvars, with = FALSE], 1, paste0, collapse = "")
        setkeyv(dt_inner, "tmp")

        # subsetting dt to those ids, we want to compute the contributing indices from
        dt <- dt[.(ids), on = "id"]

        # prepare output
        res <- vector("list", length(ids))
        names(res) <- ids

        for (i in seq_len(nrow(dt))) {
          str_id <- dt$id[i]
          if (dt$freq[i] > 0) {
            index_vec <- which(dt_inner$id == str_id)
            if (length(index_vec) > 0) {
              res[[str_id]] <- dt_inner$idx[index_vec]
            } else {
              lev_info <- vector("list", length = nr_dims)
              names(lev_info) <- dimvars
              for (dv in dimvars) {
                code <- dt[[dv]][i]
                info <- contr_codes[[dv]][[code]]
                if (!info$is_root) {
                  lev_info[[dv]] <- info$contr_codes
                } else {
                  lev_info[[dv]] <- unique_codes[[dv]]
                }
              }
              cell_indices <- sdcTable:::pasteStrVec(unlist(expand.grid(lev_info)), nr_dims)
              res[[str_id]] <- which(dt_inner$tmp %in% cell_indices)
            }
          }
        }
        return(res)
      }

      if (length(numvars) == 0) {
        ck_log("no numerical variables were defined in `ck_setup()`, ",
        "therefore no contributing indices need to be computed.")
        max_contributions <- contr_indices <- NULL
      } else {
        ck_log("compute contributing indices for (sub)-totals")
        contr_indices <- sdcTable::contributing_indices(
         prob = prob, ids = NULL
        )
        ck_log("find top k contributors for each cell and numerical variable")

        # finding top_k contributors for each cell and numerical variable
        microdat <- prob@dataObj@rawData[, c(names(dims), numvars, wvar), with = FALSE]
        microdat$.tmpid <- seq_len(nrow(x))

        # perhaps c++?
        # for each numerical variable and each cell, get the top_k contributions
        # along with its ids and its values
        .get_max_contributions <- function(indices, microdat, wvar, nv, top_k) {
          res <- vector("list", length = length(indices))
          names(res) <- names(indices)

          # we have no numerical-variables -> return empty list
          if (length(nv) == 0) {
            return(res)
          }

          for (i in seq_len(length(res))) {
            out <- vector("list", length = length(nv))
            names(out) <- nv
            xx <- subset(microdat, .tmpid %in% indices[[i]])
            top_k <- min(top_k, nrow(xx))
            for (v in nv) {
              xx$.tmpordervar <- abs(xx[[v]])
              xx$.tmpweightvar <- xx[[v]] * xx[[wvar]]
              setorderv(xx, c(".tmpordervar", wvar), order = c(-1L, -1L))
              if (nrow(xx) == 0) {
                out[[v]]$uw_vals <- out[[v]]$w_vals <- 0
                out[[v]]$uw_ids <- out[[v]]$w_ids <- NA
                out[[v]]$uw_spread <- out[[v]]$w_spread <- 0
                out[[v]]$uw_sum <- out[[v]]$w_sum <- 0
                out[[v]]$uw_mean <- out[[v]]$w_mean <- 0
              } else {
                out[[v]]$uw_vals <- xx[[v]][1:top_k]
                out[[v]]$uw_ids <- out[[v]]$w_ids <- xx$.tmpid[1:top_k]
                out[[v]]$uw_spread <- diff(range(xx[[v]], na.rm = TRUE))
                out[[v]]$uw_sum <- sum(xx[[v]], na.rm = TRUE)
                out[[v]]$uw_mean <- out[[v]]$uw_sum / nrow(xx)
                out[[v]]$w_vals <- xx$.tmpweightvar[1:top_k]
                out[[v]]$w_spread <- diff(range(xx[[v]], na.rm = TRUE))
                out[[v]]$w_sum <- sum(xx$.tmpweightvar, na.rm = TRUE)
                out[[v]]$w_mean <- out[[v]]$w_sum / sum(xx[[wvar]], na.rm = TRUE)
              }
              # we compute if the number of contributors to the cell
              # is even or odd. This information can later be used if
              # we have different ptables (parity-case)
              out[[v]]$even_contributors <- nrow(xx) %% 2 == 0
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
          top_k = 6)
      }
      res <- tab[, c("strID", "freq", dimvars, "is_bogus", nv), with = FALSE]
      cols_ck <- gen_vnames(c("total", countvars), prefix = "ck")

      if (length(countvars) > 0) {
        setnames(res, c(rkeyvar, paste0("rkey_", countvars)), cols_ck)
      } else {
        setnames(res, rkeyvar, cols_ck)
      }

      # fix cell keys
      for (i in seq_len(length(cols_ck))) {
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
      rescnts[[1]] <- res[, c("strID", dimvars, "is_bogus"), with = FALSE]

      for (i in seq_along(cv)) {
        index <- i + 1
        vv <- c(
          gen_vnames(cv[i], prefix = "ck"),
          gen_vnames(cv[i], prefix = "uwc"),
          gen_vnames(cv[i], prefix = "wc")
        )
        rescnts[[index]] <- res[, vv, with = FALSE]
      }

      nr_numvars <- length(numvars)
      if (nr_numvars > 0) {
        seq_nums <- seq_len(nr_numvars)
        resnums <- vector("list", length = length(numvars))
        names(resnums) <- numvars

        cols_uws <- gen_vnames(numvars, prefix = "uws")
        setnames(res, numvars, cols_uws)
        cols_ws <- gen_vnames(numvars, prefix = "ws")
        for (j in seq_nums) {
          cn1 <- gen_vnames("total", prefix = "ck")
          cn2 <- gen_vnames(numvars[j], prefix = "rkey_nozero")
          vv <- c(cn1, cn2, cols_uws[j], cols_ws[j])
          vv <- res[, vv, with = FALSE]
          vv[[cn2]] <- vv[[cn2]] %% 1 # normalize cellkey

          setnames(vv, cn1, gen_vnames(numvars[j], prefix = "ck"))
          setnames(vv, cn2, gen_vnames(numvars[j], prefix = "ck_nz"))
          vv[["special_protection"]] <- FALSE
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
        vname = c("total", countvars), type = "countvars"
      ))
      if (length(numvars) > 0) {
        varsdt <- rbind(varsdt, data.table(
          vname = numvars, type = "numvars"
        ))
      }
      varsdt$is_perturbed <- NA
      varsdt[type %in% c("countvars", "numvars"), is_perturbed := FALSE]

      # modifications and perturbation parameters
      mods <- pert_params <- list(
        cnts = NULL,
        nums = NULL
      )

      dupsinfo <- data.table(
        id = seq_len(nrow(rescnts[[1]])),
        is_bogus = rescnts[[1]]$is_bogus)
      rescnts[[1]]$is_bogus <- NULL

      private$.prob <- prob
      private$.max_contributions <- max_contributions
      private$.varsdt <- varsdt
      private$.is_weighted <- is_weighted
      private$.pert_params <- pert_params
      private$.modifications <- mods
      private$.results <- rescnts
      private$.dupsinfo <- dupsinfo
      private$.h_info <- h_info
      private$.is_initialized <- TRUE
      private$.validate()
      invisible(self)
    },

    #' @description Perturb a count- or magnitude variable
    #' @param v name(s) of count- or magnitude variables that should be perturbed.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. Updated data can be accessed using other exported
    #' methods like `$freqtab()` or `$numtab()`.
    perturb = function(v) {
      # important variables
      countvars <- self$cntvars()
      numvars <- self$numvars()
      avail <- c(countvars, numvars)
      .check_avail(
        v = v,
        avail = avail,
        msg = "Invalid variables specified in `v`:",
        single_v = FALSE)

      for (i in seq_len(length(v))) {
        vname <- v[i]
        if (vname %in% countvars) {
          private$.ck_perturb_cnts(v = vname)
        } else {
          private$.ck_perturb_nums(v = vname)
        }
      }
      return(invisible(self))
    },

    #' @description Extract results from already perturbed count variables as a
    #' `data.table`
    #' @param v a vector of variable names for count variables. If `NULL`
    #' (the default), the results are returned for all available count
    #' variables. For variables that have not yet perturbed, columns
    #' `puwc` and `pwc` are filled with `NA`.
    #' @param path if not `NULL`, a scalar character defining a (relative
    #' or absolute) path to which the result table should be written. A `csv`
    #' file will be generated and, if specified, `path` must have
    #' ".csv" as file-ending
    #' @return This method returns a `data.table` containing all combinations of the dimensional variables in
    #' the first n columns. Additionally, the following columns are shown:
    #' - `vname`: name of the perturbed variable
    #' - `uwc`: unweighted counts
    #' - `wc`: weighted counts
    #' - `puwc`: perturbed unweighted counts or `NA` if `vname` was not yet perturbed
    #' - `pwc`: perturbed weighted counts or `NA` if `vname` was not yet perturbed
    freqtab = function(v = NULL, path = NULL) {
      if (!is.null(path)) {
        .valid_path(path, ext = "csv", check_exists = FALSE)
      }

      avail <- self$cntvars()
      if (is.null(v)) {
        v <- avail
      } else {
        if (!is.character(v)) {
          stop("Argument `v` must be a character vector specifying variable names.", call. = FALSE)
        }
        v <- tolower(v)
        .check_avail(
          v = v,
          avail = avail,
          msg = "Invalid variables specified in `v`:",
          single_v = FALSE)
      }

      # return a table for a single perturbed count variable
      .onevar_tab_freq <- function(results, vname, is_perturbed) {
        dt1 <-  results$dims
        dt1[["strID"]] <- NULL
        dt1$vname <- vname
        dt2 <-  results[[vname]]

        cn <- names(dt2)
        setnames(dt2, gsub(paste0("_", vname), "", cn))

        if (is_perturbed) {
          dt2$pert <- NULL
          dt2$ck <- dt2$wcavg <- NULL
        } else {
          dt2$puwc <- NA_real_
          dt2$pwc <- NA_real_
          dt2$ck <- NULL
        }
        dt <- cbind(dt1, dt2)
        dt[]
      }

      results <- copy(private$.results)
      res <- rbindlist(lapply(v, function(x) {
        .onevar_tab_freq(
          results = results,
          vname = x,
          is_perturbed = private$.is_perturbed_countvar(x))
      }))

      if (!is.null(path)) {
        if (file.exists(path)) {
          file.remove(path)
        }
        res <- tryCatch(fwrite(res, file = path, sep = ";"), error = function(e) e)
        if (inherits(res, "error")) {
          warning("File ", shQuote(path), " could not be written to disk", call. = FALSE)
        }
        message("File ", shQuote(path), " successfully written to disk.")
      }
      return(res)
    },

    #' @description Extract results from already perturbed continuous variables
    #' as a `data.table`.
    #' @param v a vector of variable names of continuous variables. If `NULL`
    #' (the default), the results are returned for all available numeric variables.
    #' @param mean_before_sum (logical); if `TRUE`, the perturbed values are adjusted
    #' by a factor `((n+p))⁄n` with
    #' - `n`: the original weighted cell value
    #' - `p`: the perturbed cell value
    #'
    #' This makes sense if the the accuracy of the variable mean is considered to be
    #' more important than accuracy of sums of the variable. The default value is
    #' `FALSE` (no adjustment is done)
    #' @param path if not `NULL`, a scalar character defining a (relative or absolute)
    #' path to which the result table should be written. A `csv` file will be generated
    #' and, if specified, `path` must have ".csv" as file-ending
    #'
    #' @return This method returns a `data.table` containing all combinations of the
    #' dimensional variables in the first n columns. Additionally, the following
    #' columns are shown:
    #' - `vname`: name of the perturbed variable
    #' - `uws`: unweighted sum of the given variable
    #' - `ws`: weighted cellsum
    #' - `pws`: perturbed weighted sum of the given cell or `NA` if `vname`
    #' has not not perturbed
    numtab = function(v = NULL, mean_before_sum = FALSE, path = NULL) {
      if (!is.null(path)) {
        .valid_path(path, ext = "csv", check_exists = FALSE)
      }

      if (!rlang::is_scalar_logical(mean_before_sum)) {
        stop("Argument `mean_before_sum` needs to be a logical scalar.", call. = FALSE)
      }

      avail <- self$numvars()
      if (length(avail) == 0) {
        stop("No numerical variables found in this table-instance.", call. = FALSE)
      }

      if (is.null(v)) {
        v <- avail
      }
      if (!is.character(v)) {
        stop("Argument `v` must be a character vector specifying variable names.", call. = FALSE)
      } else {
        .check_avail(
          v = v,
          avail = avail,
          msg = "Invalid continuous variables specified in `v`:",
          single_v = FALSE)
      }

      tab <- rbindlist(lapply(v, function(x) {
        tmp <- private$.results[[x]]
        tmp$vname <- x

        if (private$.is_perturbed_numvar(x)) {
          tmp <- tmp[, c(8, 3:4, 7)]
        } else {
          tmp <- tmp[, c(6, 3:4)]
          tmp$pws <- NA_real_
        }
        setnames(tmp, c("vname", "uws", "ws", "pws"))
        cbind(private$.results$dims, tmp)[, -1]
      }))

      if (mean_before_sum) {
        ck_log("applying correction based on argument `mean_before_sum`")
        tab$pws <- tab$pws * (tab$pws / tab$ws)
      }

      if (!is.null(path)) {
        if (file.exists(path)) {
          file.remove(path)
        }
        res <- tryCatch(fwrite(tab, file = path, sep = ";"), error = function(e) e)
        if (inherits(res, "error")) {
          warning("File ", shQuote(path), " could not be written to disk", call. = FALSE)
        }
        message("File ", shQuote(path), " successfully written to disk.")
      }
      tab
    },

    #' @description Utility measures for perturbed count variables
    #' @param v name of a count variable for which utility measures
    #' should be computed.
    #' @param exclude_zeros should empty (zero) cells in the original values
    #' be excluded when computing distance measures
    #' @return This method returns a `list` containing a set of utility
    #' measures based on some distance functions. For a detailed description
    #' of the computed measures, see [ck_cnt_measures()]
    measures_cnts = function(v, exclude_zeros = TRUE) {
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
          "nor a perturbed numeric variable.\n")
        if (length(avail) > 0) {
          e <- c(e, "Possible choices are: ",
          paste(shQuote(avail), collapse = "\n"))
        } else {
          e <- c(e, "Please perturb one of the following variables first: ",
            paste(shQuote(self$cntvars()), collapse = "\n"))
        }
        stop(paste(e, collapse = ""), call. = FALSE)
      }
    },

    #' @description Utility measures for continuous variables (not yet implemented)
    #' @param v name of a continuous variable for which utility measures
    #' should be computed.
    #' @return for (now) an empty list; In future versions of the package, the
    #' Method will return utility measures for perturbed magnitude tables.
    measures_nums = function(v) {
      warning("No utility measures for continuous variables are available yet.", call. = FALSE)
      return(list())
    },

    #' @description Names of variables that can be perturbed / tabulated
    #' @return returns a `list` with the following two elements:
    #' - `cntvars`: character vector with names of available count
    #' variables for perturbation
    #' - `numvars`: character vector with names of available numerical
    #' variables for perturbation
    allvars = function() {
      list(
        cntvars = self$cntvars(),
        numvars = self$numvars()
      )
    },

    #' @description Names of count variables that can be perturbed
    #' @return a character vector containing variable names
    cntvars = function() {
      private$.ck_vars("countvars")
    },

    #' @description Names of continuous variables that can be perturbed
    #' @return a character vector containing variable names
    numvars = function() {
      private$.ck_vars("numvars")
    },

    #' @description Information about hierarchies
    #' @return a `list` (for each dimensional variable) with
    #' information on the hierarchies. This may be used to restrict output tables to
    #' specific levels or codes. Each list element is a `data.table` containing
    #' the following variables:
    #' - `code`: the name of a code within the hierarchy
    #' - `level`: number defining the level of the code; the higher the number,
    #' the lower the hierarchy with `1` being the overall total
    #' - `is_leaf`: if `TRUE`, this code is a leaf node which means no other codes
    #' contribute to it
    #' - `parent`: name of the parent code
    #'
    hierarchy_info = function() {
      return(private$.h_info)
    },

    #' @description Modifications applied to count variables
    #' @return a `data.table` containing modifications applied to count variables
    mod_cnts = function() {
      private$modifications(type = "cnts")
    },

    #' @description Modifications applied to numerical variables
    #' @return a `data.table` containing modifications applied to numerical variables
    mod_nums = function() {
      private$modifications(type = "nums")
    },

    #' @description Identify sensitive cells based on minimum frequency rule
    #' @param v a single variable name of a continuous variable (see method `numvars()`)
    #' @param n a number defining the threshold. All cells `<= n` are considered as unsafe.
    #' @param weighted if `TRUE`, the weighted number of contributors to a cell are compared to
    #' the threshold specified in `n` (default); else the unweighted number of contributors is used.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other methods (e.g `$perturb()`).
    supp_freq = function(v, n, weighted = TRUE) {
      if (!rlang::is_scalar_logical(weighted)) {
        stop("argument `weighted` is not a logical value.", call. = FALSE)
      }
      .check_avail(
        v = v,
        avail = self$numvars(),
        msg = "Invalid variable specified in `v`:",
        single_v = TRUE)

      if (!rlang::is_scalar_double(n)) {
        stop("`n` is not a number.", call. = FALSE)
      }

      if (!weighted) {
        testv <- "freq"
      } else {
        testv <- "tmpweightvarfortabulation"
      }
      res <- sdcProb2df(private$.prob, addDups = TRUE, addNumVars = TRUE)[[testv]]
      pat <- res <= n
      private$.update_supps(v = v, pat = pat, rule = "freq-rule")
      invisible(self)
    },

    #' @description Identify sensitive cells based on weighted or unweighted cell value
    #' @param v a single variable name of a continuous variable (see method `numvars()`)
    #' @param n a number defining the threshold. All cells `<= n` are considered as unsafe.
    #' @param weighted if `TRUE`, the weighted cell value of variable `v` is compared to
    #' the threshold specified in `n` (default); else the unweighted number is used.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other methods (e.g `$perturb()`).
    supp_val = function(v, n, weighted = TRUE) {
      if (!rlang::is_scalar_logical(weighted)) {
        stop("argument `weighted` is not a logical value.", call. = FALSE)
      }
      .check_avail(
        v = v,
        avail = self$numvars(),
        msg = "Invalid variable specified in `v`:",
        single_v = TRUE)

      if (!rlang::is_scalar_double(n)) {
        stop("`n` is not a number.", call. = FALSE)
      }

      if (weighted) {
        testv <- paste0("ws_", v)
      } else {
        testv <- v
      }
      res <- sdcProb2df(private$.prob, addDups = TRUE, addNumVars = TRUE)[[testv]]
      pat <- res <= n
      private$.update_supps(v = v, pat = pat, rule = "val-rule")
      invisible(self)
    },

    #' @description Identify sensitive cells based on their `names`
    #' @param v a single variable name of a continuous variable (see method `numvars()`)
    #' @param inp a `data.frame` where each colum represents a dimensional variable. Each row of
    #' this input is then used to compute the relevant cells to be identified as sensitive where
    #' `NA`-values are possible and used to match any characteristics of the dimensional variable.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other methods (e.g `$perturb()`).
    supp_cells = function(v, inp) {
      .check_avail(
        v = v,
        avail = self$numvars(),
        msg = "Invalid variable specified in `v`:",
        single_v = TRUE
      )

      if (!inherits(inp, "data.frame")) {
        stop("argument `inp` needs to be a `data.frame`", call. = FALSE)
      }

      dv <- sort(private$.ck_vars("dimvars"))
      if (!identical(sort(names(inp)), dv)) {
        e <- c(
          "invalid and/or missing columns detected in argument `inp`",
          "(needs to match names of dimensional variables)"
        )
        stop(paste(e, collapse = " "), call. = FALSE)
      }

      for (vname in names(inp)) {
        if (!all(na.omit(unique(inp[[vname]])) %in% private$.results$dims[[vname]])) {
          e <- c(
            "invalid/unknown values detected in argument `inp`",
            "in column", shQuote(vname)
          )
          stop(paste(e, collapse = " "), call. = FALSE)
        }
      }

      inp <- as.data.frame(inp)
      inp <- inp[, dv, drop = FALSE]
      res <- sdcProb2df(private$.prob, addDups = TRUE, addNumVars = TRUE, dimCodes = "original")
      res <- res[, c("strID", "freq", dv), with = FALSE]
      res$supp <- FALSE
      inp$tmpid <- seq_len(nrow(inp))
      for (i in seq_len(nrow(inp))) {
        xx <- inp[i, , drop = FALSE]
        xx <- xx[, !apply(xx, 2, is.na), drop = FALSE]
        res <- merge(res, xx, by = setdiff(names(xx), "tmpid"), all.x = TRUE)
        res[!is.na(tmpid), supp := TRUE]
        res[, tmpid := NULL]
      }
      res <- setkey(res, strID) # required because merge changes order
      private$.update_supps(v = v, pat = res$supp, rule = "cell-rule")
      invisible(self)
    },

    #' @description Identify sensitive cells based on the p%-rule rule. Please note that this rule
    #' can only be applied to positive-only variables.
    #' @param v a single variable name of a continuous variable (see method `numvars()`)
    #' @param p a number defining a percentage between `1` and `99`.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other methods (e.g `$perturb()`).
    supp_p = function(v, p) {
      .check_avail(
        v = v,
        avail = self$numvars(),
        msg = "Invalid variable specified in `v`:",
        single_v = TRUE)

      if (private$.has_negative_values(v = v)) {
        stop("dominance rule can only be computed on strictly positive variables.", call. = FALSE)
      }

      if (!rlang::is_scalar_integerish(p)) {
        stop("`p` is not an integerish number.", call. = FALSE)
      }
      if (p < 1 | p > 99) {
        stop("`p` must be >= 1 and < 100.", call. = FALSE)
      }

      res <- sdcTable::primarySuppression(
        object = private$.prob,
        type = "p",
        numVarName = v,
        p = p)
      pat <- sdcProb2df(res, addDups = TRUE)[["sdcStatus"]] == "u"
      private$.update_supps(v = v, pat = pat, rule = "p%-rule")
      invisible(self)
    },

    #' @description Identify sensitive cells based on the pq-rule. Please note that this rule
    #' can only be applied to positive-only variables.
    #' @param v a single variable name of a continuous variable (see method `numvars()`)
    #' @param p a number defining a percentage between `1` and `99`.
    #' @param q a number defining a percentage between `1` and `99`. This value must be larger than `p`.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other methods (e.g `$perturb()`).
    supp_pq = function(v, p, q) {
      .check_avail(
        v = v,
        avail = self$numvars(),
        msg = "Invalid variable specified in `v`:",
        single_v = TRUE)
      if (private$.has_negative_values(v = v)) {
        stop("dominance rule can only be computed on strictly positive variables.", call. = FALSE)
      }
      if (!rlang::is_scalar_integerish(p)) {
        stop("`p` is not an integerish number.", call. = FALSE)
      }
      if (p < 1 | p > 99) {
        stop("`p` must be >= 1 and < 100.", call. = FALSE)
      }
      if (!rlang::is_scalar_integerish(q)) {
        stop("`q` is not an integerish number.", call. = FALSE)
      }
      if (q < 1 | q > 99) {
        stop("`q` must be >= 1 and < 100.", call. = FALSE)
      }
      if (p > q) {
        stop("Argument `p` must be < than `q`", call. = FALSE)
      }
      res <- sdcTable::primarySuppression(
        object = private$.prob,
        type = "pq",
        numVarName = v,
        p = p,
        q = q)
      pat <- sdcProb2df(res, addDups = TRUE)[["sdcStatus"]] == "u"
      private$.update_supps(v = v, pat = pat, rule = "pq-rule")
      invisible(self)
    },

    #' @description Identify sensitive cells based on the nk-dominance rule. Please note that this rule
    #' can only be applied to positive-only variables.
    #' @param v a single variable name of a continuous variable (see method `numvars()`)
    #' @param n an integerish number `>= 2`
    #' @param k a number defining a percentage between `1` and `99`. All cells to which the top `n`
    #' contributers contribute more than `k%` is considered unsafe
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other methods (e.g `$perturb()`).
    supp_nk = function(v, n, k) {
      .check_avail(
        v = v,
        avail = self$numvars(),
        msg = "Invalid variable specified in `v`:",
        single_v = TRUE)
      if (private$.has_negative_values(v = v)) {
        stop("dominance rule can only be computed on strictly positive variables.", call. = FALSE)
      }
      if (!rlang::is_scalar_double(n)) {
        stop("`n` is not a number.", call. = FALSE)
      }
      if (n < 2) {
        stop("`n` must be >= 2", call. = FALSE)
      }
      if (!rlang::is_scalar_integerish(k)) {
        stop("`k` is not an integerish number.", call. = FALSE)
      }
      if (k < 1 | k > 99) {
        stop("`k` must be >= 1 and < 100.", call. = FALSE)
      }
      res <- sdcTable::primarySuppression(
        object = private$.prob,
        type = "nk",
        numVarName = v,
        n = n,
        k = k)
      pat <- sdcProb2df(res, addDups = TRUE)[["sdcStatus"]] == "u"
      private$.update_supps(v = v, pat = pat, rule = "nk-rule")
      invisible(self)
    },

    #' @description Return perturbation parameters of count variables
    #' @return a named list in which each list-element contains the
    #' active perturbation parameters for the specific count variable
    #' defined by the list-name.
    params_cnts_get = function() {
      return(private$.pert_params$cnts)
    },

    #' @description Set perturbation parameters for count variables
    #' @param val a perturbation object created with [ck_params_cnts()]
    #' @param v a character vector (or `NULL`). If `NULL` (the default),
    #' the perturbation parameters provided in `val` are set for all
    #' count variables; otherwise one may specify the names of
    #' the count variables for which the parameters should be set.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other
    #' methods (e.g `$perturb()`).
    params_cnts_set = function(val, v = NULL) {
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
        .check_avail(
          v = v,
          avail = cv,
          msg = "Invalid variables specified in `v`:",
          single_v = FALSE)
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

    #' @description reset results and parameters for already perturbed
    #' count variables
    #' @param v if `v` equals `NULL` (the default), the results are reset
    #' for all perturbed count variables; otherwise it is possible to specify
    #' the names of already perturbed count variables.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other
    #' methods (e.g `$perturb()` or `$freqtab()`).
    reset_cntvars = function(v = NULL) {
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

    #' @description reset results and parameters for already perturbed
    #' numerical variables
    #' @param v if `v` equals `NULL` (the default), the results are reset for all perturbed
    #' numerical variables; otherwise it is possible to specify the names of already
    #' perturbed continuous variables.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other
    #' methods (e.g `$perturb()` or `$numtab()`).
    reset_numvars = function(v = NULL) {
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

    #' @description reset results and parameters for all already perturbed variables.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other
    #' methods (e.g `$perturb()`, `$freqtab()` or `$numtab()`).
    reset_allvars = function() {
      avail <- c(private$.ck_perturbed_vars(what = "countvars"), private$.ck_perturbed_vars(what = "numvars"))
      if (length(avail) == 0) {
        message("No perturbed variables available.")
        return(invisible(self))
      }
      private$.reset_vars(vars = avail)
    },

    #' @description Return perturbation parameters of continuous variables
    #' @return a named list in which each list-element contains the
    #' active perturbation parameters for the specific continuous variable
    #' defined by the list-name.
    params_nums_get = function() {
      return(private$.pert_params$nums)
    },

    #' @description set perturbation parameters for continuous variables.
    #' @param val a perturbation object created with [ck_params_nums()]
    #' @param v a character vector (or `NULL`); if `NULL` (the default), the
    #' perturbation parameters provided in `val` are set for all continuous
    #' variables; otherwise one may specify the names of the numeric variables for
    #' which the parameters should be set.
    #' @return A modified `cellkey_obj` object in which private slots were
    #' updated for side-effects. These updated values are used by other
    #' methods (e.g `$perturb()`).
    params_nums_set = function(val, v = NULL) {
      nv <- self$numvars()
      if (length(nv) == 0) {
        message("no numeric variables are available!")
        return(invisible(NULL))
      }

      if (!inherits(val, "ck_params")) {
        stop("Please create the input using `ck_params_nums()`", call. = FALSE)
      }

      if (!val$type %in% c("params_m_flex", "params_m_simple")) {
        stop("Please create the input using `ck_params_nums()`", call. = FALSE)
      }

      if (is.null(v)) {
        v <- nv
      } else {
        if (!is.character(v)) {
          stop("Argument `v` needs to be a character vector.", call. = FALSE)
        }
        .check_avail(
          v = v,
          avail = nv,
          msg = "Invalid variables specified in `v`:",
          single_v = FALSE)
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

    #' @description some aggregated summary statistics about perturbed variables
    #' @return invisible `NULL`
    summary = function() {
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

      numvars_available <- TRUE
      if (numvars_available) {
        cli::cat_line(cli::boxx("Utility measures for perturbed numerical variables", padding = 0))
        nv <- private$.ck_perturbed_vars("numvars")
        if (length(nv) == 0) {
          cli::cat_line(cli::symbol$cross, " no numerical variables have been perturbed")
        } else {
          cli::cat_rule("Distribution statistics of perturbations")
          rr <- self$numtab()
          rr$pert_total <- rr$pws - rr$ws
          num_info <- rr[, as.list(get_distr_vals(pert_total)), by = "vname"]
          print(num_info)
        }
      }
      return(invisible(NULL))
      # return(invisible(
      #   list(cnt_info = cnt_info,
      #        cnt_measures = cnt_measures,
      #        num_info = num_info)
      # ))
    },

    #' @description prints information about the current table
    #' @return invisible `NULL`
    print = function() {
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
      for (i in seq_len(nrow(countvars))) {
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
        for (i in seq_len(nrow(numvars))) {
          v <- numvars$vname[i]
          if (numvars$is_perturbed[i]) {
            cli::cat_line(cli::symbol$checkbox_on, " ", shQuote(v), " (perturbed)")
          } else {
            cli::cat_line(cli::symbol$checkbox_off, " ", shQuote(v))
          }
        }
      }
      return(invisible(NULL))
    }#,
    # to delete
    #everything=function() {
    #  list(
    #    prob = private$.prob,
    #    pert_params = private$.pert_params,
    #    results = private$.results,
    #    max_contr = private$.max_contributions)
    #}
  ),
  private = list(
    .prob = NULL,
    .max_contributions = NULL,
    .varsdt = NULL,
    .is_weighted = NULL,
    .pert_params = list(),
    .modifications = list(),
    .results = list(),
    .dupsinfo = NULL,
    .h_info = NULL,
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
    .ck_perturbed_vars = function(what) {
      stopifnot(what %in% c("countvars", "numvars"))
      type <- is_perturbed <- vname <- NULL
      private$.varsdt[type == what & is_perturbed == TRUE, vname]
    },
    # actually perturb count variables
    .ck_perturb_cnts = function(v) {
      lookup_freq <- function(ptab, cellkeys, cellvals) {
        row_nr <- rep(-1, length(cellkeys))
        pert_vals <- rep(0L, length(cellkeys))
        symmetry <- max(ptab$i)
        ptab$ids <- seq_len(nrow(ptab))

        for (d in seq_len(symmetry)) {
          if (d == symmetry) {
            rkind <- cellvals >= d
          } else {
            rkind <- cellvals == d
          }
          if (sum(rkind) > 0) {
            ii <- which(ptab$i == d)
            v <- ptab$ub[ii]
            ids <- ptab$ids[ii]
            ck <- cellkeys[rkind]
            diffs <- ptab$v[ii]

            # row_ind
            rowind <- sapply(1:sum(rkind), function(x) {
              which.max(ck[x] < v)
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
        ptab = params$params$ptable,
        cellkeys = newtab[[col_ck]],
        cellvals = newtab[[col_n]])

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
      pwc <- newtab[[col_puwc]] * newtab[[col_wcavg]]
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
    .ck_perturb_nums = function(v) {
      dig <- .ck_digits()
      if (private$.is_perturbed_numvar(v)) {
        message("Variable ", shQuote(v),  " was already perturbed!")
        return(invisible(self))
      }

      # parameters
      params <- private$.pert_params$nums[[v]]$params
      if (is.null(params)) {
        e <- "Please set perturbation parameters using `$params_nums_set() first.`"
        stop(e, call. = FALSE)
      }

      ck_log("parameters:")
      ck_log("--> type: ", shQuote(params$type))
      ck_log("--> top_k: ", params$top_k)
      ck_log("--> same_key: ", params$same_key)
      ck_log("--> use_zero_rkeys: ", params$use_zero_rkeys)

      # for now; only flex-function is allowed as input!
      mult_params <- params$mult_params
      ck_log("--> epsilon: ", paste(params$mult_params$epsilon, collapse = ", "))
      ck_log("--> mu_c: ", params$mu_c)
      ck_log("--> separation: ", params$separation)
      ck_log("--> m1_sq: ", round(params$m_fixed_sq, digits = dig))
      ck_log("--> separation_point: ", round(params$zs, digit = dig))
      ck_log("--> E: ", round(params$E, digits = dig))
      ck_log("--> even_odd: ", params$even_odd)

      cl <- class(params$mult_params)
      if (cl == "params_m_simple") {
        lookup_type <- "simple"
      } else if (cl == "params_m_flex") {
        lookup_type <- "flex"
      } else {
        stop("invalid perturbation parameters received.", call. = FALSE)
      }
      ck_log("--> multiplicator: ", shQuote(lookup_type))

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
        x_vals <- lapply(max_contr, function(x) {
          list(x = na.omit(x$w_vals[1:params$top_k]), even_odd = x$even_contributors)
        })
      } else if (params$type == "mean") {
        x_vals <- lapply(max_contr, function(x) {
          list(x = x$w_mean, even_odd = x$even_contributors)
        })
      } else if (params$type == "range") {
        x_vals <- lapply(max_contr, function(x) {
          list(x = x$w_spread, even_odd = x$even_contributors)
        })
      } else if (params$type == "sum") {
        eo <- sapply(max_contr, function(x) {
          x$even_contributors
        })
        x_vals <- lapply(names(w_sums), function(x) {
          list(x = w_sums[[x]], even_odd = eo[[x]])
        })
        names(x_vals) <- names(w_sums)
      }

      if (!params$even_odd) {
        # handle the "all"-case
        x_vals <- lapply(x_vals, function(x) {
          x$even_odd <- NA
          x
        })
      }

      # the table with cell_keys (non-zero cellkeys and weighted and unweighted sums)
      newtab <- copy(private$.results[[v]])
      dim_dt <- private$.results[["dims"]]

      ck_log("selecting cell keys based on choice of argument `use_zero_rkey`")
      if (!params$use_zero_rkeys) {
        cellkeys <- newtab[[paste0("ck_nz_", v)]]
      } else {
        cellkeys <- newtab[[paste0("ck_", v)]]
      }
      cellvals <- newtab[[paste0("ws_", v)]]

      # we compute top_k-scrambled cell keys
      # required for argument `same_key`
      ck_log("scrambling cell keys depending on argument `same_key`")
      cellkeys <- parallel::mclapply(cellkeys, function(x) {
        .scramble_ckeys(
          ck = x,
          top_k = params$top_k,
          same_key = params$same_key)
      }, mc.cores = .ck_cores()) # nolint

      # we remove duplicated cells for now and add them later!
      names(cellkeys) <- dim_dt[["strID"]]
      index_nondup <- !private$.dupsinfo$is_bogus
      cellkeys <- cellkeys[index_nondup]
      cellvals <- cellvals[index_nondup]
      prot_req <- private$.results[[v]]$special_protection[index_nondup]

      lookup <- parallel::mclapply(x_vals, function(x) {
        if (is.na(x$even_odd)) {
          return(rep("all", length(x$x)))
        } else if (isTRUE(x$even_odd)) {
          return("even")
        }
        return("odd")
      }, mc.cores = .ck_cores()) # nolint

      if (lookup_type == "flex") {
        fun <- .perturb_cell_flex
      } else if (lookup_type == "simple") {
        fun <- .perturb_cell_simple
      }

      res <- parallel::mclapply(seq_len(length(x_vals)), function(y) {
        fun(
          cv = cellvals[y],
          x = x_vals[[y]]$x,
          ck = cellkeys[[y]],
          lookup = lookup[[y]],
          prot_req = prot_req[y],
          params = params)
      }, mc.cores = .ck_cores()) # nolint

      pert_result <- rbindlist(parallel::mclapply(seq_len(length(x_vals)), function(x) {
       data.table(
         cell_id = names(x_vals)[x],
         ckey  = cellkeys[[x]],
         cv = res[[x]]$cv,
         cv_pert = res[[x]]$cv_p,
         pert = sum(res[[x]]$x_hats),
         x_delta = res[[x]]$x_delta,
         lookup = res[[x]]$lookup,
         additional_protection = prot_req[x])
      }, mc.cores = .ck_cores())) # nolint

      ck_log("add rows for duplicated cells!")
      index_nondup <- !private$.dupsinfo$is_bogus
      pert_result <- merge(
        x = dim_dt[index_nondup],
        y = pert_result,
        by.x = "strID",
        by.y = "cell_id",
        all.x = TRUE)

      ck_log("updating results")
      newtab <- private$.results[[v]]

      # no duplicates
      # nolint start
      unq <- pert_result[!duplicated(strID), list(strID, pert, cv_pert)]
      # nolint end
      setnames(unq, c("strID", "pert_total", gen_vnames(v, prefix = "pws")))

      # we merge the results back; duplicated cells have identical
      # perturbation values (pws_{var} and pert_total)
      # thus, we can easily merge
      newtab[["strID"]] <- dim_dt[["strID"]]
      newtab <- merge(newtab, unq, by = "strID", all.x = TRUE)
      newtab[["strID"]] <- NULL
      private$.results[[v]] <- newtab

      ck_log("compute and add modifications")
      # note: `ckey`: the first cell key (all others were scrambled)
      mods <- copy(pert_result)
      mods[["strID"]] <- mods$cv <- mods$cv_pert <- NULL
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
    .ck_utility_cnts = function(v, exclude_zeros) {
      v <- tolower(v)
      dt <- private$.results[[v]]
      return(ck_cnt_measures(
        orig = dt$uwc,
        pert = dt$puwc,
        exclude_zeros = exclude_zeros
      ))
    },
    # utility measures for numeric variables
    .ck_utility_nums = function(v) {
      stop(".ck_utility_nums() needs to be implemented", call. = FALSE)
    },
    # is `v` a perturbed countvar?
    .is_perturbed_countvar = function(v) {
      stopifnot(is_scalar_character(v))
      v %in% private$.ck_perturbed_vars("countvars")
    },
    # is `v` a perturbed numvar?
    .is_perturbed_numvar = function(v) {
      stopifnot(is_scalar_character(v))
      v %in% private$.ck_perturbed_vars("numvars")
    },
    # resets perturbed variables
    .reset_vars = function(vars) {
      for (v in vars) {
        ck_log("removing perturbation results and parameters for ", shQuote(v))
        if (private$.is_perturbed_countvar(v)) {
          private$.results[[v]] <- private$.results[[v]][, 1:4]
          private$.modifications$cnts <- private$.modifications$cnts[countvar != v]
          private$.pert_params$cnts[[v]] <- NULL
        } else {
          private$.results[[v]] <- private$.results[[v]][, 1:5]
          private$.modifications$nums <- private$.modifications$nums[numvar != v]
          private$.pert_params$nums[[v]] <- NULL
        }
        private$.varsdt[vname == v, is_perturbed := FALSE]
      }
    },

    # updates column `special_protection` after applying primary suppression
    .update_supps = function(v, pat, rule) {
      nr_ex_supp <- sum(private$.results[[v]]$special_protection)
      index_new_supps <- which(pat)
      if (length(index_new_supps) > 0) {
        private$.results[[v]]$special_protection[index_new_supps] <- TRUE
      }
      nr_new_supps <- sum(private$.results[[v]]$special_protection)
      s <- nr_new_supps - nr_ex_supp
      message(rule, ": ", s, " new sensitive cells (incl. duplicates) found (total: ", nr_new_supps, ")")
    },
    # returns modification slot

    # returns TRUE if we have a variable containing negative values
    .has_negative_values = function(v) {
      min(private$.prob@dataObj@rawData[[v]], na.rm = TRUE) < 0
    },

    modifications = function(type = NULL) {
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
#' @return A new `cellkey_obj` object. Such objects (internally) contain the fully computed
#' statistical tables given input microdata (`x`), the hierarchical definitionals (`dims`) as
#' well as the remaining inputs. Intermediate results are stored internally and can only be
#' modified / accessed via the exported public methods described below.
#' @export
#' @examples
#' \donttest{
#' x <- ck_create_testdata()
#'
#' # create some 0/1 variables that should be perturbed later
#' x[, cnt_females := ifelse(sex == "male", 0, 1)]
#' x[, cnt_males := ifelse(sex == "male", 1, 0)]
#' x[, cnt_highincome := ifelse(income >= 9000, 1, 0)]
#' # a variable with positive and negative contributions
#' x[, mixed := sample(-10:10, nrow(x), replace = TRUE)]
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
#' numvars <- c("expend", "income", "savings", "mixed")
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
#' # information about the hierarchies
#' tab$hierarchy_info()
#'
#' # which variables have been defined?
#' tab$allvars()
#'
#' # count variables
#' tab$cntvars()
#'
#' # continuous variables
#' tab$numvars()
#'
#' # create perturbation parameters for "total" variable and
#' # write to yaml-file
#'
#' # create a ptable using functionality from the ptable-pkg
#' f_yaml <- tempfile(fileext = ".yaml")
#' p_cnts1 <- ck_params_cnts(
#'   ptab = ptable::pt_ex_cnts(),
#'   path = f_yaml)
#'
#' # read parameters from yaml-file and set them for variable `"total"`
#' p_cnts1 <- ck_read_yaml(path = f_yaml)
#'
#' tab$params_cnts_set(val = p_cnts1, v = "total")
#'
#' # create alternative perturbation parameters by specifying parameters
#' para2 <- ptable::create_cnt_ptable(
#'   D = 8, V = 3, js = 2, create = FALSE)
#'
#' p_cnts2 <- ck_params_cnts(ptab = para2)
#'
#' # use these ptable it for the remaining variables
#' tab$params_cnts_set(val = p_cnts2, v = countvars)
#'
#' # perturb a variable
#' tab$perturb(v = "total")
#'
#' # multiple variables can be perturbed as well
#' tab$perturb(v = c("cnt_males", "cnt_highincome"))
#'
#' # return weighted and unweighted results
#' tab$freqtab(v = c("total", "cnt_males"))
#'
#' # numerical variables (positive variables using flex-function)
#' # we also write the config to a yaml file
#' f_yaml <- tempfile(fileext = ".yaml")
#'
#' # create a ptable using functionality from the ptable-pkg
#' # a single ptable for all cells
#' ptab1 <- ptable::pt_ex_nums(parity = TRUE, separation = FALSE)
#'
#' # a single ptab for all cells except for very small ones
#' ptab2 <- ptable::pt_ex_nums(parity = TRUE, separation = TRUE)
#'
#' # different ptables for cells with even/odd number of contributors
#' # and very small cells
#' ptab3 <- ptable::pt_ex_nums(parity = FALSE, separation = TRUE)
#'
#' p_nums1 <- ck_params_nums(
#'   ptab = ptab1,
#'   type = "top_contr",
#'   top_k = 3,
#'   mult_params = ck_flexparams(
#'     fp = 1000,
#'     p = c(0.30, 0.03),
#'     epsilon = c(1, 0.5, 0.2),
#'     q = 3),
#'   mu_c = 2,
#'   same_key = FALSE,
#'   use_zero_rkeys = FALSE,
#'   path = f_yaml)
#'
#' # we read the parameters from the yaml-file
#' p_nums1 <- ck_read_yaml(path = f_yaml)
#'
#' # for variables with positive and negative values
#' p_nums2 <- ck_params_nums(
#'   ptab = ptab2,
#'   type = "top_contr",
#'   top_k = 3,
#'   mult_params = ck_flexparams(
#'     fp = 1000,
#'     p = c(0.15, 0.02),
#'     epsilon = c(1, 0.4, 0.15),
#'     q = 3),
#'   mu_c = 2,
#'   same_key = FALSE)
#'
#' # simple perturbation parameters (not using the flex-function approach)
#' p_nums3 <- ck_params_nums(
#'   ptab = ptab3,
#'   type = "mean",
#'   mult_params = ck_simpleparams(p = 0.25),
#'   mu_c = 2,
#'   same_key = FALSE)
#'
#' # use `p_nums1` for all variables
#' tab$params_nums_set(p_nums1, c("savings", "income", "expend"))
#'
#' # use different parameters for variable `mixed`
#' tab$params_nums_set(p_nums2, v = "mixed")
#'
#' # identify sensitive cells to which extra protection (`mu_c`) is added.
#' tab$supp_p(v = "income", p = 85)
#' tab$supp_pq(v = "income", p = 85, q = 90)
#' tab$supp_nk(v = "income", n = 2, k = 90)
#' tab$supp_freq(v = "income", n = 14, weighted = FALSE)
#' tab$supp_val(v = "income", n = 10000, weighted = TRUE)
#' tab$supp_cells(
#'   v = "income",
#'   inp = data.frame(
#'     sex = c("female", "female"),
#'     "age" = c("age_group1", "age_group3")
#'   )
#' )
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
#'
#' tab$perturb(v = "cnt_males")
#'
#' # write results to a .csv file
#' tab$freqtab(
#'   v = c("total", "cnt_males"),
#'   path = file.path(tempdir(), "outtab.csv")
#' )
#'
#' # show results containing weighted and unweighted results
#' tab$freqtab(v = c("total", "cnt_males"))
#'
#' # utility measures for a count variable
#' tab$measures_cnts(v = "total", exclude_zeros = TRUE)
#'
#' # modifications for perturbed count variables
#' tab$mod_cnts()
#'
#' # display a summary about utility measures
#' tab$summary()
#' }
ck_setup <- function(x, rkey, dims, w = NULL, countvars = NULL, numvars = NULL) {
  ck_class$new(
    x = x,
    rkey = rkey,
    dims = dims,
    w = w,
    countvars = countvars,
    numvars = numvars)
}

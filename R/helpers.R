# generate standardized variable names
gen_vnames <- function(vv, prefix) {
  return(paste0(prefix, "_", vv))
}

# simple check functions for record keys
check_rkeys <- function(rkeys) {
  if (!is.numeric(rkeys)) {
    stop("`rkeys` must be a numeric vector.", call. = FALSE)
  }
  if (!all(rkeys >= 0)) {
    stop("`rkeys` must be a >= 0.", call. = FALSE)
  }
  if (!all(rkeys <= 1)) {
    stop("`rkeys` must be a <= 1.", call. = FALSE)
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
    quantile(dd, c(seq(60, 90, by = 10), 95, 99) / 100), max(dd))
  vals <- round(vals, digits = 3)
  names(vals)[1] <- "Min"
  names(vals)[2:5] <- paste0("Q", seq(10, 40, by = 10))
  names(vals)[6] <- "Mean"
  names(vals)[7] <- "Median"
  names(vals)[8:13] <- paste0("Q", c(seq(60, 90, by = 10), 95, 99))
  names(vals)[14] <- "Max"
  vals
}

# temporary variable names
.tmpvarname <- function(type, what = "tabulation") {
  paste0("tmp", type, "for", what, collapse = "")
}

# checks if all values in `v` exist in `avail` and
# returns a nice error message in case this is not true
.check_avail <- function(v, avail, msg = "Invalid variables specified:", single_v = FALSE) {
  if (single_v) {
    if (!rlang::is_scalar_character(v)) {
      stop("Argument `v` needs to be a scalar character", call. = FALSE)
    }
  } else {
    if (!is_character(v)) {
      stop("Argument `v` needs to be a character vector.", call. = FALSE)
    }
  }

  if (!all(v %in% avail)) {
    e <- c(
      msg,
      "Possible choices are:\n",
      paste(shQuote(avail), collapse = "; ")
    )
    stop(paste(e, collapse = " "), call. = FALSE)
  }
  invisible(TRUE)
}

# debugging
.ck_debug <- function() {
  Sys.getenv("cellkey_debug") == TRUE
}
.ck_debug_on <- function() {
  Sys.setenv("cellkey_debug" = TRUE)
}
.ck_debug_off <- function() {
  Sys.setenv("cellkey_debug" = FALSE)
}
ck_log <- function(..., br = TRUE) {
  if (.ck_debug()) {
    message(...,  appendLF = br)
  }
}

# checks if file has a given extension and optionally if it does not exist
.valid_path <- function(path, ext = "yaml", check_exists = TRUE) {
  if (!rlang::is_scalar_character(path)) {
    stop("`path` is not a scalar character.", call. = FALSE)
  }

  if (tools::file_ext(path) != ext) {
    stop("file-exension of argument `path` is not ", shQuote(ext), call. = FALSE)
  }

  if (check_exists == TRUE) {
    if (file.exists(path)) {
      stop("the provided file ", shQuote(path), " already exists.", call. = FALSE)
    }
  }
  invisible(NULL)
}

# check validity of ptab-inputs from ptable-pkg
.chk_ptab <- function(ptab, type = "nums") {
  .is_ptab <- function(x, arg = "ptab") {
    if (!inherits(x, "ptable")) {
      e <- paste("argument", shQuote(arg), "was not created using ptable::pt_create_pTable()")
      stop(e, call. = FALSE)
    }
  }

  .is_cntptab <- function(x, arg = "ptab") {
    .is_ptab(x, arg = arg)
    if (slot(x, "table") != "cnts") {
      e <- paste("argument", arg, "is not a perturbation table suitable for counts.")
      stop(e, call. = FALSE)
    }
  }

  .is_numptab <- function(x, arg = "ptab") {
    .is_ptab(x, arg = arg)
    if (slot(x, "table") != "nums") {
      e <- paste("argument", arg, "is not a perturbation table suitable for numeric variables.")
      stop(e, call. = FALSE)
    }
  }

  if (is.list(ptab)) {
    cn <- names(ptab)
    if (!all(cn %in% c("all", "even", "odd", "small_cells"))) {
      stop("invalid names in input list of perturbation tables.", call. = FALSE)
    }
    if ("all" %in% cn) {
      .is_numptab(ptab$all, arg = "ptab$all")
      ptab_final <- ptab$all@pTable
    } else {
      .is_numptab(ptab$even, arg = "ptab$even")
      .is_numptab(ptab$even, arg = "ptab$odd")
      ptab_final <- rbind(ptab$even@pTable, ptab$odd@pTable)
    }

    if ("small_cells" %in% cn) {
      .is_numptab(ptab$small_cells, arg = "ptab$small_cells")
      pt_sc <- ptab$small_cells@pTable
      pt_sc$type <- "small_cells"
      ptab_final <- rbind(ptab_final, pt_sc)
    }
    ptab <- ptab_final
  } else {
    .is_ptab(ptab, arg = "ptab")
    if (slot(ptab, "table") != type) {
      if (type == "nums") {
        .is_numptab(ptab, arg = "ptab")
      } else {
        .is_cntptab(ptab, arg = "ptab")
      }
    }
    ptab <- ptab@pTable
  }

  out <- data.table::copy(ptab)
  data.table::setnames(out, "p_int_lb", "lb", skip_absent = TRUE)
  data.table::setnames(out, "p_int_ub", "ub", skip_absent = TRUE)
  out
}

# generate standardized variable names
gen_vnames <- function(vv, prefix) {
  return(paste0(prefix, "_", vv))
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
check_rkeys <- function(rkeys) {
  if (!is.numeric(rkeys)) {
    stop("`rkeys` must be a numeric vector.", call. = FALSE)
  }
  if (!all(rkeys >=0)) {
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
      paste(shQuote(avail), collapse ="; ")
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

# check validity of ptab-inputs from ptable-pkg
.chk_ptab <- function(ptab, type = "nums") {
  if (!inherits(ptab, "ptable")) {
    stop("argument `ptab` was not created using ptable::pt_create_pTable()", call. = FALSE)
  }

  if (slot(ptab, "table") != type) {
    if (type == "nums") {
      e <- "argument `ptab` is not a perturbation table suitable for numeric variables."
    } else {
      e <- "argument `ptab` is not a perturbation table suitable for counts."
    }
    stop(e, call. = FALSE)
  }

  ptab <- ptab@pTable
  if (type == "nums") {
    setnames(ptab, c("i", "j", "v", "p", "lb", "ub", "type"))
  }
  if (type == "cnts") {
    setnames(ptab, c("i", "j", "p", "v", "lb", "ub", "type"))
  }
  ptab
}

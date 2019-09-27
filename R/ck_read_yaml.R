#' Read perturbation parameters from yaml-files
#'
#' [ck_read_yaml()] allows to create perturbation parameter inputs from yaml-files
#' that were previously created using [ck_params_cnts()] or [ck_params_nums()].
#'
#' @param path a path to a yaml-input file
#'
#' @return an object object suitable as input to method `$params_nums_set()` for the perturbation
#' of continous variables in case `path` was created using [ck_params_nums()] or an object
#' suitable as input for `$params_cnts_set()` for the perturbation
#' of counts and frequencies if the input file was generated using [ck_params_cnts()].
#' @export
#' @md
#' @inherit cellkey_pkg examples
ck_read_yaml <- function(path) {
  .yaml_cnts <- function(cfg) {
    cfg$params$ptable <- as.data.table(cfg$params$ptable, stringsAsFactors = FALSE)
    cfg
  }
  .yaml_nums <- function(cfg) {
    cfg$params$ptab <- as.data.table(cfg$params$ptab, stringsAsFactors = FALSE)
    class(cfg$params$mult_params) <- cfg$type
    cfg
  }

  cfg <- yaml::yaml.load_file(input = path)
  type <- cfg$ptype
  cfg$ptype <- NULL

  cur_v <- paste0(utils::packageVersion("cellKey"), collapse = ".")
  if (cur_v != cfg$version) {
    warning("yaml-file ", shQuote(path), " was written with a different version of the package.", call. = FALSE)
  }
  cfg$version <- NULL

  if (type == "params_cnts") {
    cfg <- .yaml_cnts(cfg)
  } else {
    cfg <- .yaml_nums(cfg)
  }
  class(cfg) <- "ck_params"
  cfg
}

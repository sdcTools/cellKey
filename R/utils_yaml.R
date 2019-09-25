# checks if file to which yaml should be written is valid!
.yaml_valid_path <- function(path) {
  if (!rlang::is_scalar_character(path)) {
    stop("`path` is not a scalar character.", call. = FALSE)
  }

  if (tools::file_ext(path) != "yaml") {
    stop("file-exension of argument `path` is not ", shQuote("yaml"), call. = FALSE)
  }

  if (file.exists(path)) {
    stop("the provided file ", shQuote(path), " already exists.", call. = FALSE)
  }
  invisible(NULL)
}

.yaml_write <- function(x, path) {
  .yaml_valid_path(path)
  yaml::write_yaml(
    x = x,
    file = path,
    fileEncoding = "UTF-8",
    precision = .ck_digits() + 1)
}

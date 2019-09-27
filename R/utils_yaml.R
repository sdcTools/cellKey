.yaml_write <- function(x, path) {
  .valid_path(path = path, ext = "yaml")
  yaml::write_yaml(
    x = x,
    file = path,
    fileEncoding = "UTF-8",
    precision = .ck_digits() + 1)
}

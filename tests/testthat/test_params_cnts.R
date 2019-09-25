if (R.version$major >= 3 & R.version$minor >= 6) {
  suppressWarnings(set.seed(120, sample.kind = "Rounding"))
} else {
  set.seed(120)
}

f_yaml <- tempfile(fileext = ".yaml")

context("Testing perturbation parameters for cnt-vars")
test_that("invalid inputs are caught", {
  # path is non-character
  expect_error(
    ck_params_cnts(
      D = 5,
      V = 3,
      js = 2,
      pstay = 0.5,
      optim = 1,
      mono = TRUE,
      path = 5))

  # path is a scalar
  expect_error(
    ck_params_cnts(
      D = 5,
      V = 3,
      js = 2,
      pstay = 0.5,
      optim = 1,
      mono = TRUE,
      path = c("bla1", "bla2")))

  # file-ext is not "yaml"
  expect_error(
    ck_params_cnts(
      D = 5,
      V = 3,
      js = 2,
      pstay = 0.5,
      optim = 1,
      mono = TRUE,
      path = "bla1"))

  # file already exists
  cat("x", file = f_yaml)
  expect_error(
    ck_params_cnts(
      D = 5,
      V = 3,
      js = 2,
      pstay = 0.5,
      optim = 1,
      mono = TRUE,
      path = f_yaml))
  file.remove(f_yaml)
})

f_yaml <- tempfile(fileext = ".yaml")

suppressMessages(para <- ck_params_cnts(
  D = 5,
  V = 3,
  js = 2,
  pstay = 0.5,
  optim = 1,
  mono = TRUE,
  path = f_yaml))
para_yaml <- ck_read_yaml(f_yaml)


context("testing perturbation parameters for cnt-vars")
test_that("results from yaml match original results", {
  expect_equal(para, para_yaml, check.attributes = FALSE)
})

file.remove(f_yaml)

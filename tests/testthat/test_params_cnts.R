if (R.version$major >= 3 & R.version$minor >= 6) {
  suppressWarnings(set.seed(120, sample.kind = "Rounding"))
} else {
  set.seed(120)
}


para <- ptable::pt_create_pParams(
  D = 5,
  V = 3,
  js = 2,
  pstay = 0.5,
  optim = 1,
  mono = TRUE)
ptab <- ptable::pt_create_pTable(para)


test_that("ptable params can be used too", {
  p1 <- ptable::pt_create_pParams(D = 5, V = 2, table = "cnts")
  p2 <- ptable::pt_create_pTable(p1)
  expect_identical(
    ck_params_cnts(ptab = p1),
    ck_params_cnts(ptab = p2))
})


f_yaml <- tempfile(fileext = ".yaml")

context("Testing perturbation parameters for cnt-vars")
test_that("invalid inputs are caught", {
  # path is non-character
  expect_error(ck_params_cnts(ptab = ptab, path = 5))

  # path is not a scalar
  expect_error(ck_params_cnts(ptab = ptab, path = c("a", "b")))

  # file-ext is not "yaml"
  expect_error(ck_params_cnts(ptab = ptab, path = "bla"))

  # file already exists
  cat("x", file = f_yaml)
  expect_error(ck_params_cnts(ptab = ptab, path = f_yaml))
  file.remove(f_yaml)
})

f_yaml <- tempfile(fileext = ".yaml")

suppressMessages(para <- ck_params_cnts(
  ptab = ptab,
  path = f_yaml))
para_yaml <- ck_read_yaml(f_yaml)

context("testing perturbation parameters for cnt-vars")
test_that("results from yaml match original results", {
  expect_equal(para, para_yaml, check.attributes = FALSE)
})

file.remove(f_yaml)

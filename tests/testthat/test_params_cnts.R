set.seed(120, sample.kind = "Reject")

para <- ptable::create_cnt_ptable(
  D = 5,
  V = 3,
  js = 2,
  pstay = 0.5,
  optim = 1,
  mono = TRUE,
  create = FALSE)
ptab <- ptable::create_ptable(params = para)


test_that("ptable params can be used too", {
  p1 <- ptable::create_cnt_ptable(D = 5, V = 2, create = FALSE)
  p2 <- ptable::create_ptable(params = p1)
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

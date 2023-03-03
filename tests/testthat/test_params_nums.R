set.seed(120, sample.kind = "Reject")
f_yaml <- tempfile(fileext = ".yaml")

context("Testing perturbation parameters for num-vars")
test_that("invalid inputs are caught in ck_flexparams()", {
  expect_error(ck_flexparams(fp = "1"))
  expect_error(ck_flexparams(fp = 1:2))
  expect_error(ck_flexparams(fp = -100))
  expect_error(ck_flexparams(fp = 100, p = "0.5"))
  expect_error(ck_flexparams(fp = 100, p = 0.5))
  expect_error(ck_flexparams(fp = 100, p = c(0.2, 0.5)))
  expect_error(ck_flexparams(fp = 100, p = c(0.5, -1)))
  expect_error(ck_flexparams(fp = 100, p = c(0.25, 0.05), q = "x"))
  expect_error(ck_flexparams(fp = 100, p = c(0.25, 0.05), q = -1))
  expect_error(ck_flexparams(fp = 100, p = c(0.25, 0.05), q = 2, epsilon = "x"))
  expect_error(ck_flexparams(fp = 100, p = c(0.25, 0.05), q = 2, epsilon = -1))
  expect_error(ck_flexparams(fp = 100, p = c(0.25, 0.05), q = 2, epsilon = 2))
  expect_error(ck_flexparams(fp = 100, p = c(0.25, 0.05), q = 2, epsilon = 0.5))
  expect_error(ck_flexparams(fp = 100, p = c(0.25, 0.05), q = 2, epsilon = c(1, 0.5, 0.7)))
})

test_that("invalid inputs are caught in ck_simple()", {
  # fp must be a number
  expect_error(ck_simpleparams(p = "1"))
  expect_error(ck_simpleparams(p = -1))
  expect_error(ck_simpleparams(p = 2))

  expect_error(ck_simpleparams(p = 0.25, epsilon = "1"))
  expect_error(ck_simpleparams(p = 0.25, epsilon = -2))
  expect_error(ck_simpleparams(p = 0.25, epsilon = 2))
  expect_error(ck_simpleparams(p = 0.25, epsilon = 0.5))
  expect_error(ck_simpleparams(p = 0.25, epsilon = c(1, 0.5, 0.7)))
})

ep <- c(1, 0.5, 0.3)
flex <-  ck_flexparams(
  fp = 1000,
  p = c(0.20, 0.03),
  epsilon = ep,
  q = 2)

# perturbation parameters for continuous variables
flex_single <- flex
flex_single$epsilon <- 1

# simple
simple <- ck_simpleparams(p = 0.25, epsilon = c(1, 0.5, 0.2))

ptab <- ptable::pt_ex_nums(parity = TRUE, separation = FALSE)

test_that("invalid inputs are caught", {
  # ptab is not valid
  expect_error(ck_params_nums(ptab = 1, type = "top_contr", top_k = 3, mult_params = flex))

  # type is not an character scalar
  expect_error(ck_params_nums(ptab = ptab, type = 5, top_k = 3, mult_params = flex))

  # invalid value in `type`
  expect_error(ck_params_nums(ptab = ptab, type = "invalid", top_k = 3, mult_params = flex))

  # mu is not a number
  expect_error(ck_params_nums(ptab = ptab, type = "top_contr", top_k = 3, mult_params = flex, mu_c = "x"))

  # mu is >= 0
  expect_error(ck_params_nums(ptab = ptab, type = "top_contr", top_k = 3, mult_params = flex, mu_c = -6))

  # same_key is not logical
  expect_error(ck_params_nums(ptab = ptab, type = "top_contr", top_k = 3, mult_params = flex, same_key = 1))

  # use_zero_rkeys is not logical
  expect_error(ck_params_nums(ptab = ptab, type = "top_contr", top_k = 3, mult_params = flex, use_zero_rkeys = 1))

  # additional_checks for type == "top_contr"
  # top_k must be provided
  expect_error(ck_params_nums(ptab = ptab, type = "top_contr", mult_params = flex))

  # top_k must be >= 1
  expect_error(ck_params_nums(ptab = ptab, type = "top_contr", top_k = 0, mult_params = flex))
  # top_k must be <= 6

  flex$epsilon <- c(1, seq(0.5, 0.2, length = 6))
  expect_error(ck_params_nums(ptab = ptab, type = "top_contr", top_k = 7, mult_params = flex))
  flex$epsilon <- ep

  # argument mult_params has invalid class
  expect_error(ck_params_nums(ptab = ptab, type = "top_contr", top_k = 3, mult_params = 1))

  # argument epsilon has different length from top_k
  expect_error(ck_params_nums(ptab = ptab, type = "top_contr", top_k = 4, mult_params = flex))

  # converting top_k to 1
  expect_message(
    ck_params_nums(ptab = ptab, type = "mean", top_k = 4, mult_params = flex_single),
    "ignoring argument `top_k`")
})

set.seed(120, sample.kind = "Reject")

f_yaml <- tempfile(fileext = ".yaml")
p <- ck_params_nums(
  type = "top_contr",
  ptab = ptab,
  top_k = 3,
  mult_params = flex,
  mu_c = 2,
  same_key = FALSE,
  path = f_yaml)
p_yaml <- ck_read_yaml(path = f_yaml)

test_that("top_contr works with flex_params", {
  expect_equal(p, p_yaml, check.attributes = FALSE)
  expect_identical(p$type, "params_m_flex")
  expect_identical(p$type, "params_m_flex")
  expect_identical(p$params$type, "top_contr")
  expect_identical(p$params$top_k, 3)
  expect_identical(p$params$mu_c, 2)
  expect_equal(p$params$m_fixed_sq, NA)
  expect_equal(p$params$zs, 0)
  expect_identical(digest::sha1(unlist(p$params)), "4ae15d65ea7b68f60ea250382b1b9e6131579f0c")
})
file.remove(f_yaml)

set.seed(120, sample.kind = "Reject")
f_yaml <- tempfile(fileext = ".yaml")
p <- ck_params_nums(
  type = "top_contr",
  ptab = ptab,
  top_k = 1,
  mult_params = flex_single,
  mu_c = 2.5,
  same_key = FALSE,
  path = f_yaml)
p_yaml <- ck_read_yaml(path = f_yaml)
test_that("top_contr works with single flex", {
  expect_equal(p, p_yaml, check.attributes = FALSE)
  expect_identical(p$type, "params_m_flex")
  expect_identical(p$params$type, "top_contr")
  expect_identical(p$params$top_k, 1)
  expect_identical(p$params$mu_c, 2.5)
  expect_equal(p$params$m_fixed_sq, NA)
  expect_equal(p$params$zs, 0)
  expect_identical(digest::sha1(unlist(p$params)), "eeac40b788da70c615be49a569cb7e9515572e8d")
})
file.remove(f_yaml)

set.seed(120, sample.kind = "Reject")
simple$epsilon <- 1
f_yaml <- tempfile(fileext = ".yaml")
p <- ck_params_nums(
  type = "mean",
  ptab = ptab,
  mult_params = simple,
  mu_c = 2.5,
  same_key = FALSE,
  path = f_yaml)
p_yaml <- ck_read_yaml(path = f_yaml)
test_that("top_contr works with simple", {
  expect_equal(p, p_yaml, check.attributes = FALSE)
  expect_identical(p$type, "params_m_simple")
  expect_identical(p$params$type, "mean")
  expect_identical(p$params$top_k, 1)
  expect_identical(p$params$mu_c, 2.5)
  expect_equal(p$params$m_fixed_sq, NA)
  expect_equal(p$params$zs, 0)
  expect_identical(digest::sha1(unlist(p$params)), "998097435934b5b33cb3113bb33caee84d6e78a9")
})
file.remove(f_yaml)

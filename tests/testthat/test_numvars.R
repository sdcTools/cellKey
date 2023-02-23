context("Testing Magnitude Tables")

set.seed(120, sample.kind = "Reject")

x <- ck_create_testdata()
# create some 0/1 variables that should be perturbed later
x[, cnt_females := ifelse(sex == "male", 0, 1)]
x[, cnt_males := ifelse(sex == "male", 1, 0)]
x[, cnt_highincome := ifelse(income >= 9000, 1, 0)]
x[, mixed := sample(-20:10, nrow(x), replace = TRUE)] # mixed variable with positive and negative contributions

x$id <- seq_len(nrow(x))
x$income[1297] <- x$income[1297] * 25
x$income[3503] <- x$income[3503] * 25

x$sampling_weight <- sample(1:5, nrow(x), replace = TRUE)

dim_sex <- hier_create(root = "Total", nodes = c("male", "female"))
dim_age <- hier_create(root = "Total", paste0("age_group", 1:6))
dims <- list(sex = dim_sex, age = dim_age)

## test generation of destatis rkeys
x$rkey <- ck_generate_rkeys(dat = x, nr_digits = 8)

test_that("checking dimension and structure of generated testdata is ok", {
  expect_true(is.data.table(x))
  expect_identical(digest::sha1(x), "4c804693f7d573e9dfa24bcc312ba2b61b490ada")
  expect_identical(digest::sha1(dims), "62748837ca3246a33081dd35f50d06334caa3119")
})

## perturbation parameters for count variables
ep <- c(1, 0.5, 0.3)
flex <-  ck_flexparams(
  fp = 1000,
  p = c(0.20, 0.03),
  epsilon = ep,
  q = 2)

test_that("ptable params can be used too", {
  p1 <- ptable::create_num_ptable(
    D = 5,
    V = 2,
    step = 4,
    icat = c(1, 3, 5),
    create = FALSE
  )
  p2 <- ptable::create_ptable(params = p1)

  expect_identical(
    ck_params_nums(
      type = "top_contr",
      top_k = 3,
      ptab = p1,
      mult_params = flex),
    ck_params_nums(
      type = "top_contr",
      top_k = 3,
      ptab = p2,
      mult_params = flex))
})

set.seed(100, sample.kind = "Reject")
ptab <- ptable::pt_ex_nums(
  parity = TRUE,
  separation = FALSE)

p1 <- ck_params_nums(
  type = "top_contr",
  top_k = 3,
  ptab = ptab,
  mult_params = flex,
  mu_c = 2.5,
  same_key = FALSE,
  use_zero_rkeys = TRUE)

p2 <- ck_params_nums(
  type = "mean",
  top_k = 3,
  ptab = ptab,
  mult_params = ck_simpleparams(p = 0.1, epsilon = 1),
  mu_c = 2.5,
  same_key = TRUE,
  use_zero_rkeys = FALSE)

test_that("checking perturbation parameters", {
  expect_is(p1, "ck_params")
  expect_equal(p1$type, "params_m_flex")
  expect_identical(digest::sha1(p1), "53299f7eab1d919060e12311955ed3eb02f0b38b")
})

# set up problem
numvars <- c("savings", "income", "mixed")
tab <- ck_setup(
  x = x,
  rkey = "rkey",
  dims = dims,
  w = "sampling_weight",
  numvars = numvars)

test_that("problem was correctly generated", {
  expect_equal(tab$numvars(), numvars)
})

# set perturbation parameters for all variables
tab$params_nums_set(val = p1, v = "income")
tab$params_nums_set(val = p2, v = "savings")

test_that("parameters were correctly set", {
  expect_equal(tab$params_nums_get()$income, p1)
  expect_equal(tab$params_nums_get()$savings, p2)
})

expect_message(tab$perturb("income"), "Numeric variable 'income' was perturbed.")
expect_message(tab$perturb("savings"), "Numeric variable 'savings' was perturbed.")

test_that("variable was correctly perturbed", {
  expect_equal(digest::sha1(tab$numtab("income", mean_before_sum = FALSE)), "29eec69ec43987831d951b7e6ecd4d423b27f5e2")
  expect_equal(digest::sha1(tab$numtab("savings", mean_before_sum = FALSE)), "77abfc53132ec03b5a523764eefb52e00dca897f")
  expect_equal(digest::sha1(tab$numtab("income", mean_before_sum = TRUE)), "98f4cdbe6b4b362aa0b9ab31337fc39c4d807c1e")
  expect_equal(digest::sha1(tab$numtab("savings", mean_before_sum = TRUE)), "39a0fa3ddfd8b45778549711d9f9f58c843b0912")
})

context("Testing Magnitude Tables")

if (R.version$major >= 3 & R.version$minor >= 6) {
  suppressWarnings(set.seed(120, sample.kind = "Rounding"))
} else {
  set.seed(120)
}

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
  expect_identical(digest::digest(x), "adb7111c9f85652b0744200eb9b8b30e")
  expect_true(is.data.table(x))
})


## perturbation parameters for count variables
ep <- c(1, 0.5, 0.3)
flex <-  ck_flexparams(
  fp = 1000,
  p = c(0.20, 0.03),
  epsilon = ep,
  q = 2)

test_that("ptable params can be used too", {
  p1 <- ptable::pt_create_pParams(D=5, V=2, table="nums", step=4, icat=c(1,3,5))
  p2 <- ptable::pt_create_pTable(p1)

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

set.seed(100)
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
  expect_identical(digest::digest(p1), "0474740f1ef9141738c14496b93ffb74")
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
  expect_equal(digest::sha1(tab$numtab("income", mean_before_sum = FALSE)), "f26f2eb9f88b3f691354c7450e39b7e443aa8352")
  expect_equal(digest::sha1(tab$numtab("savings", mean_before_sum = FALSE)), "21c708fc30cd613797ac67cbb28468676b4db84c")
  expect_equal(digest::sha1(tab$numtab("income", mean_before_sum = TRUE)), "20ff3a0b92dff03a6c4cff1580ad4fcf31dcf3c7")
  expect_equal(digest::sha1(tab$numtab("savings", mean_before_sum = TRUE)), "c243c42a6cff0cb7e44bbb4af9eaed69d543a7a3")
})

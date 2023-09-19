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
  expect_identical(dim(x), c(4580L, 21L))
  expect_identical(round(mean(x$rkey), digits = 3), 0.5)
  expect_identical(round(mean(x$household_weights), digits = 3), 21.834)
  expect_identical(range(x$mixed), c(-20L, 10L))
  expect_identical(sum(x$cnt_males), 2296)
  expect_identical(sum(x$cnt_females), 2284)
  expect_identical(sum(x$cnt_highincome), 445)
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

  pp <- p1$params

  expect_identical(pp$type, "top_contr")
  expect_identical(pp$top_k, 3)
  expect_identical(dim(pp$ptab), c(96L, 7L))
  expect_identical(round(mean(pp$ptab$p), digits = 3), 0.042)

  expect_identical(pp$mu_c, 2.5)
  expect_identical(pp$m_fixed_sq, NA)
  expect_identical(pp$zs, 0)
  expect_identical(pp$E, 1.34)
  expect_identical(pp$mult_params$fp, 1000)
  expect_identical(pp$mult_params$p_small, 0.2)
  expect_identical(pp$mult_params$p_large, 0.03)
  expect_identical(pp$mult_params$epsilon, c(1, 0.5, 0.3))

  expect_identical(pp$same_key, FALSE)
  expect_identical(pp$use_zero_rkeys, TRUE)
  expect_identical(pp$even_odd, FALSE)
  expect_identical(pp$separation, FALSE)
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
  dt <- tab$numtab("income", mean_before_sum = FALSE)
  expect_identical(dim(dt), c(21L, 6L))
  expect_equal(round(mean(dt$pws), digits = 3), 13181967.1)
  expect_equal(round(mean(dt$ws), digits = 3), 13184145.52)

  dt <- tab$numtab("savings", mean_before_sum = FALSE)
  expect_identical(dim(dt), c(21L, 6L))
  expect_equal(round(mean(dt$pws), digits = 3), 1306306.6)
  expect_equal(round(mean(dt$ws), digits = 3), 1306303.048)

  dt <- tab$numtab("income", mean_before_sum = TRUE)
  expect_identical(dim(dt), c(21L, 6L))
  expect_equal(round(mean(dt$pws), digits = 3), 13179895.45)
  expect_equal(round(mean(dt$ws), digits = 3), 13184145.52)

  dt <- tab$numtab("savings", mean_before_sum = TRUE)
  expect_identical(dim(dt), c(21L, 6L))
  expect_equal(round(mean(dt$pws), digits = 3), 1306310.272)
  expect_equal(round(mean(dt$ws), digits = 3), 1306303.048)
})

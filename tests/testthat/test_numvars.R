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

x$id <- 1:nrow(x)
x$income[1297] <- x$income[1297] * 25
x$income[3503] <- x$income[3503] * 25

x$sampling_weight <- sample(1:5, nrow(x), replace = TRUE)


dim_sex <- hier_create(root = "Total", nodes = c("male", "female"))
dim_age <- hier_create(root = "Total", paste0("age_group", 1:6))
dims <- list(sex = dim_sex, age = dim_age)

## test generation of destatis rkeys
x$rkey <- ck_generate_rkeys(dat = x, nr_digits = 8)

test_that("checking dimension and structure of generated testdata is ok", {
  expect_identical(digest::digest(x), "1ad2666d5e0cd33aba567ba1ef0117fe")
  expect_true(is.data.table(x))
})

## perturbation parameters for count variables
set.seed(100)
ep <- c(1, 0.5, 0.3)
flex <-  ck_flexparams(
  fp = 1000,
  p = c(0.20, 0.03),
  epsilon = ep,
  q = 2)
p1 <- ck_params_nums(
  type = "top_contr",
  top_k = 3,
  D = 7,
  l = 0.2,
  mult_params = flex,
  mu_c = 2.5,
  same_key = FALSE,
  separation = TRUE,
  use_zero_rkeys = TRUE,
  parity = FALSE)

p2 <- ck_params_nums(
  type = "mean",
  top_k = 3,
  D = 7,
  l = 0.2,
  mult_params = ck_simpleparams(p = 0.1, epsilon = 1),
  mu_c = 2.5,
  same_key = TRUE,
  separation = TRUE,
  use_zero_rkeys = FALSE,
  parity = FALSE)


test_that("checking perturbation parameters", {
  expect_is(p1, "ck_params")
  expect_equal(p1$type, "params_m_flex")
  expect_identical(digest::digest(p1), "4c61633dadae010488b705e74f043c69")
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
  expect_equal(digest::sha1(tab$numtab("income", mean_before_sum = FALSE)), "95d346e4649fd813d1175c2e952638b6822bbb10")
  expect_equal(digest::sha1(tab$numtab("savings", mean_before_sum = FALSE)), "c96e3b0d1be6e470df6a278c4a52d45eca0c436e")
  expect_equal(digest::sha1(tab$numtab("income", mean_before_sum = TRUE)), "74432ca0c5d8faa55fdaa881513981571e3785e3")
  expect_equal(digest::sha1(tab$numtab("savings", mean_before_sum = TRUE)), "6661a7d300043ea68fb4dca7ee69aa05406a91d2")
})

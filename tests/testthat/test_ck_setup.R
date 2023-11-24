skip_on_cran()
set.seed(120, sample.kind = "Reject")
context("Testing ck_setup()")
x <- ck_create_testdata()
# create some 0/1 variables that should be perturbed later
x[, cnt_males := ifelse(sex == "male", 1, 0)]
x[, mixed := sample(-20:10, nrow(x), replace = TRUE)] # mixed variable with positive and negative contributions
x$sampling_weight <- sample(1:5, nrow(x), replace = TRUE)


dim_sex <- hier_create(root = "Total", nodes = c("male", "female"))
dim_age <- hier_create(root = "Total", paste0("age_group", 1:6))
dims <- list(sex = dim_sex, age = dim_age)

## test generation of destatis rkeys
x$rkey <- ck_generate_rkeys(dat = x, nr_digits = 8)

cv <- "cnt_males"
nv <- c("mixed", "savings")
w <- "sampling_weight"
rk <- "rkey"

test_that("invalid inputs are caught", {
  # input must be a data.frame
  expect_error(ck_setup(x = 1:5, rkey = rk, dims = dims, w = w, countvars = cv, numvars = nv))

  # check argument `rkey`
  expect_error(ck_setup(x = x, rkey = 3, dims = dims, w = w, countvars = cv, numvars = nv))
  expect_error(ck_setup(x = x, rkey = FALSE, dims = dims, w = w, countvars = cv, numvars = nv))

  # invalid record keys
  expect_error(ck_setup(x = x, rkey = "income", dims = dims, w = w, countvars = cv, numvars = nv))

  # checking `dims`
  # invalid input
  expect_error(ck_setup(x = x, rkey = rk, dims = 1:5, w = w, countvars = cv, numvars = nv))
  # variables in dims do not exist
  expect_error(ck_setup(x = x, rkey = rk, dims = list(a = 1, b = 2), w = w, countvars = cv, numvars = nv))
  # invalid input in dims
  expect_error(ck_setup(x = x, rkey = rk, dims = list(sex = 1, age = 2), w = w, countvars = cv, numvars = nv))
  # unnamed
  expect_error(ck_setup(x = x, rkey = rk, dims = list(dims$sex, dims$age), w = w, countvars = cv, numvars = nv))

  # non-existing weight variable
  expect_error(ck_setup(x = x, rkey = rk, dims = dims, w = "nonex", countvars = cv, numvars = nv))
  # non-numeric weights
  expect_error(ck_setup(x = x, rkey = rk, dims = dims, w = "sex", countvars = cv, numvars = nv))

  # countvars
  # countvars is not character
  expect_error(ck_setup(x = x, rkey = rk, dims = dims, w = w, countvars = 5, numvars = nv))
  # non-existing countvars
  expect_error(ck_setup(x = x, rkey = rk, dims = dims, w = w, countvars = "nonex", numvars = nv))
  # non-numeric countvars
  expect_error(ck_setup(x = x, rkey = rk, dims = dims, w = w, countvars = "age", numvars = nv))
  # non-binary countvars
  expect_error(ck_setup(x = x, rkey = rk, dims = dims, w = w, countvars = rk, numvars = nv))

  # numvars
  # numvars is not character
  expect_error(ck_setup(x = x, rkey = rk, dims = dims, w = w, countvars = cv, numvars = 5))
  # non-existing numvars
  expect_error(ck_setup(x = x, rkey = rk, dims = dims, w = w, countvars = cv, numvars = "nonex"))
  # non-numeric numvars
  expect_error(ck_setup(x = x, rkey = rk, dims = dims, w = w, countvars = cv, numvars = "age"))
})

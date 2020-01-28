context("Testing Frequency Tables")

if (R.version$major >= 3 & R.version$minor >= 6) {
  suppressWarnings(set.seed(120, sample.kind = "Rounding"))
} else {
  set.seed(120)
}

dat <- ck_create_testdata()

seed <- ck_create_seed_from_hash(dat)
test_that("check seed from hash generation", {
  expect_equal(seed, 1769853828)
})

dim_sex <- hier_create(root = "Total", nodes = c("male", "female"))
dim_age <- hier_create(root = "Total", paste0("age_group", 1:6))
dims <- list(sex = dim_sex, age = dim_age)

## test generation of destatis rkeys
rk1 <- ck_generate_rkeys(dat = dat, nr_digits = 5)
rk2 <- ck_generate_rkeys(dat = dat, nr_digits = 5)
test_that("check rkey generation and seed is ok", {
  expect_identical(rk1, rk2)
})
dat$rec_key <- rk1
test_that("checking dimension and structure of generated testdata is ok", {
  expect_identical(digest::digest(dat), "1bed2dc15c2b2b06bebc92f259f70c3b")
  expect_true(is.data.table(dat))
})

## perturbation parameters for count variables
set.seed(100)
suppressMessages(para <- ptable::pt_create_pParams(
  D = 5,
  V = 3,
  js = 2,
  pstay = 0.5,
  optim = 1,
  mono = TRUE))

params_cnts <- ck_params_cnts(ptab = ptable::pt_create_pTable(para))
ptab <- params_cnts$params$ptable
test_that("ck_params_cnts() is ok", {
  expect_is(params_cnts, "ck_params")
  expect_identical(params_cnts$type, "cnts")
  expect_identical(dim(ptab), c(66L, 7L))
  expect_identical(max(ptab$i), 8)
})

# we need to load a custom data set because otherwise, the
# digest based checks don't work
test_that("checking perturbation parameters for counts", {
  expect_is(params_cnts, "ck_params")
  expect_equal(params_cnts$type, "cnts")
  expect_is(params_cnts$params$ptable, "data.table")
  expect_identical(digest::digest(params_cnts), "7426d7da34a57b19b625d35413c1abdf")
})

countvars <- NULL
w <- "sampling_weight"
rkey <- "rec_key"

tab <- ck_setup(
  x = dat,
  rkey = "rec_key",
  dims = dims,
  w = w,
  countvars = countvars)

# set perturbation parameters for all variables
tab$params_cnts_set(val = params_cnts, v = NULL)

expect_message(tab$perturb("total"))
expect_message(tab$perturb("total"), "Variable 'total' was already perturbed!")

test_that("check ck_define_table() with already existing rec-keys", {
  expect_is(tab, "cellkey_obj")
  expect_identical(digest::digest(tab$freqtab("total")), "ef86182106ebaa692f91951cfb18f4b4")
})

dat$rec_key <- NULL
tab <- ck_setup(
  x = dat,
  rkey = 7,
  dims = dims,
  w = w,
  countvars = countvars)

expect_error(tab$perturb("total"))
tab$params_cnts_set(val = params_cnts, v = "total")
expect_message(tab$perturb("total"), "Count variable 'total' was perturbed.")
expect_message(tab$perturb("total"), "Variable 'total' was already perturbed!")

test_that("ck_define_table() with new record keys is ok", {
  expect_is(tab, "cellkey_obj")
  expect_identical(digest::digest(tab$freqtab("total")), "ef86182106ebaa692f91951cfb18f4b4")
})

freqtab <- tab$freqtab("total")
test_that("weighted version of ck_perturb() is ok", {
  expect_identical(digest::digest(freqtab), "ef86182106ebaa692f91951cfb18f4b4")
  expect_identical(digest::digest(tab$mod_cnts()), "d4adaecabe6d853dcc5d97780da55b97")
})

mm <- tab$measures_cnts("Total")
test_that("ck_cnt_measures() [exclude_zeros = TRUE] is ok", {
  expect_identical(digest::digest(mm), "6a4db99d1acc7d2554f1874377e84a05")
})

mm <- tab$measures_cnts("Total", exclude_zeros = FALSE)
test_that("ck_cnt_measures() [exclude_zeros = FALSE] is ok", {
  expect_identical(digest::digest(mm), "374f86f2978e1aaedd0bb26558fd75f5")
})

# no weights
tab <- ck_setup(
  x = dat,
  rkey = 7,
  dims = dims,
  w = NULL,
  countvars = countvars)

# set params
tab$params_cnts_set(params_cnts, v = NULL)

tab$perturb("total")
freqtab <- tab$freqtab("total")

test_that("checking unweighted version of perturb()", {
  expect_identical(digest::digest(freqtab), "4ab265bce8720c68acd9ee65b0a3449f")
  expect_identical(digest::digest(tab$mod_cnts()), "d4adaecabe6d853dcc5d97780da55b97")
})

context("Testing multiple countvars")
dat[, cnt_males := ifelse(sex == "male", 1, 0)]
dat[, cnt_highincome := ifelse(income >= 9000, 1, 0)]
countvars <- c("cnt_males", "cnt_highincome")

tab <- ck_setup(
  x = dat,
  rkey = 7,
  dims = dims,
  w = w,
  countvars = countvars)

# set params
tab$params_cnts_set(params_cnts, v = NULL)

tab$perturb(c("total", "cnt_males", "cnt_highincome"))

test_that("check tabulation of cnt_males is ok", {
  expect_identical(digest::digest(tab$freqtab("cnt_males")), "cb80ad1d6f657742177fe8ee87b10a22")
  expect_identical(digest::digest(tab$measures_cnts("cnt_males")), "d889e18708f07357d726b11b1c24ce1a")
})

test_that("check tabulation of cnt_highincome is ok", {
  expect_identical(digest::digest(tab$freqtab("cnt_highincome")), "633bbb4bb61d011dbab2eafb72cca3f0")
  expect_identical(digest::digest(tab$measures_cnts("cnt_highincome")), "b4f138fbe2f213e15efab44e7a82ab18")
})

test_that("check tabulation of multiple count variables is ok", {
  tt <- tab$freqtab(c("total", "cnt_males", "cnt_highincome"))
  expect_identical(digest::digest(tt), "573704a72ae88ef432ce1fbb129afc31")
})

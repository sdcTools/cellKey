context("Testing Frequency Tables")

if (R.version$major >= 3 & R.version$minor >= 6) {
  suppressWarnings(set.seed(120, sample.kind = "Rounding"))
} else {
  set.seed(120)
}

dat <- ck_create_testdata()

seed <- ck_create_seed_from_hash(dat)
test_that("check seed from hash generation", {
  expect_equal(seed, 483520795)
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
  expect_identical(digest::digest(dat), "4d1eefa7133d2532f6ff8c2d728f8374")
  expect_true(is.data.table(dat))
})

## perturbation parameters for count variables
params_cnts <- ck_params_cnts(
  D = 5,
  V = 3,
  js = 2,
  pstay = 0.5,
  optim = 1,
  mono = TRUE)

ptab <- params_cnts$params$ptable
test_that("ck_params_cnts() is ok", {
  expect_is(params_cnts, "ck_params")
  expect_identical(params_cnts$type, "cnts")
  expect_identical(dim(ptab), c(66L, 6L))
  expect_identical(max(ptab$i), 8L)
})

# we need to load a custom data set because otherwise, the
# digest based checks don't work
rm(params_cnts)
data("params_cnts", package = "cellKey")
test_that("checking perturbation parameters for counts", {
  expect_is(params_cnts, "ck_params")
  expect_equal(params_cnts$type, "cnts")
  expect_is(params_cnts$params$ptable, "data.table")
  expect_identical(digest::digest(params_cnts), "1fe9b2228b2b5be3eacfc524b8ea6104")
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
  expect_identical(digest::digest(tab$freqtab("total")), "7bbc19818c87045dc4392795e3cc6d16")
  expect_identical(digest::digest(tab$freqtab("total", type = "weighted")), "f9c380b7db5c450029f106899dd47edb")
  expect_identical(digest::digest(tab$freqtab("total", type = "unweighted")), "335afdf8d8ce54457cd9d268686ba08b")
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
  expect_identical(digest::digest(tab$freqtab("total")), "7bbc19818c87045dc4392795e3cc6d16")
})

freqtab <- tab$freqtab("total")
test_that("weighted version of ck_perturb() is ok", {
  expect_identical(digest::digest(freqtab), "7bbc19818c87045dc4392795e3cc6d16")
  expect_identical(digest::digest(tab$mod_cnts()), "4f6d21a2b4bdeb4d47f5b2d017ee81ca")
})

mm <- tab$measures_cnts("Total")
test_that("ck_cnt_measures() [exclude_zeros = TRUE] is ok", {
  expect_identical(digest::digest(mm), "71603022b3f7401d6b54ba27d47a79c0")
})

mm <- tab$measures_cnts("Total", exclude_zeros = FALSE)
test_that("ck_cnt_measures() [exclude_zeros = FALSE] is ok", {
  expect_identical(digest::digest(mm), "107665adbbaaa01786ed086ee97fb025")
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
  expect_identical(digest::digest(freqtab), "2aa7c41f2066dbab737c8e96ed46db7a")
  expect_identical(digest::digest(tab$mod_cnts()), "4f6d21a2b4bdeb4d47f5b2d017ee81ca")
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
  expect_identical(digest::digest(tab$freqtab("cnt_males")), "a99b7c88e99cdb4147d454d2ad5eee27")
  expect_identical(digest::digest(tab$measures_cnts("cnt_males")), "f25113b4d2a2ab396c7425d444d6ffc1")
})

test_that("check tabulation of cnt_highincome is ok", {
  expect_identical(digest::digest(tab$freqtab("cnt_highincome")), "ec941a87f7cf4566cc51f626f5434187")
  expect_identical(digest::digest(tab$measures_cnts("cnt_highincome")), "43be7d55fdccbcd38325bed49850c406")
  expect_identical(digest::digest(tab$freqtab("cnt_highincome", type = "weighted")), "7d29b5d5005c322ce4d826184b6fc629")
  expect_identical(digest::digest(tab$freqtab("cnt_highincome", type = "unweighted")), "23b1f9a2b517c494a7865de10c1b837d")
})

test_that("check tabulation of multiple count variables is ok", {
  tt <- tab$freqtab(c("total", "cnt_males", "cnt_highincome"))
  expect_identical(digest::digest(tt), "5001e20f3bc14991740b4f18a88c4184")
})

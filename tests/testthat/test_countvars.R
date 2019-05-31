context("Testing Frequency Tables")

set.seed(120)
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
  mono = TRUE
)

ptab <- params_cnts$params$ptable@pTable
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
  expect_is(params_cnts$params$ptable, "ptable")
  expect_identical(digest::digest(params_cnts), "3cf4a3dcde7ce420a9d9e67907678493")
})

countvars <- NULL
w <- "sampling_weight"
numvars <- NULL
params_nums <- NULL
rkey <- "rec_key"

test_that("ck_setup() fails with invalid inputs", {
  # input must be a data.frame
  expect_error(
    tab <- ck_setup(
      x = 1:5, rkey = "rec_key",
      dims = dims, w = w,
      countvars = countvars, numvars = numvars,
      params_cnts = params_cnts, params_nums = params_nums
    )
  )

  # params_cnts must be created with ck_params_*()
  expect_error(
    tab <- ck_setup(
      x = dat, rkey = "rec_key",
      dims = dims, w = w,
      countvars = countvars, numvars = numvars,
      params_cnts = 1:5, params_nums = params_nums
    )
  )
  # params_nums must be created with ck_params_*()
  expect_error(
    tab <- ck_setup(
      x = dat, rkey = "rec_key",
      dims = dims, w = w,
      countvars = countvars, numvars = numvars,
      params_cnts = params_cnts, params_nums = 1:5
    )
  )

  # numvars specified but params_nums is empty
  expect_error(
    tab <- ck_setup(
      x = dat, rkey = "rec_key",
      dims = dims, w = w,
      countvars = countvars, numvars = "income",
      params_cnts = params_cnts, params_nums = NULL
    )
  )
  # invalid record keys
  expect_error(
    tab <- ck_setup(
      x = dat, rkey = "income",
      dims = dims, w = w,
      countvars = countvars, numvars = numvars,
      params_cnts = params_cnts, params_nums = params_nums
    )
  )

  # dims is not named
  expect_error(
    tab <- ck_setup(
      x = dat, rkey = "rec_key",
      dims = 1:5, w = w,
      countvars = countvars, numvars = numvars,
      params_cnts = params_cnts, params_nums = params_nums
    )
  )

  # dims is not a named list
  expect_error(
    tab <- ck_setup(
      x = dat, rkey = "rec_key",
      dims = c("a" = 5, "b" = 3), w = w,
      countvars = countvars, numvars = numvars,
      params_cnts = params_cnts, params_nums = params_nums
    )
  )

  # invalid variable names in `dims`
  expect_error(
    tab <- ck_setup(
      x = dat, rkey = "rec_key",
      dims = list("a" = 5, "b" = 3), w = w,
      countvars = countvars, numvars = numvars,
      params_cnts = params_cnts, params_nums = params_nums
    )
  )

  # invalid input in `dims`
  expect_error(
    tab <- ck_setup(
      x = dat, rkey = "rec_key",
      dims = list("sex" = 5, "age" = 3), w = w,
      countvars = countvars, numvars = numvars,
      params_cnts = params_cnts, params_nums = params_nums
    )
  )
})

tab <- ck_setup(
  x = dat,
  rkey = "rec_key",
  dims = dims,
  w = w,
  countvars = countvars,
  numvars = numvars,
  params_cnts = params_cnts,
  params_nums = params_nums
)

expect_message(tab$perturb("total"))
expect_message(tab$perturb("total"), "Variable 'total' was already perturbed!")
expect_message(tab$freqtab())

test_that("check ck_define_table() with already existing rec-keys", {
  expect_is(tab, "cellkey_obj")
  expect_identical(digest::digest(tab$freqtab("total")), "24801744e36e9843b6affb1b1c939831")
  expect_identical(digest::digest(tab$freqtab("total", type = "weighted")), "faebd04e46fbc2e02945e0cab03bc541")
  expect_identical(digest::digest(tab$freqtab("total", type = "unweighted")), "335afdf8d8ce54457cd9d268686ba08b")
})

dat$rec_key <- NULL
tab <- ck_setup(
  x = dat,
  rkey = 7,
  dims = dims,
  w = w,
  countvars = countvars,
  numvars = numvars,
  params_cnts = params_cnts,
  params_nums = params_nums
)

expect_message(tab$perturb("total"), "Count variable 'total' was perturbed.")
expect_message(tab$perturb("total"), "Variable 'total' was already perturbed!")
expect_message(tab$freqtab())

test_that("ck_define_table() with new record keys is ok", {
  expect_is(tab, "cellkey_obj")
  expect_identical(digest::digest(tab$freqtab("total")), "24801744e36e9843b6affb1b1c939831")
})

freqtab <- tab$freqtab("total")
test_that("weighted version of ck_perturb() is ok", {
  expect_identical(digest::digest(freqtab), "24801744e36e9843b6affb1b1c939831")
  expect_identical(digest::digest(tab$mod_cnts()), "4f6d21a2b4bdeb4d47f5b2d017ee81ca")
})

mm <- tab$measures("Total")
test_that("ck_cnt_measures() is ok", {
  expect_identical(digest::digest(mm), "0449a9b8753594c85089465a29d5f592")
})

# no weights
tab <- ck_setup(
  x = dat,
  rkey = 7,
  dims = dims,
  w = NULL,
  countvars = countvars,
  numvars = numvars,
  params_cnts = params_cnts,
  params_nums = params_nums
)

tab$perturb("total")
freqtab <- tab$freqtab("total")

test_that("checking unweighted version of perturb()", {
  expect_identical(digest::digest(freqtab), "ba9e5450bb0bd31531217dd21c580ea2")
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
  countvars = countvars,
  numvars = numvars,
  params_cnts = params_cnts,
  params_nums = params_nums
)
tab$perturb(c("total", "cnt_males", "cnt_highincome"))

test_that("check tabulation of cnt_males is ok", {
  expect_identical(digest::digest(tab$freqtab("cnt_males")), "d9e7311bab76591c1155472db77ddb57")
  expect_identical(digest::digest(tab$measures("cnt_males")), "94b7636cd5e3429dee64b9e8651b1ee2")
})

test_that("check tabulation of cnt_highincome is ok", {
  expect_identical(digest::digest(tab$freqtab("cnt_highincome")), "c73f3b9d525f181bd94c45d102237eca")
  expect_identical(digest::digest(tab$measures("cnt_highincome")), "ea9b68f43524f8ec7e37b69316ef9947")
})

test_that("check tabulation of multiple count variables is ok", {
  tt <- tab$freqtab(c("total", "cnt_males", "cnt_highincome"))
  expect_identical(digest::digest(tt), "ed0b0c09ecea37ab05ee13f99a7bfa0f")
})

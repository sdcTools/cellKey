context("Testing Frequency Tables")

set.seed(120, sample.kind = "Reject")
dat <- ck_create_testdata()

seed <- ck_create_seed_from_hash(dat)
test_that("check seed from hash generation", {
  expect_equal(seed, 1474937517)
})

dim_sex <- hier_create(root = "Total", nodes = c("male", "female"))
dim_age <- hier_create(root = "Total", paste0("age_group", 1:6))
dims <- list(sex = dim_sex, age = dim_age)
test_that("dims-hash is ok", {
  expect_identical(class(dims), "list")
  expect_identical(nrow(dim_sex), 3L)
  expect_identical(max(dim_sex$level), 2)
  expect_identical(nrow(dim_age), 7L)
  expect_identical(max(dim_age$level), 2)
})

## test generation of destatis rkeys
rk1 <- ck_generate_rkeys(dat = dat, nr_digits = 5)
rk2 <- ck_generate_rkeys(dat = dat, nr_digits = 5)
test_that("check rkey generation and seed is ok", {
  expect_identical(rk1, rk2)
  expect_identical(round(mean(rk1), digits = 3), 0.501)
})
dat$rec_key <- rk1
test_that("checking dimension and structure of generated testdata is ok", {
  expect_true(is.data.table(dat))
  expect_identical(round(mean(dat$sampling_weight), digits = 3), 59.719)
  expect_identical(round(mean(dat$household_weights), digits = 3), 21.834)

  expect_identical(nrow(dat), 4580L)
  expect_identical(ncol(dat), 16L)
  expect_identical(sum(dat$sex == "male"), 2296L)
})

## perturbation parameters for count variables
set.seed(120, sample.kind = "Reject")
suppressMessages(para <- ptable::create_cnt_ptable(
  D = 5,
  V = 3,
  js = 2,
  pstay = 0.5,
  optim = 1,
  mono = TRUE,
  create = FALSE))

params_cnts <- ck_params_cnts(ptab = ptable::create_ptable(params = para))
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
  dt <- params_cnts$params$ptable
  expect_is(dt, "data.table")
  expect_identical(dim(dt), c(66L, 7L))
  expect_identical(round(mean(dt$p), digits = 3), 0.136)
  expect_identical(round(mean(dt$lb), digits = 3), 0.562)
  expect_identical(round(mean(dt$ub), digits = 3), 0.698)
  expect_identical(round(mean(dt$v), digits = 3), 1.045)
})

countvars <- NULL
w <- "sampling_weight"
rkey <- "rec_key"

tab <- ck_setup(
  x = dat,
  rkey = "rec_key",
  dims = dims,
  w = w,
  countvars = countvars
)

# set perturbation parameters for all variables
tab$params_cnts_set(val = params_cnts, v = NULL)

expect_message(tab$perturb("total"))
expect_message(tab$perturb("total"), "Variable 'total' was already perturbed!")

res_freqtab <- tab$freqtab("total")
test_that("check ck_define_table() with already existing rec-keys", {
  expect_is(tab, "cellkey_obj")
  expect_identical(dim(res_freqtab), c(21L, 7L))
  expect_identical(res_freqtab$uwc[3], 1143)
  expect_identical(res_freqtab$puwc[3], 1147)
  expect_identical(round(mean(res_freqtab$pwc), digits = 3), 52096.24)
  expect_identical(round(mean(res_freqtab$puwc), digits = 3), 872.333)
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
  dt <- tab$freqtab("total")
  expect_identical(dt$uwc[3], 1143)
  expect_identical(dt$puwc[3], 1147)
  expect_identical(round(mean(dt$pwc), digits = 3), 52096.24)
  expect_identical(round(mean(dt$puwc), digits = 3), 872.333)
})

freqtab <- tab$freqtab("total")
test_that("weighted version of ck_perturb() is ok", {
  expect_identical(freqtab$uwc[3], 1143)
  expect_identical(freqtab$puwc[3], 1147)
  expect_identical(round(mean(freqtab$pwc), digits = 3), 52096.24)
  expect_identical(round(mean(freqtab$puwc), digits = 3), 872.333)

  dt <- tab$mod_cnts()
  expect_identical(dim(dt), c(21L, 6L))
  expect_identical(round(mean(dt$ckey), digits = 3), 0.456)
  expect_identical(round(mean(dt$pert), digits = 3), -0.048)
  expect_identical(range(dt$row_nr), c(40, 65))
})

mm <- tab$measures_cnts("Total")
test_that("ck_cnt_measures() [exclude_zeros = TRUE] is ok", {
  expect_identical(range(as.numeric(mm$overview$noise)), c(1, 7))
  expect_identical(range(as.numeric(mm$overview$cnt)), c(1, 9))
  expect_identical(round(mean(as.numeric(mm$overview$pct)), digits = 3), 0.143)

  expect_identical(range(as.numeric(mm$measures$d1)), c(0, 4))
  expect_identical(range(as.numeric(mm$measures$d2)), c(0, 0.429))
  expect_identical(range(as.numeric(mm$measures$d3)), c(0, 0.517))

  expect_identical(range(mm$cumdistr_d1$cnt), c(9L, 21L))
  expect_identical(round(range(mm$cumdistr_d1$pct), digits = 3), c(0.429, 1))

  expect_identical(range(mm$cumdistr_d2$cnt), c(19L, 21L))
  expect_identical(round(range(mm$cumdistr_d2$pct), digits = 3), c(0.905, 1))

  expect_identical(range(mm$cumdistr_d3$cnt), c(12L, 21L))
  expect_identical(round(range(mm$cumdistr_d3$pct), digits = 3), c(0.571, 1))
  expect_identical(mm$false_nonzero, 0L)
  expect_identical(mm$false_zero, 0L)
  expect_identical(mm$exclude_zeros, TRUE)
})

mm <- tab$measures_cnts("Total", exclude_zeros = FALSE)
test_that("ck_cnt_measures() [exclude_zeros = FALSE] is ok", {
  expect_identical(range(as.numeric(mm$overview$noise)), c(1, 7))
  expect_identical(range(as.numeric(mm$overview$cnt)), c(1, 9))
  expect_identical(round(mean(as.numeric(mm$overview$pct)), digits = 3), 0.143)

  expect_identical(range(as.numeric(mm$measures$d1)), c(0, 4))
  expect_identical(range(as.numeric(mm$measures$d2)), c(0, 0.429))
  expect_identical(range(as.numeric(mm$measures$d3)), c(0, 0.517))

  expect_identical(range(mm$cumdistr_d1$cnt), c(9L, 21L))
  expect_identical(round(range(mm$cumdistr_d1$pct), digits = 3), c(0.429, 1))

  expect_identical(range(mm$cumdistr_d2$cnt), c(19L, 21L))
  expect_identical(round(range(mm$cumdistr_d2$pct), digits = 3), c(0.905, 1))

  expect_identical(range(mm$cumdistr_d3$cnt), c(12L, 21L))
  expect_identical(round(range(mm$cumdistr_d3$pct), digits = 3), c(0.571, 1))
  expect_identical(mm$false_nonzero, 0L)
  expect_identical(mm$false_zero, 0L)
  expect_identical(mm$exclude_zeros, FALSE)
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
  expect_identical(dim(freqtab), c(21L, 7L))
  expect_identical(freqtab$uwc[3], 1143)
  expect_identical(freqtab$puwc[3], 1147)
  expect_identical(freqtab$pwc, freqtab$puwc)
  expect_identical(round(mean(freqtab$pwc), digits = 3), 872.333)

  dt <- tab$mod_cnts()
  expect_identical(dim(dt), c(21L, 6L))
  expect_identical(round(mean(dt$ckey), digits = 3), 0.456)
  expect_identical(round(mean(dt$pert), digits = 3), -0.048)
  expect_identical(range(dt$row_nr), c(40, 65))
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
  dt <- tab$freqtab("cnt_males")
  expect_identical(dt$uwc[3], 571)
  expect_identical(dt$puwc[3], 571)
  expect_identical(round(mean(dt$pwc), digits = 3), 25757.65)
  expect_identical(round(mean(dt$puwc), digits = 3), 437.048)

  expect_identical(round(range(dt$puwc), digits = 3), c(0, 2297))
  expect_identical(round(range(dt$pwc), digits = 3), c(0, 135387.941))

  mm <- tab$measures_cnts("cnt_males")
  expect_identical(round(range(mm$overview$pct), digits = 3), c(0.095, 0.810))
  expect_identical(round(range(mm$measures$d1), digits = 3), c(0, 4))
  expect_identical(round(range(mm$measures$d2), digits = 3), c(0, 0.048))
  expect_identical(round(range(mm$measures$d3), digits = 3), c(0, 0.221))
})

test_that("check tabulation of cnt_highincome is ok", {
  dt <- tab$freqtab("cnt_highincome")
  expect_identical(dt$uwc[3], 123)
  expect_identical(dt$puwc[3], 123)
  expect_identical(round(mean(dt$pwc), digits = 3), 5063.68)
  expect_identical(round(mean(dt$puwc), digits = 3), 84.286)

  expect_identical(round(range(dt$puwc), digits = 3), c(0, 444))
  expect_identical(round(range(dt$pwc), digits = 3), c(0, 26671.928))

  mm <- tab$measures_cnts("cnt_highincome")
  expect_identical(round(range(mm$overview$pct), digits = 3), c(0.048, 0.667))
  expect_identical(round(range(mm$measures$d1), digits = 3), c(0, 3))
  expect_identical(round(range(mm$measures$d2), digits = 3), c(0, 0.158))
  expect_identical(round(range(mm$measures$d3), digits = 3), c(0, 0.359))
})

test_that("check tabulation of multiple count variables is ok", {
  dt <- tab$freqtab(c("total", "cnt_males", "cnt_highincome"))
  expect_identical(dim(dt), c(63L, 7L))
  expect_identical(dt$uwc[3], 1143)
  expect_identical(dt$puwc[3], 1148)
  expect_identical(round(mean(dt$pwc), digits = 3), 27640.21)
  expect_identical(round(mean(dt$puwc), digits = 3), 464.587)
  expect_identical(round(range(dt$puwc), digits = 3), c(0, 4582))
  expect_identical(round(range(dt$pwc), digits = 3), c(0, 273633.438))

  mm <- tab$measures_cnts("cnt_highincome")
  expect_identical(round(range(mm$overview$pct), digits = 3), c(0.048, 0.667))
  expect_identical(round(range(mm$measures$d1), digits = 3), c(0, 3))
  expect_identical(round(range(mm$measures$d2), digits = 3), c(0, 0.158))
  expect_identical(round(range(mm$measures$d3), digits = 3), c(0, 0.359))
})

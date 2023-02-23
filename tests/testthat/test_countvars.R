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
  expect_identical(digest::sha1(dim_sex), "fea2001f35be84e90b30f6773af75f03c11fbf7a")
  expect_identical(digest::sha1(dim_age), "a7648dc3f484720911f0de0e6ac563b69fd20c42")
  expect_identical(digest::sha1(dims), "62748837ca3246a33081dd35f50d06334caa3119")
})

## test generation of destatis rkeys
rk1 <- ck_generate_rkeys(dat = dat, nr_digits = 5)
rk2 <- ck_generate_rkeys(dat = dat, nr_digits = 5)
test_that("check rkey generation and seed is ok", {
  expect_identical(rk1, rk2)
  expect_identical(digest::sha1(rk1), "4de74ed6170e2142ef552ee6722921db8d091d0c")
})
dat$rec_key <- rk1
test_that("checking dimension and structure of generated testdata is ok", {
  expect_identical(digest::sha1(dat), "fb66a8be3e9044c8fecdb13c6fab5fe9ec456c25")
  expect_true(is.data.table(dat))
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
  expect_is(params_cnts$params$ptable, "data.table")
  expect_identical(digest::sha1(params_cnts), "cedf56d7064f15e55da506b1922c4cac6035765f")
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
  expect_identical(res_freqtab$uwc[3], 1143)
  expect_identical(res_freqtab$puwc[3], 1147)
  expect_identical(digest::sha1(res_freqtab), "05e71f630a385f7e428ce1fec21b5f6026bb921a")
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
  expect_identical(digest::sha1(tab$freqtab("total")), "05e71f630a385f7e428ce1fec21b5f6026bb921a")
})

freqtab <- tab$freqtab("total")
test_that("weighted version of ck_perturb() is ok", {
  expect_identical(digest::sha1(freqtab), "05e71f630a385f7e428ce1fec21b5f6026bb921a")
  expect_identical(digest::sha1(tab$mod_cnts()), "ee05433bb69fbf66094cf14e5c50320b3060eab3")
})

mm <- tab$measures_cnts("Total")
test_that("ck_cnt_measures() [exclude_zeros = TRUE] is ok", {
  expect_identical(digest::sha1(mm), "89f4aba98334930446c2fb97b0812b7ece98ef6d")
})

mm <- tab$measures_cnts("Total", exclude_zeros = FALSE)
test_that("ck_cnt_measures() [exclude_zeros = FALSE] is ok", {
  expect_identical(digest::sha1(mm), "9fb4ffe32ebc420d8ccecb6f3dab6e1431396fa0")
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
  expect_identical(digest::sha1(freqtab), "31831ab589f3bc7bdd20c4b6fd1ea1916e3edfef")
  expect_identical(digest::sha1(tab$mod_cnts()), "ee05433bb69fbf66094cf14e5c50320b3060eab3")
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
  expect_identical(digest::sha1(tab$freqtab("cnt_males")), "a7260c1f65a089d8849084a7e84c946be997e986")
  expect_identical(digest::sha1(tab$measures_cnts("cnt_males")), "b692361c523b0d37aa38bf252a360414cd941e2b")
})

test_that("check tabulation of cnt_highincome is ok", {
  expect_identical(digest::sha1(tab$freqtab("cnt_highincome")), "72c684edc0cebef0343bce02d0075050d5286a5c")
  expect_identical(digest::sha1(tab$measures_cnts("cnt_highincome")), "639f20e406ecfdd44c84e0ae04674321c2a602b1")
})

test_that("check tabulation of multiple count variables is ok", {
  tt <- tab$freqtab(c("total", "cnt_males", "cnt_highincome"))
  expect_identical(digest::sha1(tt), "1984413d2bb4ebff2d268b243b300cd8494b55fb")
})

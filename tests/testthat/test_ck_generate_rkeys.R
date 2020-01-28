if (R.version$major >= 3 & R.version$minor >= 6) {
  suppressWarnings(set.seed(120, sample.kind = "Rounding"))
} else {
  set.seed(120)
}

x <- ck_create_testdata()

context("testing generation of record keys")
test_that("invalid inputs are caught", {
  expect_error(ck_generate_rkeys(dat = 5))
  expect_error(ck_generate_rkeys(dat = x, nr_digits = "5"))
  expect_error(ck_generate_rkeys(dat = x, nr_digits = 1:2))
  expect_error(ck_generate_rkeys(dat = x, nr_digits = 4))
  expect_error(ck_generate_rkeys(dat = x, nr_digits = 25))
})

ck1 <- ck_generate_rkeys(dat = x, nr_digits = 8, seed = NULL)
ck2 <- ck_generate_rkeys(dat = x, nr_digits = 10, seed = 5)

test_that("recordkeys are correctly computed", {
  expect_equal(digest::sha1(ck1), "9462c96cea80aa0ad1eeb0a18a2c319cb2c30d29")
  expect_equal(digest::sha1(ck2), "5677f580e6ed9fa59ab670b0e7badc67a37bf14c")
})

set.seed(120, sample.kind = "Reject")
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
  expect_equal(digest::sha1(ck1), "a291a306de5db64e3fbcb94478982733a17ede20")
  expect_equal(digest::sha1(ck2), "5677f580e6ed9fa59ab670b0e7badc67a37bf14c")
})

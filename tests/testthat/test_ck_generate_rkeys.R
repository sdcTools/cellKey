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
  expect_identical(length(ck1), 4580L)
  expect_identical(ck1[1], 0.97131759)
  expect_identical(ck1[4580], 0.88019983)
  expect_identical(round(mean(ck1), digits = 3), 0.501)


  expect_identical(length(ck2), 4580L)
  expect_identical(ck2[1], 0.2002144526)
  expect_identical(ck2[4580], 0.1878750164)
  expect_identical(round(mean(ck2), digits = 3), 0.508)
})

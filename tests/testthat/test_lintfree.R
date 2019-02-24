context("Test style (lintr)")
test_that("Style should be lint-free", {
  skip_if_not(
    requireNamespace("lintr", quietly = TRUE),
    message = "Package lintr must be installed!"
  )
  lintr::expect_lint_free()
})

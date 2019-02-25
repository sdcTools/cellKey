context("Testing countVars")

set.seed(120)
dat <- ck_create_testdata()
dat[, cnt_mixed := ifelse(sex == "male", 0, 1)]
dat[, cnt_onlyzero := 0]
dat[, cnt_onlyones := 1]
dat$rkeys <-
  ck_generate_rkeys(
    dat = dat,
    max_val = 2 * nrow(dat),
    type = "abs",
    verbose = TRUE
  )

ptab <- ck_create_ptab(
  D = 5,
  V = 3,
  ptab_size = 70,
  type = "abs"
)
pert_params <- ck_create_pert_params(
  big_n = 17312941,
  small_n = 12,
  ptab = ptab
)
inp <- ck_create_input(
  dat = dat,
  def_rkey = "rkeys",
  pert_params = pert_params
)

dim_sex <- hier_create(root = "Total", nodes = c("male", "female"))
dim_age <- hier_create(root = "Total", nodes = paste0("age_group", 1:6))
dl <- list(sex = dim_sex, age = dim_age)
weightVar <- "sampling_weight"

## perturbing the table
res <- perturb_table(
  inp = inp,
  dim_list = dl,
  weightVar = weightVar,
  countVars = c("cnt_mixed", "cnt_onlyones", "cnt_onlyzero"),
  numVars = NULL
)
r0 <- ck_freq_table(res, vname = "Total")
r1 <- ck_freq_table(res, vname = "cnt_onlyones")
r2 <- ck_freq_table(res, vname = "cnt_onlyzero")

test_that("only ones equals total", {
  expect_equal(dim(r0), dim(r1))
  expect_true(all(
    sapply(1:ncol(r0), function(i) {
      identical(r0[[i]], r1[[i]])
    }
  ) == TRUE))
})

test_that("onlyzero equals total", {
  expect_equal(dim(r0), dim(r2))
  expect_true(sum(r2[, UWC_cnt_onlyzero]) == 0)
  expect_true(sum(r2[, WC_cnt_onlyzero]) == 0)
  expect_true(sum(r2[, pUWC_cnt_onlyzero]) == 0)
  expect_true(sum(r2[, pWC_cnt_onlyzero]) == 0)
  expect_true(sum(r2[, WCavg_cnt_onlyzero]) == 0)
  expect_true(sum(r2[, row_indices]) == nrow(r2))
  expect_true(sum(r2[, col_indices]) == 0)
  expect_true(sum(r2[, cellKey]) == 0)
})

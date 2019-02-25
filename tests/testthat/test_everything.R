context("Testing Frequency Tables")

set.seed(120)
big_n <- 17312941
ptab_size <- 70
small_n <- 12
dat <- ck_create_testdata()

seed <- ck_create_seed_from_hash(dat)
test_that("check seed from hash generation", {
  expect_equal(seed, 483520795)
})

dim_sex <- hier_create(root = "total", nodes = c("male", "female"))
dim_age <- hier_create(root = "total", nodes = paste0("age_group", 1:6))
dl <- list(sex = dim_sex, age = dim_age)

## test generation of destatis rkeys
rr <- ck_generate_rkeys(
  dat = dat,
  max_digits = 5,
  type = "destatis"
)
test_that("check rkey generation for destatis (and also seed)", {
  expect_equal(rr[1], 0.52107)
  expect_equal(rr[3], 0.45809)
  expect_equal(rr[5], 0.23388)
  expect_equal(min(rr), 9e-05)
  expect_equal(max(rr), 0.99986)
})

maxV <- 10 * nrow(dat)
dat$rec_key <- ck_generate_rkeys(
  dat = dat,
  max_val = maxV,
  type = "abs"
)

test_that("checking dimension and structure of generated testdata", {
  expect_equal(nrow(dat), 4580)
  expect_equal(ncol(dat), 16)
  expect_true(is.data.table(dat))
})

# perturbation input for continuous variables
mtab <- c(0.6, 0.4, 0.2)
small_c <- 12

## ptable
ptab <- ck_create_ptab(
  D = 5,
  V = 3,
  type = "abs",
  ptab_size = ptab_size
)

pert_params <- ck_create_pert_params(
  big_n = big_n,
  small_n = small_n,
  ptab = ptab,
  stab = ck_generate_stab(small_c = small_c),
  mtab = mtab
)

inp <- ck_create_input(
  dat = dat,
  def_rkey = "rec_key",
  pert_params = pert_params
)
test_that("checking output of createInput() with already existing rec-keys", {
  expect_s4_class(inp, "pert_inputdat")
  expect_s4_class(inp@pert_params, "pert_params")
  expect_equal(slot(inp@pert_params, "big_n"), 17312941)
})

dat$rec_key <- NULL
inp2 <- ck_create_input(
  dat = dat,
  def_rkey = 23984936,
  pert_params = pert_params
)
test_that("checking output of createInput() with non-existing record keys", {
  expect_s4_class(inp2, "pert_inputdat")
  expect_s4_class(inp2@pert_params, "pert_params")
  expect_equal(slot(inp2@pert_params, "big_n"), 17312941)
})

# Frequency Tables, weighted
tab_freq <- perturb_table(
  inp = inp,
  dim_list = dl,
  weightvar = "sampling_weight",
  numvars = NULL
)
res_weighted <- slot(tab_freq, "tab")

test_that("checking weighted version of perturbedFreqTable", {
  expect_equal(head(res_weighted$uwc_total, 1), 4580)
})

tab_freq_noweights <- perturb_table(
  inp = inp,
  dim_list = dl,
  weightvar = NULL,
  numvars = NULL
)
res_unweighted <- slot(tab_freq_noweights, "tab")
test_that("checking unweighted version of perturbedFreqTable", {
  expect_equal(head(res_unweighted$uwc_total, 1), 4580)
  expect_equal(head(res_unweighted$wc_total, 1), 4580)
})

context("Testing Magnitude Tables")
tab_cont <- perturb_table(
  inp = inp,
  dim_list = dl,
  weightvar = NULL,
  numvars = "savings"
)

test_that("checking unweighted version of perturbedContTable", {
  expect_s4_class(tab_cont, "pert_table")
  expect_equal(nrow(tab_cont@tab), 21)
  expect_equal(ncol(tab_cont@tab), 11)
})
context("Testing consistent record-key generation")
# for the same input-dataset and the same setting of def_rkey, the
# same record keys must be generated
inp1 <- ck_create_input(
  dat = dat,
  def_rkey = 1e5,
  pert_params = pert_params
)
inp2 <- ck_create_input(
  dat = dat,
  def_rkey = 1e5,
  pert_params = pert_params
)
test_that("check that identical record keys have been generated", {
  expect_s4_class(inp1, "pert_inputdat")
  expect_s4_class(inp2, "pert_inputdat")
  expect_identical(slot(inp1, "rkeys"), slot(inp2, "rkeys"))
})

context("Testing countvars argument")
dat[, cnt_males := ifelse(sex == "male", 1, 0)]
dat[, cnt_highincome := ifelse(income >= 9000, 1, 0)]
inp <- ck_create_input(
  dat = dat,
  def_rkey = 1e5,
  pert_params = pert_params
)

tab <- perturb_table(
  inp = inp,
  dim_list = dl,
  weightvar = "sampling_weight",
  countvars = c("cnt_males", "cnt_highincome"),
  numvars = NULL,
  by = NULL
)
tt <- ck_freq_table(tab, "cnt_males")
test_that("check tabulation of cnt_males", {
  expect_identical(tt[sex == "female", sum(uwc_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(wc_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(puwc_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(pwc_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(wc_avg_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(cellkey)], 0)

  expect_identical(max(tt[, uwc_cnt_males]), 2296)
  expect_identical(max(tt[, wc_cnt_males]), 137489)
  expect_identical(max(tt[, puwc_cnt_males]), 2294)
  expect_identical(max(tt[, pwc_cnt_males]), 137369)
  expect_identical(tt[1, cellkey], 9580150)
})

tt <- ck_freq_table(tab, "cnt_highincome")
test_that("check tabulation of cnt_males", {
  expect_identical(tt[age == "age_group6", sum(uwc_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(wc_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(puwc_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(pwc_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(wc_avg_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(cellkey)], 0)

  expect_identical(max(tt[, uwc_cnt_highincome]), 445)
  expect_identical(max(tt[, wc_cnt_highincome]), 26797)
  expect_identical(max(tt[, puwc_cnt_highincome]), 446)
  expect_identical(max(tt[, pwc_cnt_highincome]), 26857)
  expect_identical(tt[1, cellkey], 4477382)
})

context("Testing numvars with by")
tab <- perturb_table(
  inp = inp,
  dim_list = dl,
  weightvar = "sampling_weight",
  countvars = c("cnt_males"),
  numvars = c("income", "savings"),
  by = "cnt_males"
)
tt <- ck_cont_table(tab, vname = "savings", mean_before_sum = TRUE)
test_that("check tabulation of savings given cnt_males", {
  expect_identical(tt[sex == "female", sum(uw_savings)], 0)
  expect_identical(tt[sex == "female", sum(puw_savings)], 0)
  expect_identical(tt[sex == "female", sum(ws_savings)], 0)
  expect_identical(tt[sex == "female", sum(pws_savings)], 0)
  expect_identical(tt[sex == "female", sum(pwm_savings)], 0)

  expect_identical(tt[1, uw_savings], 1159816)
  expect_equal(tt[1, puw_savings], 1162403.53)
  expect_equal(tt[1, ws_savings], 69452065)
  expect_equal(tt[1, pws_savings], 69607012)
  expect_identical(round(tt[1, pwm_savings], digits = 3), 506.716)
})

pp <- attr(tt, "modifications")
test_that("check perturbation attributes of savings given cnt_males", {
  expect_identical(pp[id == "3452", magnitude], 0.4)
  expect_identical(pp[id == "3452", dir], 1)
  expect_identical(pp[id == "3452", noise], 1.37)
  expect_identical(pp[id == "3452", vals_orig], 998)
  expect_equal(pp[id == "3452", vals_pert], 546.904)
  expect_identical(pp[id == "3452", vals_mod], 1544.904)
})

context("Testing Frequency Tables (destatis)")

set.seed(120)
bigN <- 17312941
pTableSize <- 70
smallN <- 12
dat <- ck_create_testdata()

seed <- ck_create_seed_from_hash(dat)
test_that("check seed from hash generation", {
  expect_equal(seed, 483520795)
})

dim.sex <- hier_create(root = "Total", nodes = c("male", "female"))
dim.age <- hier_create(root = "Total", paste0("age_group", 1:6))
dimList <- list(sex = dim.sex, age = dim.age)

## test generation of destatis rkeys
dat$rec_key <- ck_generate_rkeys(dat = dat, max_digits = 5, type = "destatis")
test_that("check rkey generation for destatis (and also seed)", {
  expect_equal(dat$rec_key[1], 0.52107)
  expect_equal(dat$rec_key[3], 0.45809)
  expect_equal(dat$rec_key[5], 0.23388)
  expect_equal(min(dat$rec_key), 9e-05)
  expect_equal(max(dat$rec_key), 0.99986)
})

test_that("checking dimension and structure of generated testdata", {
  expect_equal(nrow(dat), 4580)
  expect_equal(ncol(dat), 16)
  expect_true(is.data.table(dat))
})

# perturbation input for continuous variables
mTable <- c(0.6, 0.4, 0.2)
smallC <- 12

## ptable
pTable <- ck_create_pTable(D = 5, V = 3, type = "destatis", pTableSize = pTableSize)
pert_params <- ck_create_pert_params(
  bigN = bigN,
  smallN = smallN,
  pTable = pTable,
  sTable = ck_generate_sTable(smallC = smallC),
  mTable = mTable)

inp <- ck_create_input(dat = dat, def_rkey = "rec_key", pert_params = pert_params)
test_that("checking output of createInput() with already existing rec-keys", {
  expect_s4_class(inp, "pert_inputdat")
  expect_s4_class(inp@pert_params, "pert_params")
  expect_equal(slot(inp@pert_params, "bigN"), 1)
  expect_equal(slot(inp@pert_params, "type"), "destatis")

})

dat$rec_key <- NULL
inp2 <- ck_create_input(dat = dat, def_rkey = 7, pert_params = pert_params)
test_that("checking output of createInput() with non-existing record keys", {
  expect_s4_class(inp2, "pert_inputdat")
  expect_s4_class(inp2@pert_params, "pert_params")
  expect_equal(slot(inp2@pert_params, "bigN"), 1)
  expect_equal(diff(range(inp2@rkeys)), 0.9997682)
})

# ## Frequency Tables
# # weighted
tab_freq <- perturbTable(inp, dimList = dimList, weightVar = "sampling_weight", numVars = NULL)
res_weighted <- slot(tab_freq, "tab")

test_that("checking weighted version of perturbedFreqTable", {
  expect_equal(head(res_weighted$UWC_Total, 1), 4580)
})

tab_freq_noweights <- perturbTable(inp = inp, dimList = dimList, weightVar = NULL, numVars = NULL)
res_unweighted <- slot(tab_freq_noweights, "tab")
test_that("checking unweighted version of perturbedFreqTable", {
  expect_equal(head(res_unweighted$UWC_Total, 1), 4580)
  expect_equal(head(res_unweighted$WC_Total, 1), 4580)
})

context("Testing Magnitude Tables (destatis)")
tab_cont <- perturbTable(inp = inp, dimList = dimList, weightVar = NULL, numVars = c("savings"))

test_that("checking unweighted version of perturbedContTable", {
  expect_s4_class(tab_cont, "pert_table")
  expect_equal(nrow(tab_cont@tab), 21)
  expect_equal(ncol(tab_cont@tab), 11)
})
context("Testing consistent record-key generation (destatis)")
# for the same input-dataset and the same setting of def_rkey, the
# same record keys must be generated
inp1 <- ck_create_input(dat = dat, def_rkey = 7, pert_params = pert_params)
inp2 <- ck_create_input(dat = dat, def_rkey = 7, pert_params = pert_params)
test_that("check that identical record keys have been generated", {
  expect_s4_class(inp1, "pert_inputdat")
  expect_s4_class(inp2, "pert_inputdat")
  expect_identical(slot(inp1, "rkeys"), slot(inp2, "rkeys"))
})

context("Testing countVars argument (destatis)")
dat[, cnt_males := ifelse(sex == "male", 1, 0)]
dat[, cnt_highincome := ifelse(income >= 9000, 1, 0)]
inp <- ck_create_input(dat = dat, def_rkey = 7, pert_params = pert_params)

tab <- perturbTable(
  inp,
  dimList = dimList,
  weightVar = "sampling_weight",
  countVars = c("cnt_males", "cnt_highincome"),
  numVars = NULL,
  by = NULL)

tt <- ck_freq_table(tab, "cnt_males")
test_that("check tabulation of cnt_males (destatis)", {
  expect_identical(tt[sex == "female", sum(UWC_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(WC_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(pUWC_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(pWC_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(WCavg_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(cellKey)], 0)

  expect_identical(max(tt[, UWC_cnt_males]), 2296)
  expect_identical(max(tt[, WC_cnt_males]), 137489)
  expect_identical(max(tt[, pUWC_cnt_males]), 2296)
  expect_identical(max(tt[, pWC_cnt_males]), 137489)
  expect_equal(tt[1, cellKey], 0.5665173)
})

tt <- ck_freq_table(tab, "cnt_highincome")
test_that("check tabulation of cnt_males (destatis)", {
  expect_identical(tt[age == "age_group6", sum(UWC_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(WC_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(pUWC_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(pWC_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(WCavg_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(cellKey)], 0)

  expect_identical(max(tt[, UWC_cnt_highincome]), 445)
  expect_identical(max(tt[, WC_cnt_highincome]), 26797)
  expect_identical(max(tt[, pUWC_cnt_highincome]), 447)
  expect_identical(max(tt[, pWC_cnt_highincome]), 26917)
  expect_equal(tt[1, cellKey], 0.9010093)
})

context("Testing numVars with by (destatis)")
tab <- perturbTable(
  inp,
  dimList = dimList,
  weightVar = "sampling_weight",
  countVars = c("cnt_males"),
  numVars = c("income", "savings"),
  by = "cnt_males")

tt <- ck_cont_table(tab, vname = "savings", meanBeforeSum = TRUE)
test_that("check tabulation of savings given cnt_males (destatis)", {
  expect_identical(tt[sex == "female", sum(UW_savings)], 0)
  expect_identical(tt[sex == "female", sum(pUW_savings)], 0)
  expect_identical(tt[sex == "female", sum(WS_savings)], 0)
  expect_identical(tt[sex == "female", sum(pWS_savings)], 0)
  expect_identical(tt[sex == "female", sum(pWM_savings)], 0)

  expect_identical(tt[1, UW_savings], 1159816)
  expect_equal(tt[1, pUW_savings], 1163298.346)
  expect_equal(tt[1, WS_savings], 69452065.3414634)
  expect_equal(tt[1, pWS_savings], 69660595.0754329)
  expect_equal(tt[1, pWM_savings], 506.663042682927)
})

pp <- attr(tt, "modifications")
test_that("check perturbation attributes of savings given cnt_males (destatis)", {
  expect_identical(pp[id == "3452", magnitude], 0.4)
  expect_identical(pp[id == "3452", dir], 1)
  expect_identical(pp[id == "3452", noise], 1.02)
  expect_identical(pp[id == "3452", vals.orig], 998)
  expect_equal(pp[id == "3452", vals.pert], 407.184)
  expect_identical(pp[id == "3452", vals.mod], 1405.184)
})

context("Testing count-measures (destatis)")
mm <- ck_cnt_measures(tab, vname = "Total")
expect_is(mm, "list")
expect_is(mm$measures, "data.table")
expect_is(mm$cumdistrA, "data.table")
expect_is(mm$cumdistrB, "data.table")
expect_equal(mm$false_zero, 0)
expect_equal(mm$false_positives, 0)
expect_equal(max(mm$measures$vals_r), 0.196)
expect_equal(mm$cumdistrA$prop_abs[1], 1 / 3)
expect_equal(mm$cumdistrA$prop_abs[2], 0.80952380952381)
expect_equal(mm$cumdistrA$prop_abs[3], 0.952380952380952)
mm <- ck_cnt_measures(tab, vname = "cnt_males")

expect_is(mm, "list")
expect_is(mm$measures, "data.table")
expect_is(mm$cumdistrA, "data.table")
expect_is(mm$cumdistrB, "data.table")
expect_equal(mm$false_zero, 0)
expect_equal(mm$false_positives, 0)
expect_equal(mm$cumdistrA$prop_abs[1], 0.619047619047619)
expect_equal(mm$cumdistrA$prop_abs[2], 0.809523809523810)
expect_equal(mm$cumdistrA$prop_abs[3], 1)

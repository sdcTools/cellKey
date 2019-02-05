context("Testing Frequency Tables (abs)")

set.seed(120)
options(digits = 20)
bigN <- 17312941
pTableSize <- 70
smallN <- 12
dat <- ck_create_testdata()

seed <- ck_create_seed_from_hash(dat)
test_that("check seed from hash generation", {
  expect_equal(seed, 483520795)
})

dim.sex <- ck_create_node(total_lab = "Total")
dim.sex <- ck_add_nodes(dim.sex, reference_node = "Total", node_labs = c("male", "female"))

dim.age <- ck_create_node(total_lab = "Total")
dim.age <- ck_add_nodes(dim.age, reference_node = "Total", node_labs = paste0("age_group", 1:6))
dimList <- list(sex = dim.sex, age = dim.age)

maxV <- 10 * nrow(dat)
dat$rec_key <- ck_generate_rkeys(dat = dat, max_val = maxV, type = "abs")

test_that("checking dimension and structure of generated testdata (abs)", {
  expect_equal(nrow(dat), 4580)
  expect_equal(ncol(dat), 16)
  expect_true(is.data.table(dat))
})

# perturbation input for continuous variables
mTable <- c(0.6, 0.4, 0.2)
smallC <- 12

## ptable
pTable <- ck_create_pTable(D = 5, V = 3, type = "abs", pTableSize = pTableSize)
pert_params <- ck_create_pert_params(
  bigN = bigN,
  smallN = smallN,
  pTable = pTable,
  sTable = ck_generate_sTable(smallC = smallC),
  mTable = mTable)

inp <- ck_create_input(dat = dat, def_rkey = "rec_key", pert_params = pert_params)
test_that("checking output of createInput() with already existing rec-keys (abs)", {
  expect_s4_class(inp, "pert_inputdat")
  expect_s4_class(inp@pert_params, "pert_params")
  expect_equal(slot(inp@pert_params, "bigN"), 17312941)
})

dat$rec_key <- NULL
inp2 <- ck_create_input(dat = dat, def_rkey = 23984936, pert_params = pert_params)
test_that("checking output of createInput() with non-existing record keys (abs)", {
  expect_s4_class(inp2, "pert_inputdat")
  expect_s4_class(inp2@pert_params, "pert_params")
  expect_equal(slot(inp2@pert_params, "bigN"), 17312941)
})

# ## Frequency Tables
# # weighted
tab_freq <- perturbTable(inp, dimList = dimList, weightVar = "sampling_weight", numVars = NULL)


test_that("print/summary of perturbed table (abs)", {
  expect_equal(print(tab_freq), NULL)
  expect_error(print(tab_freq, vname = "x"))

  ss <- summary(tab_freq)
  expect_is(ss, "list")
  expect_equal(names(ss), c("cnt_info", "cnt_measures", "num_info"))
  expect_identical(ss$cnt_info$Mean, -0.143)
  expect_identical(ss$cnt_measures$Total$measures$vals_abs[6], 1.571)
})

test_that("export of perturbed table (abs)", {
  expect_error(ck_export_table(x = tab_freq, vname="x", type = "both"))
  expect_error(ck_export_table(x = tab_freq, vname="Total", type = "x"))
  expect_equal(
    dim(ck_export_table(x = tab_freq, vname="Total", type = "both")),
    c(21, 7)
  )
  expect_equal(
    dim(ck_export_table(x = tab_freq, vname="Total", type = "weighted")),
    c(21, 5)
  )
  expect_equal(
    dim(ck_export_table(x = tab_freq, vname="Total", type = "unweighted")),
    c(21, 5)
  )
})

res_weighted <- slot(tab_freq, "tab")

test_that("checking weighted version of perturbedFreqTable (abs)", {
  expect_equal(head(res_weighted$UWC_Total, 1), 4580)
})

tab_freq_noweights <- perturbTable(inp = inp, dimList = dimList, weightVar = NULL, numVars = NULL)
res_unweighted <- slot(tab_freq_noweights, "tab")
test_that("checking unweighted version of perturbedFreqTable (abs)", {
  expect_equal(head(res_unweighted$UWC_Total, 1), 4580)
  expect_equal(head(res_unweighted$WC_Total, 1), 4580)
})

context("Testing Magnitude Tables (abs)")
tab_cont <- perturbTable(inp = inp, dimList = dimList, weightVar = NULL, numVars = c("savings"))

test_that("checking unweighted version of perturbedContTable (abs)", {
  expect_s4_class(tab_cont, "pert_table")
  expect_equal(nrow(tab_cont@tab), 21)
  expect_equal(ncol(tab_cont@tab), 11)
})
context("Testing consistent record-key generation (abs)")
# for the same input-dataset and the same setting of def_rkey, the
# same record keys must be generated
inp1 <- ck_create_input(dat = dat, def_rkey = 1e5, pert_params = pert_params)
inp2 <- ck_create_input(dat = dat, def_rkey = 1e5, pert_params = pert_params)
test_that("check that identical record keys have been generated (abs)", {
  expect_s4_class(inp1, "pert_inputdat")
  expect_s4_class(inp2, "pert_inputdat")
  expect_identical(slot(inp1, "rkeys"), slot(inp2, "rkeys"))
})


context("Testing countVars argument (abs)")
dat[, cnt_males := ifelse(sex == "male", 1, 0)]
dat[, cnt_highincome := ifelse(income >= 9000, 1, 0)]
inp <- ck_create_input(dat = dat, def_rkey = 1e5, pert_params = pert_params)

tab <- perturbTable(
  inp,
  dimList = dimList,
  weightVar = "sampling_weight",
  countVars = c("cnt_males", "cnt_highincome"),
  numVars = NULL,
  by = NULL)
tt <- ck_freq_table(tab, "cnt_males")
test_that("check tabulation of cnt_males (abs)", {
  expect_identical(tt[sex == "female", sum(UWC_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(WC_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(pUWC_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(pWC_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(WCavg_cnt_males)], 0)
  expect_identical(tt[sex == "female", sum(cellKey)], 0)

  expect_identical(max(tt[, UWC_cnt_males]), 2296)
  expect_identical(max(tt[, WC_cnt_males]), 137489)
  expect_identical(max(tt[, pUWC_cnt_males]), 2294)
  expect_identical(max(tt[, pWC_cnt_males]), 137369)
  expect_identical(tt[1, cellKey], 9580150)
})

tt <- ck_freq_table(tab, "cnt_highincome")
test_that("check tabulation of cnt_males (abs)", {
  expect_identical(tt[age == "age_group6", sum(UWC_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(WC_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(pUWC_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(pWC_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(WCavg_cnt_highincome)], 0)
  expect_identical(tt[age == "age_group6", sum(cellKey)], 0)

  expect_identical(max(tt[, UWC_cnt_highincome]), 445)
  expect_identical(max(tt[, WC_cnt_highincome]), 26797)
  expect_identical(max(tt[, pUWC_cnt_highincome]), 446)
  expect_identical(max(tt[, pWC_cnt_highincome]), 26857)
  expect_identical(tt[1, cellKey], 4477382)
})

context("Testing numVars with by (abs)")
tab <- perturbTable(
  inp,
  dimList = dimList,
  weightVar = "sampling_weight",
  countVars = c("cnt_males"),
  numVars = c("income", "savings"),
  by = "cnt_males")

tt <- ck_cont_table(tab, vname = "savings", meanBeforeSum = TRUE)
test_that("check tabulation of savings given cnt_males (abs)", {
  expect_identical(tt[sex == "female", sum(UW_savings)], 0)
  expect_identical(tt[sex == "female", sum(pUW_savings)], 0)
  expect_identical(tt[sex == "female", sum(WS_savings)], 0)
  expect_identical(tt[sex == "female", sum(pWS_savings)], 0)
  expect_identical(tt[sex == "female", sum(pWM_savings)], 0)

  expect_identical(tt[1, UW_savings], 1159816)
  expect_equal(tt[1, pUW_savings], 1162403.53)
  expect_equal(tt[1, WS_savings], 69452065)
  expect_equal(tt[1, pWS_savings], 69607012)
  expect_identical(round(tt[1, pWM_savings], digits = 3), 506.716)
})

pp <- attr(tt, "modifications")
test_that("check perturbation attributes of savings given cnt_males (abs)", {
  expect_identical(pp[id == "3452", magnitude], 0.4)
  expect_identical(pp[id == "3452", dir], 1)
  expect_identical(pp[id == "3452", noise], 1.37)
  expect_identical(pp[id == "3452", vals.orig], 998)
  expect_equal(pp[id == "3452", vals.pert], 546.904)
  expect_identical(pp[id == "3452", vals.mod], 1544.904)
})

context("Testing count-measures (abs)")
mm <- ck_cnt_measures(tab, vname = "Total")
expect_is(mm, "list")
expect_is(mm$measures, "data.table")
expect_is(mm$cumdistrA, "data.table")
expect_is(mm$cumdistrB, "data.table")
expect_equal(mm$false_zero, 0)
expect_equal(mm$false_positives, 0)
expect_equal(max(mm$measures$vals_r), 0.449)
expect_equal(mm$cumdistrA$prop_abs[1], 0.14285714285714284921)
expect_equal(mm$cumdistrA$prop_abs[2], 0.42857142857142854764)
expect_equal(mm$cumdistrA$prop_abs[3], 0.90476190476190476719)
expect_equal(mm$cumdistrA$prop_abs[4], 0.95238095238095232808)
expect_equal(mm$cumdistrA$prop_abs[5], 1)
mm <- ck_cnt_measures(tab, vname = "cnt_males")

expect_is(mm, "list")
expect_is(mm$measures, "data.table")
expect_is(mm$cumdistrA, "data.table")
expect_is(mm$cumdistrB, "data.table")
expect_equal(mm$false_zero, 0)
expect_equal(mm$false_positives, 0)
expect_equal(mm$cumdistrA$prop_abs[1], 0.52380952380952383596)
expect_equal(mm$cumdistrA$prop_abs[2], 0.61904761904761906877)
expect_equal(mm$cumdistrA$prop_abs[3], 1)

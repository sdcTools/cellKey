set.seed(100)
params_cnts <- ck_params_cnts(
  D = 5,
  V = 3,
  js = 2,
  pstay = 0.5,
  optim = 1,
  mono = TRUE
)
usethis::use_data(params_cnts, overwrite = TRUE)

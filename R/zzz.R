.onAttach <- function(lib, pkg) {
  old <- options()
  on.exit(options(old))
  options(useFancyQuotes = FALSE)
  packageStartupMessage(paste(
    "Package cellKey", utils::packageVersion("cellKey"),"has been loaded.",
    "To enable parallel processing, set `Sys.setenv('CK_RUN_PARALLEL' = TRUE)`"
  ))
}

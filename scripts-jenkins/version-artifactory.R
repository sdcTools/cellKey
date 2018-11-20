args = commandArgs(trailingOnly=TRUE)

library(desc)
desc <- description$new()
packages <- available.packages(repos = c(args[1]))
package <- subset(as.data.frame(packages, stringsAsFactors = FALSE), Package == desc$get("Package"))
if (nrow(package) == 1) {
  cat(trimws(package$Version))
} else {
  cat("NA")
}
# needed if additional packages have to be installed
# .libPaths(c("./R"))

# just to test if package is loaded from artifactory
# install.packages("brew")

library("devtools")

devtools::build(path = ".",vignettes = FALSE)
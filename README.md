
## cellKey

[![Travis build
status](https://travis-ci.org/sdcTools/cellKey.svg?branch=master)](https://travis-ci.org/sdcTools/cellKey)
[![Coverage
status](https://codecov.io/gh/sdcTools/cellKey/branch/master/graph/badge.svg)](https://codecov.io/github/sdcTools/cellKey?branch=master)
[![GitHub last
commit](https://img.shields.io/github/last-commit/sdcTools/cellKey.svg?logo=github)](https://github.com/sdcTools/cellKey/commits/master)
[![GitHub code size in
bytes](https://img.shields.io/github/languages/code-size/sdcTools/cellKey.svg?logo=github)](https://github.com/sdcTools/cellKey)

an **R** package for applying noise to statistical tables.

### Information

This package is developed within the SGA `Open source tools for
perturbative confidentiality methods`. This package is not yet tested or
optimized but already contains the core-functionality to compute
perturbed count-/magitude tables with some pre-defined hierarchies.

We have a first rough version with which interested users may play
around. Feedback (via issues) with regards to bugs or features requests
are very welcome as well as pull-requests. One the package is deemed
stable, a version will be released on CRAN too.

### Installation

The package can directly be installed from `github` using the `devtools`
package which must be installed on your system.

    if (!require("devtools")) install.packages("devtools")
    library(devtools)
    
    # update all packages
    update.packages(ask = FALSE)
    
    # finally install directly from github.com
    # install ptable
    install_github(
      repo = "sdcTools/ptable", 
      dependencies = c("Depends","Imports"), 
      force = TRUE, 
      build_opts = "--build-vignettes"
    )
      
    # install cellKey
    install_github(
      repo = "sdcTools/cellKey", 
      force = TRUE, 
      build_opts = "--build-vignettes"
    )

If you experience a timeout due to a proxy server while downloading, one
can work around this issue by specifying the proxy-server using the
`httr` package:

    httr::set_config(use_proxy(url = "xxx.xxx.xxx.xxx, port = yy))

### Usage

An example using both possible input formats for perturbation tables is
given in the main-function of the packge `perturb_table()`

    library(cellKey)
    ?perturb_table

Once finished, the package will also contain a package vignette. The
unfinished introduction vignette can be looked at using
`cellKey::ck_vignette()` or via the automatically deployed documentation
by clicking
[**here**](https://sdctools.github.io/cellKey/articles/introduction.html).
The complete [**documentation**](https://sdctools.github.io/cellKey/) is
also updated automatically and can be browes online.

### Updates

Updates/Changes are listed
[**here**](https://sdcTools.github.io/cellKey/news/index.html).

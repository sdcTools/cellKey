
## cellKey

[![R-CMD-check](https://github.com/sdcTools/cellKey/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sdcTools/cellKey/actions/workflows/R-CMD-check.yaml)
[![GitHub last
commit](https://img.shields.io/github/last-commit/sdcTools/cellKey.svg?logo=github)](https://github.com/sdcTools/cellKey/commits/master)
[![GitHub code size in
bytes](https://img.shields.io/github/languages/code-size/sdcTools/cellKey.svg?logo=github)](https://github.com/sdcTools/cellKey)

an **R** package for applying noise to statistical tables.

### Information

This package is developed within the SGA
`Open source tools for perturbative confidentiality methods`. This
package is not fully optimized for speed but already contains the core
functionality to compute perturbed count- and magnitude tables with
complex hierarchies.

Feedback (via issues) with regards to bugs or features requests are very
welcome as well as pull-requests.

## Vignette

The core-concepts of the package and its application is described in the
package vignette that can be viewed after installation of the package
using `ck_vignette()` or
[**online**](https://sdctools.github.io/cellKey/articles/introduction.html).

### Installation

The package can directly be installed from `CRAN` using the
`install.packages()` function which automatically also installs the
required dependencies (such as
[`ptable`](https://github.com/sdcTools/ptable).

    # update your R installation
    install.packages("cellKey")

### Usage

An example on how to apply the package is provided in
`?cellKey::cellkey_pkg` where also all the available methods are
described.

    ?cellKey::cellkey_pkg

The package vignette can be looked at using `cellKey::ck_vignette()` or
via the automatically deployed documentation by clicking
[**here**](https://sdctools.github.io/cellKey/articles/introduction.html).
The complete [**documentation**](https://sdctools.github.io/cellKey/) is
also updated automatically and can be viewed online.

### Updates

Updates/Changes are listed
[**here**](https://sdcTools.github.io/cellKey/news/index.html).


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
perturbative confidentiality methods`. This package is not optimized but
already contains the core functionality to compute perturbed count- and
magnitude tables with complex hierarchies.

We have a first rough version with which interested users may play
around. Feedback (via issues) with regards to bugs or features requests
are very welcome as well as pull-requests. One the package is deemed
stable, a version will be released on CRAN too.

### Important Note:

Major parts of the package were rewritten in version `>=0.17` compared
to previous versions. The following changes are now described:

  - **Definition of hierarchies** One main change is that it is now
    required to directly use functionality from
    [`sdcHierarchies`](https://cran.r-project.org/package=sdcHierarchies)
    to setup hierarchies when defining of tables. The following changes
    are possibly required in existing code.
    
      - replace `ck_create_node()` with `hier_create()`
      - replace `ck_add_nodes()` with `hier_add()`
      - replace `ck_delete_nodes()` with `hier_delete()`
      - replace `ck_rename_nodes()` with `hier_rename()`

  - **Removing “abs”-input format** For the sake of simplification, the
    differentiation between the *“abs”* and *“destatis”* format was
    removed in version `>= 0.17`. Internally, only the “destatis” format
    is used.

  - **Simplification** In order to simplify the application of the
    package, the process of defining, modifying and perturbing tables
    was modified. The implementation in versions `>=0.17` is based on
    `R6` classes and is described in detail in the new and package
    vignette that can be viewed after installation of the package using
    `ck_vignette()` or
    [**online**](https://sdctools.github.io/cellKey/articles/introduction.html).

  - **Perturbation of magnitude tables** This features was removed in
    versions `0.17.x` and added again in versions `>= 0.18.0`. For
    examples please look at `?cellkey_pkg` or at the
    [**vignette**](https://sdctools.github.io/cellKey/articles/introduction.html).

### Installation

The package can directly be installed from `github` using the `remotes`
package which is pulled in as a dependency from the `devtools` package.
The following snippet installs the package:

    # update your R installation
    update.packages(ask = FALSE)
    
    # install cellKey from github.com
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github(
      repo = "sdcTools/cellKey",
      dependencies = TRUE,
      build_opts = "--no-resave-data",
      force = TRUE)

If you experience a timeout due to a proxy server while downloading, one
can work around this issue by specifying the proxy-server using the
`httr` package:

    httr::set_config(use_proxy(url = "xxx.xxx.xxx.xxx, port = yy))

### Usage

An example on how to apply the package is provided in
`?cellKey::cellkey_pkg` where also all the available methods are
described.

    ?cellKey::cellkey_pkg

The package vignette is currently work-in-progress. It can be looked at
using `cellKey::ck_vignette()` or via the automatically deployed
documentation by clicking
[**here**](https://sdctools.github.io/cellKey/articles/introduction.html).
The complete [**documentation**](https://sdctools.github.io/cellKey/) is
also updated automatically and can be viewed online.

### Updates

Updates/Changes are listed
[**here**](https://sdcTools.github.io/cellKey/news/index.html).

## cellKey

an **R** package for applying noise to statistical tables. 

### Information
This package is developed within the SGA `Open source tools for perturbative confidentiality methods`. This package is not yet tested or optimized but already contains the core-functionality to compute perturbed count-/magitude tables with some pre-defined hierarchies.

We have a first rough version with which interested users may play around. Feedback (via issues) with regards to bugs or features requests are very welcome as well as pull-requests.

### Installation
The package can directly be installed from `github` 
```
devtools::install_github("sdcTools/cellKey", build_vignette=TRUE)
```

### Usage

A simple example is given in the main-function `perturbTable()`

```
library(cellKey)
?perturbTable
```

Once finished, the package will also contain a package vignette. The unfinished introduction vignette can be looked at using the following command:

```
vignette("introduction", pa="cellKey")
```
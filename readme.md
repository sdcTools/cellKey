## cellKey

an **R** package for applying noise to statistical tables. 

### Information
This package is developed within the SGA `Open source tools for perturbative confidentiality methods`. This package is not yet tested or optimized but already contains the core-functionality to compute perturbed count-/magitude tables with some pre-defined hierarchies.

We have a first rough version with which interested users may play around. Feedback (via issues) with regards to bugs or features requests are very welcome as well as pull-requests.

### News:

#### Version 0.9.0
- new dynamic way to specify hierarchies for tables, for an example see `?ck_manage_hierarchies`. This functionality will eventually also find its way to **sdcTable**

#### Version 0.8.2
- bugfix: rkeys need not to be integer if the "destatis"-method is used

#### Version 0.8.1
- small fixes and some exported function gained `verbose` arguments

#### Version 0.8.0
Functionality was added to specify the perturbation tables (`pTable`) in two different formats. The (default) way is to specify it as described in the original ABS-paper *Methodology for the Automatic Confidentialisation of Statistical Outputs from Remote Servers at the Australian Bureau of Statistics* (Thompson, Broadfoot, Elazar). An alternative way is to provide the perturbation tables for count tables in the *"destatis"*-format. `ck_create_pTable(type="destatis")` returns an exemplary pTable in this format. In the future, such pTables will likely be generated from another package. As the requirements regarding record keys are different in the following lookup-approach, we have already implemented some (basic) checks for validity of record keys when they are already available in the microdata used in `ck_create_input()`.

### Installation
The package can directly be installed from `github` 
```
devtools::install_github("sdcTools/cellKey", build_vignette=TRUE)
```

### Usage

An example using both possible input formats for perturbation tables is given in the main-function of the packge `perturbTable()`

```
library(cellKey)
?perturbTable
```

Once finished, the package will also contain a package vignette. The unfinished introduction vignette can be looked at using the following command:

```
vignette("introduction", pa="cellKey")
```

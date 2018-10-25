## cellKey

an **R** package for applying noise to statistical tables. 

### Information
This package is developed within the SGA `Open source tools for perturbative confidentiality methods`. This package is not yet tested or optimized but already contains the core-functionality to compute perturbed count-/magitude tables with some pre-defined hierarchies.

We have a first rough version with which interested users may play around. Feedback (via issues) with regards to bugs or features requests are very welcome as well as pull-requests.

### News:
#### Version 0.13.5
- feature: new method `ck_export_table()` that allows to save results in a simple format
- improvement: better error-message if too large values for `bigN` are specified
- improvement: better error-message if parameter `smallN` is too large in respect to the specified pTable
- improvement: display message about ignored parameters in `ck_generate_rkeys()` only if non-required parameters have been actually specified
- improvement: no warning messages that parameters are ignored in case they are irrelevant

#### Version 0.13.4
- feature: new function `ck_cnt_measures_basic()` that computes infoloss/utility measures based on two input vectors referring to original and perturbed values
- bugfix: check that record-keys in destatis-format are >= 0

#### Version 0.13.3
- feature: perturbation parameters for magnitude tables can be left empty
- new function `ck_cnt_measures()` that computes some (distance-based) information loss measures for count variables
- updated vignette and examples

#### Version 0.13.2
- feature: new method `print()` for objects returned from `perturbTable()`
- feature: new method `summary()` for objects returned from `perturbTable()`

#### Version 0.13.1
- small updates to reflect changes in [**ptable**](http://github.com/tenderle/ptable)

#### Version 0.13.0
- depend on package [**ptable**](http://github.com/sdcTools/ptable) to generate the perturbation tables by rewriting `ck_create_pTable()`; thus the package must be installed, e.g using `devtools::install_github("sdcTools/ptable", build_vignette=FALSE)`

#### Version 0.12.0
- feature: use package [**ptable**](http://github.com/sdcTools/ptable) to generate pTables in destatis format

#### Version 0.11.0
- feature: use functions to create hierarchies directly from `sdcTable` and bump version requirement of this package to `>=0.23` 

#### Version 0.10.2
- feature: if a (valid) variable is specified in argument `by` `in perturbTable()` it is automatically added to `countVars` even though not explicitely specified.

#### Version 0.10.1
- feature: `perturbTable()` gained an optional new argument `by`. In this argument one can use a variable that must also be listed in `countVars`. This variable is then used to compute the magnitute tables *by* the given 0/1 binary variable. For an example see `?perturbTable`.

#### Version 0.10.0
- feature: new argument `countVars` in `perturbTable()` which allows to additionally tabulate any number or 0/1 variables. For such variables. In such case, the record-keys of non-contribution units are set to 0 prior to the lookup in the perturbation table
- removed method `results()` and replaced it with new methods `ck_freq_table()` and `ck_cont_table()` that should be used to query specific tables from the output of `perturbTable()`
- updated examples showing new features in `perturbTable()`, `ck_freq_table()` and `ck_cont_table()`
- updated introduction vignette

#### Version 0.9.1
- bugfix for continous tables using pre-specified record keys

#### Version 0.9.0
- feature: new dynamic way to specify hierarchies for tables, for an example see `?ck_manage_hierarchies`. This functionality will eventually also find its way to **sdcTable**

#### Version 0.8.2
- bugfix: rkeys need not to be integer if the "destatis"-method is used

#### Version 0.8.1
- small fixes and some exported function gained `verbose` arguments

#### Version 0.8.0
- feature: perturbation tables (`pTable`) can now be specified in two different formats. The (default) way is to specify it as described in the original ABS-paper *Methodology for the Automatic Confidentialisation of Statistical Outputs from Remote Servers at the Australian Bureau of Statistics* (Thompson, Broadfoot, Elazar). An alternative way is to provide the perturbation tables for count tables in the *"destatis"*-format. `ck_create_pTable(type="destatis")` returns an exemplary pTable in this format. In the future, such pTables will likely be generated from another package. As the requirements regarding record keys are different in the following lookup-approach, we have already implemented some (basic) checks for validity of record keys when they are already available in the microdata used in `ck_create_input()`.

### Installation
The package can directly be installed from `github` 

```
devtools::install_github("sdcTools/ptable", build_vignette=FALSE)
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

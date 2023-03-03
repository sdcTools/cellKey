# cellKey 1.0.0
- first version on CRAN
- updated due to changes in Package `ptable`
- parallel computation is controlled using environment-variable `CK_RUN_PARALLEL`. If this is set
to `TRUE`, parallel computation is enabled, otherwise it is disabled. By default, parallel computing is
disabled.

# cellKey 0.19.2
- do not rename columns of `ptable` input object
- fix vignette and tests
- bump requirements

# cellKey 0.19.1
- performance improvement: do not compute max contributions if no numeric key variable was specified

# cellKey 0.19.0
- feature: new method `$supp_cells()` that allows to specifiy sensitive cells based on names
- bugfix: replace `sign()` with `ifelse` to enforce perturbation of cells that require additional protection

# cellKey 0.18.3
- bugfix: do not perturb cells with value 0 using flex-approach
- bugfix: use only actual number of contributors in flex-approach when perturbing numvars
- performance-improvement: speed-up computation of contributing units to cells
- bugfix: scrambling cell keys having not enough digits
- bugfix: computation separation point in flex-approach
- set default value of argument `w = NULL` in `ck_setup()`
- improve some error messages/outputs
- bugfix when looking up perturbation values; small_cells and others must not intersect
- fixing issues with empty cells for magnitude tables
- fixed computation of weighted spreads, thx @staudtlex
- do not perturb cells with 0 contributors in the flex-approach
- update tests due to updates in `digest`-pkg
- correctly compute weighted spread
- do not convert variable names to lowercase
- document R6 methods/classes via `roxygen2`
- new method `hierarchy_info()` containing some important information for each dimension
- reference new methods `create_cnt_ptable()` and `create_num_ptable()` from [ptable-pkg](https://github.com/sdcTools/ptable)
- update vignette
- improve documentation of `ck_params_nums()`
- update tests due to updates in ptable-pkg

# cellKey 0.18.2
- feature: allow objects from `ptable::pt_create_pParams` as input in `ck_params_cnts()`
- feature: allow objects from `ptable::pt_create_pParams` as input in `ck_params_nums()`

# cellKey 0.18.1
- bugfix: fallback to use a single core on windows-machines, fixing [`issue #131`](https://github.com/sdcTools/UserSupport/issues/131)
- updating dependencies and required versions
- fix vignette due to updates in `ptable`-pkg
- simplify examples by using examplary ptable from `ptable`-pkg using `ptable::pt_ex_cnts()` and `ptable::pt_ex_nums()`
- feature: allow tabulation of non-perturbed variables in `freqtab()`
- feature: allow tabulation of non-perturbed variables in `numtab()`
- code linting
- bugfix: fixing issues with "simple" approach; harmonizing and code-cleanup

# cellKey 0.18.0
- allow to save perturbation-schemes for different variables in `params_cnts_set()` and `params_nums_set()`
- allow return current active perturbation parameters for variables with `params_cnts_get()` and `params_nums_get()`
- added new convenience methods `allvars()`, `numvars()` and `cntvars()` returning variable names eligable for perturbation
- implemented the perturbation of numerical variables
  * new method `ck_params_nums()` to define perturbation parameters for continuous variables along with helper-functions `ck_flexparams()` and `ck_simpleparams()`
  * new method `numtab()` to extract numerical tables
  * new method `mod_nums()` returning modifications for numerical variables
- updated methods `print()` and `summary()` to include information about perturbed continuous variables
- new methods `reset_cntvars()`, `reset_numvars()` and `reset_allvars()` to remove perturbation results and provided perturbation parameters
- new methods to identify sensitive cells
  * `$supp_freq(v, max_n)`
  * `$supp_nk(v, max_n)`
  * `$supp_p(v, max_n)`
  * `$supp_pq(v, max_n)`
- added test-cases and improved coverage
- make use of ptables from [`ptable`](https://github.com/sdctools/ptable)-pkg
- Reproducibility:
  * allow to write perturbation parameters as yaml in `ck_params_nums()` and `ck_params_cnts()`
  * allow to import such parameters with `ck_read_yaml()`
- updated and extended package vignette

# cellKey 0.17.3
- correctly compute perturbed weighted counts

# cellKey 0.17.2
- fixed issue in lookup up values for frequency tables

# cellKey 0.17.1
- adding parameter `exclude_zero_cells` to `ck_cnt_measures()`
- updated documentation

# cellKey 0.17.0
- force usage of functionality from [`sdcHierarchies`](https://github.com/bernhard-da/sdcHierarchies) to define hierarchies
- removed features to perturb magnitude tables for now as parametrisation from `ptable`-pkg is not yet defined
- removed possiblity to specify parameters for count variables using the `ABS` definition
- removed `by`-argument in `$perturb()` method
- rewrite frequency table perturbation using `R6` classes
    * new function `ck_setup()` to define a table
    * allow multiple variables in method `$perturb()`-method
    * removed `ck_export_table()` and added arguments to method `freqtab()`
    * new method `$print()` for R6 objects
    * new method `$summary()` for R6 objects
    * new method `$mod_cnts()` returning modifications for count variables
    * new method `$params_cnts()` that allow to query and set count parameters
- updated `ck_cnt_measures()`
    * renamed `false_positives` to `false_nonzero`
    * improved documentation
    * add table of to `ck_cnt_measures` showing exact perturbations
    * harmonized output with tau-argus
- new method `$measures()` uses `ck_cnt_measures()` internally for count variables
- updated unit tests for counts based on hashes
- updated of package vignette

# cellKey 0.16.3
- fix tests due to changes in R 3.6.0

# cellKey 0.16.2
- install `ptable` from personal fork until `sdcTools/ptable` is updated

# cellKey 0.16.1
- removed placeholder for `pThreshold` in `perturbTable()`

# cellKey 0.16.0
- make use of new package `sdcHierarchies` to generate and update hierarchies
- new function `ck_rename_nodes()`
- `perturbTable()` got a new argument `pThreshold` that allows to specify a threshold above no perturbation is applied independent from the perturbation table. Currently only a placeholder and not used.

# cellKey 0.15.0
- new convenience function `ck_vignette()` that displays the package vignette in a browser
- `ck_generate_rkeys()` got a new argument `seed` that allows to overwrite the default seed computed from a hash of the input dataset.
- improvements in code-styling and readability of examples
- improvements in vignette

# cellKey 0.14.0
- feature: new method `ck_export_table()` that allows to save results in a simple format
- improvement: better error-message if too large values for `bigN` are specified
- improvement: better error-message if parameter `smallN` is too large in respect to the specified pTable
- improvement: display message about ignored parameters in `ck_generate_rkeys()` only if non-required parameters have been actually specified
- improvement: no warning messages that parameters are ignored in case they are irrelevant

# cellKey 0.13.4
- feature: new function `ck_cnt_measures_basic()` that computes infoloss/utility measures based on two input vectors referring to original and perturbed values
- bugfix: check that record-keys in destatis-format are >= 0

# cellKey 0.13.3
- feature: perturbation parameters for magnitude tables can be left empty
- new function `ck_cnt_measures()` that computes some (distance-based) information loss measures for count variables
- updated vignette and examples

# cellKey 0.13.2
- feature: new method `print()` for objects returned from `perturbTable()`
- feature: new method `summary()` for objects returned from `perturbTable()`

# cellKey 0.13.1
- small updates to reflect changes in [**ptable**](https://github.com/tenderle/ptable)

# cellKey 0.13.0
- depend on package [**ptable**](https://github.com/sdcTools/ptable) to generate the perturbation tables by rewriting `ck_create_pTable()`; thus the package must be installed, e.g using `devtools::install_github("sdcTools/ptable", build_vignette=FALSE)`

# cellKey 0.12.0
- feature: use package [**ptable**](https://github.com/sdcTools/ptable) to generate pTables in destatis format

# cellKey 0.11.0
- feature: use functions to create hierarchies directly from `sdcTable` and bump version requirement of this package to `>=0.23`

# cellKey 0.10.2
- feature: if a (valid) variable is specified in argument `by` `in perturbTable()` it is automatically added to `countVars` even though not explicitely specified.

# cellKey 0.10.1
- feature: `perturbTable()` gained an optional new argument `by`. In this argument one can use a variable that must also be listed in `countVars`. This variable is then used to compute the magnitute tables *by* the given 0/1 binary variable. For an example see `?perturbTable`.

# cellKey 0.10.0
- feature: new argument `countVars` in `perturbTable()` which allows to additionally tabulate any number or 0/1 variables. For such variables. In such case, the record-keys of non-contribution units are set to 0 prior to the lookup in the perturbation table
- removed method `results()` and replaced it with new methods `ck_freq_table()` and `ck_cont_table()` that should be used to query specific tables from the output of `perturbTable()`
- updated examples showing new features in `perturbTable()`, `ck_freq_table()` and `ck_cont_table()`
- updated introduction vignette

# Version 0.9.1
- bugfix for continous tables using pre-specified record keys

# cellKey 0.9.0
- feature: new dynamic way to specify hierarchies for tables, for an example see `?ck_manage_hierarchies`. This functionality will eventually also find its way to **sdcTable**

  # cellKey 0.8.2
  - bugfix: rkeys need not to be integer if the "destatis"-method is used

# cellKey 0.8.1
- small fixes and some exported function gained `verbose` arguments

# cellKey 0.8.0
- feature: perturbation tables (`pTable`) can now be specified in two different formats. The (default) way is to specify it as described in the original ABS-paper *Methodology for the Automatic Confidentialisation of Statistical Outputs from Remote Servers at the Australian Bureau of Statistics* (Thompson, Broadfoot, Elazar). An alternative way is to provide the perturbation tables for count tables in the *"destatis"*-format. `ck_create_pTable(type="destatis")` returns an exemplary pTable in this format. In the future, such pTables will likely be generated from another package. As the requirements regarding record keys are different in the following lookup-approach, we have already implemented some (basic) checks for validity of record keys when they are already available in the microdata used in `ck_create_input()`.

# cellKey 0.19.0
- [todo] allow grids for magnitude tables
- [todo] allow special requirements for positive cells
- [todo] magnitude_tables: 
  * [todo] implement suitable utility/risk measures
- [todo] add parameter `pos_neg_var` again  
  
# cellKey 0.18.0
- allow to save perturbation-schemes for different variables in `$params_cnts_set()` and `$params_nums_set()`
- allow return current active perturbation-schemes for variables with `params_cnts_get()` and `params_nums_get()`
- new methods `$numvars()` and `$cntvars()` returning variable names eligable for perturbation
- do not allow to change parameters once a variable has been perturbed
- new method `$params_nums()` that allow to query and set perturbation parameters for magnitude variables
- added maximum value of `6` in `params_num()` for parameter `top_k`
- compute cell keys for both variants (all units and all units contributing to a given numerical key variable)
- use environment variable `cellkey_debug` for debugging: internally we have `.ck_debug_off()`, `.ck_debug_on()` and can use `ck_log()`
- export and document (helper)functions `ck_gridparams()` and `ck_flexparams()`
- implement and document method `$numtab()` to extract numerical tables
- added parameter `mean_before_sum` in `$numtab()`
- new convinience method `$allvars()` returning a list form `$numvars()` and `$cntvars()`
- new method `ck_params_nums()` to define perturbation parameters for continuous variables
- add `mu_c` to `ck_params_nums()`
- make use of `mu_c` in `$perturb()`
- updated method `$print()` to include information about perturbed continuous variables
- new methods `$reset_cntvars()`, `$reset_numvars()` and `$reset_allvars()` to remove perturbation results and provided perturbation parameters
- implemented method `$mod_nums()` returning modifications for numerical variables (for each `top_k` values)
- added argument `pos_neg_var` to `ck_params_nums()` and use it in `$perturb()` for continuous variables
- implemented special requirements for positive variables according to section 2.5.1
- update flex-function use logical argument `scaling`
- magnitutes of perturbation depend on argument `epsilon` in `ck_flexparams()`
- implemented grid inputs as list for each `top_k` value but removed it from being accessible
- changed (internal) format of ptable for magnitude tables
    * added column `type` with values `all` or `small_cells`
- modified flex-function to implement concept of "proportional flexing"
- setting `m_fixed_sq != NA` triggers lookup in extra block of provided perturbation table (`type == "small_cells"`)
- implemented and document the following methods to identify sensitive cells
  * `$supp_freq(v, max_n)`
  * `$supp_nk(v, max_n)`
  * `$supp_p(v, max_n)`
  * `$supp_pq(v, max_n)`
- generate recordkeys with 8 digits and update tests
- internally compute parameter m1_squared:
  * new helper function `.compute_m1sq()`
  * updated `ck_params_nums()` and examples accordingly
- implement scaling = no using fixed parameters
    * removed scaling argument from `ck_flexparams()` 
    * new aux-function `ck_simpleparams()` to be used in `ck_params_nums()`
    * new aux-function `.x_delta_simle()` for computing x_delta
- combine parameters `m_small` and `m_large` to `p` in `ck_flexparams()`
- created test-cases for peter-paul
- implement `parity` (different lookups for even/odd columns) based on even/odd number of contributors to a cell
  * for small-cells: never apply parity
  * Parity=Yes does not make sense with top_k > 1. This should not be allowed.    
- removed argument `pos_neg_var` for now
- compute parameter `zs` (which was `g1`) in `ck_params_nums`
- reparametrisation and bugfixes in the lookup-procedure for numeric variables
- return a perturbation value of 0 for all cells below the separation point which are not the largest contribution
- reimplment perturbation based on flex-function according to "confused_doc"
  * no scaling with epsilons in case x_j < z_s
  * m_1 = 1 in case x_j < z_s --> x_delta = x_j in case x_j < z_s
  * v_j = 0 if x_j < z_s and j > 1 --> \hat{X}_j = 0 in this case
- `.ck_nr_digits()` returns number of digits used for rounding from env-variable `.ck_nr_digits` (default 8)
- allow to export perturbation inputs (for cnt- and numvars) as yamls in `ck_params_cnts()` and `ck_params_nums()` and read them with `ck_read_yaml()`



- [todo] check primsupp rules (not weighted?)
- [todo] implement magnitude tables
  * [todo] add tests for magitude tables
  * [todo] update method `$summary()` to include values for perturbed magnitude tables
- [todo] update vignette
- [todo] remove or re-implement `by`-argument in ck_perturb()
- [todo] remove method `$everything()` which is only here for debugging purposes




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
- small updates to reflect changes in [**ptable**](http://github.com/tenderle/ptable)

# cellKey 0.13.0
- depend on package [**ptable**](http://github.com/sdcTools/ptable) to generate the perturbation tables by rewriting `ck_create_pTable()`; thus the package must be installed, e.g using `devtools::install_github("sdcTools/ptable", build_vignette=FALSE)`

# cellKey 0.12.0
- feature: use package [**ptable**](http://github.com/sdcTools/ptable) to generate pTables in destatis format

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

#' A real-world data set on persons
#'
#' 820000 obervations in 5 Variables without sampling weights.
#'
#' @name ck_dat_hc92
#' @docType data
#' @format ck_dat_hc92: a data frame with 820000 observations on the following 6 variables.
#' - `id`: a numeric identifier
#' - `geo_m`: a character vector defining regions
#' - `sex a`: character vector defining gender
#' - `age_m`: a character vector containing age groups
#' - `yae_h`: a character vector
#' - `rkey`: a numeric vector holding record keys
#' @references https://ec.europa.eu/eurostat/cros/content/3-random-noise-cell-key-method_en
#' @keywords datasets
#' @md
#' @examples
#' data(ck_dat_hc92)
#' head(ck_dat_hc92)
NULL

#' A real-world data set on household income and expenditures
#'
#' 4580 Obervations in 15 Variables; This dataset also contains sampling weights!
#'
#' @name testdata
#' @docType data
#' @format testdata: a data frame with 4580 observations on the following 15 variables.
#' - `urbrur`: a numeric vector
#' - `roof`: a numeric vector
#' - `walls`: a numeric vector
#' - `water`: a numeric vector
#' - `electcon`: a numeric vector
#' - `relat`: a numeric vector
#' - `sex`: a numeric vector
#' - `age`: a numeric vector
#' - `hhcivil`: a numeric vector
#' - `expend`: a numeric vector
#' - `income`: a numeric vector
#' - `savings`: a numeric vector
#' - `ori_hid`: a numeric vector
#' - `sampling_weight`: a numeric vector
#' - `household_weights`: a numeric vector
#' @references The International Household Survey Network, www.ihsn.org
#' @keywords datasets
#' @md
#' @examples
#' data(testdata)
#' head(testdata)
NULL

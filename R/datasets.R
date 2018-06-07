#' A real-world data set on persons
#'
#' 820000 obervations in 5 Variables without sampling weights.
#'
#' @name mdata
#' @docType data
#' @format mdata: a data frame with 820000 observations on the following 5 variables.
#' \describe{
#' \item{geo_m}{a character vector defining regions}
#' \item{sex}{a character vector defining gender}
#' \item{age_m}{a character vector containing age groups}
#' \item{yae_h}{a character vector}
#' \item{rkey}{a numeric vector holding record keys}}
#' @references https://ec.europa.eu/eurostat/cros/content/3-random-noise-cell-key-method_en
#' @keywords datasets
#' @examples
#' data(mdata)
#' head(mdata)
NULL

#' A real-world data set on household income and expenditures
#'
#' 4580 Obervations in 15 Variables; This dataset also contains sampling weights!
#'
#' @name testdata
#' @docType data
#' @format testdata: a data frame with 4580 observations on the following 15 variables.
#' \describe{
#' \item{urbrur}{a numeric vector}
#' \item{roof}{a numeric vector}
#' \item{walls}{a numeric vector}
#' \item{water}{a numeric vector}
#' \item{electcon}{a numeric vector}
#' \item{relat}{a numeric vector}
#' \item{sex}{a numeric vector}
#' \item{age}{a numeric vector}
#' \item{hhcivil}{a numeric vector}
#' \item{expend}{a numeric vector}
#' \item{income}{a numeric vector}
#' \item{savings}{a numeric vector}
#' \item{ori_hid}{a numeric vector}
#' \item{sampling_weight}{a numeric vector}
#' \item{household_weights}{a numeric vector}}
#' @references The International Household Survey Network, www.ihsn.org
#' @keywords datasets
#' @examples
#' data(testdata)
#' head(testdata)
NULL

#' A perturbation table based on DESTATIS format
#'
#' 66 rows in 6 Variables
#'
#' @name ptable_destatis
#' @docType data
#' @format ptable_destatis: a data frame with 88 rows on the following 6 variables.
#' \describe{
#' \item{i}{a numeric vector}
#' \item{j}{a numeric vector}
#' \item{p}{a numeric vector}
#' \item{kum_p_u}{a numeric vector}
#' \item{kum_p_o}{a numeric vector}
#' \item{diff}{a numeric vector}}
#' @keywords datasets
#' @examples
#' data(ptable_destatis)
#' head(ptable_destatis)
NULL
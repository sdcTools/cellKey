#' ck_vignette
#'
#' starts the package vignette that gets you started with the package
#'
#' @return a browser windows/tab with showing the vignette
#' @export
#' @examples
#' \dontrun{
#' ck_vignette()
#' }
ck_vignette <- function() {
  RShowDoc("introduction", package = "cellKey")
}

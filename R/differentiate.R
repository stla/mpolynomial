#' Title
#'
#' @param poly 
#' @param coefficients 
#' @param powers 
#' @param diff 
#'
#' @return
#' @export
#'
#' @examples
differentiate <- function(poly, coefficients=NULL, powers=NULL, diff){
  if(!missing(poly)){
    burkardt_polynomial_dif(poly$value, poly$index, diff)
  }
}
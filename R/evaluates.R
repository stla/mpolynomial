#' Evaluate a polynomial
#' @description Evaluation of a multivariate polynomial at some values of the 
#' variables.
#' 
#' @param poly a polynomial given as a list of two elements: the powers in an 
#' integer matrix and the coefficients of the monomials composing the polynomial, 
#' like a \code{\link[spray]{spray}} object
#' @param coefficients,powers an alternative way to input the polynomial, to use 
#' with missing \code{poly} 
#' @param X a matrix whose rows are the evaluation points
#'
#' @return A numeric vector, the values of the polynomial at the points given by 
#' \code{X}.
#' @export
#'
#' @examples
evalPolynomial <- function(poly, coefficients=NULL, powers=NULL, X){
  if(!is.matrix(X)){
    X <- rbind(X)
  }
  if(!missing(poly)){
    if(is_valid_spray(poly) && ncol(X) == ncol(poly[[1]])){
      burkardt_polynomial_value(poly[[2]], poly[[1]], X)
    }
  }else{
    if(is_valid_spray(list(powers, coefficients)) && ncol(X) == ncol(powers)){
      burkardt_polynomial_value(coefficients, powers, X)
    }
  }
}
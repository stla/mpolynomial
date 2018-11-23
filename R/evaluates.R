#' Evaluate a polynomial
#' @description Evaluation of a multivariate polynomial at some values of the 
#' variables.
#' 
#' @param poly a polynomial given as a list of two elements: the powers in an 
#' integer matrix and the coefficients of the monomials composing the polynomial, 
#' like a \code{\link[spray]{spray}} object
#' @param coefficients,powers an alternative way to input the polynomial, to use 
#' with missing \code{poly} 
#'
#' @return A function with one argument, a matrix whose rows are the evaluation 
#' points, and which returns the values of the polynomial at these points.
#' @export
#'
#' @examples
#' f0 <- function(x,y) (x+y)^2 - 2*x*y
#' P <- f0(lone(1,2),lone(2,2))
#' f <- evalPolynomial(P)
#' f(rbind(c(1,1),c(2,2),c(3,3)))
evalPolynomial <- function(poly, coefficients=NULL, powers=NULL){
  if(missing(poly)){
    poly <- list(powers, coefficients)
  }
  function(X){
    if(!is.matrix(X)){
      X <- rbind(X)
    }
    if(is_valid_spray(poly) && ncol(X) == ncol(poly[[1]])){
      burkardt_polynomial_value(poly[[2]], poly[[1]], X)
    }else{
      stop(sprintf("X should have %d columns", ncol(poly[[1]])))
    }
  }
  # if(!missing(poly)){
  #   if(is_valid_spray(poly) && ncol(X) == ncol(poly[[1]])){
  #     burkardt_polynomial_value(poly[[2]], poly[[1]], X)
  #   }
  # }else{
  #   if(is_valid_spray(list(powers, coefficients)) && ncol(X) == ncol(powers)){
  #     burkardt_polynomial_value(coefficients, powers, X)
  #   }
  # }
}
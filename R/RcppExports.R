# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title A realization of stochastic gradient descent(SGD) algorithm using Rcpp
#' @description A SGD realization using Rcpp
#' @param matrix the features of training data
#' @param result the estimation of the interesting feature
#' @param theta initialization of paramaters
#' @param loss initialization of flag variable to judge the convergency
#' @return the optimized paramaters \code{n}
#' @examples
#' \dontrun{
#' matrix <- matrix(c(1,4,2,5,5,1,4,2), nrow = 4, ncol = 2)
#' result <- c(19,26,19,20)
#' theta <- c(2,5)
#' loss <- 1000.0
#' res <- SGD_C(matrix, result, theta, loss)
#' print(res)
#' }
#' @import Rcpp
#' @useDynLib StatComp22004
#' @export
SGD_C <- function(matrix, result, theta, loss) {
    .Call('_StatComp22004_SGD_C', PACKAGE = 'StatComp22004', matrix, result, theta, loss)
}


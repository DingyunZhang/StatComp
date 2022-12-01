#' @title Variational Method in 'Sparse Reduced-Rank Regression for Simultaneous Dimension Reduction and Variable Selection'
#' @description The algorithm is described in "https://www.tandfonline.com/doi/full/10.1080/01621459.2012.734178".
#' @param X predictor matrix (matrix)
#' @param Y response matrix (matrix)
#' @param lambda eigenvalue vector (vector)
#' @param r rank of Bt(A)
#' @return the estimated A and B \code{n}
#' @examples
#' \dontrun{
#' r <- 2
#' p <- 6
#' q <- 6
#' lambda <- rep(0.5, p)
#' X <- matrix(0, nrow = n, ncol = p)
#' Y <- matrix(0, nrow = n, ncol = q)
#' res <- SRRR.variational(X, Y, lambda, r)
#' }
#' @import MASS
#' @export

SRRR.variational <- function(X, Y, lambda, r){
  obj <- function(Y, X, B, A, lambda){
    norm(Y-X%*%B%*%t(A),type = "F")^2 + sum(lambda*apply(B,1,function(x) sqrt(sum(x^2))))
  }
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  Sxx <- t(X)%*%X/n
  Sxy <- t(X)%*%Y/n
  Syx <- t(Y)%*%X/n
  V <- eigen(Syx%*%solve(Sxx)%*%Sxy,symmetric = TRUE)$vector[,1:r]
  A <- V
  B <- solve(Sxx)%*%Sxy%*%V
  A.hat <- matrix(numeric(q*r), nrow = q, ncol = r)
  B.hat <- matrix(numeric(p*r), nrow = p, ncol = r)
  
  mu <- numeric(p)
  
  while(abs(obj(Y,X,B.hat,A.hat,lambda)-obj(Y,X,B,A,lambda))>1e-5){
    A.hat <- A
    svd.yxb <- svd(t(Y)%*%X%*%B)
    u <- svd.yxb$u
    v <- svd.yxb$v
    A <- u%*%t(v)
    
    while(norm(B-B.hat,type = "F")>1e-3){
      B.hat <- B
      mu <- apply(B,1,function(x) 1/sqrt(sum(x^2)))
      B <- ginv(t(X)%*%X+diag(lambda*mu)/2)%*%t(X)%*%Y%*%A
    }
  }
  
  return(list(A=A,B=B))
}
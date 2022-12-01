## -----------------------------------------------------------------------------
# Data Initialization

knitr::opts_chunk$set(echo = TRUE)
library(MASS)

set.seed(22004)
n = 20
p = 6
q = 6
p0 = 1
r = 2
rho.x = 0.1
rho.e = 0.1
lambda <- rep(0.5, p)

sample.generator <- function(n, p, p0, q, r, rho.x, rho.e){
  X <- matrix(0, nrow = n, ncol = p)
  Y <- E <- matrix(0, nrow = n, ncol = q)
  B <- matrix(c(rnorm(p0*r),rep(0,(p-p0)*r)), nrow = p, ncol = r, byrow = TRUE)
  A <- matrix(rnorm(q*r), nrow = q, ncol = r)
  C <- B%*%t(A)
  
  sigma.x <- matrix(rho.x, nrow = p, ncol = p)
  sigma.x <- sigma.x+diag(1-rho.x, p)
  X <- mvrnorm(n, mu = rep(0,p), Sigma = sigma.x)
  
  sigma.e <- matrix(rho.e, nrow = q, ncol = q)
  sigma.e <- sigma.e+diag(1-rho.e, q)
  sgm2 <- sum(diag(t(C)%*%sigma.x%*%C))/q
  E <- mvrnorm(n, mu = rep(0,q), Sigma = sgm2*sigma.e)
  
  Y <- X%*%C+E
  
  return(list(X=X,Y=Y, E=E, C=C, B=B, A=A))
}

X <- sample.generator(n, p, p0, q, r, rho.x, rho.e)$X
Y <- sample.generator(n, p, p0, q, r, rho.x, rho.e)$Y
E <- sample.generator(n, p, p0, q, r, rho.x, rho.e)$E
C <- sample.generator(n, p, p0, q, r, rho.x, rho.e)$C
A <- sample.generator(n, p, p0, q, r, rho.x, rho.e)$A

## -----------------------------------------------------------------------------
# Variational Method
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

## -----------------------------------------------------------------------------
# Test
SRRR.variational(X, Y, lambda, r)

## ----warning=FALSE------------------------------------------------------------
# Test
library(Rcpp)

# Data Initialization
matrix <- matrix(c(1,4,2,5,5,1,4,2), nrow = 4, ncol = 2)
result <- c(19,26,19,20)
theta <- c(2,5)
loss <- 1000.0

sourceCpp('C:/Users/10313/Desktop/StatisticalComputing/HW/HW11/SGD_C.cpp')
SGD_C(matrix, result, theta, loss)



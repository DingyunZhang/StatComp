#include <Rcpp.h>
using namespace Rcpp;

//' @title A realization of stochastic gradient descent(SGD) algorithm using Rcpp
//' @description A SGD realization using Rcpp
//' @param matrix the features of training data
//' @param result the estimation of the interesting feature
//' @param theta initialization of paramaters
//' @param loss initialization of flag variable to judge the convergency
//' @return the optimized paramaters \code{n}
//' @examples
//' \dontrun{
//' matrix <- matrix(c(1,4,2,5,5,1,4,2), nrow = 4, ncol = 2)
//' result <- c(19,26,19,20)
//' theta <- c(2,5)
//' loss <- 1000.0
//' res <- SGD_C(matrix, result, theta, loss)
//' print(res)
//' }
//' @export
// [[Rcpp::export]]
NumericVector SGD_C(IntegerMatrix matrix, NumericVector result, NumericVector theta, float loss){
  for(int i =0 ;i<100 && loss>0.001; ++i){
    float error_sum = 0.0;
    int j=i%4;//j=0,1,2,3
    float h = 0.0;
    for(int k=0; k<2; ++k){
      h += matrix(j,k)*theta[k];
    }
    error_sum = result[j]-h;
    for(int k=0;k<2;++k){
      theta[k] += 0.01*(error_sum)*matrix(j,k);
    }
    printf("parameters:%f,%f\n",theta[0],theta[1]);
    float loss = 0.0;
    for(int j = 0;j<4;++j){
      float sum=0.0;
      for(int k = 0;k<2;++k){
        sum += matrix(j,k)*theta[k];
      }
      loss += (sum-result[j])*(sum-result[j]);
    }
    printf("iter:%d\n",i);
    printf("loss:%f\n",loss);
  }
  return theta;
}
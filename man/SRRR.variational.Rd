% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SRRR.variational.R
\name{SRRR.variational}
\alias{SRRR.variational}
\title{Variational Method in 'Sparse Reduced-Rank Regression for Simultaneous Dimension Reduction and Variable Selection'}
\usage{
SRRR.variational(X, Y, lambda, r)
}
\arguments{
\item{X}{predictor matrix (matrix)}

\item{Y}{response matrix (matrix)}

\item{lambda}{eigenvalue vector (vector)}

\item{r}{rank of Bt(A)}
}
\value{
the estimated A and B \code{n}
}
\description{
The algorithm is described in "https://www.tandfonline.com/doi/full/10.1080/01621459.2012.734178".
}
\examples{
\dontrun{
r <- 2
p <- 2
q <- 2
n <- 2
lambda <- rep(0.5, p)
X <- matrix(c(1.0604068, -0.454802, 0.7287505,-2.014123), nrow = n, ncol = p)
Y <- matrix(c(-0.5825189, 0.3734626, -0.4302581, -0.1078283), nrow = n, ncol = q)
res <- SRRR.variational(X, Y, lambda, r)
}
}

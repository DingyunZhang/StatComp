## -----------------------------------------------------------------------------
# 加载包
knitr::opts_chunk$set(echo = TRUE)
library(snowfall)
library(ggplot2)
library(tibble)
library(bootstrap)   
library(corrplot)     
library(boot)     
library(DAAG)
library(Rcpp)
library(microbenchmark)

## -----------------------------------------------------------------------------
x <- c(rep(1, 11),rep(2, 7),rep(3, 9),rep(4, 4),rep(5, 6))
b <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
a <- c("A", "B", "C", "D", "E")
d <- terrain.colors(5)
hist(x, breaks = b, labels = a, col = d)

## -----------------------------------------------------------------------------
set.seed(22004)
pareto <- function(a,b) {
  # Generate random numbers with cdf F(x)
  u <- runif(10000)
  x <- b*(1-u)^(-1/a)
  
  # Draw the histogram of random numbers generated
  hist(x, prob = TRUE, main = paste('Pareto(',a,',',b,')'))
  
  # Draw the density function f(x)
  y <- seq(0, max(x), 0.1)
  lines(y, a*b^a/(y^(a+1)))
}

pareto(2, 2)

## -----------------------------------------------------------------------------
set.seed(22004)
beta <- function(a,b) {
  # Calculate constant c
  x0 <- (a-1)/(a+b-2)
  c <- x0^(a-1)*(1-x0)^(b-1)  # constant in pdf can be ignored
  
  # Generate random numbers with pdf f(x)
  n <- 10000
  k <- 0
  y <- numeric(n)
  while (k < n) {
    u <- runif(1)
    x <- runif(1) # random variate from g(x)
    if (x^(a-1)*(1-x)^(b-1) / c > u) {
      # accept x
      k <- k + 1
      y[k] <- x
    }
  }
  
  # Draw the histogram of random numbers generated
  hist(y, prob = TRUE, main = paste('Beta(',a,',',b,')'), xlab = "x")
  
  # Draw the density function f(x)
  z <- seq(0, 1, 0.01)
  lines(z, z^(a-1)*(1-z)^(b-1)*gamma(a+b)/gamma(a)/gamma(b))
}

beta(3, 2)

## -----------------------------------------------------------------------------
set.seed(22004)
expgamma <- function(r, beta) {
  # Generate random numbers from the mixture
  n <- 1000
  x <- rgamma(n, r, beta)
  y <- rexp(n, x)
  return(y)
}

r <- 4; beta <- 2
rnd = expgamma(r, beta)

## -----------------------------------------------------------------------------
# Draw the histogram of random numbers generated
hist(rnd, prob = TRUE, main = paste('Pareto(',r,',',beta,')'), xlab = "y")

# Draw the density function f(y)
y <- seq(0, max(rnd), 0.01)
lines(y, r*beta^r/(beta+y)^(r+1))

## -----------------------------------------------------------------------------
set.seed(0)
# This part is copied from bb
quick_sort <- function(x){
  num <- length(x)
  if(num==0||num==1){return(x)
  }else{
    a <- x[1]
    y <- x[-1]
    lower <- y[y<a]
    upper <- y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}#form a loop
}


test<-sample(1:1e4)
system.time(quick_sort(test))[1]
test <- quick_sort(test)
# show the result of fast sort algorithm
test[1:10]
test[9991:10000]

## -----------------------------------------------------------------------------
set.seed(0)
n <- c(1e4, 2e4, 4e4, 6e4, 8e4)
computation_time <- function(n){
  t <- numeric(100)
  set.seed(0)
  for(i in 1:100){
    test <- sample(1:n)
    t[i] <- system.time(quick_sort(test))[1]
  }
  t_mean <- mean(t)
  return(t_mean)
}


an <- c(computation_time(n[1]),computation_time(n[2]),computation_time(n[3]),
       computation_time(n[4]),computation_time(n[5]))
an

## -----------------------------------------------------------------------------
tn <- n*log(n)
mylm <- lm(an~tn)
x <- seq(0,1e6,length.out=100)
b <- coefficients(mylm)
plot(tn, an, main="Regression line")
lines(x, b[1]+b[2]*x, col="red")

## -----------------------------------------------------------------------------
set.seed(0)
m <- 1e4
U <- runif(m)
theta1 <- mean(exp(U))                # simple MC estimator
theta2 <- mean((exp(U)+exp(1-U))/2)   # antithetic variables estimator
var1 <- var(exp(U))                   # sample variance of simple MC
var2 <- var((exp(U)+exp(1-U))/2)      # sample variance of antithetic variables
theta1
theta2
100*(var1-var2)/var1      # empirical estimator of percent reduction of variance

## -----------------------------------------------------------------------------
g <-function(x) x^2/sqrt(2*pi)*exp(-x^2/2)*(x>1)

f1 <- function(x) 1/sqrt(2*pi)*exp(-x^2/2)

f2 <- function(x) x^2*exp(-x)/2

## -----------------------------------------------------------------------------
m <- 10000
theta.hat <- se<- numeric(2)
fg<-matrix(0,nrow=2,ncol=m)

# using f1
set.seed(15)
x <- rnorm(m)
fg[1,] <- g(x)/f1(x)
theta.hat[1] <-mean(fg[1,])
se[1] <- sd(fg[1,])

# using f2
set.seed(15)
x <- rgamma(m,3,1)
fg[2,] <- g(x)/f2(x)
theta.hat[2] <-mean(fg[2,])
se[2] <- sd(fg[2,])

## -----------------------------------------------------------------------------
rbind(theta.hat,se)

## -----------------------------------------------------------------------------
summary(fg[1,])
summary(fg[2,])

## -----------------------------------------------------------------------------
m <- 10000
g <- function(x) exp(-x)/(1+x^2)*(x>0)*(x<1)
f <- function(x){ exp(-x)/(1-exp(-1))*(x>0)*(x<1)}

# f3, inverse transform method
set.seed(15)
u <- runif(m)
x <- -log(1-u*(1-exp(-1)))
fg <- g(x)/f(x)
theta.im <- mean(fg)
se.im <-sd(fg)

## -----------------------------------------------------------------------------
set.seed(15)
k<-5
n<-m/k
theta_s <- var_s <-numeric(k)
for(i in 1:k){
  u <- runif(n,(i-1)/5,i/5)
  x <- -log(1-(1-exp(-1))*u)
  fg <- g(x)/k/f(x)
  theta_s[i]<-mean(fg)
  var_s[i]<-var(fg)
}

## -----------------------------------------------------------------------------
sum(theta_s)

## -----------------------------------------------------------------------------
sqrt(sum(var_s))

## -----------------------------------------------------------------------------
# sample generation function
sample_gen <- function(n, mu=0, sigma=1){
  x <- rlnorm(n=n, meanlog = mu, sdlog = sigma)
  return(x)
}

# data analysis function (constuct a confidence interval with level alpha)
CI <- function(x, alpha=0.05){
  n <- length(x)
  y <- log(x)
  mu.hat <- mean(y)
  sigma2.hat <- var(y)
  lower <- mu.hat+qt(alpha/2,df=n-1)*sqrt(sigma2.hat/n)
  upper <- mu.hat+qt(1-alpha/2,df=n-1)*sqrt(sigma2.hat/n)
  return(c("lower.bound"=lower,"upper.bound"=upper))
}

## -----------------------------------------------------------------------------
set.seed(0)
m <- 1e4
lower <- upper <- numeric(m)

for(i in 1:m){
  Sample <- sample_gen(n=10, mu=0, sigma=1)
  lower[i] <- CI(x=Sample)[1]
  upper[i] <- CI(x=Sample)[2]
}

CP <- mean((lower<0)&(upper>0))
cat("CP =",CP)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
# The functions of "Count Five" test is copied from the book
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}

count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

F.test <- function(x, y, alpha=0.05){
  S1 <- var(x)
  S2 <- var(y)
  m <- length(x)
  n <- length(y)
  f <- S2/S1
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(f>qf(1-alpha/2,df1 = n-1,df2 = m-1)||
                           f<qf(alpha/2,df1 = n-1,df2 = m-1)))
}

## -----------------------------------------------------------------------------
power_count5test <- function(m, n1, n2, sigma1, sigma2){
  mean(replicate(m, expr={
    x <- rnorm(n1, 0, sigma1)
    y <- rnorm(n2, 0, sigma2)
    count5test(x, y)
  }))
}

power_F.test <- function(m, n1, n2, sigma1, sigma2){
  mean(replicate(m, expr = {
    x <- rnorm(n1, 0, sigma1)
    y <- rnorm(n2, 0, sigma2)
    F.test(x, y, alpha = 0.055)
  }))
}

## -----------------------------------------------------------------------------
set.seed(0)
m <- 1e4
# generate samples under H1 to estimate power
sigma1 <- 1
sigma2 <- 1.5
result1 <- numeric(3)
result2 <- numeric(3)
n <- c(20,100,1000)

for(i in 1:3){
  result1[i] <- power_count5test(m, n1=n[i], n2=n[i], sigma1, sigma2)
  result2[i] <- power_F.test(m, n1=n[i], n2=n[i], sigma1, sigma2)
}


pander::pander(data.frame("size"=c(20,100,200),"count five test"=result1,
                          "F test"=result2))

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
mat <-
  matrix(c(6510, 3490, 10000, 6760, 3240, 10000, 13270, 6730, 20000), 3, 3,
         dimnames = list(
           c("Rejected", "Accepted", "total"),
           c("method A", "method B", "total")
         ))
mat

## -----------------------------------------------------------------------------
rm(list = ls())

library(boot)

Sample1 <- function(x){
  samp <- aircondit[x]
  samp
}

Rate1 <- function(samp, i) {
  rat <- 1/mean(as.matrix(samp[i, ]))
  rat
}

Result1 <- function(samp,func,Rr){
  bo <- boot(samp, statistic = func, R = Rr)
  print(bo)
}

set.seed(1234)
samp <- Sample1(1)
resu <- Result1(samp,Rate1,2000)

detach(package:boot)

rm(list = ls())

## -----------------------------------------------------------------------------
rm(list = ls())

library(boot)

Sample2 <- function(x){
  samp <- aircondit[x]
  samp
}

Meant2 <- function(x, i) {
  mea <- mean(as.matrix(x[i, ]))
  mea
}

Result2 <- function(samp,func,Rr){
  bo <- boot(samp, statistic = func, R = Rr)
  re <- boot.ci(bo, type = c("norm", "perc", "basic", "bca"))
  print(bo)
  print(re)
  hist(bo$t, prob = TRUE, main = " ")
  points(bo$t0, 0, cex = 2, pch = 16)
  bo
}

set.seed(1234)
samp <- Sample2(1)
resu <- Result2(samp,Meant2,2000)

detach(package:boot)

rm(list = ls())


## -----------------------------------------------------------------------------
rm(list = ls())

skewness <- function(x,i) {
  #computes the sample skewness coeff.
  x_bar <- mean(x[i])
  x_bar
}

Sample3 <- function(n, mea, sd){
  samp <- rnorm(n, mea, sd)
  samp
}

Analysis3 <- function(m, func, Rr, n, mea, sd){
  library(boot)
  nornorm <- matrix(0, m, 2)
  norbasi <- matrix(0, m, 2)
  norperc <- matrix(0, m, 2)
  for (i in 1:m) {
    Samp <- Sample3(n, mea, sd)
    Skew <- boot(Samp, statistic = func, R=Rr)
    Nor <- boot.ci(Skew, type=c("norm","basic","perc"))
    nornorm[i,] <- Nor$norm[2:3]
    norbasi[i,] <- Nor$basic[4:5]
    norperc[i,] <- Nor$percent[4:5]
  }
  #Calculate the coverage probability of a normal distribution
  norm <- mean(nornorm[,1] <= s & nornorm[,2] >= s)
  basi <- mean(norbasi[,1] <= s & norbasi[,2] >= s)
  perc <- mean(norperc[,1] <= s & norperc[,2] >= s)
  #Calculate the probability of the left side of the normal distribution
  normleft <- mean(nornorm[,1] >= s )
  basileft <- mean(norbasi[,1] >= s )
  percleft <- mean(norperc[,1] >= s )
  #Calculate the right side probability of a normal distribution
  normright <- mean(nornorm[,2] <= s )
  basiright <- mean(norbasi[,2] <= s )
  percright <- mean(norperc[,2] <= s )
  analyresu <- c(norm, basi, perc, normleft, basileft, percleft, normright, basiright, percright)
  analyresu
}

Result3 <- function(sd, analyresu){
  dnam <- paste("N ( 0 ,", as.character(sd^2),")",seq="")
  Distribution <- c(dnam)
  Type <- c("basic", "norm", "perc")
  Left <- analyresu[4:6]
  Right <- analyresu[7:9]
  P.coverage <- analyresu[1:3]
  result <- data.frame(Distribution, Type, Left, Right, P.coverage)
  result
}

s <- 0
n <- 20
m <- 1000
R <- 1000

mea <- 0
sd <- 3 

# We can set n, m, R, mea, sd any way we want.

set.seed(1234)
library(boot)

Analyresu <- Analysis3(m, skewness, R, n, mea, sd)
Resu <- Result3(sd, Analyresu)

knitr::kable (Resu, align="c")

rm(list = ls())

## -----------------------------------------------------------------------------
rm(list=ls())
invisible(gc())

library(bootstrap)
set.seed(22004)

bias_se.jack <- function(scor){
  #'*Compute original theta.hat*
  scor.cov <- cov(scor)
  ev <- eigen(scor.cov)$values
  theta.hat <- max(ev)/sum(ev)
  
  #'*Define function for each jackknife estimate*
  jack.scor <- function(scor,i){
    d.scor <- scor[-i,]
    d.scor.cov <- cov(d.scor)
    d.ev <- eigen(d.scor.cov)$values
    max(d.ev)/sum(d.ev)
  }
  
  #'*Iteration*
  n <- nrow(scor)
  theta.jack<- sapply(1:n,jack.scor,scor=scor)
  
  #'*Return list containing jackknife bias & se*
  list(
  bias.jack = (n-1)*(mean(theta.jack)-theta.hat),
  se.jack = sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2)))
}

print(bias_se.jack(scor))

detach(package:bootstrap)
rm(list=ls())
gc()


## -----------------------------------------------------------------------------
rm(list=ls())
invisible(gc())

library(DAAG)
attach(ironslag)
set.seed(22004)

model.validation <- function(y,x,valid.num){
  
  #'*Step 1:Define a function that computes the sum of error square for different model & validation points(could be leave-k-out)*
  validerror <- function(y,x,index,fun){
    #'*Delete validation points from sample*
    y_d <- y[-index]
    x_d <- x[-index]
    #'*Choose model according to fun*
    if(fun=="linear"){
      J <- lm(y_d ~ x_d)
      yhat <- J$coef[1] + J$coef[2] * x[index]
    }
    else if(fun=="quadratic"){
      J <- lm(y_d ~ x_d + I(x_d^2))
      yhat <- J$coef[1] + J$coef[2] * x[index] + J$coef[3] * x[index]^2
    }
    else if(fun=="exponential"){
      J <- lm(log(y_d) ~ x_d)
      yhat <- exp(J$coef[1] + J$coef[2] * x[index])
    }
    else if(fun=="log-log"){
      J <- lm(log(y_d) ~ log(x_d))
      yhat <- exp(J$coef[1] + J$coef[2] * log(x[index]))
    }
    else{
      print("Error!Invalid argument.")
    }
    e <- sum((y[index] - yhat)^2)
  }
  
  #'*Step 2:leave k out*
  n <- length(y)
  num <- c(1:n)
  
  #'*Obtain all combinations of 2 from n*
  valid.index <- combn(num,valid.num)
  
  #'*The following sapply function uses dataframe instead of matrix*
  valid.index <- as.data.frame(valid.index)
  
  fun_names <- c("linear","quadratic","exponential","log-log")
  
  for(name in fun_names){
    #'*For each column of dataframe apply validerror function resulting a vector*
    er <- sapply(valid.index,validerror,y=y,x=x,fun=name)
    #'*For each name iteration,assign a different name to er*
    assign(paste0(name,".e"),er)
  }
  
  list(
    `linear validation error`=mean(linear.e),
    `quadratic validation error`=mean(quadratic.e),
    `exponential validation error`=mean(exponential.e),
    `log-log validation error`=mean(`log-log.e`)
  )
}

y <- magnetic
x <- chemical

print(model.validation(y,x,2))

## -----------------------------------------------------------------------------
detach(package:DAAG,ironslag)
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
rm(list=ls())
invisible(gc())

library(boot)
set.seed(22004)

spcor.compare<- function(z){
  #'*Function for each permutation statistic*
  p.spcor <- function(z,ix){
    x <- z[,1]
    y <- z[ix,2]
    return(cor.test(x,y,method="spearman",exact = FALSE)$estimate)
  }
  #'*Use boot to sample permutation statistics*
  obj <- boot(data=z,statistic = p.spcor,R=10000,sim="permutation")
  ts <- c(obj$t0,obj$t)
  
  list(
    #'*p-value of spearman's test using cor.test*
    pvalue = cor.test(x, y, method = "spearman", exact = FALSE)$p.value,
    #'*p-value of spearman's test using permutation*
    p.pvalue = mean(abs(ts)>=abs(ts[1]))
  )
}

x <- iris[1:50,1]
y <- iris[1:50,3]
z <- cbind(x,y)

print(spcor.compare(z))

detach(package:boot)
rm(list=ls())
gc()


## -----------------------------------------------------------------------------
# clear memory and set seed
rm(list = ls())
set.seed(22034)

rl.metropolis <- function(sigma, x0, N) {
  # sigma: sd of proposal distribution N(xt,sigma^2)
  # x0: initial value
  # N: length of chain
  
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0  # to calculate acceptance rate
  for (t in 2:N) {
    y <- rnorm(1, x[t-1], sigma)
    if (u[t] <= exp(abs(x[t-1]) - abs(y))) { x[t] <- y; k <- k + 1 }
    else { x[t] <- x[t-1] }
  }
  return(list(mc = x, acc.prob = k / N))
}

N <- 10000
b <- 1000
k <- 4
sigma <- c(0.5, 1, 4, 16)
x0 <- c(-5, -2, 2, 5)
X <- matrix(nrow = k, ncol = N)
acc.prob <- numeric(k)
for (i in 1:k) {
  rl <- rl.metropolis(sigma[i], x0[i], N)
  X[i, ] <- rl$mc
  acc.prob[i] <- rl$acc.prob
}
acc.prob

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2, 2))
for (i in 1:k) {
  plot(X[i,], type = "l", xlab = bquote(sigma == .(sigma[i])),
       ylab = "X", ylim = range(X[i,]))
}

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2, 2))
x <- seq(-6, 6, 0.01)
fx <- exp(-abs(x)) / 2
for (i in 1:k) {
  hist(X[i, -(1:b)], breaks = "Scott", freq = FALSE, main = "",
       xlab = bquote(sigma == .(sigma[i])), xlim = c(-6, 6), ylim = c(0, 0.5),)
  lines(x, fx, col = 2, lty = 2)
}

## -----------------------------------------------------------------------------
z <- rexp(100, 1)
z <- c(-rev(z), z) # generate laplace random numbers
p <- c(0.05, seq(0.1, 0.9, 0.1), 0.95)
Q <- quantile(z, p)
mc <- X[, -(1:b)]
Qmc <- apply(mc, 1, function(x) quantile(x, p))
QQ <- data.frame(round(cbind(Q, Qmc), 3))
names(QQ) <- c('True', 'sigma=0.5', 'sigma=1', 'sigma=4', 'sigma=16')
knitr::kable(QQ)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(phi) {
  phi <- as.matrix(phi)
  k <- nrow(phi); n <- ncol(phi)
  phi.means <- rowMeans(phi)
  B <- n * var(phi.means)
  phi.w <- apply(phi, 1, var)
  W <- mean(phi.w)
  v.hat <- W * (n - 1) / n + B / n
  r.hat <- v.hat / W
  return(r.hat)
}

# ergodic mean plot
phi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(phi)) {
  phi[i,] <- phi[i,] / (1:ncol(phi))
}
for (i in 1:k) {
  if (i == 1) {
    plot((b+1):N, phi[i, (b+1):N], ylim = c(-0.5, 0.5),
         type = "l", xlab = 'Index', ylab = bquote(phi))
  } else { lines(phi[i, (b+1):N], col = i) }
}

# plot of R_hat
rhat <- rep(0, N)
for (j in (b+1):N) {
  rhat[j] <- Gelman.Rubin(phi[, 1:j])
}
plot(rhat[(b+1):N], type = "l", xlab = "", ylab = "R")
abline(h = 1.2, lty = 2)

## -----------------------------------------------------------------------------
# clear memory and set seed
#rm(list = ls())
set.seed(22004)

rbn.metropolis <- function(mu, sigma, rho, initial, N) {
  # mu, sigma, rho: parameter of bivariate normal distribution.
  # initial: initial value
  # N: length of chain
  
  X <- Y <- numeric(N)
  s <- sqrt(1 - rho^2) * sigma
  X[1] <- initial[1]; Y[1] <- initial[2]
  for (i in 2:N) {
    y <- Y[i-1]
    m1 <- mu[1] + rho * (y - mu[2]) * sigma[1] / sigma[2]
    X[i] <- rnorm(1, m1, s[1])
    x <- X[i]
    m2 <- mu[2] + rho * (x - mu[1]) * sigma[2] / sigma[1]
    Y[i] <- rnorm(1, m2, s[2])
  }
  return(list(X = X, Y = Y))
}

N <- 10000
b <- 1000
rho <- 0.9
mu <- c(0, 0)
sigma <- c(1, 1)
XY <- rbn.metropolis(mu, sigma, rho, mu, N)
X <- XY$X[-(1:b)]; Y <- XY$Y[-(1:b)]
plot(X, Y, xlab = bquote(X[t]), ylab = bquote(Y[t]),
     main = "", cex = 0.5, ylim = range(Y))
cov(cbind(X, Y))

## ----fig.height=4, fig.width=9------------------------------------------------
k <- 4
x0 <- matrix(c(2,2,-2,-2,4,-4,-4,4), nrow = 2, ncol = k)
Xmc <- Ymc <- XYmc <- matrix(0, nrow = k, ncol = N)
for (i in 1:k) {
  XY <- rbn.metropolis(mu, sigma, rho, x0[,i], N)
  Xmc[i,] <- XY$X; Ymc[i,] <- XY$Y
  XYmc[i,] <- Xmc[i,] * Ymc[i,]
}

# ergodic mean plot
cal_phi <- function(X) {
  phi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(phi)) {
    phi[i,] <- phi[i,] / (1:ncol(phi))
  }
  return(phi)
}
phiX <- cal_phi(Xmc)
phiY <- cal_phi(Ymc)
phiXY <- cal_phi(XYmc)

plot_erg_mean <- function(phi, rg) {
  for (i in 1:k) {
    if (i == 1) {
      plot((b+1):N, phi[i, (b+1):N], type = "l", ylim = rg,
           xlab = "Index", ylab = bquote(phi))
    }
    else { lines(phi[i, (b+1):N], col = i) }
  }
}
par(mfrow = c(1, 3))
plot_erg_mean(phiX, rg = c(-0.5, 0.5))
plot_erg_mean(phiY, rg = c(-0.5, 0.5))
plot_erg_mean(phiXY, rg = c(0.7, 1.1))

## ----fig.height=4, fig.width=9------------------------------------------------
Gelman.Rubin <- function(phi) {
  phi <- as.matrix(phi)
  k <- nrow(phi); n <- ncol(phi)
  phi.means <- rowMeans(phi)
  B <- n * var(phi.means)
  phi.w <- apply(phi, 1, var)
  W <- mean(phi.w)
  v.hat <- W * (n - 1) / n + B / n
  r.hat <- v.hat / W
  return(r.hat)
}

# plot of R_hat
plot_R_hat <- function(phi) {
  rhat <- rep(0, N)
  for (j in (b+1):N) {
    rhat[j] <- Gelman.Rubin(phi[, 1:j])
  }
  plot(rhat[(b+1):N], type = "l", xlab = "", ylab = "R", ylim = c(1, 1.25))
  abline(h = 1.2, lty = 2)
}
par(mfrow = c(1, 3))
plot_R_hat(phiX)
plot_R_hat(phiY)
plot_R_hat(phiXY)

## ----comment = ''-------------------------------------------------------------
lm.fit <- lm(Y ~ X)
summary(lm.fit)

## ----fig.height=5, fig.width=8------------------------------------------------
par(mfrow = c(1, 2))
e <- lm.fit$residuals
qx <- seq(-2, 2, 0.01)
hist(e, breaks = "Scott", freq = FALSE, main = "", xlim = c(-2, 2), ylim = c(0, 1))
lines(qx, dnorm(qx, 0, sqrt(0.19)), col = 2, lwd = 1.5)
qqnorm(e)
qqline(e, col = 2, lwd = 2, lty = 2)

## -----------------------------------------------------------------------------
set.seed(123)

# The function to generate the random sample
RSample <- function(n,alpha,beta){
  X <- runif(n,10,20)
  gamma <- 1;aM <- 0.5;aY <- 1
  M <- aM+alpha*X+rnorm(n)
  Y <- aY+beta*M+gamma*X+rnorm(n)
  return(list(X,M,Y))
}

# The function of test statistics computation
Ttest <- function(X,M,Y){
  fit1 <- summary(lm(M~X))
  fit2 <- summary(lm(Y~X+M))
  a <- fit1$coefficients[2,1]
  sea <- fit1$coefficients[2,2]
  b <- fit2$coefficients[3,1]
  seb <- fit2$coefficients[3,2]
  return(a*b/((a*seb)^2+(b*sea)^2)^0.5)
}

# The function to implement the test hypothesis
Imptest <- function(N,n,X,M,Y,T0){
  T1 <- T2 <- T3 <- numeric(N)
  # Condition 1
  for(i in 1:N){
    n1 <- sample(1:n, size=n, replace=FALSE)
    n2 <- sample(1:n, size=n, replace=FALSE)
    X1 <- X[n1];M1 <- M[n2];Y1 <- Y[n2]
    T1[i] <- Ttest(X1,M1,Y1)
  }
  # Condition 2
  for(i in 1:N){
    n1 <- sample(1:n, size = n, replace = FALSE)
    n2 <- sample(1:n, size = n, replace = FALSE)
    X2 <- X[n1];M2 <- M[n1];Y2 <- Y[n2]
    T2[i] <- Ttest(X2,M2,Y2)
  }
  # Condition 3
  for(i in 1:N){
    n1 <- sample(1:n, size = n, replace = FALSE)
    n2 <- sample(1:n, size = n, replace = FALSE)
    M3 <- M[n1];X3 <- X[n2];Y3 <- Y[n2]
    T3[i] <- Ttest(X3,M3,Y3)
  }
  # The p-value of Condition1
  p1 <- mean(abs(c(T0,T1))>abs(T0))
  # The p-value of Condition2
  p2 <- mean(abs(c(T0,T2))>abs(T0))
  # The p-value of Condition3
  p3 <- mean(abs(c(T0,T3))>abs(T0))
  return(c(p1,p2,p3))
}

N <- 1000 # The number of simulation
n <- 100 # The number of random sample
T0 <- numeric(3)
p <- matrix(0,3,3)
# The real values of parameters
alpha <- c(0,0,1);beta <- c(0,1,0)

for(i in 1:3){
  result <- RSample(n,alpha[i],beta[i])
  X <- result[[1]]
  M <- result[[2]]
  Y <- result[[3]]
  # The original value of test statistics
  T0[i] <- Ttest(X,M,Y)
  p[i,] <- Imptest(N,n,X,M,Y,T0[i])
}

## ----setup12, fig.height=10, fig.width=10, echo=T, eval=T---------------------
# Result reporting
colnames(p) <- c("Condition 1","Condition 2","Condition 3")
rownames(p) <- c("alpha=0,beta=0","alpha=0,beta=1","alpha=1,beta=0")
p

# Clean the memory of the variables
rm(list=ls())

## -----------------------------------------------------------------------------
set.seed(22004)

# The function to solve the alpha
solve <- function(N,b1,b2,b3,f0){
  x1 <- rpois(N,1)
  x2 <- rexp(N,1)
  x3 <- rbinom(N,1,0.5)
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
  }
  solution <- uniroot(g,c(-50,0))
  return(round(unlist(solution),5)[1])
}

N <- 1e6;b1 <- 0;b2 <- 1;b3 <- -1
f0 <- c(0.1,0.01,0.001,0.0001)
alpha <- numeric(length(f0))
for(i in 1:length(f0)){
  alpha[i] <- solve(N,b1,b2,b3,f0[i])
}
result <- rbind(f0,alpha)
rownames(result) <- c("f0","alpha")
result

par(mfrow=c(1,2))
# Draw the scatter plot of f0 and alpha
plot(f0,alpha,main="The scatter plot of f0 and alpha")
# Draw the scatter plot of log(f0) and alpha
plot(log(f0),alpha,main="The scatter plot of log(f0) and alpha")
par(mfrow=c(1,1))

# Clean the memory of the variables
rm(list=ls())

## -----------------------------------------------------------------------------
rm(list=ls())
invisible(gc())


## -----------------------------------------------------------------------------
#'*lower bound*
u <- c(11,8,27,13,16,0,23,10,24,2)
#'*upper bound*
v <- c(12,9,28,14,17,1,24,11,25,3)


#'*Observed data likelihood*
o.likelihood <- function(lambda){
  sum((v*exp(-lambda*v)-u*exp(-lambda*u))/(exp(-lambda*u)-exp(-lambda*v)))
}



solution <- uniroot(o.likelihood,interval = c(0.05,0.1),extendInt = "yes")

k <- round(unlist(solution),5)[1:3]


MLE <- k[1]

#'*EM algorithm*

lambda.old <- 0.0000000001
N <- 1e5

tol <- .Machine$double.eps

options(digits=10)
for(j in 1:N) {
  
  lambda <- length(u)/(sum((u*exp(-lambda.old*u)-v*exp(-lambda.old*v))/(exp(-lambda.old*u)-exp(-lambda.old*v)))+length(u)/lambda.old)

  if ((abs(lambda - lambda.old)/lambda.old) < tol) break
  lambda.old <- lambda
}




## -----------------------------------------------------------------------------
rm(list=ls())
invisible(gc())


## -----------------------------------------------------------------------------
v <- c(1,2,3)

dim(v)


## -----------------------------------------------------------------------------

rm(list=ls())
invisible(gc())


## -----------------------------------------------------------------------------

m <- matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T)

is.matrix(m)

is.array(m)


## -----------------------------------------------------------------------------

rm(list=ls())
invisible(gc())


## -----------------------------------------------------------------------------
#'*atrributes of dataframe*
d <- data.frame(a=c(1,2,3),b=c(4,5,6))

attributes(d)

#'*atrributes of matrix*
x <- cbind(a = 1:3, pi = pi)

attributes(x)


## -----------------------------------------------------------------------------

rm(list=ls())
invisible(gc())


## -----------------------------------------------------------------------------

df <- data.frame(x = 1:3, y = I(list(1:2, c("a","b","c"), c(F,T,F,T))))

df

as.matrix(df)


## -----------------------------------------------------------------------------

rm(list=ls())
invisible(gc())


## -----------------------------------------------------------------------------

#'*dataframe with 0 rows and 2 cols*

d1 <- data.frame(a=matrix(nrow = 0,ncol=1),b=matrix(nrow = 0,ncol=1))

dim(d1)

#'*dataframe with 2 rows and 0 cols*

d2 <- data.frame(row.names=c("a","b"))

dim(d2)


## -----------------------------------------------------------------------------

rm(list=ls())
invisible(gc())


## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
x <- data.frame(x1 = c(1.5, 2.5, 3.5, 4.5), x2 = rnorm(4, 4, 4))
str(x)

## -----------------------------------------------------------------------------
res1 <- data.frame(lapply(x, scale01))
res1

## -----------------------------------------------------------------------------
# add a non-numeric column
x$x3 = c(rep("A", 2), rep("B", 2))

res2 <- data.frame(lapply(x, function(x) if (is.numeric(x)) scale01(x) else x))
res2

## -----------------------------------------------------------------------------
rm(list = ls())
x <- data.frame(x1 = c(0.6, 1.3, 7.6, 2.4), x2 = rnorm(4, 2, 2))
str(x)

## -----------------------------------------------------------------------------
res1 <- vapply(x, sd, 1)
res1

## -----------------------------------------------------------------------------
# add a non-numeric column
x$x3 = c(rep("A", 2), rep("B", 2))

res2 <- vapply(x[vapply(x, is.numeric, TRUE)], sd, 1)
res2

## -----------------------------------------------------------------------------
# clear memory and set seed
rm(list = ls())
set.seed(22004)

gibbsR <- function(mu, sigma, rho, initial, N) {
  # mu, sigma, rho: parameter of bivariate normal distribution.
  # initial: initial value
  # N: length of chain
  
  X <- Y <- numeric(N)
  s <- sqrt(1 - rho^2) * sigma
  X[1] <- initial[1]; Y[1] <- initial[2]
  for (i in 2:N) {
    y <- Y[i-1]
    m1 <- mu[1] + rho * (y - mu[2]) * sigma[1] / sigma[2]
    X[i] <- rnorm(1, m1, s[1])
    x <- X[i]
    m2 <- mu[2] + rho * (x - mu[1]) * sigma[2] / sigma[1]
    Y[i] <- rnorm(1, m2, s[2])
  }
  return(list(X = X, Y = Y))
}

## ----fig.height=4, fig.width=8------------------------------------------------
library(Rcpp)
# generate chains
N <- 10000
b <- 1000
rho <- 0.9
mu <- c(0, 0)
sigma <- c(1, 1)
XYR <- gibbsR(mu, sigma, rho, mu, N)
XR <- XYR$X[-(1:b)]; YR <- XYR$Y[-(1:b)]
sourceCpp('C:/Users/10313/Desktop/StatisticalComputing/HW/HW11/gibbsC.cpp')
XYC <- gibbsC(mu, sigma, rho, mu, N)
XC <- XYC[-(1:b), 1]; YC <- XYC[-(1:b), 2]

par(mfrow = c(1, 2))
qqplot(XR, XC, plot.it = TRUE)
abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
qqplot(YR, YC, plot.it = TRUE)
abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)

## ----fig.height=4, fig.width=8------------------------------------------------
par(mfrow = c(1, 2))
plot(XR, YR, cex = 0.5)
plot(XC, YC, cex = 0.5)

## -----------------------------------------------------------------------------
# import package
ts <- microbenchmark(gibbsR = gibbsR(mu, sigma, rho, mu, N),
                     gibbsC = gibbsC(mu, sigma, rho, mu, N))
summary(ts)[, c(1, 3, 5, 6)]


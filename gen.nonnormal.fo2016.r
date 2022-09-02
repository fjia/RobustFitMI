#### Foldnes, N., & Olsson, U. H. (2016). 
#### A Simple Simulation Technique for Nonnormal Data with Prespecified Skewness, Kurtosis, and Covariance Matrix. 
#### Multivariate Behavioral Research, 51(2–3), 207–219.

# install.packages("nleqslv")
# install.packages("PearsonDS")

library(lavaan)
library(nleqslv)
library(PearsonDS)
library(psych)


gen.nonnormal.f02016 <- function(npop, obs.sigma, sk, ku, k){
#### npop = pop size 
#### obs.sigma = targe covariance matrix
#### sk = vector of target skewness
#### ku = vector of target kurtosis
#### k = number of variables, must equal length(sk) and length(ku)

## STEP 1
#target matrix:
sigma <- obs.sigma

#target skewness and kurtosis
observed.skew <- sk
observed.excesskurt <- ku

### STEP 2
A <- t(chol(sigma))


### STEP 3
function.skew <- function(IGvalues, k){
fval <- numeric(k)
for (i in 1:k)
fval[i] <- A[i, ]^3 %*% IGvalues/(sum(A[i,]^2)^(3/2))
fval-observed.skew
}
IGskew <- nleqslv(x=observed.skew, function.skew, k = k)$x

function.kurt <- function(IGvalues, k){
fval <- numeric(k)
for (i in 1:k)
fval[i] <- A[i, ]^4 %*% IGvalues/(sum(A[i,]^2)^(2))
fval-observed.excesskurt
}
IGkurt.excess <- nleqslv(x=observed.excesskurt, function.kurt, k = k)$x

### STEP 4
parlist <-  list()
for (i in 1:k)
parlist[[i]] <- pearsonFitM(moments=c(mean=0,variance=1,skewness=IGskew[i],3+IGkurt.excess[i]))

### STEP 5
N <- npop
IGdata <- matrix(ncol=k, nrow=N)
for (i in 1:k)
IGdata[,i] <- rpearson(N, parlist[[i]])
simulated.sample <- IGdata %*% t(A)
#colnames(simulated.sample) <- colnames(sigma)

return(simulated.sample)

# testing
#round(cov(simulated.sample)-sigma,2) #target covariance
#round(skew(simulated.sample)- observed.skew,2)# target skew
#round(kurtosi(simulated.sample)-observed.excesskurt,2)#target kurtosis
}

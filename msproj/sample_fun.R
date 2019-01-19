## load library
library("knitr")
library("mvtnorm")
library("SpatialExtremes")
library("rmutil")



#### Function for generating the random sample

#  The arguments of the function sample.fun are:
#   n  = sample size
#   dist.type = type of distribution to by used for random sample 

sample_fun <- function(n, dist.type){
  # random sample 
  if(dist.type == "Normal"){
    grp1 <- rnorm(n, mean = 0, sd = 1)
    grp2 <- rnorm(n, mean = 1, sd = 1)
    grp3 <- rnorm(n, mean = 2, sd = 1)
    dat  <- cbind(grp1, grp2, grp3)
  } else if (dist.type == "Exponential"){
    grp1 <- rexp(n, rate = 1)
    grp2 <- rexp(n, rate = 2)
    grp3 <- rexp(n, rate = 3)
    dat  <- cbind(grp1, grp2, grp3)
  } else if (dist.type == "Cauchy"){
    grp1 <- rcauchy(n, location = 0, scale = 1)
    grp2 <- rcauchy(n, location = 1, scale = 1)
    grp3 <- rcauchy(n, location = 2, scale = 1)
    dat  <- cbind(grp1, grp2, grp3)
  } else if (dist.type == "Laplace"){
    grp1 <- rlaplace(n, m = 0, s = 1)
    grp2 <- rlaplace(n, m = 1, s = 1)
    grp3 <- rlaplace(n, m = 2, s = 1)
    dat  <- cbind(grp1, grp2, grp3)
  } else if (dist.type == "GEV"){
    grp1 <- rgev(n, loc = 0, scale = 1, shape = 0)
    grp2 <- rgev(n, loc = 1, scale = 1, shape = 0)
    grp3 <- rgev(n, loc = 2, scale = 1, shape = 0)
    dat  <- cbind(grp1, grp2, grp3)
  } else if(dist.type == "NormalMixture"){
    # group 1
    dists <- runif(n)          
    grp1  <- vector(length = n)
    for(i in 1: n){if(dists[i] < 0.5){
      grp1[i] <-rnorm(1, mean=0, sd=1)
    } else {grp1[i] <- rnorm(1, mean=1, sd=1)}}
    # group 2
    dists <- runif(n)          
    grp2  <- vector(length = n)
    for(i in 1: n){if(dists[i] < 0.5){
      grp2[i] <- rnorm(1, mean=1, sd=1)
    } else {grp2[i] <- rnorm(1, mean=2, sd=1)}}
    # group 3
    dists <- runif(n)          
    grp3  <- vector(length = n)
    for(i in 1: n){if(dists[i] < 0.5){
      grp3[i] <- rnorm(1, mean=2, sd=1)
    } else {grp3[i] <- rnorm(1, mean=3, sd=1)}}
    dat  <- cbind(grp1, grp2, grp3)
  } else if(dist.type == "Lognormal"){
    grp1 <- rlnorm(n, meanlog = 0, sdlog = 1)
    grp2 <- rlnorm(n, meanlog = 1, sdlog = 1)
    grp3 <- rlnorm(n, meanlog = 2, sdlog = 1)
    dat  <- cbind(grp1,grp2, grp3)
  }
  return(dat)
}

## clearing the memory
rm(list = ls())

## load library
library("knitr")
library("mvtnorm")
library("SpatialExtremes")
library("rmutil")

###
### Universal Parameters
###


###  Number of Replications
R <- 10000

###  Sample sizes
N <- c(10, 30, 100, 500)

### Confidence level
ci_level <- 0.95

### quantile
p <- 0.5




####################################################
####                                            #### 
####     Ratio and difference of Quantiles      ####
####                                            ####
####################################################



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


## kernel density function
kd_fun <- function(i, dat, md){
  # treatment
  den.grp <- density(dat[,i])
  datfram.grp <- data.frame(x = den.grp$x, y = den.grp$y)
  ix = which.min(abs(den.grp$x - md[i]))
  datfram.grp[ix,"y"]
}

Coverage_Prob = function(n,R,p,ci_level,comp.type, dist.type, den.type, equi.q){
  
  if(comp.type == "ap"){
        ## get the design matrix
        design.matrix <- matrix(c(1,-1, 0,1,0,-1,0,1,-1), nrow = 3, byrow = T)
         
        ## True DELTA -- for difference
        if(dist.type == "Normal"){ 
          true.DELTA = design.matrix%*%c(0,1,2)
        } else if (dist.type == "Exponential"){
          true.DELTA = design.matrix%*%c(log(2), 0.5*log(2), (1/3)*log(2))
        } else if (dist.type == "Cauchy"){
          true.DELTA = design.matrix%*%c(0,1,2)
        } else if (dist.type == "Laplace"){
          true.DELTA = design.matrix%*%c(0,1,2)
        } else if (dist.type == "GEV"){
          true.DELTA = design.matrix%*%c(-log(log(2)), 1- log(log(2)), 2 -log(log(2)))
        } else if (dist.type == "NormalMixture"){
          true.DELTA = design.matrix%*%c((0.5),(3/2),(5/2))
        } else if (dist.type == "Lognormal"){
          true.DELTA = design.matrix%*%c(exp(0),exp(1),exp(2))
        }
        
        ## True density at the median
        if (den.type == "tde"){
          # True density at the median is based on the distn in question
          if(dist.type == "Normal"){
            est.f.md <- c(dnorm(0,mean = 0, sd = 1),
                          dnorm(1,mean = 1, sd = 1),
                          dnorm(2,mean = 2, sd = 1))
          } else if (dist.type == "Exponential"){
            est.f.md <- c(dexp(      log(2), rate = 1),
                          dexp(  0.5*log(2), rate = 2),
                          dexp((1/3)*log(2), rate = 3))
          } else if (dist.type == "Cauchy"){
            est.f.md <- c(dcauchy(0,location = 0, scale = 1),
                          dcauchy(1,location = 1, scale = 1),
                          dcauchy(2,location = 2, scale = 1))
          } else if (dist.type == "Laplace"){
            est.f.md <- c(dlaplace(0, m = 0, s = 1),
                          dlaplace(1, m = 1, s = 1),
                          dlaplace(2, m = 2, s = 1))
          } else if (dist.type == "GEV"){
            est.f.md <- c(dgev(  (-log(log(2))), loc = 0, scale = 1, shape = 0),
                          dgev((1- log(log(2))), loc = 1, scale = 1, shape = 0),
                          dgev((2 -log(log(2))), loc = 2, scale = 1, shape = 0))
          } else if (dist.type == "NormalMixture"){
            est.f.md <- c(0.5*dnorm(0.5,mean = 0, sd = 1) + 0.5*dnorm(0.5,mean = 1, sd = 1),
                          0.5*dnorm(1.5,mean = 1, sd = 1) + 0.5*dnorm(1.5,mean = 2, sd = 1),
                          0.5*dnorm(2.5,mean = 2, sd = 1) + 0.5*dnorm(2.5,mean = 3, sd = 1))
          }  else if (dist.type == "Lognormal"){
            est.f.md <- c(dlnorm(exp(0), meanlog = 0, sdlog = 1),
                          dlnorm(exp(1), meanlog = 1, sdlog = 1),
                          dlnorm(exp(2), meanlog = 2, sdlog = 1))
          }
        }
        
        CI_replication <- replicate(R,  {
            # get sample for the three groups
            dat = sample_fun(n, dist.type)
            # get the median of each group
            md = apply(dat, 2, median)
            ## density estimates at the median
            if(den.type == "kde"){est.f.md <- sapply(1:3, kd_fun, dat, md)}
            # vector of variance 1 for each group
            sig2     <- (p*(1 - p))/(n*(est.f.md^2))
            # get the SIGMA
            SIGMA    <- diag(sig2)
            # covariance matrix
            cov.mat  <- design.matrix%*%SIGMA%*%t(design.matrix)
            # get correlation matrix
            D.inv    <- solve(sqrt(diag(diag(cov.mat))))
            corr.mat <- D.inv%*%cov.mat%*%D.inv
        
            # equicoordinate quantile
            if (equi.q == "MVNq"){
              q <- qmvnorm(ci_level, corr = corr.mat, tail = "both")$quantile
            } else if (equi.q == "MVTq"){
              q <- qmvt(ci_level, corr = corr.mat, df = 3*(n -1), tail = "both")$quantile
            }
            
            # get delta
            DELTA = design.matrix%*%md
            # get std.err and quantile
            C <- q*sqrt(diag(cov.mat))
            # confidence interval
            CI = cbind(DELTA - C, DELTA + C)
            
            ## get the average length for difference
            Length_diff    <- CI[,2] - CI[,1]
            L_tmp_diff     <- mean((Length_diff))
          
          
            ci_diff     <- ((CI[1,1] < true.DELTA[1,1]) & (true.DELTA[1,1] < CI[1,2] ) &
                            (CI[2,1] < true.DELTA[2,1]) & (true.DELTA[2,1] < CI[2,2] ) &
                            (CI[3,1] < true.DELTA[3,1]) & (true.DELTA[3,1] < CI[3,2] ) )
          
            c(ci_diff, round(L_tmp_diff,1))
        }) # end of the replication
        
        res <- rowMeans( CI_replication )
        CP_d = res[1]
        L_d  = round(res[2],1)
        
        ap.out <- c(CP_d, L_d)
        ap.out
  } 
  
  if(comp.type == "mc"){
    
      # container for CI of ratio
      ci <- matrix(NA, 2,2)
      ## get the design matrix
      design.matrix <- matrix(c(1,0,-1,0,1,-1), nrow = 2, byrow = T)
      
      ## True DELTA -- for difference
      if(dist.type == "Normal"){ 
        true.DELTA = design.matrix%*%c(0,1,2)
      } else if (dist.type == "Exponential"){
        true.DELTA = design.matrix%*%c(log(2), 0.5*log(2), (1/3)*log(2))
      } else if (dist.type == "Cauchy"){
        true.DELTA = design.matrix%*%c(0,1,2)
      } else if (dist.type == "Laplace"){
        true.DELTA = design.matrix%*%c(0,1,2)
      } else if (dist.type == "GEV"){
        true.DELTA = design.matrix%*%c(-log(log(2)), 1- log(log(2)), 2 -log(log(2)))
      } else if (dist.type == "NormalMixture"){
        true.DELTA = design.matrix%*%c((0.5),(3/2),(5/2))
      } else if (dist.type == "Lognormal"){
        true.DELTA = design.matrix%*%c(exp(0),exp(1),exp(2))
      }
      
      
        ## True rho -- for ratio
        if(dist.type == "Normal"){
          true.RHO = c(0/2, 1/2)
        } else if (dist.type == "Exponential"){
          true.RHO = c(    (log(2))/((1/3)*log(2)), 
                           (0.5*log(2))/((1/3)*log(2)))
        } else if (dist.type == "Cauchy"){
          true.RHO = c(0/2, 1/2)
        } else if (dist.type == "Laplace"){
          true.RHO = c(0/2, 1/2)
        } else if (dist.type == "GEV"){
          true.RHO = c( (-log(log(2)))/(2-log(log(2))),
                        (1-log(log(2)))/(2-log(log(2))))
        } else if (dist.type == "NormalMixture"){
          true.RHO = c((0.5)/(5/2),(3/2)/(5/2))
        } else if (dist.type == "Lognormal"){
          true.RHO = c(exp(0)/exp(2), exp(1)/exp(2))
        }
       
      
      ## True density at the median
      if (den.type == "tde"){
        # True density at the median is based on the distn in question
        if(dist.type == "Normal"){
          est.f.md <- c(dnorm(0,mean = 0, sd = 1),
                        dnorm(1,mean = 1, sd = 1),
                        dnorm(2,mean = 2, sd = 1))
        } else if (dist.type == "Exponential"){
          est.f.md <- c(dexp(      log(2), rate = 1),
                        dexp(  0.5*log(2), rate = 2),
                        dexp((1/3)*log(2), rate = 3))
        } else if (dist.type == "Cauchy"){
          est.f.md <- c(dcauchy(0,location = 0, scale = 1),
                        dcauchy(1,location = 1, scale = 1),
                        dcauchy(2,location = 2, scale = 1))
        } else if (dist.type == "Laplace"){
          est.f.md <- c(dlaplace(0, m = 0, s = 1),
                        dlaplace(1, m = 1, s = 1),
                        dlaplace(2, m = 2, s = 1))
        } else if (dist.type == "GEV"){
          est.f.md <- c(dgev(  (-log(log(2))), loc = 0, scale = 1, shape = 0),
                        dgev((1- log(log(2))), loc = 1, scale = 1, shape = 0),
                        dgev((2 -log(log(2))), loc = 2, scale = 1, shape = 0))
        } else if (dist.type == "NormalMixture"){
          est.f.md <- c(0.5*dnorm(0.5,mean = 0, sd = 1) + 0.5*dnorm(0.5,mean = 1, sd = 1),
                        0.5*dnorm(1.5,mean = 1, sd = 1) + 0.5*dnorm(1.5,mean = 2, sd = 1),
                        0.5*dnorm(2.5,mean = 2, sd = 1) + 0.5*dnorm(2.5,mean = 3, sd = 1))
        }  else if (dist.type == "Lognormal"){
          est.f.md <- c(dlnorm(exp(0), meanlog = 0, sdlog = 1),
                        dlnorm(exp(1), meanlog = 1, sdlog = 1),
                        dlnorm(exp(2), meanlog = 2, sdlog = 1))
        }
      }
      
      CI_replication <- replicate(R,  {
        
        g <- 2
        
        while (g >= 1 ) { ## to avoid sqrt of negative and zero
          
          # get sample for the three groups
          dat = sample_fun(n, dist.type)
          # get the median of each group
          md = apply(dat, 2, median)
          
          ## density estimates at the median
          if(den.type == "kde"){est.f.md <- sapply(1:3, kd_fun, dat, md)}
          
          # vector of variance 1 for each group
          sig2     <- (p*(1 - p))/(n*(est.f.md^2))
          # get the SIGMA
          SIGMA    <- diag(sig2)
          
          ##
          ## For ratio ci
          ##
          
          # estimated rho 
          rho_13 <- md[1]/md[3]
          rho_23 <- md[2]/md[3]
          rho    <- c(rho_13, rho_23)
          
          # design matrix for ratio using the plug-in approach
          # Dilba et al, 2006
          design.matrixx <- matrix(c(1,0,-rho_13,0,1,-rho_23), nrow = 2, byrow = T)
          # covariance matrix
          cov.matt  <- design.matrixx%*%SIGMA%*%t(design.matrixx)
          # get correlation matrix 
          D.invv    <- solve(sqrt(diag(diag(cov.matt))))
          corr.matt <- D.invv%*%cov.matt%*%D.invv
          
          # equicoordinate quantile for MVN and MVT
          if (equi.q == "MVNq"){
            Q <- qmvnorm(ci_level, corr = corr.matt, tail = "both")$quantile
          } else if (equi.q == "MVTq"){
            Q <- qmvt(ci_level, corr = corr.matt, df = 3*(n -1), tail = "both")$quantile
          }
          
          
          # get the g
          g <- ((Q*sig2[3])/md[3])^2
          
          
        } # end of the while loop
        
        
          ## Confidence Interval for ratio
        ci[1,] <- (1/(1-g))*(rho[1]*c(1,1) + c(-1,1)*((Q/md[3])*sqrt(sig2[1]*(1-g) + (rho[1]^2)*sig2[3])))
        ci[2,] <- (1/(1-g))*(rho[2]*c(1,1) + c(-1,1)*((Q/md[3])*sqrt(sig2[2]*(1-g) + (rho[2]^2)*sig2[3])))
        
        ## get the average length for ratio
        Length_ratio    <- ci[,2] - ci[,1]
        L_tmp_ratio     <- mean((Length_ratio))
        
        ci_ratio    <- ((ci[1,1] < true.RHO[1]) & (true.RHO[1] < ci[1,2]) &
                        (ci[2,1] < true.RHO[2]) & (true.RHO[2] < ci[2,2]) )
        ##
        ## For difference
        ##
        
        # covariance matrix
        cov.mat  <- design.matrix%*%SIGMA%*%t(design.matrix)
        # get correlation matrix 
        D.inv    <- solve(sqrt(diag(diag(cov.mat))))
        corr.mat <- D.inv%*%cov.mat%*%D.inv
        
        # equicoordinate quantile
        if (equi.q == "MVNq"){
          q <- qmvnorm(ci_level, corr = corr.mat, tail = "both")$quantile
        } else if (equi.q == "MVTq"){
          q <- qmvt(ci_level, corr = corr.mat, df = 3*(n -1), tail = "both")$quantile
        }
        
        # get delta
        DELTA = design.matrix%*%md
        # get std.err and quantile
        C <- q*sqrt(diag(cov.mat))
        # confidence interval
        CI = cbind(DELTA - C, DELTA + C)
        
        ## get the average length for difference
        Length_diff    <- CI[,2] - CI[,1]
        L_tmp_diff     <- mean((Length_diff))
        
        ci_diff = ((CI[1,1] < true.DELTA[1,1]) & ( true.DELTA[1,1] < CI[1,2] ) &
                   (CI[2,1] < true.DELTA[2,1]) & ( true.DELTA[2,1] < CI[2,2] ) )
        
        c(ci_diff, round(L_tmp_diff,1), ci_ratio, round(L_tmp_ratio,1))
       
      }) # end of the replication
      
      res <- rowMeans( CI_replication )
      
      CP_d = res[1]
      L_d  = round(res[2],1)
      CP_r = res[3]
      L_r  = round(res[4],1)
      
      mc.out <- c(CP_d,L_d,CP_r,L_r)
      mc.out
  }
  
  if(comp.type == "ap"){
    return(ap.out)
  } else {
    return(mc.out)
  }
  
}


system.time({ res <- Coverage_Prob(n             = 10      , 
                                   R             = 10    ,
                                   p             = 0.5     ,
                                   ci_level      = 0.95    ,
                                   comp.type     = "mc"    , 
                                   dist.type     = "NormalMixture", 
                                   den.type      = "kde"   , 
                                   equi.q        = "MVNq"  )
}); res





##############################################################################################
###                                                                                        ###
###      "Multiple comparison to control - kde and tde - difference of quantiles - MVNq"   ###
###                                                                                        ###
##############################################################################################

startTime1 <- Sys.time()
system.time({
# kernel density estimate at the median
CP.n_mc_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Normal",        "kde", "MVNq")
CP.e_mc_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Exponential",   "kde", "MVNq")
CP.c_mc_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Cauchy",        "kde", "MVNq")
CP.l_mc_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Laplace",       "kde", "MVNq")
CP.g_mc_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "GEV",           "kde", "MVNq")
CP.nm_mc_kde <- sapply(N,Coverage_Prob,R,p, ci_level, "mc", "NormalMixture", "kde", "MVNq")
# True Density estimate at the median
CP.n_mc_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Normal",        "tde", "MVNq")
CP.e_mc_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Exponential",   "tde", "MVNq")
CP.c_mc_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Cauchy",        "tde", "MVNq")
CP.l_mc_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Laplace",       "tde", "MVNq")
CP.g_mc_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "GEV",           "tde", "MVNq")
CP.nm_mc_tde <- sapply(N,Coverage_Prob,R,p, ci_level, "mc", "NormalMixture", "tde", "MVNq")
})
stopTime1 <- Sys.time()
stopTime1 - startTime1
# Time difference of 6.339286 hours

Cov.Prob_mc_diff = data.frame( Density= c(rep("kde",length(N)), rep("tde", length(N))),
                               n      = c(N, N),
                               Norm   = c(CP.n_mc_kde[1,], CP.n_mc_tde[1,]),
                               L_Norm = c(CP.n_mc_kde[2,], CP.n_mc_tde[2,]),
                               Exp    = c(CP.e_mc_kde[1,], CP.e_mc_tde[1,]),
                               L_Exp  = c(CP.e_mc_kde[2,], CP.e_mc_tde[2,]),
                               Cau    = c(CP.c_mc_kde[1,], CP.c_mc_tde[1,]),
                               L_Cau  = c(CP.c_mc_kde[2,], CP.c_mc_tde[2,]),
                               Lap    = c(CP.l_mc_kde[1,], CP.l_mc_tde[1,]),
                               L_Lap  = c(CP.l_mc_kde[2,], CP.l_mc_tde[2,]),
                               GEV    = c(CP.g_mc_kde[1,], CP.g_mc_tde[1,]),
                               L_GEV  = c(CP.g_mc_kde[2,], CP.g_mc_tde[2,]),
                               NMix   = c(CP.nm_mc_kde[1,], CP.nm_mc_tde[1,]),
                               L_NMix = c(CP.nm_mc_kde[2,], CP.nm_mc_tde[2,]))

kable(Cov.Prob_mc_diff,
      caption = "Multiple comparison to control - kde and tde - difference of quantiles - MVNq",
      format = "pandoc")
write.csv(Cov.Prob_mc_diff, "Cov.Prob_mc_diff.csv")

##############################################################################################
###                                                                                        ###
###        "Multiple comparison to control - kde and tde - ratio of quantiles - MVNq"      ###
###                                                                                        ###
##############################################################################################



Cov.Prob_mc_ratio = data.frame( Density= c(rep("kde",length(N)), rep("tde", length(N))),
                                n      = c(N, N),
                                Norm   = c(CP.n_mc_kde[3,], CP.n_mc_tde[3,]),
                                L_Norm = c(CP.n_mc_kde[4,], CP.n_mc_tde[4,]),
                                Exp    = c(CP.e_mc_kde[3,], CP.e_mc_tde[3,]),
                                L_Exp  = c(CP.e_mc_kde[4,], CP.e_mc_tde[4,]),
                                Cau    = c(CP.c_mc_kde[3,], CP.c_mc_tde[3,]),
                                L_Cau  = c(CP.c_mc_kde[4,], CP.c_mc_tde[4,]),
                                Lap    = c(CP.l_mc_kde[3,], CP.l_mc_tde[3,]),
                                L_Lap  = c(CP.l_mc_kde[4,], CP.l_mc_tde[4,]),
                                GEV    = c(CP.g_mc_kde[3,], CP.g_mc_tde[3,]),
                                L_GEV  = c(CP.g_mc_kde[4,], CP.g_mc_tde[4,]),
                                NMix   = c(CP.nm_mc_kde[3,], CP.nm_mc_tde[3,]),
                                L_NMix = c(CP.nm_mc_kde[4,], CP.nm_mc_tde[4,]))

kable(Cov.Prob_mc_ratio,
      caption = "Multiple comparison to control - kde and tde - ratios of quantiles - MVNq",
      format = "pandoc")
write.csv(Cov.Prob_mc_ratio, "Cov.Prob_mc_ratio.csv")




##############################################################################################
###                                                                                        ###
###      "Multiple comparison to control - kde and tde - difference of quantiles - MVTq"   ###
###                                                                                        ###
##############################################################################################


startTime2 <- Sys.time()
system.time({
# kernel density estimate at the median
CP.n_mc_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Normal",        "kde", "MVTq")
CP.e_mc_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Exponential",   "kde", "MVTq")
CP.c_mc_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Cauchy",        "kde", "MVTq")
CP.l_mc_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Laplace",       "kde", "MVTq")
CP.g_mc_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "GEV",           "kde", "MVTq")
CP.nm_mc_kde_T <- sapply(N,Coverage_Prob,R,p, ci_level, "mc", "NormalMixture", "kde", "MVTq")
# True Density estimate at the median
CP.n_mc_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Normal",        "tde", "MVTq")
CP.e_mc_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Exponential",   "tde", "MVTq")
CP.c_mc_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Cauchy",        "tde", "MVTq")
CP.l_mc_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "Laplace",       "tde", "MVTq")
CP.g_mc_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "mc", "GEV",           "tde", "MVTq")
CP.nm_mc_tde_T <- sapply(N,Coverage_Prob,R,p, ci_level, "mc", "NormalMixture", "tde", "MVTq")
})
stopTime2 <- Sys.time()
stopTime2 - startTime2
# Time difference of 7.307918 hours



# user   system  elapsed 
# 17160.83     5.92 17270.80
# > 17271/3600
# [1] 4.7975 -- this is approximately 5 hours

Cov.Prob_mc_diff_T = data.frame( Density= c(rep("kde",length(N)), rep("tde", length(N))),
                                 n      = c(N, N),
                                 Norm   = c(CP.n_mc_kde_T[1,], CP.n_mc_tde_T[1,]),
                                 L_Norm = c(CP.n_mc_kde_T[2,], CP.n_mc_tde_T[2,]),
                                 Exp    = c(CP.e_mc_kde_T[1,], CP.e_mc_tde_T[1,]),
                                 L_Exp  = c(CP.e_mc_kde_T[2,], CP.e_mc_tde_T[2,]),
                                 Cau    = c(CP.c_mc_kde_T[1,], CP.c_mc_tde_T[1,]),
                                 L_Cau  = c(CP.c_mc_kde_T[2,], CP.c_mc_tde_T[2,]),
                                 Lap    = c(CP.l_mc_kde_T[1,], CP.l_mc_tde_T[1,]),
                                 L_Lap  = c(CP.l_mc_kde_T[2,], CP.l_mc_tde_T[2,]),
                                 GEV    = c(CP.g_mc_kde_T[1,], CP.g_mc_tde_T[1,]),
                                 L_GEV  = c(CP.g_mc_kde_T[2,], CP.g_mc_tde_T[2,]),
                                 NMix   = c(CP.nm_mc_kde_T[1,], CP.nm_mc_tde_T[1,]),
                                 L_NMix = c(CP.nm_mc_kde_T[2,], CP.nm_mc_tde_T[2,]))

kable(Cov.Prob_mc_diff_T,
      caption = "Multiple comparison to control - kde and tde - difference of quantiles - MVTq",
      format = "pandoc")
write.csv(Cov.Prob_mc_diff_T, "Cov.Prob_mc_diff_T.csv")

##############################################################################################
###                                                                                        ###
###        "Multiple comparison to control - kde and tde - ratio of quantiles - MVNq"      ###
###                                                                                        ###
##############################################################################################



Cov.Prob_mc_ratio_T=data.frame( Density= c(rep("kde",length(N)), rep("tde", length(N))),
                                n      = c(N, N),
                                Norm   = c(CP.n_mc_kde_T[3,], CP.n_mc_tde_T[3,]),
                                L_Norm = c(CP.n_mc_kde_T[4,], CP.n_mc_tde_T[4,]),
                                Exp    = c(CP.e_mc_kde_T[3,], CP.e_mc_tde_T[3,]),
                                L_Exp  = c(CP.e_mc_kde_T[4,], CP.e_mc_tde_T[4,]),
                                Cau    = c(CP.c_mc_kde_T[3,], CP.c_mc_tde_T[3,]),
                                L_Cau  = c(CP.c_mc_kde_T[4,], CP.c_mc_tde_T[4,]),
                                Lap    = c(CP.l_mc_kde_T[3,], CP.l_mc_tde_T[3,]),
                                L_Lap  = c(CP.l_mc_kde_T[4,], CP.l_mc_tde_T[4,]),
                                GEV    = c(CP.g_mc_kde_T[3,], CP.g_mc_tde_T[3,]),
                                L_GEV  = c(CP.g_mc_kde_T[4,], CP.g_mc_tde_T[4,]),
                                NMix   = c(CP.nm_mc_kde_T[3,], CP.nm_mc_tde_T[3,]),
                                L_NMix = c(CP.nm_mc_kde_T[4,], CP.nm_mc_tde_T[4,]))

kable(Cov.Prob_mc_ratio_T,
      caption = "Multiple comparison to control - kde and tde - ratios of quantiles - MVTq",
      format = "pandoc")
write.csv(Cov.Prob_mc_ratio_T, "Cov.Prob_mc_ratio_T.csv")





##############################################################################################
###                                                                                        ###
###         "All pairwise comparison - kde and tde - difference of quantiles - MVNq"       ###
###                                                                                        ###
##############################################################################################


startTime3 <- Sys.time()
system.time({
# kernel density estimate at the median
CP.n_ap_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Normal",        "kde", "MVNq")
CP.e_ap_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Exponential",   "kde", "MVNq")
CP.c_ap_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Cauchy",        "kde", "MVNq")
CP.l_ap_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Laplace",       "kde", "MVNq")
CP.g_ap_kde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "GEV",           "kde", "MVNq")
CP.nm_ap_kde <- sapply(N,Coverage_Prob,R,p, ci_level, "ap", "NormalMixture", "kde", "MVNq")
# True Density estimate at the median
CP.n_ap_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Normal",        "tde", "MVNq")
CP.e_ap_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Exponential",   "tde", "MVNq")
CP.c_ap_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Cauchy",        "tde", "MVNq")
CP.l_ap_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Laplace",       "tde", "MVNq")
CP.g_ap_tde <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "GEV",           "tde", "MVNq")
CP.nm_ap_tde <- sapply(N,Coverage_Prob,R,p, ci_level, "ap", "NormalMixture", "tde", "MVNq")
})
stopTime3 <- Sys.time()
stopTime3 - startTime3
# Time difference of 6.918426 hours

# user   system  elapsed 
# 29535.14     7.33 29668.41
# > 29668/3600
# [1] 8.241111 -- this approximately 8 and half hours

Cov.Prob_ap_diff=data.frame( Density= c(rep("kde",length(N)), rep("tde", length(N))),
                             n      = c(N, N),
                             Norm   = c(CP.n_ap_kde[1,], CP.n_ap_tde[1,]),
                             L_Norm = c(CP.n_ap_kde[2,], CP.n_ap_tde[2,]),
                             Exp    = c(CP.e_ap_kde[1,], CP.e_ap_tde[1,]),
                             L_Exp  = c(CP.e_ap_kde[2,], CP.e_ap_tde[2,]),
                             Cau    = c(CP.c_ap_kde[1,], CP.c_ap_tde[1,]),
                             L_Cau  = c(CP.c_ap_kde[2,], CP.c_ap_tde[2,]),
                             Lap    = c(CP.l_ap_kde[1,], CP.l_ap_tde[1,]),
                             L_Lap  = c(CP.l_ap_kde[2,], CP.l_ap_tde[2,]),
                             GEV    = c(CP.g_ap_kde[1,], CP.g_ap_tde[1,]),
                             L_GEV  = c(CP.g_ap_kde[2,], CP.g_ap_tde[2,]),
                             NMix   = c(CP.nm_ap_kde[1,], CP.nm_ap_tde[1,]),
                             L_NMix = c(CP.nm_ap_kde[2,], CP.nm_ap_tde[2,]))
kable(Cov.Prob_ap_diff,
      caption = "All pairwise comparison - kde and tde - difference of quantiles - MVNq",
      format = "pandoc")
write.csv(Cov.Prob_ap_diff, "Cov.Prob_ap_diff.csv")


##############################################################################################
###                                                                                        ###
###            "All pairwise - kde and tde - difference of quantiles - MVTq"               ###
###                                                                                        ###
##############################################################################################

startTime4 <- Sys.time()
system.time({
# kernel density estimate at the median
CP.n_ap_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Normal",        "kde", "MVTq")
CP.e_ap_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Exponential",   "kde", "MVTq")
CP.c_ap_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Cauchy",        "kde", "MVTq")
CP.l_ap_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Laplace",       "kde", "MVTq")
CP.g_ap_kde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "GEV",           "kde", "MVTq")
CP.nm_ap_kde_T <- sapply(N,Coverage_Prob,R,p, ci_level, "ap", "NormalMixture", "kde", "MVTq")
# True Density estimate at the median
CP.n_ap_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Normal",        "tde", "MVTq")
CP.e_ap_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Exponential",   "tde", "MVTq")
CP.c_ap_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Cauchy",        "tde", "MVTq")
CP.l_ap_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "Laplace",       "tde", "MVTq")
CP.g_ap_tde_T <- sapply(N,Coverage_Prob,R, p, ci_level, "ap", "GEV",           "tde", "MVTq")
CP.nm_ap_tde_T <- sapply(N,Coverage_Prob,R,p, ci_level, "ap", "NormalMixture", "tde", "MVTq")
})
stopTime4 <- Sys.time()
stopTime4 - startTime4
# Time difference of 17.60888 hours

Cov.Prob_ap_diff_T=data.frame( Density= c(rep("kde",length(N)), rep("tde", length(N))),
                               n      = c(N, N),
                               Norm   = c(CP.n_ap_kde_T[1,], CP.n_ap_tde_T[1,]),
                               L_Norm = c(CP.n_ap_kde_T[2,], CP.n_ap_tde_T[2,]),
                               Exp    = c(CP.e_ap_kde_T[1,], CP.e_ap_tde_T[1,]),
                               L_Exp  = c(CP.e_ap_kde_T[2,], CP.e_ap_tde_T[2,]),
                               Cau    = c(CP.c_ap_kde_T[1,], CP.c_ap_tde_T[1,]),
                               L_Cau  = c(CP.c_ap_kde_T[2,], CP.c_ap_tde_T[2,]),
                               Lap    = c(CP.l_ap_kde_T[1,], CP.l_ap_tde_T[1,]),
                               L_Lap  = c(CP.l_ap_kde_T[2,], CP.l_ap_tde_T[2,]),
                               GEV    = c(CP.g_ap_kde_T[1,], CP.g_ap_tde_T[1,]),
                               L_GEV  = c(CP.g_ap_kde_T[2,], CP.g_ap_tde_T[2,]),
                               NMix   = c(CP.nm_ap_kde_T[1,], CP.nm_ap_tde_T[1,]),
                               L_NMix = c(CP.nm_ap_kde_T[2,], CP.nm_ap_tde_T[2,]))

kable(Cov.Prob_ap_diff_T,
      caption = "All pairwise comparison - kde and tde - difference of quantiles - MVTq",
      format = "pandoc")
write.csv(Cov.Prob_ap_diff_T, "Cov.Prob_ap_diff_T.csv")








####
#### Putting the tables together for difference of quantiles
####

equi.quan = c(rep("MVNq",length(N)*2), rep("MVT-q",length(N)*2))
comp_mc_diff   = rep("mc_diff", length(N)*4)
comp_mc_ratio  = rep("mc_ratio", length(N)*4)
comp_ap_diff   = rep("ap_diff"   , length(N)*4)

## For multipe comparison to control
Multi_Comp_control_diff = cbind(Comp = comp_mc_diff, 
                                equi.q = equi.quan,
                                (rbind(Cov.Prob_mc_diff,Cov.Prob_mc_diff_T)))
write.csv(Multi_Comp_control_diff, "Multi_Comp_control_diff.csv")
## For all pairwise comparison
All_pairwise_diff       = cbind(Comp = comp_ap_diff, 
                                equi.q = equi.quan,
                                (rbind(Cov.Prob_ap_diff,Cov.Prob_ap_diff_T)))
write.csv(All_pairwise_diff, "All_pairwise_diff.csv")
## For multiple compariosn to control
Multi_Comp_control_ratio = cbind(Comp = comp_mc_ratio, 
                                equi.q = equi.quan,
                                (rbind(Cov.Prob_mc_ratio,Cov.Prob_mc_ratio_T)))
write.csv(Multi_Comp_control_ratio, "Multi_Comp_control_ratio.csv")

## get six tables in one
six_tables = rbind(Multi_Comp_control_diff, 
                   All_pairwise_diff, 
                   Multi_Comp_control_ratio)
kable(six_tables, "rst")
write.csv(six_tables, "sixtables.csv")
###
### Plot for Coverage Probability
###


par(mfrow = c(2,3))
for(i in 1:6){
  plot(six_tables[, ((i+2)*2)-1], type = "l" ,ylim = c(0.75,1), lwd = 1.5,
       main = paste0(colnames(six_tables)[((i+2)*2)-1]), las = 1,
       ylab = "Coverage Probability")
  abline(h = 0.95,lty = 2, lwd = 2, col = "magenta")
}
par(mfrow = c(1,1))



###
### Plot for Coverage Probability
###

par(mfrow = c(2,3))
for(i in 1:6){
  plot(six_tables[, ((i+2)*2)], type = "l" , lwd = 1.5,
       main = paste0(colnames(six_tables)[((i+2)*2)]), las = 1,
       ylab = "Average Length")
  # abline(h = 0.95,lty = 2, lwd = 2, col = "magenta")
}
par(mfrow = c(1,1))



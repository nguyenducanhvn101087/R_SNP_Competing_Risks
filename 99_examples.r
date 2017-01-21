# Created:  30 Dec 2016
# Modified: 21 Jan 2017
# Anh Nguyen Duc - OUCRU
#================================================================================

# Contents:
# Simulating Competing Risks Data from Underlying SNP-based Model
# - demo of func. simttesnp

# Estimating Competing Risks Data Using SNP-based Model
# - demo of func. snp.crreg

# CIF Comparision using SNP-based IWD
# - demo of func. iwd.snp.onecif

# CLEAR DEATH EXAMPLE
# - demo of func. snp.cif

#================================================================================
rm(list = ls())
source("1_snpcr.r")
# library("timereg")
# require(graphics)
# require(etm)

graphics.off()
#================================================================================



############### Simulating Competing Risks Data from Underlying SNP-based Model ###############
######################### The Important Function is << simttesnp >> ###########################
seed <- 21
set.seed(seed=seed)

J <- 3 # number of competing risks (CRs)
n.dat <- 3e2 # size of data set
marprobs=c(.32, .28, .4) # marginal probabilities i.e. probs of eventually having each of the events

# parameter vectors for each competing risk, 
# for each vector the last 2 elements are mu and sigma,
# the other elements are the spherical coordinates of the resp. snp-polynomial
par.list = list(c(-.1, .1), 
                c(pi/9, -.5, .5), 
                c(pi/9, pi, -.4, .5)) 

# SNP- or error base distribution i.e. 
# dist of Z in log(T) = mu + sig * Z for stdnorm case or 
# dist of exp(Z) for exp case
errdists=c("snp.stdnorm", "snp.exp", "snp.stdnorm")

# initialize the vector of event status and time, 
# where event status is generate according to the pre-specified marginal probabilities
dat.evstat <- sample(1:J, size=n.dat, replace=T, prob=marprobs)
dat.tte    <- rep(NA, n.dat)  

for(j in 1:J) { # loop through each CR   
  
  if(errdists[j]=="snp.stdnorm") { # if the base density is that of a stdnorm
    
    if(length(par.list[[j]])==2) { # check if the resp. snp-polynomial is degenerate or not i.e. if 
      tmp.phi <- NULL              # the resp. spherical coordinates are not NULL
    } else {
      tmp.phi <- head(par.list[[j]], -2)
    }                                          
    # generate time to event type j for those destined to have this event
    dat.tte[dat.evstat==j] <- simttesnp(n   = length(dat.tte[dat.evstat==j]), 
                                        mu  = tail(par.list[[j]],2)[1], 
                                        sig = tail(par.list[[j]],2)[2],
                                        phi = tmp.phi, 
                                        base.dens="stdnorm")
  } else if(errdists[j]=="snp.exp") { # if the base density is that of a stdexp
    
    if(length(par.list[[j]])==2) {
      tmp.phi <- NULL
    } else {
      tmp.phi <- head(par.list[[j]], -2)
    }      
    dat.tte[dat.evstat==j] <- simttesnp(n   = length(dat.tte[dat.evstat==j]), 
                                        mu  = tail(par.list[[j]],2)[1], 
                                        sig = tail(par.list[[j]],2)[2],
                                        phi = tmp.phi,
                                        base.dens="exp")
  }    
  
} # end of for(j in 1:J)

t_m  <- 1.94 # admin. censoring at prefixed time point
rate <- 0.42 # non-admin. independent right censoring following exp dist.

# generate right-censored only data --------------------------------------------------------
dat.rct                        <- pmin(rexp(n=n.dat, rate=rate), t_m)
dat.evstat[dat.tte >= dat.rct] <- 0
dat.tte[dat.evstat==0]         <- dat.rct[dat.evstat==0]

# final right-censored only data set with irrelevant covariates added
dat.rc <- data.frame(tte = dat.tte, evstat = dat.evstat, 
                      x1 = rnorm(length(dat.tte)), 
                      x2 = rexp(length(dat.tte)))


# generate interval censored data --------------------------------------------------------
nints <- 10 # number of intervals
ints  <- seq(from=0, to=t_m, length.out = nints + 1) # set the intervals

# individual-level intervals
ind.ints <- matrix(ints, nrow=1)[rep(1, length(dat.tte)),]

# add some random noise to reflect random visit windows
Noise.mat<- matrix(rnorm(n=(nints-1)*length(dat.tte), mean=0, sd=ints[2]/5), nrow=length(dat.tte))
ind.ints[, 2:(ncol(ind.ints)-1)] <- ind.ints[, 2:(ncol(ind.ints)-1)] + Noise.mat

ind.ints <- t(apply(ind.ints, MARGIN=1, FUN=sort))

L <- R <- rep(NA, length(dat.tte))
R[dat.evstat==0] <- Inf
L[dat.evstat==0] <- dat.tte[dat.evstat==0]

for(r in 1:nrow(ind.ints)) {
  if(dat.evstat[r]!=0) {      
    R[r] <- ind.ints[r,][findInterval(dat.tte[r], ind.ints[r,])+1]
    L[r] <- ind.ints[r,][findInterval(dat.tte[r], ind.ints[r,])]
  } # end ofif(dat.evstat[r]!=0)
} # end of for(r in nrow(ind.ints))

# final interval ensored data set with irrelevant covariates added
# of note, a final column specifying late entry time or left truncation must be given!
dat.ic <- data.frame(L = L, R = R, evstat = dat.evstat, LT=rep(0, n.dat), 
                     x1= rnorm(length(dat.tte)), 
                     x2= rexp(length(dat.tte)))   
dat.ic$L[dat.ic$L==0] <- dat.ic$L[dat.ic$L==0] + .001 # otherwise optimization crashes!


############### Estimating Competing Risks Data Using SNP-based Model ###############
#################### The Important Function is << snp.crreg >> ######################

# use the dat data set above
formula <- Hist(tte, evstat, cens.code=0) ~ 1 # regression model, in our exp there's no covariate
# In general:,this should be 
# A list of formulae given as ouputs from function Hist in package prodlim, 
# the first one is for P(T|D=1), the second one is for
# P(T|D=2). The last one is for the multinomial logistic model for all 
# P(D=j).There must be at least one formula which can be given as a one element
# list or an object of class formula, in that case the models for all P(D) and 
# P(T|D) will follow  that formula. If there are more than one formulae but 
# still fewer fomulae than needed, a model with only intercept is assumed for
# the multinomial logistic model and a NULL model for any remaining AFT models.
# Only the LHS of the first formula will be used to create the time response data.
# To directly specify NULL model for the AFT models a model with only -1 on the
# right-hand side can also be used. The formula for the multinomial logistic component
# can take any variable in the dataset on its left hand side except for `.` and blank.


criterion <- "HQCn" # info criterion used in the greedy algorithm to choose the best model
                    # supported values are: AIC, HQCn, HQCd, BICn and BICd

check.hess <- FALSE # if T (default) estimated model whose hessian matrix has negative eigenvalues will be removed.

parallelism <- TRUE # if to use parallel package
ncores      <- 4    # number of cores used for prallelism if any

kms <- c(3, 3, 3) # A vector having the maximum polynomial degrees in fitting the snp distribution of P(T|D)

anal.grad <- snpcr.grad.loglik # function implementing the analytical gradient for the loglikelihood,. 
                               # This has already been implemented in snpcr.grad.loglik, hence should always be used!


# fit snp model for exp base density --------------------------------------------------------
snpe.mod <- snp.crreg(formula   = formula, data = dat.rc, 
                      base.dens = c("exp", "exp", 'exp'), kms = kms, 
                      criterion = criterion, anal.grad = anal.grad, check.hess = check.hess,
                      parallelism=parallelism, ncores=ncores, 
                      control=list(maxiter=100))
# show result
snpe.mod

# compare the resulting snp polynomial with the real one
polycofsnp(base.dens = "stdnorm", phi = pi/9) # for the 2nd CR
polycofsnp(base.dens = "stdnorm", phi = tail(snpe.mod$final.outs, 1)[[1]]$phis[[2]])


# fit snp model for stdnorm base density --------------------------------------------------------
snpn.mod <- snp.crreg(formula = formula, data = dat.rc, 
                      base.dens = c("stdnorm", "stdnorm", 'stdnorm'), kms = kms, 
                      criterion = criterion, anal.grad = anal.grad, check.hess = check.hess,
                      parallelism=parallelism, ncores=ncores, 
                      control=list(maxiter=100))
# snpn.mod

# fit snp model with mixed base densities matching the true one plus having x1 and x2 -----------

snpm.mod <- snp.crreg(formula = update(formula, ~ . + x1 + x2), data = dat.rc, 
                      base.dens = c("stdnorm", "exp", 'stdnorm'), kms = kms, 
                      criterion = criterion, anal.grad = anal.grad, check.hess = check.hess,
                      parallelism=parallelism, ncores=ncores, 
                      control=list(maxiter=100))
# snpm.mod 

# for interval censored data --------------------------------------------------------
# fit snp model with mixed base densities matching the true one plus having x1 and x2
formula.ic <- Hist(list(L,R), evstat, cens.code=0) ~ 1
snpm.mod.ic <- snp.crreg(formula = update(formula.ic, ~ . + x1 + x2), data = dat.ic, 
                         base.dens = c("stdnorm", "exp", 'stdnorm'), kms = kms, 
                         criterion = criterion, anal.grad = anal.grad, check.hess = check.hess,
                         parallelism=parallelism, ncores=ncores, 
                         control=list(maxiter=100))
# snpm.mod.ic


#################### CIF Comparision using SNP-based IWD ############################
############### The Important Function is << iwd.snp.onecif >> ######################
# ONLY WORKS FOR PURE CIF ESTIMATION I.E. NO COVARIATE INVOLVMENT!

parlist1 <- tail(snpe.mod$final.outs, 1)[[1]]
parlist1$base.dens <- 'exp'

parlist2 <- tail(snpn.mod$final.outs, 1)[[1]]
parlist2$base.dens <- 'stdnorm'

iwd.out <- iwd.snp.onecif(parlist1 = parlist1, parlist2 = parlist2, 
                          j = 2, wfun = function(x) 1, wfix = 1, lower= 0, upper=t_m)
iwd.out


############################## CLEAR DEATH EXAMPLE ################################
#################### Fitting data with interval censoring #########################

dat <- read.csv('cleardeath.csv', header = T)

dat <- subset(dat, subset = r.arm.long != 'Ampho plus fluconazole') # not use this arm
dim(dat) # 175

# Pure CIF estimation for each treatment arm
dat.mono <- subset(dat, subset = r.arm.long == 'Ampho mono') # mono  therapy
dat.comb <- subset(dat, subset = r.arm.long != 'Ampho mono') # combo therapy

formula.app <- Hist(list(tt.clearDeath.lower, tt.clearDeath.upper), ev.clearDeath, cens.code=0) ~ 1

# SNP-based estimation
snp.mod.mono <- snp.crreg(formula   = formula.app, data = dat.mono, 
                          base.dens = c("exp", "exp"), kms = c(3, 3), 
                          criterion = 'HQCn', anal.grad = anal.grad, check.hess = F,
                          parallelism=T, ncores=4, 
                          control=list(maxiter=100))

snp.mod.comb <- snp.crreg(formula   = formula.app, data = dat.comb, 
                          base.dens = c("stdnorm", "stdnorm"), kms = c(3, 3), 
                          criterion = 'HQCn', anal.grad = anal.grad, check.hess = F,
                          parallelism=T, ncores=4, 
                          control=list(maxiter=100))

## Un-comment below to try diff. base densities and compare HQCn
# snp.mod.mono2 <- snp.crreg(formula   = formula.app, data = dat.mono, 
#                           base.dens = c("stdnorm", "stdnorm"), kms = c(3, 3), 
#                           criterion = 'HQCn', anal.grad = anal.grad, check.hess = F,
#                           parallelism=T, ncores=4, 
#                           control=list(maxiter=100))
# 
# snp.mod.comb2 <- snp.crreg(formula   = formula.app, data = dat.comb, 
#                           base.dens = c("exp", "exp"), kms = c(3, 3), 
#                           criterion = 'HQCn', anal.grad = anal.grad, check.hess = F,
#                           parallelism=T, ncores=4, 
#                           control=list(maxiter=100))


# Plot results - Demo of funciton snp.cif --------------------------------------------------
times <- seq(0, 28, length.out = 100)

# get best result from each model 
best.snp.mono <- tail(snp.mod.mono$final.outs, 1)[[1]]
best.snp.comb <- tail(snp.mod.comb$final.outs, 1)[[1]]

# event type 1 is fungal clearance
cif.mono.clear <- snp.cif(tt = times, theta = c(best.snp.mono$phis[[1]], best.snp.mono$mu.sigs[1,]), 
                          base.dens = snp.mod.mono$base.dens[1], bet = NULL, y = NULL, P = best.snp.mono$P[1,1])

cif.comb.clear <- snp.cif(tt = times, theta = c(best.snp.comb$phis[[1]], best.snp.comb$mu.sigs[1,]), 
                          base.dens = snp.mod.comb$base.dens[1], bet = NULL, y = NULL, P = best.snp.comb$P[1,1])

# event type 2 is prior death
cif.mono.death <- snp.cif(tt = times, theta = c(best.snp.mono$phis[[2]], best.snp.mono$mu.sigs[2,]), 
                          base.dens = snp.mod.mono$base.dens[2], bet = NULL, y = NULL, P = best.snp.mono$P[1,2])

cif.comb.death <- snp.cif(tt = times, theta = c(best.snp.comb$phis[[2]], best.snp.comb$mu.sigs[2,]), 
                          base.dens = snp.mod.comb$base.dens[2], bet = NULL, y = NULL, P = best.snp.comb$P[1,2])

x11()
plot(times, cif.mono.clear, type='l', bty= 'n', xlab='Days since randomization', ylab = 'Probability', 
     xlim=c(0,28), ylim=c(0,1), lwd=2, main='CIF of Fungale  Clearance and 1 - CIF of Prior Death')
lines(times, cif.comb.clear, col=gray(.8), lwd=2)

lines(times, 1-cif.mono.death, lwd=2)
lines(times, 1-cif.comb.death, col=gray(.8), lwd=2)
text(x=20, y=0.2, labels = 'Fungal Clearance', cex = 2)
text(x=20, y=0.9, labels = 'Prior Death', cex = 2)
legend(x = 0, y=.7, lwd=2, col = c('black', gray(0.8)), legend = c('Ampho mono', 'Ampho + flucytosine'), bty = 'n')

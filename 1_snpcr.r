
# BASED ON VERSION 2.7.9.3 MAIN - USING MAXLIK

# ANH NGUYEN DUC - OXFORD UNIVERSITY CLINICAL RESEARCH UNIT

# SUPERVISED BY:  MARCEL WOLBERS
# CREATED:        APR    10  2012
# LAST MODIFIED:  OCT    13  2016
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# FUNCTIONS INCLUDED:
# print.SNPMixtureCIFMod   
# print.PredSNPMixtureCIFMod
# predict.SNPMixtureCIFMod (working on)
# snp.crreg           

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Best version, does not have iterlim=12, can allow any K_max, working on LT and ITC
require(cmprsk)
require(eha)
require(parallel)
# require(BB)
# require(Rcgmin)
# require(rgenoud)
require(methods)
require(maxLik)
require(MLEcens)
require(Matrix)
require(matrixcalc)
require(varComp)
require(prodlim)

require(MyRcpp)

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Set up GLOBAL variables , remember to export them to parallel nodes

source.dir <- getwd()

# Cholesky decomposition of moment matrix for stdnorm ---------
k <- 50

xn <- seq(from = 0, to = 2*k+2, by = 1) # Also for function snp.meanvar

# Create a moment sequence from 0 to 2 times k
moment.n <- (1 - abs(2 * round(xn/2) - xn)) * factorial(xn) / 
            (2^(xn/2) * factorial(xn/2))

# Create the moment matrix based on the above moment sequence
A.n <- outer(1:(k+1), 1:(k+1), function(x,y) moment.n[x+y-1])

# Cholesky decomposition
Lc.n <- lapply(1:ncol(A.n), function(i) {
  tmp <- mychol(A.n[1:i, 1:i, drop=F]) 
  
  list(L=tmp$L, Inv=tmp$Inv, InvXL=tmp$Inv %*% tmp$L)
})


# Cholesky decomposition of moment matrix for exp ---------
xe <- seq(from = 0, to = 2*k+2, by = 1) # Also for function snp.meanvar

# Create a moment sequence from 0 to 2 times k
moment.e <- factorial(xe)

# Create the moment matrix based on the above moment sequence
A.e <- outer(1:(k+1), 1:(k+1), function(x,y) moment.e[x+y-1])

# Cholesky decomposition
Lc.e <- lapply(1:ncol(A.e), function(i) {
  tmp <- mychol(A.e[1:i, 1:i, drop=F]) 
  
  list(L=tmp$L, Inv=tmp$Inv, InvXL=tmp$Inv %*% tmp$L)
})


# This is for function snp.meanvar
II <- HH <- rep(NA, 2*k+1)  
II[1] <- 0.5772   # mean of the Gumbel dist with mu=0 and bet=1 (wiki)
HH[1] <- 1.978094 # 2nd central moment (EX^2) of the same dist

for (i in 2:length(II)) {
  II[i] <- (i-1)*II[i-1] - factorial(i-2)
  HH[i] <- -2*II[i-1]+(i-1)*HH[i-1]
} # end of for (i in...)


# Prevent unwanted effect
remove("k", "xn", "A.n", "xe", "A.e", "i") 


#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
source("2_snpcr.uti.r")
source("3_snpcr.grad.r") # This now also has the log-likelihood
source("4_snp.aft.r") 
source("5_snp.r") 


#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# # This hasn't yet considered the case where global.cl exists but stopCluster(global.cl) was called
# if (!exists("global.cl")) {
# 
#   global.cl <- makeCluster(detectCores(), methods=F)
#     
# } else if (exists("global.cl") & class(global.cl)[1] != 'SOCKcluster') {
#   
#   global.cl <- makeCluster(detectCores(), methods=F)
#   
# }# else if (exists("global.cl") & class(global.cl)[1] == 'SOCKcluster') {
#   
# #   try(stopCluster(global.cl), T)
# #   remove(global.cl)
# # }
# 
# # Remember to do the same for cl inside the main function
# clusterExport(cl=global.cl, list("maxLik", "sumKeepAttr", "snp.crlogLik", "mychol",
#                                  "moment.n", "Lc.n", "Lc.e", "moment.e", "II", "HH"))
# clusterEvalQ(global.cl, {
#   source("2_snpcr.uti.r")
#   source("3_snpcr.grad.r")
#   source("4_snp.aft.r") 
#   source("5_snp.r") 
# })  
# 
# 
# 
# # end of if (!exists(global.cl) | (exists(global.cl) & class(global.cl)[1] != 'SOCKcluster')) else ...



#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Create class and the default print function  
# setClass("SNPMixtureCIFMod", representation(final.outs="list", formula="list", base.dens="character",
#                                             criterion="character", dim.gridps="numeric", call="call"))

setClass("SNPMixtureCIFMod", representation(obj="list")) # This is more user-friendly

# The print function should mimic function CSC in package riskRegression for the conditional component
# and function multinom in package nnet for the marginal results
print.SNPMixtureCIFMod <- function(x, digits=2, ...) {
  
  nc <- nrow(x$final.outs[[1]]$mu.sigs)
  bo <- tail(x$final.outs, 1)[[1]]
  
  cat("Call:\n")
  
  print(x$call)
  
  cat("\n")
  cat(paste("Original no.Observations incl. rows with >= 1 NA (if any):", x$n.wna))
  cat("\n")
  
  print(x$response, right=F)
  
  cat("\n")
  
  print(data.frame(Cause=1:nc, Name=getStates(x$response)), row.names=F, right=F)
  
  # cat("\n")
  
  cat("----------------------------------------------------------------------------")
  
  cat("\n\n")
  cat("Multinomial logistic model for the marginal components: P(Cause)")
  cat("\n")
  cat("Model: ")
  print(update(x$formula[[nc+1]], Cause ~ .)) # better replace Cause with the name of the 
  # variable used in Hist(...) that indicates event type
  
  cat("\n")
  
  cat("Coefficients:\n")
  tgam <- bo$gammas  
  
  if (ncol(tgam) == 1) { # if there's only intercept -> give the P
    
    tgam <- cbind(tgam, bo$P[1, 1:(nc-1)])
    colnames(tgam)[length(colnames(tgam))] <- "P(Cause)"
    
  } # end of if (ncol(tgam) == 1)
  
  rownames(tgam) <- unlist(lapply(1:(nc-1), function(j) paste("Cause", j)))
  
  print(tgam, digits=digits)
  
  cat("\n")
  
  cat("Std. Errors:\n")
  tgam.se <- bo$gammas.se
  colnames(tgam.se) <- colnames(bo$gammas)
  rownames(tgam.se) <- rownames(tgam)
  print(tgam.se, digits=digits)
  
  # cat("\n")
  cat("----------------------------------------------------------------------------")      
  
  cat("\n\n")
  cat("AFT models for conditional components: P(Time | Cause)\n")    
  
  for (j in 1:nc) {
    
    cat("\n")
    cat(paste("# Cause ", j, " (", x$base.dens[j], ")", " - Model: ", sep=""))
    print(x$formula[[j]])
    
    cat("\n")
    
    if (length(bo$bets[[j]]) > 0) {
      
      cat("Regression coefficients:\n")
      
      tbet <- t(t(bo$bets[[j]]))
      tbet <- cbind(tbet, t(t(bo$bets.se[[j]])))
      tbet <- cbind(tbet, t(bo$bets.ci[[j]])[, 2:1, drop=F])      
      
      colnames(tbet) <- c("coef", "se(coef)", "lower .95", "upper .95")
      
      print(tbet, digits=digits)
      
    } else {
      
      cat("No covariates specified for this component.")
      cat("\n")
      
    } # end of if (length(bo$bets[[j]]) > 0)        
    
    cat("\n")
    cat("SNP parameters:\n")
    
    tmsi <- t(t(bo$mu.sigs[j,]))
    tmsi <- cbind(tmsi, t(t(bo$mu.sigs.se[j,])))
    rownames(tmsi) <- c(expression(mu), expression(sigma))
    colnames(tmsi) <- c("param", "se(param)")
  
    print(tmsi, digits=digits)
    
    cat("\n")
    cat("SNP squared polynomial: ")
    
    if (length(bo$phis[[j]]) > 0) {
      
      phi <- bo$phis[[j]]
      print(polynomial(polycofsnp(base.dens = x$base.dens[j], phi = phi))^2, digits=digits)
      
    } else {
      
      print(polynomial(c(1,0))^2)
      
    }# end of if (length(bo$phis[[j]]) > 0) else ...    
    
    if (j < nc) cat("............................................................................")    
    else cat("----------------------------------------------------------------------------")    
    
    cat("\n")
  } # end of for (j in 1:length(x$final.outs))
  
  cat("\n")
  cat(paste(x$criterion, ": ", sep=""))
  cat(bo[[x$criterion]])
  
  cat("\t\t\t\t")
  cat(paste("Num. of parameters:", bo$lp))
  
  
} # end of print.SNPMixtureCIFMod
#===============================================================================

# Need a new class PredSNPMixtureCIFMod
setClass("PredSNPMixtureCIFMod", representation(obj="list"))

print.PredSNPMixtureCIFMod <- function(x, digits=2, ...) {
  
  no <- nrow(x$newdata)
  nt <- length(x$time)
  nc <- ncol(x$ps)
  
  cat(paste("Original no.Data rows:", x$n.wna))
  cat("\n\n")
  cat(paste("No.Data rows with >= 1 used covar. being NA:", x$n.wna - no))
  cat("\n")
  cat("............................................................................")
  cat("\n\n")
  
  if (no > nt) { # More rows than time points
    
    cat("Fully complete data:\n")
    print(x$newdata)  
    
    cat("\n\n")
    
    cat("P(Cause):\n")
    
    print(x$ps, digits=digits)
    
    cat("\n\n")
    
    cat("Cumulative incidence function:\n\n")
    
    print(x$cifs, digits=digits)
    
    if (length(x$se.cifs) > 0) {
      
      cat("\n\n")
      
      cat("Standard error:\n\n")
      
      print(x$se.cifs, digits=digits)    
      
    } # end of if (length(x$se.cifs) > 0)

    
  } else { # More time points than rows
    
    for (k in 1:no) {
      
      cat("Input data:\n")
      print(x$newdata[k,])
      
      cat("\n\n")
      
      cat("P(Cause):\n")
      
      print(x$ps[k,], digits=digits)
      
      cat("\n\n")
      
      cat("Cumulative incidence function:\n")

      for (j in 1:nc) {

        cat("\n")
        cat(paste("Cause", j))
        cat("\n")
        
        tab <- matrix(NA, nrow=nt, ncol=3)
        colnames(tab) <- c("CIF Est.", "Lower. 95", "Upper. 95")   
        tab[,1] <- t(x$cifs[[j]][k,,drop=F])
        
        if (length(x$se.cifs) > 0)  {
          tab[,2] <- t(x$lo.cifs[[j]][k,,drop=F])
          tab[,3] <- t(x$up.cifs[[j]][k,,drop=F])
     
        } # end of if (length(x$se.cifs) > 0)

        rownames(tab) <- colnames(x$cifs[[j]])
        
        print(tab, digits=digits)
               
      } # end of for (j in 1:nc)      
      
      cat("............................................................................")
      cat("\n\n") 
    } # end of for (k in 1:no)
    
  } # end of if (no > nt) else ...
  
} # end of print.PredSNPMixtureCIFMod
#===============================================================================

predict.SNPMixtureCIFMod <- function(x, time, newdata, se=F, ci.method="cloglog") { 

  bo <- tail(x$final.outs, 1)[[1]]
  
  nt <- length(time)
  no <- nrow(newdata)
  
  formula <- x$formula
  nc      <- length(formula) - 1
  
  if (is.null(rownames(newdata))) rownames(newdata) <- 1:no

  # remove rows with NA in at least 1 of the relevant covariate(s)
  terms <- attr( delete.response(terms(formula[[1]])), "term.labels")
  allfrm <- paste("~ 1")
  n.wna    <- no
  newdata$id  <- as.integer(rownames(newdata))
  
  if (length(terms) > 0) {
    
    for (i in 2:length(formula)) {
      terms <- c(terms, attr( delete.response(terms(formula[[i]])), "term.labels"))
    } # end of for (i in 2:length(formula))
    terms <- unique(terms)
    
    for (i in 1:length(terms)) {
      allfrm <- paste(allfrm, "+", terms[i])
    } # end of for (i in 1:length(terms))
    
    alldat   <- model.matrix(as.formula(allfrm), data = newdata)
    
    notna.id <- as.integer(rownames(alldat))
    newdata     <- subset(newdata, subset = id %in% notna.id)    
    
    no <- nrow(newdata)
  } # end of if (length(terms) > 0)
  
  gam.dat <- model.matrix(tail(formula, 1)[[1]], data = newdata)
  
  bet.dats<- replicate(nc, matrix(numeric(0), nrow=no))

  for (j in 1:nc) {
    
    tmp <- model.matrix(formula[[j]], data = newdata)
    
    if (ncol(tmp) != 0) {
      # For AFT model even though an intercept is not expected but we can't use the models with -1.
      # Instead we must use the model with an intercept and remove the intercept column from the design matrix. But if 
      # only CIF estimation is desired then we have to use the model with -1 so that the matrix has 0 column which triggers
      # a "null detection" mechanism later on.        
      tmpf <- update(formula[[j]], . ~ . + 1)
      bet.dats[[j]] <- model.matrix(tmpf, data = newdata)[, -1, drop=F]
      
    } # if covariate model is supplied to the current AFT model
    
  } # end of for (j in 1:nc)
  
  ps   <- exp(logprob(X = gam.dat, Gamma = bo$gammas))
  colnames(ps) <- unlist(lapply(1:nc, function(i) paste("Cause", i)))
  
  cifs <- list()
  lo.cifs <- up.cifs <- se.cifs <- list() 

  for (j in 1:nc) {
    
    theta <- c(bo$phis[[j]], bo$mu.sigs[j,])
    cifs[[j]] <- matrix(NA, nrow = no, ncol = nt)
    
    rownames(cifs[[j]]) <- rownames(newdata)
    colnames(cifs[[j]]) <- unlist(lapply(time, function(tt) paste("t=", tt, sep="")))
    
    if (!se) {
      
      for (k in 1:no) {
        cifs[[j]][k,] <- snp.cif(tt = time, theta = theta, base.dens = x$base.dens[j], 
                                 y = bet.dats[[j]][k,,drop=F], bet = bo$bets[[j]], P = ps[k,j,drop=F])
      } # end of for (k in 1:no)
      
      
    } else {      

      lo.cifs[[j]] <- up.cifs[[j]] <- se.cifs[[j]] <- matrix(NA, nrow = no, ncol = nt)
      
      parinfo <- list(nc=nc, gammas=ncol(bo$gammas), phis=unlist(lapply(bo$phis, length)),
                      bets=unlist(lapply(bo$bets, length)))
      
      params  <- c(as.vector(t(bo$gammas)), 
                   as.vector(t(bo$mu.sigs)), 
                   as.vector(t(unlist(bo$phis))), 
                   unlist(bo$bets))           
      
      for (k in 1:no) {
        tmp <- snp.cif.cb(time, params, x$base.dens[j], bet.dats[[j]][k,,drop=F], gam.dat[k,,drop=F],
                          ps[k,,drop=F], cimethod=ci.method, parinfo=parinfo, hessian=bo$hessian,
                          j=j)
        
        cifs[[j]][k,] <- tmp$snp.CIF
        lo.cifs[[j]][k,] <- tmp$low
        up.cifs[[j]][k,] <- tmp$upp
        se.cifs[[j]][k,] <- tmp$se
      } # end of for (k in 1:no)

      rownames(lo.cifs[[j]]) <- rownames(up.cifs[[j]]) <- rownames(se.cifs[[j]]) <- rownames(cifs[[j]])
      colnames(lo.cifs[[j]]) <- colnames(up.cifs[[j]]) <- colnames(se.cifs[[j]]) <- colnames(cifs[[j]])
    } # end of if (!se) else ...
    

    
  } # end of for (j in 1:nc)

  names(cifs) <- colnames(ps)
  
  if (length(lo.cifs) > 0) names(lo.cifs) <- names(up.cifs) <- names(se.cifs) <- names(cifs)
  
  res <- list(time=time, ps=ps, cifs=cifs, lo.cifs=lo.cifs, up.cifs=up.cifs,
              se.cifs=se.cifs, newdata=newdata, n.wna=n.wna, ci.method=ci.method)
  
  class(res) <- "PredSNPMixtureCIFMod"

  res
  
} # end of predict.SNPMixtureCIFMod
#===============================================================================

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# The interface of the main function should mimic that of function CSC in package riskRegression
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# subfunctions using ... : survreg, maxLik, optim.    
snp.crreg <- function(formula, data, base.dens="stdnorm", kms=3, criterion="HQd", 
                      anal.grad=snpcr.grad.loglik, check.hess=T, init.params=NULL, dim.gridps=NULL, 
                      method="BFGS", maxstep=NULL, maxit=100, parallelism=F, ncores=NULL, cl=NULL, ...) {                                                       
  
  # Description: 
  
  # - Find the snp mixture factorization competing risks model. It is advised that user should at least center the covariates
  # in all models.
  
  # - Note that only time independent covariates are dealt with. Currently only right-, interval-censoring
  # and left truncation are covered.
  
  #------------------------------
  # Arguments:

  # - formula                     A list of formulae, the first one is for P(T|D=1), the second one is for
  #                               P(T|D=2). The last one is for the multinomial logistic model for all 
  #                               P(D=j).There must be at least one formula which can be given as a one element
  #                               list or an object of class formula, in that case the models for all P(D) and 
  #                               P(T|D) will follow  that formula. If there are more than one formulae but 
  #                               still fewer fomulae than needed, a model with only intercept is assumed for
  #                               the multinomial logistic model and a NULL model for any remaining AFT models.
  #                               Only the LHS of the first formula will be used to create the time response data.
  #                               To directly specify NULL model for the AFT models a model with only -1 on the
  #                               RHS can also be used. The formula for the multinomial logistic component
  #                               can take any variable in the dataset on its left hand side except for `.` and
  #                               blank. 
  
  # - data (Data)                 A data for fitting the models.
  
  # - base.dens                   A vector of base density used in the models of P(T|D) for each event type.
  #                               Default to stdnorm. Possible values are 'stdnorm' and 'exp'. If there are 
  #                               fewer elements in base.dens than event types then the last element will be
  #                               used for the remaining event types.
  
  # - kms (Kmax)                  A vector having the maximum polynomial degrees in fitting the snp 
  #                               distribution of P(T|D). Default to 3. If there are fewer elements in kms 
  #                               than event types then 3 will be used for the remaining event types. 
  #                               (check positive integer).
  #                               As the moment matrices for stdnorm and exp dist. are currently hardcoded to
  #                               speed up estimation, only k <= 50 is supported.
  
  # - criterion                   Model selection criterion for stepwise selection: AIC, BICn, HQn, BICd or HQd.
  
  # - anal.grad                   If snpcr.grad.loglik (default), analytic gradient of -loglik is used. Other
  #                               implementation can be used but must conform with the data. NULL means using
  #                               numerical optimization.
  
  # - check.hess                  If T (default) estimated parameters whose hessian matrix has negative 
  #                               eigenvalues will be removed.
  
  # - init.params                 (Optional) a list having user selected initial parameters withe elements:
  #
  # $ gammas                      A matrix with rows storing the parameters for the linear predictors in the
  #                               multinomial logistic model for all P(D=j), j=1, ..., nc-1. nc is the # of
  #                               event types.
  #
  # $ mu.sigs                     A matrix with rows storing mu and sigma of the AFT model for P(T|D=j), 
  #                               j=1, ..., nc.
  #
  # $ phis                        A list whose jth element is the vector of spherical coordinates for the 
  #                               snp polynomial of the jth event type.
  #
  # $ bets                        A list whose jth element is the regression coefficients of P(T|D=j).
  #
  #                               Note that if this is set to NOT NULL and the jth snp polynomial degree
  #                               derived from $phis[[j]] is larger than kms[j] then kms[j] is set to 
  #                               length(init.params$phis[[j]]).
  #
  #                               For phis and bets, a NULL member of these lists is theoretically possible.
  #                               However use must supply numeric(0) in this case.
  
  # - dim.gridps                  Vector having the numbers of gridpoints for estimating each element of phi
  #                               when the number of elements is 1, ...., kms. Default to NULL means a grid of 
  #                               16 points when phi is a scalar, 4x4 grids for 2-dimensional phi and 3^n
  #                               grid for phi having n>2 dimension.
  
  # - method                      The methods used by R's maxLik. Default to BFGS.
  
  # - maxstep                     The maximum number of iterations in the forward stepwise procedure.
  #                               BOTH maxstep and kms act as hard stopping criteria simultaneously i.e.
  #                               the procedure stops when either of these is met. DEFAULT to NULL and then will be set to J * kmax
  
  # - parallelism                 A flag indicating whether or not to use parallelism. Default to F.
  
  # - ncores                      The number of cores to be used by parallelism if parallelism is T.
  
  # - cl                          The cluster for parallel (if any). Default value is global.cl which saves time.
  
  # - ...                         Extra arguments for survreg and maxLik.
  
  
  #------------------------------
  # Output:
  
  # A list having sublists for the best fits for each step in the estimation procedure with:
  
  # - call                    the call.
  
  # - models                  a list whose elements are also lists with 2 subelements: a model for P(D) and one for P(T|D).
  
  # - response                the event history response.
  
  # - eventTimes              the sorted (unique) event times.
  
  # - causes                  the event types.
  
  # - base.dens               the vector of base densities.
  
  # - kms                     the vector having the maximum polynomial degrees used to fit the snp distribution of P(T|D).
  
  # - ps                      a matrix having the estimated P(D=j) for each individual in the data.
  
  # - mloglik                 negative log-likelihood.
  
  # - AIC, BICn, BICd, HQCn and HQCd.
      
  # - p                       the number of parameters.
  
  # - params                  a vector having all the parameters with names.
  
  # - outputs from maxLik.
  
  # MORE TO BE ADDED.
  
  # print the result should gives a nice output similar to that of CSC.
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008.
  
  # - ZhangDavidian_SmoothRegression_Appendix.
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data.
  
  # - Upcoming Anh's publication.
  
  #------------------------------
  # Created:       May 28 2012
  # Last modified: Oct 05 2016
  
  tryCatch({ # To catch unprecedented errors   

    if (!is.data.frame(data)) {
      stop("Data must be a data.frame object!")
    } # end of if (!is.data.frame(data))
    
    rownames(data) <- 1:nrow(data)
    
    threedots <- names(list(...))
    
    arg.survreg <- names(formals(survreg))
    arg.maxLik  <- names(formals(maxLik))
    
    id.survreg <- threedots %in% arg.survreg
    id.maxLik  <- threedots %in% arg.maxLik
    
    larg.survreg <- list(...)[id.survreg]
    larg.maxLik  <- list(...)[id.maxLik]      
    
    if (class(formula) == "formula") {
      
      formula <- list(formula)
      
    } else if (class(formula)=="list") {
      
      lapply(formula, 
             function(x) 
               if (class(x)!="formula") 
                 stop("Formula can either be a single formula to be used for all components for all CIFs or a list of formulae")
             )
      
    } else {
      
      stop("Formula can either be a single formula to be used for all components for all CIFs or a list of formulae")
      
    } # end of if (class(formula) == "formula") else ...
    
    call <- match.call()

    res.formula <- reformulate("1", formula[[1]][[2]]) # take the left hand side of the 1st formula. '~' is formula[[1]][[1]]
    mf          <- model.frame(res.formula, data = data, na.action = na.omit)
    response    <- model.response(mf)
 
    unique.types <- getStates(response) # unique event types as passed in by user (character, factor or numeric)

    nc <- length(unique.types) # number of uniques event types or competing risks, old Js already changed to nc
    
    if (is.null(maxstep)) maxstep <- nc * max(kms)
    
    if (nc < 2) stop("There must be at least two event types!")    
        
    if (length(base.dens) < nc) base.dens[(length(base.dens)+1) : nc] <- tail(base.dens, 1)
    
    if (is.null(dim.gridps)) {
      
      dim.gridps <- c(16, 4, 3)
      
      if (length(dim.gridps) < nc) dim.gridps[(length(dim.gridps)+1) : nc] <- rep(3, length.out = nc - length(dim.gridps))
      
    } # end of if (is.null(dim.gridps))
    
    if (length(formula) == 1) {
      
      formula[2 : (nc+1)] <- replicate(nc, tail(formula, 1))
      
    } else if (length(formula) < nc+1) {
      
      formula[(length(formula)+1) : (nc+1)] <- replicate(nc+1-length(formula), update(head(formula,1)[[1]], . ~ -1))
      
      formula[[nc + 1]] <- update(tail(formula,1)[[1]], . ~ +1) 
      # unlike AFT models, the multinomial logistic model still needs an intercept
      
    } # end of if (length(formula) == 1) else ...

    for (i in 2:length(formula)) {
      
      if (length(formula[[i]]) < 3) {
        formula[[i]][[3]] <- formula[[i]][[2]] 
      } # end of if (length(formula[[i]]) < 3)
      
      formula[[i]][[2]] <- formula[[1]][[2]]
      
    } # end of for (i in 1:length(formula))
    
    formula[1:nc] <- lapply(formula[1:nc], function(x) update(x, . ~ . - 1))    
       
    # Remove rows with any NA
    terms <- attr( delete.response(terms(formula[[1]])), "term.labels")
    allfrm <- paste("~ 1")
    n.wna    <- nrow(data)
    data$id  <- as.integer(rownames(data))
    
    if (length(terms) > 0) { # If not pure CIF estimation
      
      for (i in 2:length(formula)) {
        terms <- c(terms, attr( delete.response(terms(formula[[i]])), "term.labels"))
      } # end of for (i in 2:length(formula))
      terms <- unique(terms)
      
      for (i in 1:length(terms)) {
        allfrm <- paste(allfrm, "+", terms[i])
      } # end of for (i in 1:length(terms))
      
      alldat   <- model.matrix(as.formula(allfrm), data=data)
      notna.id <- as.integer(rownames(alldat))
      data     <- subset(data, subset = id %in% notna.id)
      
    } # end of if (length(terms) > 0)

    # data     <- data[notna.id, , drop=F] # this doesn't work
    
    # Recalculate reponse with new data set as the below line did only remove NA cases notified by the first formula
    mf        <- model.frame(res.formula, data = data, na.action = na.omit)
    
    response  <- model.response(mf)
    
    status012 <- response[, "status"] # only 1-event (regardless of cause) and 0-censored, actually 2 for interval censored
    
    ev.facts  <- getEvent(response) # same length as times with factor values in 1,2,..,J and "unknown"
    
    ev.nums   <- as.numeric(ev.facts) # same as ev.facts but members are numeric and "unknown" becomes 0
    ev.nums[ev.nums == max(ev.nums)] <- 0
    
    entry <- rep(0, length(ev.facts))
    
    if (attr(response, "entry.type") == "left-truncated") {
      
      entry <- response[, "entry"]
      
    } # end of if (attr(response, "entry.type") == "left-truncated")
    
    # Check if we're dealing with interval-censoring only  
    if (attr(response, "cens.type") == "intervalCensored") {
      
      is.ic <- T
      L <- response[, "L"] # checked, no matter how we name L and R in the dataset, R gives L and R
      R <- response[, "R"]
      res.dat <- data.frame(L, R, ev.nums, entry, status012)
      
    } else {
      
      is.ic <- F
      times <- response[, "time"] # the time data, what happen if there's interval censoring
      res.dat <- data.frame(times, ev.nums, entry, status012)
      
    } # end of if (attr(response, "cens.type") == "intervalCensored") else ...          
    
    gam.dat <- model.matrix(tail(formula, 1)[[1]], data = data) # already has an intercept
    
    bet.dats<- replicate(nc, matrix(numeric(0), nrow=1))
    
    for (j in 1:nc) {
      
      tmp <- model.matrix(formula[[j]], data = data)
      
      if (ncol(tmp) != 0) {
        # For AFT model even though an intercept is not expected but we can't use the models with -1.
        # Instead we must use the model with an intercept and remove the intercept column from the design matrix. But if 
        # only CIF estimation is desired then we have to use the model with -1 so that the matrix has 0 column which triggers
        # a "null detection" mechanism later on.        
        tmpf <- update(formula[[j]], . ~ . + 1)
        bet.dats[[j]] <- model.matrix(tmpf, data=data)[, -1, drop=F]
        
      } # if covariate model is supplied to the current AFT model
      
    } # end of for (j in 1:nc)
    
    if (length(kms) < nc) kms[(length(kms)+1) : nc] <- 3
      
    if (!is.null(init.params)) { # When the user doesn pass in initial pars
      
      init.gammas <- init.params$gammas
      init.mu.sigs<- init.params$mu.sigs
      init.bets   <- init.params$bets
      init.phis   <- init.params$phis          
      
      if (ncol(init.gammas) != ncol(gam.dat)) {
        stop("Invalid number of variables (columns) in init.params$gammas!")
      }
      
      if (nrow(init.gammas) != nc-1) {
        stop("Invalid number of rows in init.params$gammas!")
      }
      
      if (nrow(init.mu.sigs) != nc) {
        stop("Invalid number of rows in init.params$mu.sigs!")
      }
      
      if (ncol(init.mu.sigs) != 2) {
        stop("Invalid number of columns in init.params$mu.sigs!")
      }   
      
      if (length(init.phis) != nc) {
        stop("Invalid number of members in init.params$phis!")
      }   
      
      if (length(init.bets) != nc) {
        stop("Invalid number of members in init.params$bets!")
      }   
      
      lapply(1:nc, function(i) if (length(init.bets[[i]]) != ncol(bet.dats[[i]])) 
        stop(paste("Invalid length in init.params$bets[", i, "]", sep=""))
      )
      
      init.kms <- unlist(lapply(init.phis, length))
      kms      <- pmax(kms, init.kms)
      
    } else { # If the user does not put in initial parameters

      cat("Initiating parameter values ...\n")
      
      if (!is.ic) { # RC and LT case
        ftime   <- res.dat$times
        fstatus <- res.dat$ev.nums
        
      } else { # only IC case
        ftime   <- (res.dat$L + res.dat$R) / 2
        fstatus <- res.dat$ev.nums    
        ftime[res.dat$R == Inf] <- res.dat$L[res.dat$R == Inf]
        
      } # end of if (!is.ic) else ...
      
      # results of this will be used to get initial values for mus, sigmas and gammas
      tmp.cum <- cuminc(ftime, fstatus, cencode=0)
      
      # In this variant of getting initial values for the intercepts of Gammas,
      # for each set of initial values for Mus and Sigmas, we optimize the 
      # loglikelihood with only the intercepts of Gammas as variable. This
      # will be done later in the first iteration of the algorithm.
      init.gammas <- matrix(0, nrow=nc-1, ncol=ncol(gam.dat))
      
      # However we still need starting values for optim, and this will be taken
      # from the current method of getting initial values for gammas at the first
      # step
      init.ps <- rep(NA, nc) # will also be used for getting the weights for survreg later
      
      # Potential problems will occur if one or more CR doesn't have any
      # corresponding observed event time (not yet addressed, actually user are not 
      # allowed to put in such data)
      for (j in 1:length(init.ps)) {
        init.ps[j] <- max(tmp.cum[[j]]$est)
      }                 
        
      init.ps <- init.ps / sum(init.ps)    
      
      # Get the starting values for mus and sigmas when ks = 0s
      init.mu.sigs <- matrix(NA, nrow=nc, ncol=2)
      init.bets    <- replicate(nc, numeric(0))
      
      if (!is.ic) { # If RC and LT case
        
        # create the weight vectors     
        pDT <- matrix(NA, nrow=nrow(res.dat), ncol=nc)
        for (j in 1:nc) {
          id        <- res.dat$ev.nums==0
          pDT[id,j] <- init.ps[j] - my.timepoints(tmp.cum, res.dat$times[id])$est[j,]
          # for very small and negative number due to numerical rounding and cases
          # where last obs is event
          pDT[id,j] <- ifelse(pDT[id,j] <= 0, .Machine$double.eps, pDT[id,j])
        } # end of for (j in 1:nc)

        for (j in 1:nc) {
          id      <- res.dat$ev.nums==j | res.dat$ev.nums==0
          tmp.dat <- res.dat[id, c("times", "status012", "entry")] 
          
          # create the weight vectors      
          weights <- pDT[id,j] / apply(pDT[id,], 1, sum)
          weights[is.na(weights)] <- 1
          
          if (length(bet.dats[[j]])>0) {
            tmpX <- bet.dats[[j]][id,,drop=F]
            colnames(tmpX) <- paste("X", 1:ncol(tmpX), se="")
            tmp.dat <- cbind(tmp.dat, tmpX)
          } # end of if (length(bet.dats[[j]])>0)
          
          tmp.formula <- Surv(times, status012) ~ . -times -status012 -entry
          larg.survreg$formula <- tmp.formula
          larg.survreg$data    <- data.frame(tmp.dat)
          
          if (base.dens[j]=="stdnorm") { # In case the base density is stdnorm 
            
            larg.survreg$dist <- "lognormal"
            
            if (any(is.na(weights))) {
              
              tmp <- do.call(survreg, larg.survreg)
              
            } else {            
              
              larg.survreg$weights <- weights
              tmp <- do.call(survreg, larg.survreg)
                
            } # end of if (any(is.na(weights))) else ...
            
          } else { # In case the base density is weibull
            
            larg.survreg$dist <- "weibull"
            
            if (any(is.na(weights))) {
              
              tmp <- do.call(survreg, larg.survreg)
              tmp <- survreg(formula, data.frame(tmp.dat), dist="weibull")  
              
            } else {            
              
              larg.survreg$weights <- weights
              tmp <- do.call(survreg, larg.survreg)

            } # end of if (any(is.na(weights))) else ...              
            
          } # end of if (base.dens[j]=="stdnorm")
          
          if (any(is.na(coef(tmp)))) {
            
            larg.survreg$weights <- NULL
            tmp <- do.call(survreg, larg.survreg)
            
          } # end of if (any(is.na(coef(tmp))))
          
          mu                <- coef(tmp)[1]
          attr(mu, "names") <- NULL
          sig               <- tmp$scale  
          bet               <- coef(tmp)[-1]
          init.mu.sigs[j,]  <- c(mu, sig)
          init.bets[[j]]    <- bet
          
        } # end of for (j in 1:nc)
        
      } else { # If only IC case

        # create the weight vectors     
        pDT <- matrix(NA, nrow=nrow(res.dat), ncol=nc)
        for (j in 1:nc) {
          id        <- res.dat$ev.nums==0
          pDT[id,j] <- init.ps[j] - my.timepoints(tmp.cum, res.dat$L[id])$est[j,]
          # for very small and negative number due to numerical rounding and cases
          # where last obs is event
          pDT[id,j] <- ifelse(pDT[id,j] <= 0, .Machine$double.eps, pDT[id,j])
        } # end of for (j in 1:nc)
        
        for (j in 1:nc) {
          
          id      <- res.dat$ev.nums==j | res.dat$ev.nums==0
          tmp.dat <- res.dat[id, c("L", "R", "status012", "entry")] 
          tmp.dat$status012[tmp.dat$status012 == 2] <- 3
          
          # create the weight vectors
          weights <- pDT[id,j] / apply(pDT[id,], 1, sum)
          weights[is.na(weights)] <- 1
          
          if (length(bet.dats[[j]])>0) {
            tmpX <- bet.dats[[j]][id,,drop=F]
            colnames(tmpX) <- paste("X", 1:ncol(tmpX), se="")
            tmp.dat <- cbind(tmp.dat, tmpX)
          }
          
          tmp.formula <- Surv(time=L, time2=R, event=status012, type='interval') ~. -L -R -status012 - entry
          
          # Right now survreg won't work if there is 0 in the left-points
          idl0 <- tmp.dat$L == 0
          tmp.dat$L[idl0] <- tmp.dat$L[idl0] + min(tmp.dat$L[tmp.dat$L != 0]) / 10000
          tmp.dat$R[idl0] <- tmp.dat$R[idl0] + min(tmp.dat$R[tmp.dat$L != 0]) / 10000
          
          larg.survreg$formula <- tmp.formula
          larg.survreg$data    <- data.frame(tmp.dat)
  
          if (base.dens[j]=="stdnorm") { # In case the base density is stdnorm    
            
            larg.survreg$dist <- "lognormal"
            
            if (any(is.na(weights))) {
              
              tmp <- do.call(survreg, larg.survreg)

            } else {
              
              larg.survreg$weights <- weights
              tmp <- do.call(survreg, larg.survreg)

            } # end of if (any(is.na(weights)))            
            
          } else { # In case the base density is weibull
            
            larg.survreg$dist <- "weibull"
            
            if (any(is.na(weights))) {
              
              tmp <- do.call(survreg, larg.survreg)

            } else {
              
              larg.survreg$weights <- weights
              tmp <- do.call(survreg, larg.survreg)

            } # end of if (any(is.na(weights)))
            
          } # End of ifelse
  
          if (any(is.na(coef(tmp)))) {
            
            larg.survreg$weights <- NULL
            tmp <- do.call(survreg, larg.survreg)
            
          } # end of if (any(is.na(coef(tmp))))
          
          mu                <- coef(tmp)[1]
          attr(mu, "names") <- NULL
          sig               <- tmp$scale
          bet               <- coef(tmp)[-1]
          init.mu.sigs[j,]  <- c(mu, sig)
          init.bets[[j]]    <- bet
        } # end of forloop 
        
      } # end of ifelse      
      
      init.gammas[,1] <- log(init.ps[-length(init.ps)] / init.ps[length(init.ps)]) 
      
    } # end of if (!is.null(init.params)) else ...        
    
    # End of setting up initial values for mus, sigmas, betas, gammas
    
    # Set up the initial parameters' info list   
    opt.parinfo        <- list()
    opt.parinfo$nc     <- nc
    opt.parinfo$gammas <- ncol(gam.dat) # correct, the length of each gamma and each
    # gamma has the same length, not the number of gammas
    
    if (is.null(init.params)) {    
      opt.parinfo$phis <- rep(0, nc)      
    } else {
      opt.parinfo$phis <- sapply(init.params$phis, length)  
    }
    
    # Here, if any of the element of Y is NULL or numeric(0) i.e. not a decent
    # matrix, ncol w.r.t to that element will result in NULL!!!
    tmpinfobet  <- sapply(bet.dats, ncol)
    opt.parinfo$bets <- rep(0, nc)

    if (length(tmpinfobet) != 0) {
      for (k in 1:length(tmpinfobet)) {
        if (is.null(tmpinfobet[[k]])) opt.parinfo$bets[k] <- 0
        else opt.parinfo$bets[k] <- tmpinfobet[[k]]
      }    
    }    
    # End of setting up the initial parameters' info list  
 
    # The list storing parameters info w.r.t all optimal output at each iteration
    parinfos <- list()
    
    # The list storing all the best result of each loop
    outs <- list()  
    
    # Flag to stop forward procedure
    halt <- F

    # Running indicator
    step <- 1 # Being used
    
    # The latest best result
    opt.out <- NULL

    use.global.cl <- T
    
    # Initiating parallelism
    if (parallelism & is.null(cl)) {      
      cl <- makeCluster(ifelse(length(ncores)==0, detectCores(), ncores), methods=F) #, outfile="parallelerror.txt")
      clusterExport(cl=cl, list("maxLik", "sumKeepAttr", "snp.crlogLik", "mychol",
                                "moment.n", "Lc.n", "Lc.e", "moment.e", "II", "HH", 
                                "source.dir"))
      # clusterSetRNGStream(cl = cl, iseed = 12345)
      clusterEvalQ(cl, {
        setwd(source.dir)
        source("2_snpcr.uti.r")
        source("3_snpcr.grad.r")
        source("4_snp.aft.r") 
        source("5_snp.r") 
      })    

      use.global.cl <- F
    } # End of if (parallelism & is.null(cl))

    # Get the number of all events
    noe <- length(res.dat$status012[res.dat$status012 != 0])
    
    # Get the number of each event
    evcol <- res.dat$ev.nums[res.dat$ev.nums != 0]
    noes  <- sapply(sort(unique(evcol)), function(i) {length(evcol[evcol==i])})
    
    # Get the number of observations
    N <- nrow(res.dat)
    
    # This variable is a vector whose length is the number of CR (nc), each 
    # element tells if the algorithm should stop increasing the K for certain
    # CR(s). Currently the condition is based on if one K is larger than 2 or not
    kstop   <- rep(F, nc)
    
    # message to say if any warning occurs or maximum # of iterations reached
    message <- ""    
    
    grad.fun<- anal.grad   
    
    obs.res <- method=="BHHH"      
    
    onedims <- list()
    onedims[[1]] <- seq(-1.5, 1.5, length.out=dim.gridps[1])
    for (i in 2:length(dim.gridps)) {
      onedims[[i]] <- seq(-1.5, 1.5, length.out=dim.gridps[i])
    } # end of for (i in 2:length(dim.gridps))
        
    # The forward like procedure  
    # Now we start at phis=NULL for the first step  
    cat("Optimizing ...\n")
    
    while(!halt) {  
      
      cat(paste("Step:", step))
      cat("\n")
      
      phis.grid   <- list()
      
      # Temp list storing all outputs
      tmp.outs    <- list()
      
      # The parameters' info of the previous step will be used for the current step
      tmp.parinfo <- list() 
      
      # Prepare the phis.grid for the first step 
      if (step == 1) { 
        
        tmp.parinfo[[1]] <- opt.parinfo
        
        if (is.null(init.params)) { # If init.params is null
          
          # There is no phis to optimize over when Ks=0s
          phis.grid[[1]] <- list()
          
          for (j in 1:nc) {
            phis.grid[[1]][[j]] <- numeric(0)
          }
          
        } else { # If init.params is NOT null
          
          phis.grid[[1]] <- list()
          
          for (j in 1:nc) {
            phis.grid[[1]][[j]] <- init.params$phis[[j]]
          }
          
        } # End if (is.null(init.params)) else      
        
      } # End of if (step == 1)
       
      # Prepare the phis.grid for an intermediate step 
      if (step > 1) {# 2nd step onward, at this step, for a CR j, its K could 
        # be either 0 or 1. If it is 0 then it will be increased to
        # 1 then we will do as in step 2. If it is 1 then it will be
        # increased to 2 then we will do as Davidian suggested. Now we allow for
        # upto any K_max
        
        tmpphiold <- list() # first get the best fits of spherical coordinates 
        # from the previous step for all CRs    
        
        for (j in 1:nc) {
          # Find the starting position of the previous phis of CR # j
          stapos <- opt.parinfo$gammas * (nc-1) + 2*nc + ifelse(j>1, sum(opt.parinfo$phis[1 : (j-1)]) + 1, 1)          
          
          if (opt.parinfo$phis[j]!=0) {
            tmpphiold[[j]] <- outs[[step-1]]$estimate[stapos:(stapos + opt.parinfo$phis[j] - 1)]
          } else {
            tmpphiold[[j]] <- numeric(0)
          }
          
        } # end of for (j in 1:nc)
        
        for (j in 1:nc) {
          if (!kstop[j]) {# if K_j can still be increased
            lphiold <- length(tmpphiold[[j]])
            
            if (lphiold==0) { # The current CRj has Kj=0
              
              for (i in 1:dim.gridps[1]) {
                
                tmp <- list()
                tmp[[j]] <- onedims[[1]][i]
                
                for (k in (1:nc)[(1:nc)!=j]) {
                  tmp[[k]] <- tmpphiold[[k]]
                }
                
                phis.grid[[length(phis.grid)+1]] <- tmp
                
              }  # end of for (i in 1:len1)
              
            } # end of if (length(tmpphiold[[j]])==0)
            
            if (lphiold>0) { # The current CRj has Kj>=1
              tmp0 <- rep( list( onedims[[lphiold+1]] ), lphiold+1 )
              tmp0 <- do.call(expand.grid, tmp0)
              
              for (i in 1:nrow(tmp0)) {
                
                tmp      <- list() 
                tmp[[j]] <- unlist(tmp0[i,])
                
                for (k in (1:nc)[(1:nc)!=j]) {
                  tmp[[k]] <- tmpphiold[[k]]
                }
                
                phis.grid[[length(phis.grid)+1]] <- tmp                
                
              } # end of for (i in 1:len2)
              
            } # end of if (length(tmpphiold[[j]])==1)
          } # end of if (!kstop[j])
        } # end of for (j in 1:nc)
        
      } # end of if (step > 1)
      
      # Get the parameter list corresponding to each combination of phis
      params <- list()
      
      id.par <- 1
      
      for (i in 1:length(phis.grid)) {  
        
        # For each CR, get the mean and variance of the base line t_0 w.r.t. phis
        tmp.musig     <- matrix(NA, nrow=nc, ncol=2)
        
        # Adjustment for phis and sigmas
        for (j in 1:nc) {
          if (step > 1) { 
            # Actually parlist will be created multiple times, however, just leave
            # it this way for the time being
            parlist <- parvect2parlist(opt.out$estimate, opt.parinfo, exp.sig=T)
            
            tmp.meanvar.Z.pre <- snp.meanvar(parlist$Phi[[j]], base.dens[j])
            
            tmp.mu.logT.pre   <- parlist$Mu_Sig[j,1] + parlist$Mu_Sig[j,2] * tmp.meanvar.Z.pre[1] 
            
            tmp.sig.logT.pre  <- parlist$Mu_Sig[j,2] * sqrt(tmp.meanvar.Z.pre[2])
            
            tmp.musig[j,] <- get.meanvar.on.phi(phis.grid[[i]][[j]], tmp.mu.logT.pre, tmp.sig.logT.pre, base.dens[j])        
            
          } else { # if the first step       
            tmp.musig[j,]     <- init.mu.sigs[j,]
            
          } # end of if (step > 1) else ...        
          
          tmp.musig[j,2] <- log(tmp.musig[j,2]) 
          
        } # end of for (j in 1:nc)
        
        # Create initial values for an intermediate step
        if (step > 1) { # an intermediate step
          # Reset the pooled parameter vector by adding new spherical coordinates
          # Here we must use opt.parinfo (from previous step) instead of tmp.parinfo
          # As we just wanna extract the previous phiss without adding the new ones
          
          parlist <- parvect2parlist(opt.out$estimate, opt.parinfo) 
          # Here it doesn't matter if we set exp.sig=T or not, as we just use gammas
          # and bets
          
          # The following is only needed when there's RC
          if (!is.ic) {
            
            # To know which CRj has Kj increased, we rely on onedimgridp, nc and the current i
            j <- which(sapply(phis.grid[[i]], length) - opt.parinfo$phis == 1) 
            
            # polycofsnp can't deal with phi=numeric(0), but we never meet that
            a <- c(polycofsnp(base.dens=base.dens[j], phi=phis.grid[[i]][[j]])) 
            
            tmp.sig <- tmp.musig[j,2]
            len.mus <- 1 # must be at least 1
            musigj.grid <- matrix(NA, ncol=2, nrow=len.mus*length(tmp.sig))
            musigj.grid[,2] <- rep(tmp.sig, each=len.mus)
            
            for (h in 1:length(tmp.sig)) { 
              tmp.sigj <- exp(tmp.sig[h])                             
              
              # get the corresponding mu and check if it leads to a root of the poly
              tmp.mujs <- rep(NA, len.mus)
              
              if (h != 1) {
                parlist <- parvect2parlist(opt.out$estimate, opt.parinfo, exp.sig=T)
                
                tmp.meanvar.Z.pre <- snp.meanvar(parlist$Phi[[j]], base.dens[j])
                
                tmp.mu.logT.pre   <- parlist$Mu_Sig[j,1] + parlist$Mu_Sig[j,2] * tmp.meanvar.Z.pre[1] 
                
                tmp.mujs[1] <- tmp.mu.logT.pre - tmp.sigj * snp.meanvar(phis.grid[[i]][[j]], base.dens[j])[1]
              } else {
                tmp.mujs[1] <- tmp.musig[j,1]
                
              } # end of if (h != 1) else ...
              
              if (len.mus > 1) {
                for (kk in 2:len.mus) {
                  tmp.mujs[kk] <- tmp.mujs[1] + tmp.sigj * (-1)^kk * ifelse( kk==2 | kk==3, 1/2, 1 ) 
                  
                } # end of for (kk in 2:len.mus)
                
              } # end of if (len.mus > 1)
              
              musigj.grid[ ( (h-1)*len.mus + 1 ):( (h-1)*len.mus + len.mus ), 1] <- tmp.mujs
            } # end of for (h in 1:length(tmp.sig))
            
            for (h in 1 : nrow(musigj.grid)) {                
              tmp.parinfo[[id.par]]     <- opt.parinfo
              tmp.parinfo[[id.par]]$phis<- sapply(phis.grid[[i]], length)
              
              musig.mod       <- tmp.musig
              musig.mod[j,]   <- musigj.grid[h]
              params[[id.par]]<- c(as.vector(t(parlist$Gamma)), as.vector(t(musig.mod)), 
                                   as.vector(t(unlist(phis.grid[[i]]))), 
                                   unlist(parlist$Bet) )
              
              id.par <- id.par + 1
              
            } # end of for (h in 1 : nrow(musigj.grid))
            
          } else { # if we only have interval censoring => no need to do all above
            params[[id.par]]<- c(as.vector(t(parlist$Gamma)), as.vector(t(tmp.musig)), 
                                 as.vector(t(unlist(phis.grid[[i]]))), 
                                 unlist(parlist$Bet))    
            
            tmp.parinfo[[id.par]] <- opt.parinfo
            tmp.parinfo[[id.par]]$phis <- sapply(phis.grid[[i]], length)
            
            id.par <- id.par + 1
            
          } # end of if (!is.ic) else ...
          
        } else { # If the first step
          
          if (is.null(init.params)) { # If this is the first step and there is no
            # initial values input, we diversify the set of starting
            # values by adding more values to mus and sigmas according to
            # Davidian. Also we diversify the starting values for gammass as
            # mentioned earlier
            tmp <- list()
            
            for (j in 1:nc) {
              tmp[[j]]    <- rep(NA, 3)
              tmp[[j]][1] <- tmp.musig[j,1]
              
              tmp[[j]][2] <- tmp[[j]][1] - exp(tmp.musig[j,2])/2
              tmp[[j]][3] <- tmp[[j]][1] + exp(tmp.musig[j,2])/2   
              
            } # end of for (j in 1:nc)
            
            mu.grid <- as.matrix(expand.grid(tmp))          
            
            foo.gamma <- function(gamma, musig, phi, bet, parinfo) {

              -snp.crlogLik(c(as.vector(t(gamma)), musig, phi, bet), res.dat, parinfo, Z=list(X=gam.dat, Y=bet.dats), 
                            base.dens, obs.res=F)
              
            } # end of foo.gamma

            for (k in 1:nrow(mu.grid)) {

              tmp.parinfo[[id.par]] <- opt.parinfo          
              
              npar       <- opt.parinfo$gammas * (nc-1)
              
              control    <- list(ndeps=rep(1e-3, npar), maxit=maxit)              
              
              gamma.outs <- try(optim(init.gammas, foo.gamma, gr=NULL, 
                                      musig=as.vector(rbind(mu.grid[k,], t(tmp.musig)[2,])), 
                                      phi=as.vector(t(unlist(phis.grid[[i]]))),
                                      bet=rep(0,sum(unlist(opt.parinfo$bets))),
                                      parinfo=tmp.parinfo[[id.par]], 
                                      control=control, method="BFGS", hessian=F), silent=T)
              if (class(gamma.outs)[1] != "try-error") {
                tmp.gammas <- gamma.outs$par
              } else {
                tmp.gammas <- init.gammas
              }
              
              params[[id.par]] <- c(as.vector(t(tmp.gammas)), 
                                    as.vector(rbind(mu.grid[k,], t(tmp.musig)[2,])), 
                                    as.vector(t(unlist(phis.grid[[i]]))),
                                    unlist(init.bets))
              
              id.par <- id.par + 1
              
            } # end of for (k in 1:nrow(mu.grid)) 
            
          } else { # if there is initial values
            # need to take the log of the sigma
            tmp.mu.sigs    <- init.params$mu.sigs
            tmp.mu.sigs[,2]<- log(tmp.mu.sigs[,2])
            
            params[[i]]   <- c(as.vector(t(init.params$gammas)), 
                               as.vector(t(tmp.mu.sigs)), 
                               as.vector(t(unlist(phis.grid[[i]]))), 
                               unlist(init.params$bets))            
            # Here we don't need to care about tmp.parinfo as it's already dealt with earlier
            
          } # end of if (is.null(init.params)) else...         
          
        } # end of if (step>1) else ...        
        
      } # End of for (i in 1:length(phis.grid))
      
      # Using parallelism
      # But first if there is only one param or parallel=F we use the conventional method  
      if (length(params) == 1 | parallelism == F) {
        
        AIC <- BICn <- HQCn <- BICd <- HQCd <- rep(NA, length(params))      
        
        for (i in 1:length(params)) {
          larg.maxLik$logLik <- snp.crlogLik
          larg.maxLik$grad   <- grad.fun
          larg.maxLik$start  <- params[[i]]          
          larg.maxLik$method <- method
          larg.maxLik$Data   <- res.dat
          larg.maxLik$parinfo<- tmp.parinfo[[i]]
          larg.maxLik$Z      <- list(X=gam.dat, Y=bet.dats)
          larg.maxLik$Base.dens <- base.dens
          # larg.maxLik$optonsig  <- F
          larg.maxLik$obs.res   <- obs.res
                    
          # at least grad.fun can return NA during optimization
          tmp.outs[[i]] <- try(do.call(maxLik, larg.maxLik), T)                  
          
          mlk <- q <- numeric(0)
          
          if (class(tmp.outs[[i]])[1]!="try-error") {
            mlk <- -1 * tmp.outs[[i]]$maximum
            
            q   <- length(tmp.outs[[i]]$estimate)
          }
          
          if (length(mlk)==0) {
            AIC[i]  <- tmp.outs[[i]]$AIC  <- NA
            BICd[i] <- tmp.outs[[i]]$BICd <- NA
            HQCd[i] <- tmp.outs[[i]]$HQCd <- NA
            BICn[i] <- tmp.outs[[i]]$BICn <- NA
            HQCn[i] <- tmp.outs[[i]]$HQCn <- NA
          } else {
            AIC[i]  <- tmp.outs[[i]]$AIC  <- 2*(mlk + q)
            BICd[i] <- tmp.outs[[i]]$BICd <- 2*mlk + q*log(noe)         
            HQCd[i] <- tmp.outs[[i]]$HQCd <- 2*(mlk + q*log(log(noe)))       
            BICn[i] <- tmp.outs[[i]]$BICn <- 2*mlk + q*log(N)
            HQCn[i] <- tmp.outs[[i]]$HQCn <- 2*(mlk + q*log(log(N))) 
          } # end of if (length(mlk)==0) else ...
          
        } # end of for (i in 1:length(params))
        
      } else { # Now if parallelism is desired and there are more than one params
        
        tmp.outs  <- parLapplyLB(cl, 1:length(params),
                               function(i, res.dat, params, tmp.parinfo, Z, base.dens, method, 
                                        maxit, check.hess, grad.fun, obs.res, larg.maxLik, mychol) {
                                   mychol <- mychol
                                         
                                   larg.maxLik$logLik <- snp.crlogLik
                                   larg.maxLik$grad   <- grad.fun
                                   larg.maxLik$start  <- params[[i]]
                                   larg.maxLik$method <- method
                                   larg.maxLik$Data   <- res.dat
                                   larg.maxLik$parinfo<- tmp.parinfo[[i]]
                                   larg.maxLik$Z      <- list(X=gam.dat, Y=bet.dats)
                                   larg.maxLik$Base.dens <- base.dens
                                   # larg.maxLik$optonsig  <- F
                                   larg.maxLik$obs.res   <- obs.res
                                   
                                   re <- try(do.call(maxLik, larg.maxLik), T)
                                         
                                   mlk <- q <- numeric(0)
                                   
                                   if (class(re)[1]!="try-error") {
                                     mlk <- -1 * re$maximum
                                     
                                     q   <- length(re$estimate)
                                   }
                                   
                                   if (length(mlk)==0) {
                                     re$AIC <- NA
                                     re$BICd<- NA
                                     re$HQCd<- NA
                                     re$BICn<- NA
                                     re$HQCn<- NA
                                   } else {
                                     re$AIC <- 2*(mlk + q)   
                                     re$BICd<- 2*mlk + q*log(noe)
                                     re$HQCd<- 2*(mlk + q*log(log(noe)))
                                     re$BICn<- 2*mlk + q*log(N)
                                     re$HQCn<- 2*(mlk + q*log(log(N)))
                                   } # end of if (length(mlk)==0) else ...
                                   
                                   re
                                },
                                res.dat=res.dat, params=params, tmp.parinfo=tmp.parinfo, 
                                Z=list(X=gam.dat, Y=bet.dats), base.dens=base.dens, method=method, 
                                maxit=maxit, check.hess=check.hess, grad.fun=grad.fun,
                                obs.res=obs.res, larg.maxLik, mychol)  # still not work mychol

        AIC <- unlist(sapply(tmp.outs, function(re) re$AIC))
        BICn<- unlist(sapply(tmp.outs, function(re) re$BICn))
        HQCn<- unlist(sapply(tmp.outs, function(re) re$HQCn))
        BICd<- unlist(sapply(tmp.outs, function(re) re$BICd))
        HQCd<- unlist(sapply(tmp.outs, function(re) re$HQCd))
        
      } # end of if (length(params)==1 | parallel=F)       
      
      # End of using parallelism

      # remove results with non-semidefinite hessian
      idomit <- NULL
      
      for (l in 1:length(tmp.outs)) {
        
        if (length(tmp.outs[[l]]$estimate)>0) { # sometimes the result contains only fields with numeric(0)!       
          
          if (is.na(tmp.outs[[l]]$maximum)) { # some times convergence is 0 or par is not NULL but value is NaN
            
            idomit <- c(idomit, l)
            
          } else {
            
            if (tmp.outs[[l]]$maximum == -Inf) idomit <- c(idomit, l)  
            
          } # end of if (is.na(tmp.outs[[l]]$value)) else ...          
          
          if (check.hess) {
            eig.vals <- try(eigen(tmp.outs[[l]]$hessian)$values, T)
            
            if (class(eig.vals)[1]!="try-error") {
              
              if (any(eig.vals >= 0)) idomit <- c(idomit, l)
              
            } else {
              idomit <- c(idomit, l)
              
            } # end of if (class(eig.vals)...) else...  
            
          } # end of if (check.hess)
          
          if (!is.element(el=tmp.outs[[l]]$code, set=c(0, 1, 2, 8, 4))) idomit <- c(idomit, l)
          
        } else {
          idomit <- c(idomit, l)
          
        } # end of if (length(tmp.outs[[l]]$estimate>0)) else ...
        
      } # end of for (l in 1:....)
      
      idomit <- unique(idomit)
      
      if (length(idomit) == length(tmp.outs)) {
        cat("All initial values at the latest step lead to non-convergent optimization problem!\n")
        return()
      } # end of if (length(idomit) == length(tmp.outs))
      
      if (!is.null(idomit)) {
        tmp.outs <- tmp.outs[-idomit] # this is syntaxtically correct
        AIC  <- AIC[-idomit]
        BICn <- BICn[-idomit]
        HQCn <- HQCn[-idomit]
        BICd <- BICd[-idomit]
        HQCd <- HQCd[-idomit]
        tmp.parinfo <- tmp.parinfo[-idomit]
      }
      
      if (length(AIC) == 0) break()
      
      stop.loop <- F
      
      while(!stop.loop) {
        
        outs[[step]] <- list()
        
        # actually, within the same step it doesn't matter which criterion is used
        id.op <- which.min(AIC)
        
        outs[[step]] <- opt.out <- tmp.outs[[id.op]]
        
        outs[[step]]$AIC  <- opt.out$AIC <- AIC[[id.op]]
        outs[[step]]$BICn <- opt.out$BICn<- BICn[[id.op]]
        outs[[step]]$HQCn <- opt.out$HQCn<- HQCn[[id.op]]
        outs[[step]]$BICd <- opt.out$BICd<- BICd[[id.op]]
        outs[[step]]$HQCd <- opt.out$HQCd<- HQCd[[id.op]]
        
        parinfos[[step]]  <- opt.parinfo <- tmp.parinfo[[id.op]]                      
        
        if (check.hess) {
          
          if (class(outs[[step]]$hessian)[1]=="try-error") {
            outs[[step]]$hessian <- NULL
            
          } else {
            if (any(abs(outs[[step]]$hessian)==Inf) | any(is.na(outs[[step]]$hessian))) { 
              outs[[step]]$hessian <- NULL
              
            } else if (any(eigen(outs[[step]]$hessian)$values >= 0)) {
              outs[[step]]$hessian <- NULL
              
            } # end of if (any(abs(outs[[step]]$hessian)==Inf) | any(is.na(outs[[step]]$hessian)) else...
            
          } # end of if (class(outs[[step]]$hessian)[1]=="try-error") else ...
          
        } # end of if (!check.hess)
        
        stop.loop <- T # always do this here
        
        if (step > 1) { # if not the first step
          halt <- (step > maxstep)
          
          # if stop due to  worse information criterion or if hessian ill-behaved=> not record the last element
          if (outs[[step]][[criterion]] >= outs[[step-1]][[criterion]]) { # if information of the best fit
            # at current step is worse than previous
            halt         <- T
            outs[[step]] <- NULL
            stop.loop    <- T
            
          } else { # if information of currest best fit is better than previous step
            
            if (is.null(outs[[step]]$hessian)) { # if hessian is singular thus NULL
              
              id.op2nd <- which.min(AIC[-id.op]) # find the second best fit (according to AIC or mloglik)
              
              diff.1st2nd <- (AIC[-id.op][id.op2nd]-AIC[id.op]) / AIC[id.op]
              
              if (length(diff.1st2nd)==0)  diff.1st2nd  <- NA # to avoid the case AIC[-id.op]=NULL
              
              if (diff.1st2nd<=1e-2 & !is.na(diff.1st2nd)) { # if the AIC from the 2nd best is not much 
                # worse than the current best fit (both at same step)
                
                tmp.outs <- tmp.outs[-id.op] # remore the current best
                AIC  <- AIC[-id.op]
                BICn <- BICn[-id.op]
                HQCn <- HQCn[-id.op]
                BICd <- BICd[-id.op]
                HQCd <- HQCd[-id.op]
                tmp.parinfo <- tmp.parinfo[-id.op]
                
                outs[[step]] <- NULL # remove the current best
                stop.loop <- F # continue to loop
                
              } else { # if the AIC from next best is much worse than the current best (at the same step)
                
                outs[[step]] <- NULL # remove current best
                halt         <- stop.loop <- T # stop loop and stop algorithm (big loop)
                message      <- "More complex model with better information criterion found but also with singular hessian"
                cat(message)
                cat("\n")
              } # end of if ((tmp.outs[[id.op2nd]]$AIC-tmp.outs[[id.op]]$AIC)<=1e-4) else ...
              
            } else { # if hessian is not singular
              stop.loop <- T # it is time to move to next step or atleast get out of this loop
              
            } # end of if (is.null(outs[[step]]$hessian)) else ...
            
          } # end of if (outs[[step]][[criterion]] >= outs[[step-1]][[criterion]]) else ....
          
        } # end of if (step > 1)
        
      } # end of while(!stop.loop)  
      
      # check if any of the CRs has the corresponding K > the resp. element in kms
      kstop <- opt.parinfo$phis >= kms
      
      step  <- step + 1
      
      if (!halt) {# no need to consider below if already have to stop due to singular hessian and not able 
        # to find "next best candidate" at the current step => go back to previous step
        
        if (all(kstop)) {
          message <- paste(message, "All Ks equal K_maxs!", sep="")
          cat("All Ks equal K_maxs!\n")  
          halt <- T
        } # end of if (all(kstop))
        
        if (step > maxstep) {
          message <- paste(message, "Maximum number of iteration reached!", 
                           sep="")
          cat("Maximum number of iteration reached!\n")
          halt <- T
        } # end of if (step > maxstep)
        
      } # end of if (!halt)

    } # end of while (forward procedure)  
    
    # Always close any created cluster
    if (parallelism & !use.global.cl) {
      stopCluster(cl)
      remove(cl)
    } # end of if (parallelism & !use.global.cl)
    

    # Polishing the outputs    
    final.outs <- list()
    
    for (i in 1:length(outs)) {
      parinfo <- parinfos[[i]]
      params  <- outs[[i]]$estimate
      # Extract parameters from params. Note as we assumed gammas_nc = 0,
      # the number of rows is nc-1 
      parlist <- parvect2parlist(params, parinfo, exp.sig=T)
      
      gammas  <- parlist$Gamma
      
      # This is a nc X 2 matrix
      mu.sigs <- parlist$Mu_Sig
      
      # This is a list with nc elements
      phis    <- parlist$Phi
      
      # Also get back the polynomials coefficients
      pcoef   <- list()
      
      # And the squared version
      spcoef  <- list()
      
      for (j in 1:nc) {
        if (parinfo$phis[j] != 0) {            
          pcoef[[j]]  <- polycofsnp(phi=phis[[j]], base.dens=base.dens[j])
          spcoef[[j]] <- coef(polynomial(coef=pcoef[[j]])^2)
          
        } else {
          pcoef[[j]]  <- numeric(0)
          spcoef[[j]] <- numeric(0)
        }
      }

      # This is a list with nc elements
      bets <- parlist$Bet
      
      P    <- exp(logprob(gam.dat, gammas))
      
      final.outs[[i]]         <- outs[[i]]
      
      colnames(gammas)        <- colnames(gam.dat)
      final.outs[[i]]$gammas  <- gammas
      
      final.outs[[i]]$mu.sigs <- mu.sigs        
      
      final.outs[[i]]$phis    <- phis
      final.outs[[i]]$pcoef   <- pcoef
      final.outs[[i]]$spcoef  <- spcoef
      
      for (j in 1:nc) {
        names(bets[[j]]) <- colnames(bet.dats[[j]])
      } # end of for (j in 1:nc)
      final.outs[[i]]$bets    <- bets    
      
      final.outs[[i]]$P       <- P
      
      final.outs[[i]]$mloglik <- -outs[[i]]$maximum
      
      final.outs[[i]]$AIC  <- outs[[i]]$AIC
      final.outs[[i]]$BICn <- outs[[i]]$BICn
      final.outs[[i]]$HQCn <- outs[[i]]$HQCn
      final.outs[[i]]$BICd <- outs[[i]]$BICd
      final.outs[[i]]$HQCd <- outs[[i]]$HQCd
      
      final.outs[[i]]$noe  <- noe
      
      final.outs[[i]]$hessian     <- outs[[i]]$hessian
      final.outs[[i]]$gradient.mat<- outs[[i]]$gradientObs
      
      if (method!="BHHH") 
        final.outs[[i]]$gradient.mat <- try(snpcr.grad.loglik(params=params, Data=res.dat, parinfo=parinfo,
                                                              Z=list(X=gam.dat, Y=bet.dats), Base.dens=base.dens, 
                                                              obs.res=T), T)  
      
      final.outs[[i]]$gradient <- outs[[i]]$gradient      
      
      K <- 0
      for (h in 1:N) {
        K <- K + final.outs[[i]]$gradient.mat[h,] %*% t(final.outs[[i]]$gradient.mat[h,])
      }
      
      final.outs[[i]]$K           <- K
      
      final.outs[[i]]$iterations  <- outs[[i]]$iterations
      final.outs[[i]]$code        <- outs[[i]]$code    
      final.outs[[i]]$lp          <- length(params)
      final.outs[[i]]$params      <- params # params always has the centered mus and gamma_0s        
      
      if (!is.null(outs[[i]]$hessian)) {
        # As part of the hessian matrix is calculated w.r.t. sigmas on logscale,  
        # the Jacobian is not really an identity matrix        
        foo <- function(x, parinfo) {

          pl <- parvect2parlist(x, parinfo, exp.sig=T)        
          
          x <- c(as.vector(t(pl$Gamma)), 
                 as.vector(t(pl$Mu_Sig)), 
                 as.vector(t(unlist(pl$Phi))), 
                 unlist(pl$Bet))
        }
        
        # invhes <- try(svd.inverse(final.outs[[i]]$hessian), T)
        invhes <- try(-mychol(-final.outs[[i]]$hessian)$Inv, T)
        # For negative definite matrix like this hessian we need the double - trick otherwise mychol's result
        # differs substantially from solve or svd.inverse
        
        jcb   <- jacobian(foo, final.outs[[i]]$params, parinfo=parinfo) # delta method

        covar <- try(jcb %*% (-invhes)  %*% t(jcb), T) # correct
        
        final.outs[[i]]$cov <- covar
        
        if (class(covar)[1]!="try-error") {
          
          se <- ifelse(diag(covar)<0, Inf, sqrt(diag(covar)))
          
          gammas.se <- t( matrix(se[1:( parinfo$gammas * (nc-1) )], nrow=parinfo$gammas) )
          
          alpha     <- 0.05
          
          # Also get the 95% CIs for gammass
          gammas.up <- gammas + qnorm(1-alpha/2) * gammas.se
          gammas.lw <- gammas - qnorm(1-alpha/2) * gammas.se
          
          gammas.ci <- NULL
          for (j in 1:(nc-1)) {
            gammas.ci <- rbind(gammas.ci, gammas.up[j,], gammas.lw[j,])
          }
          
          # Reset the counter
          cnt <- parinfo$gammas * (nc-1)
          
          mu.sigs.se <- t( matrix(se[( cnt + 1 ):( cnt + 2*nc)], nrow=2) )
          # Also get the 95% CIs for mu.sigss
          mu.sigs.up <- mu.sigs + qnorm(1-alpha/2) * mu.sigs.se
          mu.sigs.lw <- mu.sigs - qnorm(1-alpha/2) * mu.sigs.se
          
          mu.sigs.ci   <- NULL
          for (j in 1:nc) {
            mu.sigs.ci <- rbind(mu.sigs.ci, mu.sigs.up[j,], mu.sigs.lw[j,])
          }
          
          cnt <- cnt + 2 * nc + sum(unlist(sapply(1:nc, function(j)
                                                        {parinfo$phis[j]})))
          
          # This is a list with nc elements
          bets.se <- list()
          # Also get the 95% CI for bets
          bets.ci <- list() 
          
          if (length(parinfo$bets)!=0) {
            for (j in 1:nc) {
              bets.ci[[j]] <- numeric(0)
              
              if (parinfo$bets[j]!=0) {
                bets.se[[j]] <- se[(cnt + 1):(cnt + parinfo$bets[j])]
                bets.up <- bets[[j]] + qnorm(1-alpha/2) * bets.se[[j]]
                bets.lw <- bets[[j]] - qnorm(1-alpha/2) * bets.se[[j]]
                bets.ci[[j]] <- rbind(bets.up, bets.lw)
              } else {
                bets.se[[j]] <- NULL
              }
              cnt <- cnt + parinfo$bets[j]
            }      
          }
          
          final.outs[[i]]$gammas.se    <- gammas.se
          final.outs[[i]]$gammas.ci    <- gammas.ci
          
          final.outs[[i]]$mu.sigs.se   <- mu.sigs.se
          final.outs[[i]]$mu.sigs.ci   <- mu.sigs.ci
          
          final.outs[[i]]$bets.se      <- bets.se
          final.outs[[i]]$bets.ci      <- bets.ci
          
        } # end of if (class(covar)[1]!="try-error") 
        
        
        # For sandwich SE
        covar.sw<- try(jcb %*% invhes %*% final.outs[[i]]$K %*% invhes %*% t(jcb), silent=T) # all correct
        
        final.outs[[i]]$cov.sw <- covar.sw
        
        if (class(covar.sw)[1] != "try-error") {
          
          se.sw  <- ifelse(diag(covar.sw)<0, Inf, sqrt(diag(covar.sw)))
          
          gammas.se.sw   <- t(matrix(se.sw[1:( parinfo$gammas * (nc-1) )], 
                                     nrow=parinfo$gammas) )
          
          alpha      <- 0.05
          
          # Also get the 95% CIs for gammas
          gammas.up.sw   <- gammas + qnorm(1-alpha/2) * gammas.se.sw
          gammas.lw.sw   <- gammas - qnorm(1-alpha/2) * gammas.se.sw
          
          gammas.ci.sw   <- NULL
          for (j in 1:(nc-1)) {
            gammas.ci.sw <- rbind(gammas.ci.sw, gammas.up.sw[j,], gammas.lw.sw[j,])
          }
          
          # Reset the counter
          cnt <- parinfo$gammas * (nc-1)
          
          mu.sigs.se.sw <- t( matrix(se.sw[( cnt + 1 ):( cnt + 2*nc)], nrow=2) )
          # Also get the 95% CIs for mu.sigss
          mu.sigs.up.sw <- mu.sigs + qnorm(1-alpha/2) * mu.sigs.se.sw
          mu.sigs.lw.sw <- mu.sigs - qnorm(1-alpha/2) * mu.sigs.se.sw
          
          mu.sigs.ci.sw <- NULL
          for (j in 1:nc) {
            mu.sigs.ci.sw <- rbind(mu.sigs.ci.sw, mu.sigs.up.sw[j,], mu.sigs.lw.sw[j,])
          }
          
          cnt <- cnt + 2 * nc + sum(unlist(sapply(1:nc, function(j)
                                                        {parinfo$phis[j]})))
          
          # This is a list with nc elements
          bets.se.sw <- list()
          # Also get the 95% CI for bets
          bets.ci.sw <- list() 
          
          if (length(parinfo$bets)!=0) {
            for (j in 1:nc) {
              bets.ci.sw[[j]] <- numeric(0)
              
              if (parinfo$bets[j]!=0) {
                bets.se.sw[[j]] <- se.sw[(cnt + 1):(cnt + parinfo$bets[j])]
                bets.up.sw <- bets[[j]] + qnorm(1-alpha/2) * bets.se.sw[[j]]
                bets.lw.sw <- bets[[j]] - qnorm(1-alpha/2) * bets.se.sw[[j]]
                bets.ci.sw[[j]] <- rbind(bets.up.sw, bets.lw.sw)
              } else {
                bets.se.sw[[j]] <- NULL
              }
              cnt <- cnt + parinfo$bets[j]
            }      
          }
          
          final.outs[[i]]$gammas.se.sw  <- gammas.se.sw
          final.outs[[i]]$gammas.ci.sw  <- gammas.ci.sw
          
          final.outs[[i]]$mu.sigs.se.sw <- mu.sigs.se.sw
          final.outs[[i]]$mu.sigs.ci.sw <- mu.sigs.ci.sw
          
          final.outs[[i]]$bets.se.sw    <- bets.se.sw
          final.outs[[i]]$bets.ci.sw    <- bets.ci.sw
          
        } # end of if (class(covar.sw)[1]!="try-error")        
        
      } # end of if (!is.null(Outs[[i]]$hessian))
      
      if (i == length(outs)) final.outs[[i]]$message <- message           
      
    } # end of for (i in 1:length(outs))
    
    res <- list(final.outs=final.outs, formula=formula, response=response,
                base.dens=base.dens, criterion=criterion, dim.gridps=dim.gridps, call=call,
                n.wna=n.wna)
    
    class(res) <- "SNPMixtureCIFMod"
    
    res

  }, error=function(ee) {cat("From function: snp.crreg"); ee; browser()})
  
} # End of snp.crreg

#===============================================================================

# snp.crloglik moved to gradcr.loglik file
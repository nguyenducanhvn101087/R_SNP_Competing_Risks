
# VERSION 2.5 - ONLY AFT

# ANH NGUYEN DUC - OXFORD UNIVERSITY CLINICAL RESEARCH UNIT

# SUPERVISED BY:  MARCEL WOLBERS
# CREATED:        May      31  2013
# LAST MODIFIED:  Jun      04  2015
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# FUNCTIONS INCLUDED:

# logf_K.aft                          (tested 27 Jul 2012)
# logS_Kdif.aft
# logS_K.aft                          (tested 30 Jul 2012)
# logF_K.aft
# S_K.aft                             (tested 30 Jul 2012)
# F_K.aft
# f_K.aft                             (tested 27 Jul 2012)
# h_K.aft

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# require(polynom)
# require(numDeriv)
# require(survival)
# require(Rlab)

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



#===============================================================================

logf_K.aft <- function(tt=NULL, theta=NULL, base.dens=NULL, bet=NULL, x=NULL, 
                       bx=NULL, e=NULL, logf_0tt=NULL, logf_0z=NULL, 
                       S_0tt=NULL, cc=NULL, a=NULL) {
  # The estimated LOG density function of t with respect to the base density 
  # specified in base.dens for any of the following regression models: PH, AFT and 
  # PO
  
  # Note that only time independent covariates are dealt with
  #------------------------------
  # Caller:
  
  # - snp.logLik 
  
  #------------------------------
  # Input:
  
  # - t           the time
  
  # - theta       the parameter vectors containing
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - bet         the vector of regression coefficients (beta)
  
  # - x           the design matrix where each row is the covariate vector of the
  #               corresponding individuals
  
  # To improve speed performance of snp.logLik, this function will take the 
  # following extra arguments which should be computed inside snp.logLik. Otherwise
  # , they will be recomputed inside this function.
  
  # - bx          bet %*% x
  
  # - e           exp(bx)
  
  # - logf_0tt    logf_0(tt, theta, base.dens)
  
  # - logf_0z     logf_0(tt * (1/e), theta, base.dens) (only needed when mod is 
  #               "AFT"
  
  # - S_0tt       S_0(tt, theta, base.dens) (only needed when mod is "PH"
  #               or "PO"
  
  # - cc          only needed for AFT mod. When base.dens is stdnorm, cc is 
  #               log(tt*exp(-XB)-mu)/sig = log(tt-XB-mu)/sig,
  #               and when base.dens is exp, z is 
  #               (tt*exp(-XB)/exp(mu))^(1/sig) = (tt*/exp(mu+XB))^(1/sig)
  
  # - a           the polynomial coefficients
  
  #------------------------------
  # Output:
  
  # - The value of the estimated LOG density at time t (>0) given the regression 
  # coefficients and predictors
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Dec 30 2010
  # Last modified: Sep 10 2013
  
  
  tryCatch({ # To catch unprecedented errors
    
    if ((length(bet)==0 | is.null(x) | prod(dim(x))==0) & (is.null(bx) & is.null(e) & is.null(logf_0z))) { 
      if (is.null(logf_0tt)) logf_0tt <- logf_0(tt, theta, base.dens, a)
      
      logf_K <- logf_0tt           
      
    } else {   
      # Accelerated failure time model, for this we don't case f_0, but directly calculate as
      # this is used for the CR models
      if (is.null(bx)) bx<- x %*% bet              
      if (is.null(e)) e  <- exp(bx)
      
      # In case the polynomial is degenerate, or theta contains only mu and sigma
      phi <- NULL
      
      l <- length(theta)
      if (l>2) {
        phi <- theta[1:(l-2)]     
      } # end of if (l>2)
      
      mu  <- theta[l-1]
      sig <- theta[l]
      # if (sig < 0) stop("sigma is negative!")
                  
      # Change of variable for each base-line case  
      if (base.dens=="stdnorm") {          
        if (is.null(cc)) cc   <- (log(tt)-bx-mu)/sig
        cons <- -(log(tt) - bx) - log(sig)
        
      } else {
        if (is.null(cc)) cc   <- (tt/exp(mu+bx))^(1/sig)
        cons <- -mu/sig - log(sig) + 
          (1/sig - 1) * (log(tt) - bx) 
        
      } # end of if (base.dens=="stdnorm") else..        
      
      if (is.null(logf_0z)) { 

        logf_0z <- cons + logdsnp(cc, phi, base.dens, a)          
      }        
      
      logf_K    <- -bx + logf_0z
      
    } # End of ifelse
    
    #c(logf_K)
    logf_K
    
  }, error=function(ee) {print("From function: logf_K.aft"); print(ee); browser()})
  
} # End of logf_K.aft

#===============================================================================
# Old name: logreg.surdif.est
logS_Kdif.aft <- function(L, R, theta, base.dens=NULL, bet=NULL, x=NULL, 
                      bx=NULL, e=NULL, S_0L=NULL, S_0Lz=NULL, ccL=NULL, ccR=NULL,
                      a=NULL) {#, roundoff=F, digits=10) {
  # The estimated LOG of the difference in survival function between L(left) and 
  # R(right) with respect to the base density specified in base.dens for any of 
  # the following regression models: PH, AFT and PO
  
  # Note that only time independent covariates are dealt with
  #------------------------------
  # Caller:
  
  # - snp.logLik
  
  #------------------------------
  # Input:
  
  # - L           the left time point
  
  # - R           the right time point
  
  # - theta       the parameter vectors containing
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - bet         the vector of regression coefficients (beta)
  
  # - x           the design matrix where each row is the covariate vector of the
  #               corresponding individuals

  # To improve speed performance of snp.logLik, this function will take the 
  # following extra arguments which should be computed inside snp.logLik. Otherwise
  # , they will be recomputed inside this function.
  
  # - bx          bet %*% x
  
  # - e           exp(bx)
  
  # - S_0Lz        S_0(L * (1/e), theta, base.dens) (only needed when mod is 
  #               "AFT"
  
  # - S_0L       S_0(L, theta, base.dens) (only needed when mod is "PH"
  #               or "PO"
  
  # - ccL         only needed for AFT mod. When base.dens is stdnorm, ccL is 
  #               log(L*exp(-XB)-mu)/sig = log(L-XB-mu)/sig,
  #               and when base.dens is exp, z is 
  #               (L*exp(-XB)/exp(mu))^(1/sig) = (L*/exp(mu+XB))^(1/sig)
  
  # - ccR         same as ccL but for R
  
  # - a           the polynomial coefficients

  #------------------------------
  # Output:
  
  # - The value of the estimated LOG of the difference in survival between the 
  # left time point (L) and right time point (R) given the regression coefficients 
  # and predictors
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Dec 31 2010
  # Last modified: Jun 04 2015
  
  tryCatch({

    logS_Kdif <- ifelse(R!=Inf,
                        log(-F_K.aft(L, theta, base.dens, bet, x, 
                                 bx, e, S_0L, S_0Lz, ccL, a)
                            + F_K.aft(R, theta, base.dens, bet, x, bx, e, NULL, NULL, ccR, a)),
                        logS_K.aft(L, theta, base.dens, bet, x, 
                               bx, e, S_0L, S_0Lz, ccL, a)
    ) 
    
    # c(logS_Kdif)
    logS_Kdif
  }, error=function(ee) {print("From function: logS_Kdif.aft"); ee; browser()})
  
} # End of logS_Kdif.aft

#===============================================================================

logS_K.aft <- function(tt, theta, base.dens=NULL, bet=NULL, x=NULL, 
                   bx=NULL, e=NULL, S_0tt=NULL, S_0z=NULL, cc=NULL, a=NULL) {
                   #roundoff=F, digits=10) {
  # The estimated LOG survival function of t with respect to the base density
  # specified in base.dens for any of the following regression models: PH, AFT and 
  # PO
  
  # Note that only time independent covariates are dealt with
  #------------------------------
  # Caller:
  
  # - snp.logLik
  
  #------------------------------
  # Input:
  
  # - tt          the time
  
  # - theta       the parameter vectors containing
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - bet         the vector of regression coefficients (beta)
  
  # - x           the design matrix where each row is the covariate vector of the
  #               corresponding individuals  
  
  # To improve speed performance of snp.logLik, this function will take the 
  # following extra arguments which should be computed inside snp.logLik. Otherwise
  # , they will be recomputed inside this function.
  
  # - bx          bet %*% x
  
  # - e           exp(bx)
  
  # - S_0tt       S_0(tt, theta, base.dens) (only needed when mod is "PH"
  #               or "PO" or when bet or x is NULL
  
  # - S_0z        S_0(tt * (1/e), theta, base.dens) (only needed when mod is 
  #               "AFT"
  
  # - cc          only needed for AFT mod. When base.dens is stdnorm, cc is 
  #               log(tt*exp(-XB)-mu)/sig = log(tt-XB-mu)/sig,
  #               and when base.dens is exp, z is 
  #               (tt*exp(-XB)/exp(mu))^(1/sig) = (tt*/exp(mu+XB))^(1/sig)
  
  # - a           the polynomial coefficients

  #------------------------------
  # Output:
  
  # - The value of the estimated LOG survival at time t (>0) given the regression 
  # coefficients and predictors
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Dec 30 2010
  # Last modified: Jun 04 2015
  
  tryCatch({ # To catch unprecedented errors
      
    logS_K <- rep(NaN, length(tt))    
    
    if ((length(bet)==0 | is.null(x) | prod(dim(x))==0) & (is.null(bx) & is.null(e) & is.null(S_0z))) {
      if (is.null(S_0tt)) S_0tt <- S_0(tt, theta, base.dens, a)#, roundoff, digits)
      
      logS_K <- log(S_0tt) 
      
    } else { 
      
      # Accelerated failure time model
      if (is.null(S_0z)) {
        
        l   <- length(theta)
        
        # In case the polynomial is degenerate, or theta contains only mu and sigma
        phi <- NULL
        
        if (length(theta)!=2) {
          phi <- theta[1:(l-2)]           
        }
        
        if (is.null(cc)) {
          
          if (is.null(bx)) bx <- x %*% bet              
          if (is.null(e)) e   <- exp(bx)            
          
          mu  <- theta[l-1]
          sig <- theta[l]
          # if (sig < 0) stop("sigma is negative!")
          
          if (base.dens == "stdnorm") {
            if (is.null(cc)) {
              cc   <- (log(tt)-bx-mu)/sig
              S_0z <- sur.snp(cc, phi, base.dens, a)#, roundoff, digits) 
            } 
          } else {
            if (is.null(cc)) {
              cc   <- (tt/exp(mu+bx))^(1/sig) # taking log then exp here doesn't help
              S_0z <- sur.snp(cc, phi, base.dens, a)#, roundoff, digits) 
            } 
          } # end of if (base.dens=="stdnorm") else..
          
        } else {
          S_0z <- sur.snp(cc, phi, base.dens, a)#, roundoff, digits) 
        }# end of if (is.null(cc)) else ...
        
      }# end of if (is.null(S_0z))
      
      logS_K    <- log(S_0z) #my.log(S_0z)             
      
    } # End of ifelse      
    
    # c(logS_K)
    logS_K
    
  }, error=function(ee) {print("From function: logS_K.aft"); print(ee); browser()})
  
} # End of logS_K.aft

#===============================================================================
logF_K.aft <- function(tt, theta, base.dens=NULL, bet=NULL, x=NULL, 
                       bx=NULL, e=NULL, F_0tt=NULL, F_0z=NULL, cc=NULL, a=NULL) {
                       #roundoff=F, digits=10) {
  # The estimated LOG cdf of t with respect to the base density
  # specified in base.dens for any of the following regression models: PH, AFT and 
  # PO
  
  # Note that only time independent covariates are dealt with
  #------------------------------
  # Caller:
  
  # - snp.logLik
  
  #------------------------------
  # Input:
  
  # - tt          the time
  
  # - theta       the parameter vectors containing
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - bet         the vector of regression coefficients (beta)
  
  # - x           the design matrix where each row is the covariate vector of the
  #               corresponding individuals
  
  # To improve speed performance of snp.logLik, this function will take the 
  # following extra arguments which should be computed inside snp.logLik. Otherwise
  # , they will be recomputed inside this function.
  
  # - bx          bet %*% x
  
  # - e           exp(bx)
  
  # - F_0tt       S_0(tt, theta, base.dens) (only needed when mod is "PH"
  #               or "PO" or when bet or x is NULL
  
  # - F_0z        S_0(tt * (1/e), theta, base.dens) (only needed when mod is 
  #               "AFT"
  
  # - cc          only needed for AFT mod. When base.dens is stdnorm, cc is 
  #               log(tt*exp(-XB)-mu)/sig = log(tt-XB-mu)/sig,
  #               and when base.dens is exp, z is 
  #               (tt*exp(-XB)/exp(mu))^(1/sig) = (tt*/exp(mu+XB))^(1/sig)
  
  # - a           the polynomial coefficients

  #------------------------------
  # Output:
  
  # - The value of the estimated LOG cdf at time t (>0) given the regression 
  # coefficients and predictors
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Feb 19 2013
  # Last modified: Jun 04 2015
  tryCatch({ # To catch unprecedented errors        
      
    logF_K <- rep(NaN, length(tt))

    if ((length(bet)==0 | is.null(x) | prod(dim(x))==0) & (is.null(bx) & is.null(e) & is.null(F_0z))) {
      if (is.null(F_0tt)) F_0tt <- F_0(tt, theta, base.dens, a)#, roundoff, digits)
      
      logF_K <- log(F_0tt) # my.log(F_0tt)
      
    } else { 
      
      if (is.null(F_0z)) {
        
        l   <- length(theta)
        
        # In case the polynomial is degenerate, or theta contains only mu and sigma
        phi <- NULL
        
        if (length(theta)!=2) {
          phi <- theta[1:(l-2)]           
        }
        
        if (is.null(cc)) {
          
          # if (is.null(bx)) bx       <- x[id,,drop=F] %*% bet 
          if (is.null(bx)) bx <- x %*% bet 
          if (is.null(e)) e   <- exp(bx)            
          
          mu  <- theta[l-1]
          sig <- theta[l]
          # if (sig < 0) stop("sigma is negative!")
          
          if (base.dens == "stdnorm") {
            if (is.null(cc)) {
              cc   <- (log(tt)-bx-mu)/sig
              F_0z <- psnp(cc, phi, base.dens, a)#, roundoff, digits) 
            } 
          } else {
            if (is.null(cc)) {
              cc   <- (tt/exp(mu+bx))^(1/sig) # taking log then exp here doesn't help
              F_0z <- psnp(cc, phi, base.dens, a)#, roundoff, digits) 
            } 
          } # end of if (base.dens=="stdnorm") else..
          
        } else {
          F_0z <- psnp(cc, phi, base.dens, a)#, roundoff, digits) 
        }# end of if (is.null(cc)) else ...
        
      }# end of if (is.null(F_0z))
      
      logF_K    <- log(F_0z) 
      
    } # End of ifelse      
        
    # c(logF_K)
    logF_K
    
  }, error=function(ee) {print("From function: logF_K.aft"); print(ee); browser()})
  
} # End of logF_K.aft

#===============================================================================
# Old name: reg.sur.est
S_K.aft <- function(tt, theta, base.dens=c("exp","stdnorm"), bet=NULL, x=NULL, 
                    bx=NULL, e=NULL, S_0tt=NULL, S_0z=NULL, cc=NULL,
                    a=NULL) {#, roundoff=F, digits=10) {
  # The estimated survival function of t with respect to the base density specified in
  # base.dens for any of the following regression models: PH, AFT and PO
  
  # Note that only time independent covariates are dealt with
  #------------------------------
  # Caller:
  
  # - h_K, snp.logLik, snp.sur.CI
  
  #------------------------------
  # Input:
  
  # - tt          the time
  
  # - theta       the parameter vectors containing
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - bet         the vector of regression coefficients (beta)
  
  # - x           the design matrix where each row is the covariate vector of the
  #               corresponding individuals

  #------------------------------
  # Output:
  
  # - The value of the estimated survival at time t (>0) given the regression 
  # coefficients and predictors
  
  # To improve speed performance of snp.logLik, this function will take the 
  # following extra arguments which should be computed inside snp.logLik. Otherwise
  # , they will be recomputed inside this function.
  
  # - bx          bet %*% x
  
  # - e           exp(bx)
  
  # - S_0z        S_0(tt * (1/e), theta, base.dens) (only needed when mod is 
  #               "AFT"
  
  # - S_0tt       S_0(tt, theta, base.dens) (only needed when mod is "PH"
  #               or "PO"
  
  # - cc          only needed for AFT mod. When base.dens is stdnorm, cc is 
  #               log(tt*exp(-XB)-mu)/sig = log(tt-XB-mu)/sig,
  #               and when base.dens is exp, z is 
  #               (tt*exp(-XB)/exp(mu))^(1/sig) = (tt*/exp(mu+XB))^(1/sig)
  
  # - a           the polynomial coefficients

  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Oct 11 2010
  # Last modified: Jun 04 2015
  
  tryCatch({ # To catch unprecedented errors
    
    S_K          <- rep(NaN, length(tt))
    
    if ((length(bet)==0 | is.null(x) | prod(dim(x))==0) & (is.null(bx) & is.null(e) & is.null(S_0z))) {
      if (is.null(S_0tt)) S_0tt <- S_0(tt, theta, base.dens, a)#, roundoff, digits)
      
      S_K <- S_0tt
      
    } else { 
      # Accelerated failure time model, only need to mod here 1214 31 Jan 2013
      if (is.null(S_0z)) {
        
        l   <- length(theta)
        
        # In case the polynomial is degenerate, or theta contains only mu and sigma
        phi <- NULL
        
        if (length(theta)!=2) {
          phi <- theta[1:(l-2)]           
        }
        
        if (is.null(cc)) {
          
          # if (is.null(bx)) bx       <- x[id,,drop=F] %*% bet           
          if (is.null(bx)) bx <- x %*% bet           
          if (is.null(e)) e   <- exp(bx)            
          
          mu  <- theta[l-1]
          sig <- theta[l]
          # if (sig < 0) stop("sigma is negative!")
          
          if (base.dens == "stdnorm") {
            if (is.null(cc)) {
              cc   <- (log(tt)-bx-mu)/sig
              S_0z <- sur.snp(cc, phi, base.dens, a)#, roundoff, digits) 
            } 
          } else {
            if (is.null(cc)) {
              cc   <- (tt/exp(mu+bx))^(1/sig) # taking log then exp here doesn't help
              S_0z <- sur.snp(cc, phi, base.dens, a)#, roundoff, digits) 
            } 
          } # end of if (base.dens=="stdnorm") else..
          
        } else {
          S_0z <- sur.snp(cc, phi, base.dens, a)#, roundoff, digits) 
        }# end of if (is.null(cc)) else ...
        
      }# end of if (is.null(S_0z))
      
      S_K    <- S_0z
      
    } # End of ifelse
      
    # c(S_K)
    S_K
    
  }, error=function(ee) {print("From function: S_K.aft"); print(ee); browser()})
  
} # End of S_K.aft

#===============================================================================

F_K.aft <- function(tt, theta, base.dens=c("exp","stdnorm"), bet=NULL, x=NULL, 
                    bx=NULL, e=NULL, F_0tt=NULL, F_0z=NULL, cc=NULL,
                    a=NULL) {#, roundoff=F, digits=10) {
  # The estimated cdf of t with respect to the base density specified in
  # base.dens for any of the following regression models: PH, AFT and PO
  
  # Note that only time independent covariates are dealt with
  #------------------------------
  # Caller:
  
  # - h_K, snp.logLik, snp.sur.CI
  
  #------------------------------
  # Input:
  
  # - tt          the time
  
  # - theta       the parameter vectors containing
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - bet         the vector of regression coefficients (beta)
  
  # - x           the design matrix where each row is the covariate vector of the
  #               corresponding individuals

  #------------------------------
  # Output:
  
  # - The value of the estimated cdf at time t (>0) given the regression 
  # coefficients and predictors
  
  # To improve speed performance of snp.logLik, this function will take the 
  # following extra arguments which should be computed inside snp.logLik. Otherwise
  # , they will be recomputed inside this function.
  
  # - bx          bet %*% x
  
  # - e           exp(bx)
  
  # - F_0z        F_0(tt * (1/e), theta, base.dens) (only needed when mod is 
  #               "AFT"
  
  # - F_0tt       F_0(tt, theta, base.dens) (only needed when mod is "PH"
  #               or "PO"
  
  # - cc          only needed for AFT mod. When base.dens is stdnorm, cc is 
  #               log(tt*exp(-XB)-mu)/sig = log(tt-XB-mu)/sig,
  #               and when base.dens is exp, z is 
  #               (tt*exp(-XB)/exp(mu))^(1/sig) = (tt*/exp(mu+XB))^(1/sig)
  
  # - a           the polynomial coefficients

  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Feb 19 2013
  # Last modified: Jun 04 2015
  
  tryCatch({ # To catch unprecedented errors
    # browser()
    F_K          <- rep(NaN, length(tt))
    
    if ((length(bet)==0 | is.null(x) | prod(dim(x))==0) & (is.null(bx) & is.null(e) & is.null(F_0z))) {
      if (is.null(F_0tt)) F_0tt <- F_0(tt, theta, base.dens, a)#, roundoff, digits)
      
      F_K <- F_0tt
      
    } else { 
      
      if (is.null(F_0z)) {
        
        l   <- length(theta)
        
        # In case the polynomial is degenerate, or theta contains only mu and sigma
        phi <- NULL
        
        if (length(theta)!=2) {
          phi <- theta[1:(l-2)]           
        }
        
        if (is.null(cc)) {
          
          # if (is.null(bx)) bx       <- x[id,,drop=F] %*% bet 
          if (is.null(bx)) bx <- x%*% bet 
          if (is.null(e)) e   <- exp(bx)            
          
          mu  <- theta[l-1]
          sig <- theta[l]
          # if (sig < 0) stop("sigma is negative!")
          
          if (base.dens == "stdnorm") {
            if (is.null(cc)) {
              cc   <- (log(tt)-bx-mu)/sig
              F_0z <- psnp(cc, phi, base.dens, a)#, roundoff, digits) 
            } 
          } else {
            if (is.null(cc)) {
              cc   <- (tt/exp(mu+bx))^(1/sig) # taking log then exp here doesn't help
              # if (any(is.na(cc))) browser()
              F_0z <- psnp(cc, phi, base.dens, a)#, roundoff, digits) 
            } 
          } # end of if (base.dens=="stdnorm") else..
          
        } else {
          F_0z <- psnp(cc, phi, base.dens, a)#, roundoff, digits) 
        }# end of if (is.null(cc)) else ...
        
      }# end of if (is.null(F_0z))
      
      F_K    <- F_0z            
      
    } # End of ifelse      
          
    # c(F_K)
    F_K
    
  }, error=function(ee) {print("From function: F_K.aft"); print(ee); browser()})
  
} # End of F_K.aft

#===============================================================================

f_K.aft <- function(tt, theta, base.dens=c("exp","stdnorm"), bet=NULL, x=NULL, cc=NULL, a=NULL) {
  # The estimated density function of t with respect to the base density specified
  # in base.dens for any of the following regression models: PH, AFT and PO
  
  # Note that only time independent covariates are dealt with
  #------------------------------
  # Caller:
  
  # - h_K, snp.logLik 
  
  #------------------------------
  # Input:
  
  # - t           the time
  
  # - theta       the parameter vectors containing
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - bet         the vector of regression coefficients (beta)
  
  # - x           the design matrix where each row is the covariate vector of the
  #               corresponding individuals
  
  # - cc          only needed for AFT mod. When base.dens is stdnorm, cc is 
  #               log(tt*exp(-XB)-mu)/sig = log(tt-XB-mu)/sig,
  #               and when base.dens is exp, z is 
  #               (tt*exp(-XB)/exp(mu))^(1/sig) = (tt*/exp(mu+XB))^(1/sig)
  
  # - a           the polynomial coefficients
  
  #------------------------------
  # Output:
  
  # - The value of the estimated density at time t (>0) given the regression 
  # coefficients and predictors
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Oct 11 2010
  # Last modified: Sep 10 2013
  
  return(c(exp(logf_K.aft(tt, theta, base.dens, bet, x, cc=cc, a=a))))  
  
}

#===============================================================================
# Old name: reg.haz.est  
h_K.aft <- function(tt, theta, base.dens=c("exp","stdnorm"), bet=NULL, x=NULL) {#, roundoff=F, digits=10)  {
  # The estimated survival function of t with respect to the base density specified
  # in base.dens for any of the following regression models: PH, AFT and PO
  
  #------------------------------
  # Caller:
  
  # - snp.haz.CI
  
  #------------------------------
  # Input:
  
  # - t           the time
  
  # - theta       the parameter vectors containing
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - bet         the vector of regression coefficients (beta)
  
  # - x           the design matrix where each row is the covariate vector of the
  #               corresponding individuals

  #------------------------------
  # Output:
  
  # - The value of the estimated survival function of snp at t (>0) given the 
  # regression coefficients and predictors
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Oct 11 2010
  # Last modified: Jun 04 2015  
  
  h_K <- f_K.aft(tt, theta, base.dens, bet, x)  / 
         S_K.aft(tt, theta, base.dens, bet, x)#, roundoff=roundoff, digits=digits)
  
  # c(h_K)
  h_K
  
} # end of h_K.aft

#===============================================================================
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# DIFFERENCES FROM THE ORIGINAL PAPER:

# - instead of calculating the survival function we evaluate the CDF directly
# for any x in the real line.

# - the iteration works as follows (a bit differnt from those of the paper):
# + I(0,x) = pnorm(x)
# + I(1,x) = - dnorm(x)
# + I(k,x) = - x^(k-1) * dnorm(x) + (k-1) * I(k-2,x) for k >=2

#===============================================================================
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# REMAINING ISSUES:

# - Don't know how to use switch case statement  (UNCURED !!!)
#-------------------------------------------------------------------------------

# - Abnormality appears whenever phi contains more than 22 elements !!!
#-------------------------------------------------------------------------------

# - In snp.opt.init need a more sophisticated way to calculate the sample mean
# and variance based on T or logT. Marcel suggested that a simple model could be
# use such as exponential -> E(T) = 1/lambda , lambda = number of events / total
# living time. Then use delta method to get E(logT)=log(ET)
#-------------------------------------------------------------------------------

# - At t=0, both the estimated survival and hazard function as well as their
# corresponding confidence intervals, standard errors are NaN !!! This is
# because of the nature of the snp (taking log {0})
#-------------------------------------------------------------------------------

# - THE METHOD MENTIONED IN THE PAPER ISN'T SAFE !!!
#-------------------------------------------------------------------------------

# Oct 15th 2010
# - Can't find a  better way to solve the issue in dsnp
# Also ometimes have non-finite finite-difference value [4]. This shouldn't have 
# anything to do with the dimension issue above. This is because of R's optim
# (DONE) (SEE DISCUSSION BELOW)
#-------------------------------------------------------------------------------

# Oct 15th 2010
# - In snp.survreg, receive message: 
# Warning message:
# In sqrt(Var) : NaNs produced
# This may be due to taking sqrt of negative number, the function where the 
# really takes place is var.delta.rule 

# - The answer is that the optimization didn't numerically converge (DONE)
#-------------------------------------------------------------------------------

# Oct 15th 2010
# - In var.delta.rule
# Sometimes when under AFT model and PH model we see
# Error in solve.default(par.cov) : 
# system is computationally singular: reciprocal condition number = 2.29577e-18
# Theoretically, we're trying to find a minimal -loglik hence maximal loglik. 
# Therefore, the hessian (information matrix) should be positive definite hence
# invertible. This could be numerical error or even worse i.e something is wrong
# with the log-likelihood. Or that could be due to the error below (DONE)
#-------------------------------------------------------------------------------

# Oct 21st 2010
# - During the optimization process, the snp.logLik function may return NaN this 
# could be due to one of the following case: Inf / Inf, 0/0 , 0 * Inf or similar
# cases where Inf is replaced by -Inf. The first attempt to go around this is to
# use pmax, pmin and the constants .Machine$double.eps and .Machine$double.xmax
# howver, this somehow geoperdises the optimization results. The underlying 
# reasons for those abnormalities my be when the input times to events are large
# and after a lot of exponentiating (especially in PH model), we end up with Inf

# - Marcel thinks that scaling the covariates (x) by : (x - mean(x))/var(x) and
# also the times to events by : t/t.max can solve the problem (DONE) 
# (TEMPORARILY)
#-------------------------------------------------------------------------------

# Dec 23rd 2010
# - The result given by snp.logLik doesn't adhere to that given by Zhang's code
# (DONE), (SEE BELOW)
#-------------------------------------------------------------------------------

# Dec 28th 2010
# - The log-likelihood given by this implementation is systematically diffrent 
# from that produced by Zhang's code. I suspect the deviance is caused by the 
# function iter.ints as compared to, for example 
# tart N_2_ph(b) global(pi,x,v,invbn2,d,n). Specifically, we still don't know 
# the reason for adding the term ez (stdnorm) or elogz (exp)

# - The answer is that this might be due to a different in the parametrization
# , we still don't know for sure if this is neccessary (CHECKED)
#-------------------------------------------------------------------------------

# Dec 31st 2010
# - Can't deal with factor design (SOLVED)

# - The regression don't give consistent estimate when we use exponential  base
# even when the real data were generated from a snp with either normal or 
# exponential base (SOLVED) (SEE BELOW)

# - The above issue is due to the fact that the case of exponential base was 
# misimplemented. Actually, it's z* = exp(z) that has a snp density with a 
# standard exponential density so that z itself has an extreme distribution!

# - Also when one of the covariates is categorical, the regression doesn't work
# very well.  (UNSOLVED)

# - Even when all covariates are normally distributed, and when the number of 
# them is 3, the estimates of mu and  the last coefficient are very inconsistent
# given that the real value is small. In general it's not robust again ZERO 
# coefficient, eventhough when the real value is 0, the estimate is say 0.001!
#-------------------------------------------------------------------------------

# Jan 18th 2011
# - The QQ plots of estimated vs real model under PH and PO also AFT usually
# looks not so great at the upper tail. However, this could be due to the fact 
# that the SNP is not unimodal. Indeed, MC carried out a simulation where 
# data generated from a mixture distribution which is bimodal and even real
# vs real gives the same phenomenon (TEMPORARILY SOVED)

# - Also try to use R's coxph to fit the PH data generated from sim.ph and see
# if even R can estimate it well. Go to R task view - survival analsyis website
# or see timereg.pdf for more details
#-------------------------------------------------------------------------------
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# ON-GOING WORKS:

# - Develop hazard rate function (ain't found a simpler formula, so use the ratio
# h(t) = f(t)/S(t) DONE !!
#-------------------------------------------------------------------------------

# - Return the ACI and HQ and similar information (need knowledge on those)
#   + AIC = 2*(-l(theta) + p)
#   + BIC = -2*l(theta) +plog(n)
#   + HQ  = 2*(-l(theta) + plog{log(n)})
#-------------------------------- THIS WAS DONE --------------------------------

# - Accomodate also :
#   + left-truncated observations - f (x)/ S (Y_L)
# This was DONE via modifying snp.logLik.RightC and snp.logLik.IntC to account 
# for left truncated data if any
#-------------------------------------------------------------------------------

# - If possible, include every of the following cases
#   + exact lifetimes     = f(x)                          DONE
#   + right censored      = S(C_r)                        DONE
#   + left censored       = 1 - S(C_l)
#   + interval censored   = S(C_l)-S(C_r)                 DONE
#   + left truncated      = f(x) / S(Y_l)                 DONE
#   + right truncated     = f(x) / [1 - S^Y_r]
#   + interval truncated  = f(x) / [S(Y_l) - S(Y_r)]
#-------------------------------------------------------------------------------

# - Implementation of Censored Data Regression Analysis based on SNP for :
#   + PH model
#   + AFT model
#   + PO model
# This was DONE via creating NEW S_K, f_K and h_K .
#-------------------------------------------------------------------------------

# - Marcel thinks that we should also get the starting values for beta by using 
# regression functions in R's survival packages. => Read Rodriguez's materials &
# modifying snp.opt  DONE
#-------------------------------------------------------------------------------

# - Need to create imaginary data to see how the outputs of R's look like for 
#all a possible cases.  DONE
#-------------------------------------------------------------------------------

# - Need to OOP the code UNDONE (No need for now) 
#-------------------------------------------------------------------------------

# Dec 23rd 2010
# - Allow user to set initial values for sigma and mu. DONE
#-------------------------------------------------------------------------------

# Dec 23rd 2010
# - Deal with case where k (degree of polynomial) = 0. DONE
#-------------------------------------------------------------------------------

# Dec 23rd 2010
# - Report all the best estimates corresponding to each initial values. DONE
#-------------------------------------------------------------------------------

# Jan 04th 2011
# - Allow user to put in a vector called baseline whose elements are the values
# corresponding to which the baseline survival time is i.e each column in the 
# design matrix will be centered (shifted) w.r.t each element in this vector. 
# If this is left null, the df values are the column-wise means for all entries
# (UNDONE)

# - Besides, each column will also be scaled by divided by their empirical 
# standard deviations for continuous covariates and by 2*sd for binary ones.
# For categorical covariates with more than 2 categories, the answer is still 
# open !. (UNDONE)
#-------------------------------------------------------------------------------

# Oct 10th 2011
# Sometimes function my.log receives x as NAs. One of the reason is that in 
# function logdsnp, particularly in the chunk as.fucntion...., when x has Inf, 
# the evaluated value will be NA !!!
#-------------------------------------------------------------------------------

# Curretnly working on ...
#-------------------------------------------------------------------------------

#===============================================================================
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

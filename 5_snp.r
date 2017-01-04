
# VERSION 1.0 - ONLY SNP. ALL IDIOT PROOFING REMOVED TO SAVE TIME

# ANH NGUYEN DUC - OXFORD UNIVERSITY CLINICAL RESEARCH UNIT

# SUPERVISED BY:  MARCEL WOLBERS
# CREATED:        May      31  2013
# LAST MODIFIED:  Jun      02  2015
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# FUNCTIONS INCLUDED:
# moment                              (tested 26 Jul 2012)
# polycofsnp                          (tested 27 Jul 2012)
# logf_0                              (tested 27 Jul 2012)
# my.log                              (tested 26 Jul 2012)
# logdsnp                             (tested 27 Jul 2012)
# logdbase                            (tested 27 Jul 2012)
# S_0                                 (tested 27 Jul 2012)
# F_0
# sur.snp                             (tested 27 Jul 2012)
# psnp                                (tested 27 Jul 2012)
# iter.ints                           (tested 27 Jul 2012) (blocked)
# iter.ints1                           
# snp.sur.CI
# snp.haz.CI
# var.delta.rule                      (tested 31 Jul 2012)  
# dbase                               (tested 27 Jul 2012)
# dsnp                                (tested 27 Jul 2012)
# f_0                                 (tested 27 Jul 2012)
# h_0
# rsurv.snp                           (tested 01 Aug 2012)
# rsnp                                (tested 31 Jul 2012)
# sub.sim.stdnorm.snp                 (tested 31 Jul 2012)
# sub.sim.stdnorm.snp.old             (tested 31 Jul 2012)
# sub.sim.exp.snp    
# sub.sim.exp.snp.old                 (tested 31 Jul 2012)

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# This must be kept at lest here or in 4_snp.aft.r to pass the required 
# libraries into the global parallel cluster
require(polynom)
require(numDeriv)
require(survival)
require(Rlab)

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



#===============================================================================

moment <- function(base.dens=NULL, x=NULL) { # still used by snp.meanvar
  # Output the moment sequence of the orders contained in x based on the density
  # given by base.dens
  # This function mainly serves DSNP and PSNP
  
  #------------------------------
  # Input:
  
  # - base.dens   name of the density function
  
  # - x           vector storing the moment-orders
  
  #------------------------------
  # Output:
  
  # - A sequence of moments
  
  #------------------------------
  # ref 
  # http://en.wikipedia.org/wiki/Normal_distribution#Moments
  # http://en.wikipedia.org/wiki/Exponential_distribution
  
  # Last modified: 14 Jan 2013
  
  tryCatch({ # To catch unprecedented errors    
    # Idiot proofing    

    if (base.dens=="stdnorm") { # This is hopefully faster but can only support
      # up to 10
      
      # moment.Seq <- (1 - abs(2 * round(x/2) - x)) * factorial(x) / 
      #   (2^(x/2) * factorial(x/2))
      # This is syntactically correct !
      moment.Seq <- c( 1, 0, 1, 0, 3, 0, 15, 0, 105, 0, 945)[x+1]      
      
    } else if (base.dens=="exp") {
      # moment.Seq <- factorial(x)
      moment.Seq <- c(1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800)[x+1]
      
    } else {
      stop("\n base.dens must be either exp or stdnorm!")
    }     
    
  }, error=function(ee) {print("From function: moment"); ee; browser()})
  
} # End of moment

#===============================================================================

polycofsnp <- function(base.dens=NULL, phi=NULL) {
  # Name:
  
  # - Polynomial Coefficients of SNP
  
  #------------------------------
  # Output the vector of coefficents a s/t (a_0 + ... + a_k * z^ k) ^2 = 1, where
  # z has the density function base.dens
  # This function mainly serves DSNP and PSNP
  
  #------------------------------
  # Input:
  
  # - base.dens   name of the density function
  
  # - phi         the coordinates of a in spherical system
  
  #------------------------------
  # Output:
  
  # - A sequence of coefficients
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  #------------------------------
  # Created:        Oct 11 2010
  # Last modified:  Jul 15 2013
  
  tryCatch({ # To catch unprecedented errors
    
    # Get the length of phi
    k          <- length(phi)
    
    # Step 1: Calculate co from phi using spherical transformation ---------------
    cos.phi    <- cos(phi) # my.cos(phi)
    sin.phi    <- sin(phi)
    
    co         <- cumprod(c(1, cos.phi))
    co[1L:k]    <- co[1L:k] * sin.phi

    # Step 2: Create the moment matrix A based on the base density ---------------
    # tmp.Seq    <- seq(from = 0, to = 2*k, by = 1)
    
    # # Create a moment sequence from 0 to 2 times k
    # moment.Seq <- moment(base.dens, tmp.Seq)
    
    # # Create the moment matrix based on the above moment sequence
    # A          <- outer(1:(k+1), 1:(k+1), function(x,y) moment.Seq[x+y-1])
    
    # if (base.dens=="stdnorm") {
    #   A <- rbind(c(1, 0, 1, 0, 3, 0),
    #              c(0, 1, 0, 3, 0, 15),
    #              c(1, 0, 3, 0, 15, 0),
    #              c(0, 3, 0, 15, 0, 105),
    #              c(3, 0, 15, 0, 105, 0),
    #              c(0, 15, 0, 105, 0, 945))
    # } else {
    #   A <- rbind(c(1, 1, 2, 6, 24, 120),
    #              c(1, 2, 6, 24, 120, 720),
    #              c(2, 6, 24, 120, 720, 5040),
    #              c(6, 24, 120, 720, 5040, 40320),
    #              c(24, 120, 720, 5040, 40320, 362880),
    #              c(120, 720, 5040, 40320, 362880, 3628800))
    # } # end of if (base.dens=="stdnorm") else ...
    
    # A <- A[1:(k+1), 1:(k+1)]
    
    # # Step 3: Decompose A using Cholesky's method --------------------------------
    # B          <- chol(A)
    # Step 4: Create a the polynomial's coefficenct vector -----------------------
    # a          <- chol2inv(B) %*% (t(B) %*% co) # checked, correct, (...) is for time reduction
    #     } else {
    #       a          <- c(1, rep(0, length(phi)))
    #       
    #     } # end of if (round(phi...)) else...
    
    # # Using Rcpp mychol is a bit faster
    # Lc <- mychol(A)    
    # a <- Lc$Inv %*% (Lc$L %*% co)    
    
    # A faster and dirtier option is, besides hard coding A, to first declare A as a global variable and
    # then precalculate Lc using mychol for all degrees available in A. This should be done in 1_snpcr.r
    
    if (base.dens=="stdnorm") {
      
      # a <- Lc.n[[k+1L]]$InvXL %*% co # For some reason this leads to results different from old ones
      a <- Lc.n[[k+1L]]$Inv %*% (Lc.n[[k+1L]]$L %*% co)
      
    } else {
      
      # a <- Lc.e[[k+1L]]$InvXL %*% co
      a <- Lc.e[[k+1L]]$Inv %*% (Lc.e[[k+1L]]$L %*% co)
      
    } # end of if (base.dens=="stdnorm") else ...
    
    a
    
  }, error=function(ee) {print("From function: polycofsnp"); ee; browser()})
  
} # End of polycofsnp

#===============================================================================

logf_0 <- function(tt, theta=NULL, base.dens=NULL, a=NULL) {
  # The estimated log density function of t with respect to the base density 
  # specified in base.dens
  
  #------------------------------
  # Caller:
  
  # - snp.logLik, logf_K
  
  #------------------------------
  # Input:
  
  # - tt        the time
  
  # - theta     the parameter vectors containing
  
  # + phi       the spherical coordinates
  
  # + mu        MU
  
  # + sig       Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - a         the polynomial coefficients
  
  #------------------------------
  # Output:
  
  # - The value of the estimated LOG density at time t (>0)
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Dec 30 2010
  # Last modified: Sep 10 2013
  
  tryCatch({ # To catch unprecedented errors    
    
    l   <- length(theta)
    
    # In case the polynomial is degenerate, or theta contains only mu and sigma
    phi <- NULL
    
    if (length(theta)>2) {
      phi <- theta[1L:(l-2L)]     
    }
    
    mu  <- theta[l-1L]
    sig <- theta[l]
    # if (sig < 0) stop("sigma is negative!")
    
    # Change of variable for each base-line case
    if (base.dens == "stdnorm") {
      x    <- (log(tt)- mu) / sig    
      cons <- -log(tt) - log(sig)     
      
    } else {

      x    <- ( tt/exp(mu) ) ^ (1/sig)

      cons <- -mu/sig - log(sig) + (1/sig - 1)*log(tt)           
      
    }

    logf_0 <- cons + logdsnp(x, phi, base.dens, a)    
    
    # c(logf_0)
    
  }, error=function(ee) {print("From function:logf_0"); print(ee); browser()})
  
  
} # End of logf_0

#===============================================================================

my.log <- function(x) {
  # To deal with situations in which x is 0 or negative (small in magnitude) due to
  # numerical error
  
  #------------------------------
  # Caller:
  
  # - snp.logLik
  
  #------------------------------
  # Created:       Mar 19th 2010
  # Last modified: Feb 02nd 2012
  
  tryCatch({ # To catch unprecedented errors
    
    # Idiot proofing
    # if (is.nanull(x)) stop("\nx must be specified!")
    if (length(x) > 0) {
      re <- log(pmax(x,.Machine$double.eps))  
    } else {
      re <- NA
    }
    
  }, error=function(ee) {print("From function: my.log"); ee; browser()})
  
} # End of my.log

#===============================================================================

logdsnp <- function(x, phi=NULL, base.dens=NULL, a=NULL) {
  # The log density function of snp with respect to the base density specified in
  # base.dens
  
  #------------------------------
  # Caller:
  
  # - logf_0
  
  #------------------------------
  # Input:
  
  # - x     the data
  
  # - phi   the parameter vectors
  
  # - a     the polynomial coefficients
  
  #------------------------------
  # Output:
  
  # - The value of the log density snp at x
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  #------------------------------
  # Created:       Dec 30 2010
  # Last modified: Jul 04 2013
  
  tryCatch({ # To catch unprecedented errors
    
    # Get the values at the base density
    logd   <- logdbase(base.dens, x)
    
    # In case the polynomial isn't degenerate i.e phi isn't NULL
    if (!is.null(phi) & length(phi)!=0) {
      
      # Step 1 , 2 ,3 ref polycofsnp
      
      # Step 4: Create the polynomial's coefficient vector -----------------------
      if (is.null(a)) a <- polycofsnp(base.dens, phi) # Note a is a COLUMN vector !!!
      
      
      # Step 5: Get the density value
      
      if (length(x)==1) {
        logd <- 2 * log(abs(sum(sapply(0:(length(a)-1), function(i) a[i+1] * x^i)))) + logd
      } else {

        # logd <- 2 * log(abs(apply(sapply(0:(length(a)-1), function(i) a[i+1] * x^i), 1, sum))) + logd
        # Below is 10% faster
        tmp <- matrix(NA, ncol = length(a), nrow = length(x))
        
        for (i in 0L:(length(a)-1L)) {
          tmp[,i+1L] <- a[i+1L, 1L] * x^i
        }
        
        logd <- 2 * log(abs(apply(tmp, 1L, sum))) + logd
      }
      
    }
    
    #c(logd)    
    logd
    
  }, error=function(ee) {print("From function: logdsnp"); ee; browser()})
  
} # End of logdsnp

#===============================================================================

logdbase <- function(base.dens=NULL, x) {
  # Output the log density values at inputs contained in x based on the density
  # given by base.dens
  # This function mainly serves LOGDSNP
  
  #------------------------------
  # Input:
  
  # - base.dens   name of the density function
  
  # - x           vector storing the intputs
  
  #------------------------------
  # Output:
  
  # - A sequence of log density value
  
  #------------------------------
  # Created:        Dec 30th 2010
  # Last modified:  Sep 10th 2013
  
  tryCatch({ # To catch unprecedented errors      
    
    if (base.dens=="stdnorm") {
      d <- dnorm(x=x, mean=0, sd=1, log=T)
    } else if (base.dens=="exp") {
      d <- -x
    } 
    
  }, error=function(ee) {print("From function: logdbase"); ee; browser()})
  
} # End of logdbase

#===============================================================================

S_0 <- function(tt, theta, base.dens=NULL, a=NULL) {#, roundoff=F, digits=10) { #, bet, x)  {
  # The estimated survival function of t with respect to the base density specified in
  # base.dens
  
  #------------------------------
  # Caller:
  
  # - h_0, S_K, f_K, snp.logLik 
  
  #------------------------------
  # Input:
  
  # - tt        the time
  
  # - theta     the parameter vectors containing
  
  # + phi       the spherical coordinates
  
  # + mu        MU
  
  # + sig       Sigma
  
  # - base.dens either "exp" or "stdnorm" specifying the base density
  
  # - a         the polynomial coefficients

  #------------------------------
  # Output:
  
  # - The value of the estimated survival function of snp at t (>0)
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Mar 18 2010
  # Last modified: Jun 04 2015
  
  tryCatch({ # To catch unprecedented errors     
    
    # Handle extreme cases (t=0 and Inf) first  
    S_0         <- rep(NaN, length(tt))    
    id          <- (tt!=Inf & tt!=0)
    S_0[tt==Inf]<- 0
    S_0[tt==0]  <- 1
    l   <- length(theta)
    
    # In case the polynomial is degenerate, or theta contains only mu and sigma
    phi <- NULL
    
    if (length(theta)!=2) {
      phi <- theta[1L:(l-2L)]           
    }
    
    mu  <- theta[l-1L]
    sig <- theta[l]
    # if (sig < 0) stop("sigma is negative!")
    
    if (base.dens == "stdnorm") {
      x <- (log(tt[id]) - mu) / sig
    } else {
      x <- ( tt[id]/exp(mu) ) ^ (1/sig) # taking log then exp here doesn't help
      
    }
    
    S_0[id] <- sur.snp(x, phi, base.dens, a)#, roundoff, digits)          
    
    # c(S_0) 
    S_0
    
    # browser()
  }, error=function(ee) {print("From function: S_0"); print(ee); browser()})
  
} # End of S_0

#===============================================================================

F_0 <- function(tt, theta, base.dens=NULL, a=NULL) {#, roundoff=F, digits=10) { 
  # The estimated cdf of t with respect to the base density specified in
  # base.dens
  
  #------------------------------
  # Caller:
  
  # - F_K, snp.logLik 
  
  #------------------------------
  # Input:
  
  # - tt        the time
  
  # - theta     the parameter vectors containing
  
  # + phi       the spherical coordinates
  
  # + mu        MU
  
  # + sig       Sigma
  
  # - base.dens either "exp" or "stdnorm" specifying the base density
  
  # - a         the polynomial coefficients

  #------------------------------
  # Output:
  
  # - The value of the estimated cdf of snp at t (>0)
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Feb 19 2013
  # Last modified: Jun 04 2015
  
  tryCatch({ # To catch unprecedented errors    

    F_0         <- rep(NaN, length(tt))
    id          <- (tt!=Inf & tt!=0)
    F_0[tt==Inf]<- 1
    F_0[tt==0]  <- 0
    l   <- length(theta)
    
    # In case the polynomial is degenerate, or theta contains only mu and sigma
    phi <- NULL
    
    if (length(theta)!=2) {
      phi <- theta[1L:(l-2L)]           
    }
    
    mu  <- theta[l-1L]
    sig <- theta[l]
    # if (sig < 0) stop("sigma is negative!")
    
    if (base.dens == "stdnorm") {
      x <- (log(tt[id]) - mu) / sig
    } else {
      x <- ( tt[id]/exp(mu) ) ^ (1/sig) # taking log then exp here doesn't help      
    }
    
    F_0[id] <- psnp(x, phi, base.dens, a)#, roundoff, digits)      

    # c(F_0)     
    F_0
    
  }, error=function(ee) {print("From function: F_0"); print(ee); browser()})
  
} # End of S_0

#===============================================================================

psnp <- function(x, phi, base.dens=NULL, a=NULL) {#, roundoff=F, digits=10) {
  # The cumulative distribution function of snp with respect to the base density specified in base.dens
  
  #------------------------------
  # Caller:
  
  # - S_0
  
  #------------------------------
  # Input:
  
  # - x           the data
  
  # - phi         the parameter vectors
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density  
  
  # - a           the polynomial coefficients
  
  #------------------------------
  # Output:
  
  # - The value of the survival function of snp at x
  
  #-------------------------------
  # Note:
  
  # - For the exponential case, x must be nonnegative!!!
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  #------------------------------
  # Last modified: Jun 04 2015
  
  # In case the polynomial isn't degenerate i.e phi isn't NULL
  tryCatch({
  if (!is.nanull(phi) & length(phi)!=0) {
    cdf <- rep(NaN, length(x))
    
    if (base.dens=="stdnorm") {
      id  <- (x!=Inf & x!=-Inf)
      cdf[x==Inf] <- 1
      cdf[x==-Inf]<- 0
    } else if (base.dens=="exp") {
      id  <- (x!=Inf & x!=0)
      cdf[x==Inf] <- 1
      cdf[x==0]   <- 0
    }
    
    # Get the length of phi
    k <- length(phi)
    
    # Step 1: Create the coefficenct vector of the polynomial inside the square --
    if (is.null(a)) a <- polycofsnp(base.dens, phi)  # Note a is a (k+1)x1 matrix
    
    # Step 2: Create the coefficient vector of the quadratic form ----------------
    d <- coef(polynomial(a)^2)
        
    # Step 3: Evaluate the values of each term-integral and put 'em in a vector -- 
    Ints <- iter.ints(x[id], base.dens, k)
    
    # Step 4: Evaluate the final cdf
    #cdf <- t(d) %*% Ints
    
    # The line below protects our code from the case where the last elements in 
    # a are 0's. But it will increase computation time dramatically and the 
    # final result is amazingly the same !!    
    if (length(d) < dim(Ints)[2L]) d <- c(d, rep(0, dim(Ints)[2L]-length(d)))
    cdf[id] <- Ints %*% d
    
    # If we want this function to serve purposes other than survival analysis, we
    # shouldn't uncomment the following lines
    
  } else {
    if (base.dens=="stdnorm") {
      cdf <- pnorm(x)
      
    } else if (base.dens=="exp") {
      cdf <- pexp(x)
      
    }
    
  } # end of if (!is.nanull(phi) & length(phi)!=0) else...
  
  # if (roundoff) cdf <- round(cdf, digits=digits)
  # c(cdf)
  cdf
  
  }, error=function(ee) {print("From function: psnp"); print(ee); browser()})
} # End of psnp

#===============================================================================

sur.snp <- function(x, phi=NULL, base.dens=NULL, a=NULL) {#, roundoff=F, digits=10) {
  # The survival function of snp with respect to the base density
  # specified in base.dens
  
  #------------------------------
  # Caller:
  
  # - sur.snp
  
  #------------------------------
  # Input:
  
  # - x           the data
  
  # - phi         the parameter vectors
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - a           the polynomial coefficients

  #------------------------------
  # Output:
  
  # - The value of the cdf of snp at x
  
  #-------------------------------
  # Note:
  
  # - For the exponential case, x must be nonnegative!!!
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  #------------------------------
  # Last modified: Jun 04 2014
  
  tryCatch({ # To catch unprecedented errors
    
    # In case the polynomial isn't degenerate i.e phi isn't NULL
    if (!is.nanull(phi) & length(phi)!=0) {
      sur <- rep(NaN, length(x))

      if (base.dens=="stdnorm") {
        id  <- (x!=Inf & x!=-Inf)
        sur[x==Inf] <- 0
        sur[x==-Inf]<- 1
      } else if (base.dens=="exp") {
        id  <- (x!=Inf & x!=0)
        sur[x==Inf] <- 0
        sur[x==0]   <- 1
      }
      
      # Get the length of phi
      k <- length(phi)
      
      # Step 1: Create the coefficenct vector of the polynomial inside the square --
      if (is.null(a)) a <- polycofsnp(base.dens, phi)  # Note a is a (k+1)x1 matrix
      
      # Step 2: Create the coefficient vector of the quadratic form ----------------
      d <- coef(polynomial(a)^2)
      
      
      # Step 3: Evaluate the values of each term-integral and put 'em in a vector -- 
      Ints <- iter.ints1(x[id], base.dens, k)
      
      # Step 4: Evaluate the final cdf
      #cdf <- t(d) %*% Ints
      
      # The line below protects our code from the case where the last elements in 
      # a are 0's. But it will increase computation time dramatically and the 
      # final result is amazingly the same !!    
      if (length(d) < dim(Ints)[2L]) d <- c(d, rep(0, dim(Ints)[2L]-length(d)))
      sur[id]  <- Ints %*% d
      
      # If we want this function to serve purposes other than survival analysis, we
      # shouldn't uncomment the following lines
      
    } else {
      if (base.dens=="stdnorm") {
        sur <- 1-pnorm(x)
        
      } else if (base.dens=="exp") {
        sur <- 1-pexp(x)
        
      }
      
    }

    
    # if (roundoff) sur <- round(sur, digits=digits)
    # c(sur)
    
    sur

  }, error=function(ee) {print("From function: sur.snp"); print(ee); browser()})
  
} # end of sur.snp

#===============================================================================

iter.ints <- function(x, base.dens=c("stdnorm","exp"), k) {
  # Name:
  
  # - Iterated integrals
  
  #------------------------------
  # Outputs a series of recursive integrals based on the base density base.dens
  # This function mainly serves PSNP
  
  #------------------------------
  # Input:
  
  # - x           the data
  
  # - base.dens   the base density
  
  # - k           the number of iterations
  
  #------------------------------
  # Output:
  
  # - a matrix of integrals
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008 page 569
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  #------------------------------
  # Last modified: Jan 08 2013
  
  
  num.x <- length(x)
  
  Ints  <- matrix(nrow=num.x, ncol=2*k+1)
  
  if (base.dens=="stdnorm") {
    Ints[,1] <- pnorm(x)
    Ints[,2] <- -dnorm(x)
    
    for (i in 3L:(2L*k+1L)) {
      tmp      <- -((x)^(i-2L)) * (-Ints[,2L])
      # Ints[,i] <- ifelse(is.na(tmp), 0, tmp) + (i-2) * Ints[,i-2]      
      Ints[,i] <- tmp + (i-2L) * Ints[,i-2L]      
      # Must be (i-2) in the power and the multiplier to conform with the paper bc
      # R counts from 1 (not 0)
    }
    
  } else if (base.dens=="exp") {
    dex        <- dexp(x)
    Ints[,1]   <- 1-dex    

    for (i in 2L:(2L*k+1L)) {
      tmp      <- -(x^(i-1)) * dex
      # Ints[,i] <- ifelse(is.na(tmp), 0, tmp) + (i-1) * Ints[,i-1]
      Ints[,i] <- tmp + (i-1) * Ints[,i-1L]
      
    }
    # Must be (i-1) in the power and the multiplier to conform with the paper bc
    # R counts from 1 (not 0)
  }

  Ints
}

#===============================================================================

iter.ints1 <- function(x, base.dens=c("stdnorm","exp"), k) {
  # Name:
  
  # - Iterated integrals corresponding to calculate sur.snp directly
  
  #------------------------------
  # Outputs a series of recursive integrals based on the base density base.dens
  # This function mainly serves PSNP
  
  #------------------------------
  # Input:
  
  # - x           the data
  
  # - base.dens   the base density
  
  # - k           the number of iterations
  
  #------------------------------
  # Output:
  
  # - a matrix of integrals
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008 page 569
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  #------------------------------
  # Last modified: Jan 08 2013
  
  
  num.x <- length(x)
  
  Ints  <- matrix(nrow=num.x, ncol=2*k+1)
  
  if (base.dens=="stdnorm") {
    Ints[,1] <- 1-pnorm(x)
    Ints[,2] <- dnorm(x)
    
    for (i in 3L:(2L*k+1L)) {
      tmp      <- ((x)^(i-2)) * Ints[,2L]      
      Ints[,i] <- tmp + (i-2) * Ints[,i-2L]      
      # Must be (i-2) in the power and the multiplier to conform with the paper bc
      # R counts from 1 (not 0)
    }
    
  } else if (base.dens=="exp") {
    dex        <- dexp(x)
    Ints[,1L]   <- dex    
    
    for (i in 2L:(2L*k+1L)) {
      tmp      <- (x^(i-1)) * dex
      Ints[,i] <- tmp + (i-1) * Ints[,i-1L]
      
    }
    # Must be (i-1) in the power and the multiplier to conform with the paper bc
    # R counts from 1 (not 0)
  }
  
  Ints
}

#===============================================================================

# obsolete
snp.sur.CI <- function(tte, alpha, theta, base.dens=c("exp","stdnorm"), info.mat, method=c("Richardson","simple"), method.args=list(), bet=NULL, x=NULL, mod=c("PH", "AFT","PO")) {
  # Calculate the approximated survival CI based on delta rule
  
  #------------------------------
  # Caller:
  # - snp.surv.haz.fit 
  
  #------------------------------
  # Input:
  
  # - tte           The time points
  
  # - alpha         The significance level
  
  # - theta         =(phi, Mu, Sig)
  
  # - base.dens     either "exp" or "stdnorm" specifying the base density
  
  # - info.mat      The estimated covariance matrix mentioned above
  
  # - method        "Richardson" or "simple" used to get the first derivative
  #                 vector based on the function jacobian in package numDeriv
  
  # - method.args   A vector of argument to be passed to the above methods
  
  # - bet           the vector of regression coefficients (beta)
  
  # - x             the design matrix where each row is the covariate vector of 
  #                 the corresponding individuals
  
  # - mod           either "PH", "AFT" or "PO" referring to the Proportional Hazards
  #                 model, Accelerated Failure Time model and the Porportional Odds 
  #                 model repsectively
  
  #------------------------------
  # Output:
  
  # - A list structure with the following fields:
  
  #   + sur.est     estimated survival function at given time points
  
  #   + upp.band    the upperband of the confidence curve
  
  #   + low.band    the lowerband of the confidence curve
  
  #   + stdErr      the estimated standard error at given time points
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  # - Statistical Models A.C.Davison, p122.
  
  # - jacobian (numDeriv)
  
  #------------------------------
  # Created:        Apr 22nd 2010
  # Last modified:  May 17th 2011
  
  
  # Get the dimension of theta and bet
  l.the <- length(theta)
  l.bet <- length(bet)
  
  # Create a dummy function
  foo   <- function(params, tte, x, mod) {    
    S_K(tte, params[1:l.the], base.dens, 
        params[(l.the+1):(l.the+l.bet)], x, mod)
  }
  
  # Get the total number of time points
  l     <- length(tte)
  
  # Get the estimated variance then the estimated standard errors
  stdErr<- rep(0,l)
  
  for (i in 1L:l) {
    Var       <- var.delta.rule(c(theta, bet), info.mat, foo, method, method.args, 
                                tte[i], x[i,], mod)
    
    stdErr[i] <- sqrt(Var)
    
  }
  
  # Get the quantiles
  Z.alpha     <- qnorm(1-alpha)
  
  # Get the estimated survival curve and its upper and lower bands
  sur.est     <- S_K(tte, theta, base.dens, bet, x, mod)
  upp.band    <- sur.est + Z.alpha*stdErr
  low.band    <- sur.est - Z.alpha*stdErr
  
  res         <- list(sur.est = sur.est, upp.band = upp.band, 
                      low.band = low.band, stdErr = stdErr)
}

#===============================================================================
# obsolete
snp.haz.CI <- function(tte, alpha, theta, base.dens=c("exp","stdnorm"), info.mat, method=c("Richardson","simple"), method.args=list(), bet=NULL, x=NULL, mod=c("PH", "AFT","PO")) {
  # Calculate the approximated hazard CI based on delta rule
  
  #------------------------------
  # Caller:
  # - snp.surv.haz.fit
  
  #------------------------------
  # Input:
  
  # - t             The time points
  
  # - alpha         The significance level
  
  # - theta         =(phi, Mu, Sig)
  
  # - base.dens     either "exp" or "stdnorm" specifying the base density
  
  # - info.mat      The estimated covariance matrix mentioned above
  
  # - method        "Richardson" or "simple" used to get the first derivative
  #                 vector based on the function jacobian in package numDeriv
  
  # - method.args   A vector of argument to be passed to the above methods
  
  # - bet           the vector of regression coefficients (beta)
  
  # - x             the design matrix where each row is the covariate vector of 
  #                 the corresponding individuals
  
  # - mod           either "PH", "AFT" or "PO" referring to the Proportional Hazards
  #                 model, Accelerated Failure Time model and the Porportional Odds 
  #                 model repsectively
  
  #------------------------------
  # Output:
  
  # - A list structure with the following fields:
  
  #   + haz.est       estimated hazard function at given time points
  
  #   + upp.band      the upperband of the confidence curve
  
  #   + low.band      the lowerband of the confidence curve
  
  #   + stdErr        the estimated standard error at given time points
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  # - Statistical Models A.C.Davison, p122.
  
  # - jacobian (numDeriv)
  
  #------------------------------
  # Created:        Apr 22nd 2010
  # Last modified:  May 17th 2011
  
  # Get the dimension of theta and bet
  l.the   <- length(theta)
  l.bet   <- length(bet)
  
  # Create a dummy function
  foo     <- function(params, tte, x, mod) {    
    h_K(tte, params[1L:l.the], base.dens, params[(l.the+1L):(l.the+l.bet)], 
        x, mod)
  }
  
  # Get the total number of time points
  l       <- length(tte)
  
  # Get the estimated variance then the estimated standard errors
  stdErr  <- rep(0,l)
  
  for (i in 1L:l) {
    Var       <- var.delta.rule(c(theta, bet), info.mat, foo, method, 
                                method.args, tte[i], x[i,], mod)
    
    stdErr[i] <- sqrt(Var)
    
  }
  
  # Get the quantiles
  Z.alpha <- qnorm(1-alpha)
  
  # Get the estimated survival curve and its upper and lower bands
  haz.est <- h_K(tte, theta, base.dens, bet, x, mod)
  upp.band<- haz.est + Z.alpha*stdErr
  low.band<- haz.est - Z.alpha*stdErr
  
  res     <- list(haz.est = haz.est, upp.band = upp.band, low.band = low.band, 
                  stdErr = stdErr)
}

#===============================================================================

var.delta.rule <- function(par.est, par.hess, fun, method="Richardson",
                           method.args=list(), ...) {
  # Calculate the approximated variance based on delta rule
  
  #------------------------------
  # Caller(s):
  # - snp.sur.CI
  # - snp.haz.CI
  
  #------------------------------
  # Input:
  
  # - par.est       The estimated parameters whose estimated covariance matrix is
  #                 the inverse of par.hess
  
  # - par.hess      The hessian matrix mentioned above
  
  # - fun           A function of par.est that we want to estimate its variance at
  #                 par.est
  
  # - method        "Richardson" or "simple" used to get the first derivative
  #                 vector based on the function jacobian in package numDeriv
  
  # - method.args   A vector of argument to be passed to the above methods
  
  # - ...           Additional arguments to be passed to fun
  
  #------------------------------
  # Output:
  
  # - A scalar which is the estimated of the variance of a function of
  # par.est using delta method.
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  # - Statistical Models A.C.Davison, p122.
  
  # - jacobian (numDeriv)
  
  #------------------------------
  # Created:        Apr 22nd 2010
  # Last modified:  Jun 20th 2012
  
  # browser()
  J       <- jacobian(fun, par.est, method, method.args, ...)
  
  # To overcome this, use svd-based pseudo-inverse instead
  p       <- svd(par.hess)
  p$d     <- ifelse(p$d == 0, 0, 1/p$d)
  var.fun <- J %*% p$v %*% diag(p$d, ncol=length(p$d), 
                                nrow=length(p$d)) %*% t(p$u) %*% t(J)
}

#===============================================================================

dbase <- function(base.dens=c("exp","stdnorm"), x) {
  # Output the density values at inputs contained in x based on the density
  # given by base.dens
  # This function mainly serves DSNP
  
  #------------------------------
  # Input:
  
  # - base.dens   name of the density function
  
  # - x           vector storing the intputs
  
  #------------------------------
  # Output:
  
  # - A sequence of density values
  
  #------------------------------
  # Created: March 18th 2010
  
  return(exp(logdbase(base.dens, x)))
  
}

#===============================================================================

dsnp <- function(x, phi, base.dens=c("exp","stdnorm"), a=NULL) {
  # The density function of snp with respect to the base density specified in
  # base.dens
  
  #------------------------------
  # Caller:
  
  # - f_0
  
  #------------------------------
  # Input:
  
  # - x     the data
  
  # - phi   the parameter vectors
  
  # - a     the polynomial coefficients
  
  #------------------------------
  # Output:
  
  # - The value of the density snp at x
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  #------------------------------
  # Created:       18 Mar 2010
  # Last modified: 31 Jan 2013
  
  return(exp(logdsnp(x, phi, base.dens, a)))  
  
}

#===============================================================================
# Old name: dsnp.est
f_0 <- function(tt, theta, base.dens=c("exp","stdnorm")) {#, bet, x) {
  # The estimated density function of t with respect to the base density specified in
  # base.dens
  
  #------------------------------
  # Caller:
  
  # - snp.logLik, h_0, f_K
  
  #------------------------------
  # Input:
  
  # - tt        the time
  
  # - theta     the parameter vectors containing
  
  # + phi       the spherical coordinates
  
  # + mu        MU
  
  # + sig       Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  #------------------------------
  # Output:
  
  # - The value of the estimated density at time t (>0)
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Mar 19th 2010
  # Last modified: Jan 23rd 2011
  
  return(exp(logf_0(tt, theta, base.dens)))
  
}

#===============================================================================
# Old name: haz.snp.est
h_0 <- function(tt, theta, base.dens=c("exp","stdnorm"))  {
  # The estimated survival function of t with respect to the base density specified in
  # base.dens
  
  #------------------------------
  # Caller:
  
  # - 
  
  #------------------------------
  # Input:
  
  # - tt        the time
  
  # - theta     the parameter vectors containing
  
  # + phi       the spherical coordinates
  
  # + mu        MU
  
  # + sig       Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  #------------------------------
  # Output:
  
  # - The value of the estimated survival function of snp at t (>0)
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Mar 18th 2010
  # Last modified: Jan 23rd 2011
  
  h_0 <- f_0(tt, theta, base.dens) / S_0(tt, theta, base.dens)
  
  # c(h_0)
  h_0
  
}

#===============================================================================

rsurv.snp <- function(n, phi, sig=1, mu=0, base.dens=c("stdnorm", "exp")) {
  # Simulates a random vector of times to event (survival times) based on a snp
  # density
  
  #------------------------------
  # Input:
  
  # - n               the number of disired numbers
  
  # - phi             is a matrix whose rows are parameter vectors (phis). Thus,
  #                   the number of columns represents the dimension of the
  #                   parameter space.
  
  # - sig             Sigma
  
  # - mu              Mu
  
  # - base.dens       the base density which is either "stdnorm" or "exp"
  
  #------------------------------
  # Output:
  
  # - A simulated times to event vector based on a snp density
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       May 06th 2010
  # Last modified: May 17th 2011
  
  Z   <- rsnp(n, phi, base.dens)
  
  tte <- exp(mu + sig*Z)
}

#===============================================================================

rsnp <- function(n, phi=NULL, base.dens=c("stdnorm", "exp")) {
  # Simulates a random vector having a snp density (Not the time to event)
  
  #------------------------------
  # Input:
  
  # - n               the number of disired numbers
  
  # - phi             is a matrix whose rows are parameter vectors (phis). Thus,
  #                   the number of columns represents the dimension of the
  #                   parameter space. If phi==NULL then simply rnorm or log(rexp)
  #                   will be used depending on base.dens
  
  # - base.dens       the base density which is either "stdnorm" or "exp"
  
  #------------------------------
  # Output:
  
  # - A simulated snp random vector
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       May 06 2010
  # Last modified: Mar 27 2014
  
  re <- rep(NaN, n)
  
  if (is.null(phi) | length(phi)==0) { # if there's no phi=> simply parametric
    if (base.dens == "stdnorm") {
      re <- rnorm(n)
    } else if (base.dens == "exp") {
      re <- log(rexp(n))
    }
    
  } else { # if there's really phi
    
    # Call the appropriate routine based on the specified base density
    if (base.dens == "stdnorm") {

      if (n==1) {
        re <- sub.sim.stdnorm.snp.old(phi)
      } else {
        re <- sub.sim.stdnorm.snp(phi, N=n)
      } # end of if (n==1) else
      
    } else if (base.dens == "exp") {
      
      if (n==1) {
        re <- sub.sim.exp.snp.old(phi)
      } else {
        re <- sub.sim.exp.snp(phi, N=n)
      } # end of if (n==1) else
      
    } # end of if (base.dens==...)
    
  } # end of ifelse
  
  re
}

#===============================================================================

sub.sim.stdnorm.snp <- function(phi, N=1) { # not sure if works with N=1 => not work, use .old version
  # Simulates a random number having a snp density (Not the time to event) with a
  # standard normal base density. So far this function hasn't been vectorized yet!
  # Actually, this is is way faster than using a vectorized function.
  
  #------------------------------
  # Input:
  
  # - phi             is a matrix whose rows are parameter vectors (phis). Thus,
  #                   the number of columns represents the dimension of the
  #                   parameter space.
  
  #------------------------------
  # Output:
  
  # - A simulated snp random number
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  # - Gallant Tauchen 1990
  
  #------------------------------
  # Created:       May 06 2010
  # Last modified: Mar 27 2014
  
  tryCatch({ # To catch unprecedented errors
    
    # Create the coefficients from phi and base.dens
    a  <- polycofsnp("stdnorm", phi)
    
    # Take the absolute values of a and re-assign
    a  <- abs(a)
    
    p  <- polynomial(a)
    
    # Get the squared coefficients
    p2 <- p^2
    a2 <- coef(p2)
    n  <- length(a2)
    k  <- length(phi)
    if (n < (2*k+1)) {
      a2   <- c(a2, rep(0, 2*k+1-n))
      n    <- length(a2)
    }
    
    # Define the envelope function
    b  <- function(z) {
      foo <- function(i) {
        abs(z)^i
      }
      
      Z   <- sapply(0:(n-1), foo)
      
      re  <- a2%*%t(Z) * dnorm(z)
    } # End of function b
    
    # Get the vector of gamma evaluated at (i+1)/2  , i from 0
    get.gamma <- function(i) {
      re  <- gamma((i+1)/2)
    } # End of function get.gamma
    
    # Get the vector of powers of 2 i.e 2^((i-1)/2)  or  1 / 2^((1-i)/2), i from 0
    get.2.pow <- function(i) {
      re  <- 2^((i-1)/2)
    } # End of function get.2.pow
    
    # Get the unnormalized weights
    uw <- a2 * sapply(0:(n-1), get.gamma) * sapply(0:(n-1), get.2.pow)
    
    # Normalize the weights
    w  <- uw / sum(uw)
    
    # Initial value
    u  <- Inf # Later u will be generated from a uniform distribution on [0,1]
    
    hb <- -Inf # later hb will be h(v)/b(v),
    # where h is the snp density and b is the envelope function defined above.
    # V is generated from density g(v) = b(v) / int{b(s)}ds
    
    # Get the length of w
    lw <- length(w)
    
    re <- NULL
    
    m  <- 0 # number of successes in the rejection

    # Start the rejection method
    while (m < N) {
      # Use sample to specify the degree of freedom of the chi
      i   <- sample(x=1:lw, size=N-m, replace=TRUE, prob=w)      
      
      # Get the chi number by taking square root of a chisq number
      # This should be i instead of i+1 since the power starts from 0
      s.2 <- sapply(1:(N-m), function(id) rchisq(1, df=i[id]))
      
      # Randomly get the sign
      Sign<- rbern(N-m, 0.5)
      
      z   <- sqrt(s.2) * (-1)^Sign
      
      u   <- runif (N-m)
      
      hb  <- dsnp(z, phi, "stdnorm") / as.vector(b(z))
      
      re  <- c(re, z[u<=hb])
      re  <- na.omit(re)
      
      m   <- length(re)
    } # End of while
    
    re  # The final result
    
  }, error=function(ee) {print("From function: sub.sim.stdnorm.snp"); print(ee); browser()})
  
}

#===============================================================================

sub.sim.stdnorm.snp.old <- function(phi) {
  # Simulates a random number having a snp density (Not the time to event) with a
  # standard normal base density. So far this function hasn't been vectorized yet!
  # Actually, this is is way faster than using a vectorized function.
  
  #------------------------------
  # Input:
  
  # - phi             is a matrix whose rows are parameter vectors (phis). Thus,
  #                   the number of columns represents the dimension of the
  #                   parameter space.
  
  #------------------------------
  # Output:
  
  # - A simulated snp random number
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  # - Gallant Tauchen 1990
  
  #------------------------------
  # Created:       May 06 2010
  # Last modified: Jul 31 2012
  
  tryCatch({ # To catch unprecedented errors
    
    # Create the coefficients from phi and base.dens
    a  <- polycofsnp("stdnorm", phi)
    
    # Take the absolute values of a and re-assign
    a  <- abs(a)
    
    p  <- polynomial(a)
    
    # Get the squared coefficients
    p2 <- p^2
    a2 <- coef(p2)
    n  <- length(a2)
    k  <- length(phi)
    if (n < (2*k+1)) {
      a2   <- c(a2, rep(0, 2*k+1-n))
      n    <- length(a2)
    }
    
    # Define the envelope function
    b  <- function(z) {
      foo <- function(i) {
        abs(z)^i
      }
      
      Z   <- sapply(0:(n-1),foo)
      
      re  <- a2%*%Z * dnorm(z)
    } # End of function b
    
    # Get the vector of gamma evaluated at (i+1)/2  , i from 0
    get.gamma <- function(i) {
      re  <- gamma((i+1)/2)
    } # End of function get.gamma
    
    # Get the vector of powers of 2 i.e 2^((i-1)/2)  or  1 / 2^((1-i)/2), i from 0
    get.2.pow <- function(i) {
      re  <- 2^((i-1)/2)
    } # End of function get.2.pow
    
    # Get the unnormalized weights
    uw <- a2 * sapply(0:(n-1), get.gamma) * sapply(0:(n-1), get.2.pow)
    
    # Normalize the weights
    w  <- uw / sum(uw)
    
    # Initial value
    u  <- Inf # Later u will be generated from a uniform distribution on [0,1]
    
    hb <- -Inf # later hb will be h(v)/b(v),
    # where h is the snp density and b is the envelope function defined above.
    # V is generated from density g(v) = b(v) / int{b(s)}ds
    
    # Get the length of w
    lw <- length(w)
    
    z  <- NaN
    
    # Start the rejection method
    while (u > hb) {
      # Use sample to specify the degree of freedom of the chi
      i   <- sample(x=1:lw, size=1, replace=TRUE, prob=w)      
      
      # Get the chi number by taking square root of a chisq number
      # This should be i instead of i+1 since the power starts from 0
      s.2 <- rchisq(1, df=i)
      
      # Randomly get the sign
      Sign<- rbern(1, 0.5)
      
      z   <- sqrt(s.2) * (-1)^Sign
      
      u   <- runif(1)
      
      hb  <- dsnp(z, phi, "stdnorm") / b(z)
    } # End of while
    
    z  # The final result
    
  }, error=function(ee) {print("From function: sub.sim.stdnorm.snp"); print(ee); browser()})
  
}

#===============================================================================

sub.sim.exp.snp <- function(phi, N=1) { # not sure if works with N=1 => not work, use .old version
  # Simulates a random number having a log snp density (Not the time to event) with 
  # a standard exponential base density
  
  #------------------------------
  # Input:
  
  # - phi             is a matrix whose rows are parameter vectors (phis). Thus,
  #                   the number of columns represents the dimension of the
  #                   parameter space.
  
  #------------------------------
  # Output:
  
  # - A simulated snp random number
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       May 07 2010
  # Last modified: Mar 27 2014
  
  tryCatch({ # To catch unprecedented errors
    
    # Create the coefficients from phi and base.dens
    a  <- polycofsnp("exp", phi)
    
    # Take the absolute values of a and re-assign
    a  <- abs(a)
    
    p  <- polynomial(a)
    
    # Get the squared coefficients
    p2 <- p^2
    a2 <- coef(p2)
    n  <- length(a2)
    k  <- length(phi)
    if (n < (2*k+1)) {
      a2   <- c(a2, rep(0, 2*k+1-n))
      n    <- length(a2)
    }
    
    # Define the envelope function
    b  <- function(z) {
      
      foo <- function(i){z^i}
      Z   <- sapply(0:(n-1),foo)
      re  <- a2%*%t(Z) * dexp(z)
    } # End of function b
    
    # Get the vector of gamma evaluated at (i+1)/2
    get.gamma <- function(i) {
      re  <- gamma(i+1)
    } # End of get.gamma
    
    # Get the unnormalized weights
    uw <- a2 * sapply(0:(n-1), get.gamma)
    
    # Normalize the weights
    w  <- uw / sum(uw)
    
    # Initial value
    u  <- Inf # Later u will be generated from a uniform distribution on [0,1]
    
    hb <- -Inf # later v will be h(v)/b(v),
    # where h is the snp density and b is the envelope function defined above.
    # V is generated from density g(v) = b(v) / int{b(s)}ds
    
    # Get the length of w
    lw <- length(w)
    
    re <- NULL
    
    m  <- 0 # number of successes in the rejection
    # browser()
    # Start the rejection method
    while (m < N) {
      # Use sample to specify the shape parameter of the gamma distribution
      i  <- sample(x=1:lw, size=N-m, replace=TRUE, prob=w)
      
      # This should be i instead of i+1 as the power starts from 0
      z  <- sapply(1:(N-m), function(id) rgamma(1, shape=i[id], rate=1))
      
      u  <- runif(N-m)
      
      hb <- dsnp(z, phi, "exp") / as.vector(b(z))
      
      re  <- c(re, z[u<=hb])
      re  <- na.omit(re)
      
      m   <- length(re)
    } # End of while
    
    log(re)  # The final result
    
  }, error=function(ee) {print("From function: sub.sim.exp.snp"); ee; browser()})
}

#===============================================================================

sub.sim.exp.snp.old <- function(phi) {
  # Simulates a random number having a log snp density (Not the time to event) with 
  # a standard exponential base density
  
  #------------------------------
  # Input:
  
  # - phi             is a matrix whose rows are parameter vectors (phis). Thus,
  #                   the number of columns represents the dimension of the
  #                   parameter space.
  
  #------------------------------
  # Output:
  
  # - A simulated snp random number
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       May 07 2010
  # Last modified: Jul 31 2012
  
  tryCatch({ # To catch unprecedented errors
    
    # Create the coefficients from phi and base.dens
    a  <- polycofsnp("exp", phi)
    
    # Take the absolute values of a and re-assign
    a  <- abs(a)
    
    p  <- polynomial(a)
    
    # Get the squared coefficients
    p2 <- p^2
    a2 <- coef(p2)
    n  <- length(a2)
    k  <- length(phi)
    if (n < (2*k+1)) {
      a2   <- c(a2, rep(0, 2*k+1-n))
      n    <- length(a2)
    }
    
    # Define the envelope function
    b  <- function(z) {
      
      foo <- function(i){z^i}
      Z   <- sapply(0:(n-1),foo)
      re  <- a2%*%Z * dexp(z)
    } # End of function b
    
    # Get the vector of gamma evaluated at (i+1)/2
    get.gamma <- function(i) {
      re  <- gamma(i+1)
    } # End of get.gamma
    
    # Get the unnormalized weights
    uw <- a2 * sapply(0:(n-1), get.gamma)
    
    # Normalize the weights
    w  <- uw / sum(uw)
    
    # Initial value
    u  <- Inf # Later u will be generated from a uniform distribution on [0,1]
    
    hb <- -Inf # later v will be h(v)/b(v),
    # where h is the snp density and b is the envelope function defined above.
    # V is generated from density g(v) = b(v) / int{b(s)}ds
    
    # Get the length of w
    lw <- length(w)
    
    z  <- NaN
    
    # Start the rejection method
    while (u > hb) {
      # Use sample to specify the shape parameter of the gamma distribution
      i  <- sample(x=1:lw, size=1, replace=TRUE, prob=w)
      
      # This should be i instead of i+1 as the power starts from 0
      z  <- rgamma(1, shape=i, rate=1)
      
      u  <- runif(1)
      
      hb <- dsnp(z, phi, "exp") / b(z)
    } # End of while
    
    log(z)  # The final result
    
  }, error=function(ee) {print("From function: sub.sim.exp.snp"); ee; browser()})
}

#===============================================================================

is.nanull <- function(x) {
  # Check if x is NA or NULL and return TRUE in such case. Otherwise, return FALSE
  # Note: for variable of type list, data.frame, array, matrix which is NOT NULL
  # this function will check if there is any NA value inside.
  
  #------------------------------
  # Input:
  
  # - x     any variable of type numeric and character
  
  #------------------------------
  # Output:
  
  # - A logical TRUE or FALSE
  
  #------------------------------
  # Created:       Oct 10th 2011
  # Last modified: Oct 10th 2011    
  
  re <- is.null(unlist(x))
  
  re <- ifelse(re, re, re | any(is.na(unlist(x))))
  
  # c(re)
  re
}
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

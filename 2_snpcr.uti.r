# VERSION 2.3 UTI

# ANH NGUYEN DUC - OXFORD UNIVERSITY CLINICAL RESEARCH UNIT

# SUPERVISED BY:  MARCEL WOLBERS
# CREATED:        APR    12  2012
# LAST MODIFIED:  DEC    28  2016
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# FUNCTIONS INCLUDED:      
# get.meanvar.on.phi
# snp.meanvar         (tested 06 Aug 2012)
# parvect2parlist     (tested 06 Aug 2012)
# snp.cif.cb
# snp.cif             (tested 06 Aug 2012)
# logprob             (tested 06 Aug 2012)
# normalize           (tested 06 Aug 2012)
# simttesnp           (tested 06 Aug 2012)
# my.cos              (tested 06 Aug 2012) (blocked)
# snp.RMuSig
# slog
# dslog
# iwd.snp.onecif.old
# iwd.snp.onecif
# iwd.npmle.onecif
# my.timepoints
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

get.meanvar.on.phi <- function(phi, mu.in, sig.in, base.dens) {
  # Sub function used to get back mu and sigma from mean and variance of t_0
  # and those of the SNP varialbe assuming bet=0 w.r.t. phis
  
  meanvar <- snp.meanvar(phi, base.dens) # mean and variance of Z w.r.t phi 
  # base.dens
  sig.out <- sig.in / sqrt(meanvar[2]) # checked correct
  mu.out  <- mu.in - sig.out * meanvar[1]  # checked correct
  
  c(mu.out, sig.out)
  
} # End of get.meanvar.on.phi

#===============================================================================

snp.meanvar <- function(phi, base.dens) {
  # Find the mean and variance of a snp random variable w.r.t. phi and base.dens.
  # Note that when base.dens is "exp", it's the mean and variance of log(Z) that
  # are returned
  
  #------------------------------
  # Caller(s):
  # - snp.cropt
  
  #------------------------------
  # Input:  
  
  # - phi         Spherical coordinates
  
  # - base.dens   base density, either "exp" or "stdnorm"
  #------------------------------
  # Output:
  
  # mean and variance of the snp random variable
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  # - http://en.wikipedia.org/wiki/Gumbel_distribution
  
  #------------------------------
  # Created:       Aug 09 2011
  # Last modified: Jan 03 2013
  
  tryCatch({ # To catch unprecedented errors
    if (length(phi)!=0) {
      
      a <- polycofsnp(base.dens, phi)
      
      p <- polynomial(a)
      
      p2 <- p^2
      a2 <- coef(p2)
      n  <- length(a2)    
      k  <- length(phi)
      
      # To deal with situations where the last coefficients of p2 are 0s
      if (n < (2*k+1)) {
        a2   <- c(a2, rep(0, 2*k+1-n))
        n    <- length(a2)
      }
      
      if (base.dens=="stdnorm") {
        # E.Z    <- 0      
        # Get the moments, using 2*k+2 instead of n as we also need to get higher 
        # moments to compute the variance using EZ^2 - (EZ)^2
        # moments <- moment(base.dens, 1:(2*k+2)) # This is correct
        moments <- moment.n[(1:(2*k+2))+1] # using this global var. is faster
        
        # # This is the expansion of the expectation of X whose density is 
        # # {a_0+a_1*x+...+a_k*x^k}^2 * phi(x)         
        # foo <- function(i) {     
        #
        #   re      <- 0
        #   
        #   for (j in 1:n) {
        #     re    <- re + a2[j]*moments[j+i]
        #   }
        #   
        #   re    <- re - i*E.Z^2
        # }
        
        # E.Z   <- foo(0)
        # Var.Z <- foo(1)
        
        E.Z   <- sum(a2 * moments[1:n])
        Var.Z <- sum(a2 * moments[2:(n+1)]) - E.Z^2
        
      } else if (base.dens=="exp") {
        # In case of exponential base density. Note that it's Z*=e^Z that has a
        # density approximated by P_k(z)^2 * dexp(z). Thus we need to use delta
        # method to convert the expectation and variance of Z* to those of Z=logZ*
        # Now we nolonger use delta method but an iterative way
        
        # I <- H <- rep(NA, n)  
        # I[1] <- 0.5772   # mean of the Gumbel dist with mu=0 and bet=1 (wiki)
        # H[1] <- 1.978094 # 2nd central moment (EX^2) of the same dist
        
        # for (i in 2:n) {
        #   I[i] <- (i-1)*I[i-1] - factorial(i-2)
        #   H[i] <- -2*I[i-1]+(i-1)*H[i-1]
        # } # end of for (i in...)

        # Hard coded II and HH and make them global var. in 1_snpcr.r to gain speed
                
        E.Z  <- -sum(a2*II[1:n])
        Var.Z<- sum(a2*HH[1:n]) - (E.Z)^2
      } # end of if (base.dens=="stdnorm")
      
    } else { # When there is no phis i.e. only the base density
      if (base.dens=="stdnorm") {
        E.Z  <- 0
        Var.Z<- 1
      } else if (base.dens=="exp") {
        E.Z  <- -0.5772
        Var.Z<- pi^2/6
      }
    } # End of if (!is.null(phi) & length(NULL)==0)
    
    c(E.Z, Var.Z)   
    
  }, error=function(ee) {print("From function: snp.meanvar"); ee; browser()})
  
} # End of snp.meanvar

#===============================================================================

parvect2parlist <- function(params, parinfo, exp.sig=F) {#, t_m=Inf) {
  # Return a list of Gammas, Mu_Sigs, Phis and Bets given its vetorized version
  # and parinfo
  
  #------------------------------
  # Caller(s):
  # - snp.cropt
  
  #------------------------------
  # Input:  
  
  # - params     
  
  # - parinfo
  
  # - exp.sig     TRUE/FALSE whether we should take the exponentiation of Sigmas  
  #------------------------------
  # Output:
  
  # A list of Gammas, Mu_Sigs, Phis and Bets
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Aug 10 2011
  # Last modified: Jun 04 2015
  
  tryCatch({ # To catch unprecedented errors
    
    J          <- parinfo$nc
    
    Gamma      <- t( matrix(params[1:( parinfo$gammas * (J-1) )], 
                            nrow=parinfo$gammas) )
    
    # A counter
    cnt        <- parinfo$gammas * (J-1)
    
    # This is a J X 2 matrix
    Mu_Sig     <- t( matrix(params[( cnt + 1 ):( cnt + 2*J)], nrow=2) )
    if (exp.sig) {
      Mu_Sig[,2] <- exp(Mu_Sig[,2])
    }
    
    cnt        <- cnt + 2 * J
    
    # This is a list with J elements
    Phi    <- list()
    
    for (j in 1:J) {
      if (parinfo$phis[j]!=0) {
        Phi[[j]]   <- params[(cnt + 1):(cnt + parinfo$phis[j])]
        
      } else {
        Phi[[j]]   <- numeric(0) # set to NULL causes the last null emelent removed
        
      }
      cnt      <- cnt + parinfo$phis[j]
    }
    
    # This is a list with J elements
    Bet <- list()
    
    if (length(parinfo$bets)!=0) {
      for (j in 1:J) {  
        if (parinfo$bets[j]!=0) {
          Bet[[j]] <- params[(cnt + 1):(cnt + parinfo$bets[j])]
        } else {
          Bet[[j]] <- numeric(0) # NULL (old version)
        }
        cnt     <- cnt + parinfo$bets[j]
      }      
    }
    
    list(Gamma=Gamma, Mu_Sig=Mu_Sig, Phi=Phi, Bet=Bet)
    
  }, error=function(ee) {print("From function: parvect2parlist"); ee; browser()})
  
} # End of parvect2parlist

#===============================================================================

snp.cif.cb <- function(tt, params, base.dens, y=NULL, x=NULL, P, mod="AFT", cimethod="cloglog", 
                       parinfo=NULL, hessian=NULL, K=NULL, exp.sig=F, j=NULL,
                       quan.type="stdnorm", quan.df=NULL) {
  # The estimated cif of tt for one competing risk togethe with the se and the 
  # lower and upper band with respect to the base 
  # density specified in base.dens for any of the following regression models: 
  # PH, AFT and PO
  
  # It seems that right now this function only works with t_m=Inf
  #------------------------------
  # Caller:
  
  # - snp.haz.CI
  
  #------------------------------
  # Input:
  
  # - tt          the time
  
  # - params      parameters of all CRs, note that if sigs is on logscale then
  #               exp.sig must be T
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - y           the design matrix with respect to bet 
  
  # - x           the design matrix with respect to gamma
  
  # - P           the vector of marginal probabilities based on gamma and x. 

  # - mod         either "PH", "AFT" or "PO" referring to the Proportional Hazards
  #               model, Accelerated Failure Time model and the Porportional Odds 
  #               model repsectively

  # - cimethod    method used to calculate the confidence ban for the cif. 
  #               "linear" means wald confi band for cif, "cloglog" means using
  #               delta method to calculate the confi band for the cloglog 
  #               transform of the cif and then transform the confi band back
  
  # - parinfo     information on the parameters (dimension)
  
  # - hessian     hessian matrix of all parameters of all crs     
  
  # - K           the K matrix used in the calculation of the sandwich variance
  #               together with hessian. If we want to use only the hessian =>
  #               set this to NULL
  
  # - exp.sig     if the sigma parameters in params are actually log(sigma)
  
  # - j           the competing risk for which we want to get the confi band
  
  # - quan.type   the type of distribution used to calculate the se: norm or t
  
  # - quant.df    only needed when quan.type is t
  
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
  # Created:       Apr 16 2012
  # Last modified: Jun 04 2015
  tryCatch({  
    # browser()
    foo1 <- function(p) {
      # p =(Gammas (for all CRs), Phi_j, Mu_j, Sig_j, Bet_j)    
      
      pl <- parvect2parlist(p, parinfo, exp.sig=exp.sig)
    
      re <- snp.cif(tt, theta=c(pl$Phi[[j]], pl$Mu_Sig[j,]),
                    base.dens, pl$Bet[[j]], y,
                    exp(logprob(X=x, pl$Gamma))[1,j], 
                    mod="AFT")
      # when re is already small, close to 0, this transform can reslut in -Inf which causes 
      # the jacobian to be Nan. But just happen at the earliest timepoints
      if (cimethod=="cloglog") re <- log(-log(1-re))
      
      re
    } # end of foo1
    
    parlist <- parvect2parlist(params, parinfo, exp.sig=exp.sig)
    if (is.null(x)) x <- matrix(rep(1, ncol(parlist$Gamma)), nrow=1)
    # browser()
    # note that here P is just a single number of CR_j
    snp.CIF <- snp.cif(tt, theta=c(parlist$Phi[[j]], parlist$Mu_Sig[j,]), 
                       base.dens, parlist$Bet[[j]], y, P[j], mod)
    
    res     <- snp.CIF
    if (cimethod=="cloglog") res <- log(-log(1-res))
    
    # get the covariance matrix from the hessian matrix or the sandwich version 
    # if K is not NULL, the hessian should correspond to the NEGATIVE loglik
    if (length(K)==0) { 
      covmat  <- try(solve(-hessian), silent=T)
    } else {
      covmat  <- try(solve(hessian) %*% K %*% solve(hessian), silent=T)
    } # end of if (is.null(K)) else...   
    
    if (class(covmat)[1] != "try-error") {
      # if (j==2) browser()
      parlist  <- parvect2parlist(params, parinfo, exp.sig=F)
      jcb1     <- jacobian(foo1, params) # delta method
      
      #ids      <- match(c(parlist$Gamma, parlist$Mu_Sig[j,], parlist$Phi[[j]], 
      #                    parlist$Bet[[j]]), params)
      #browser()
      se1      <- sqrt(diag(jcb1 %*% covmat %*% t(jcb1))) 
      #se1      <- sqrt(diag(jcb1[,ids] %*% covmat[ids,ids] %*% t(jcb1[,ids]))) 
      
      if (quan.type=="stdnorm") {    
        low<- res - qnorm(1-.05/2)*se1
        upp<- res + qnorm(1-.05/2)*se1
      } else {
        low<- res - qt(1-.05/2, df=quan.df)*se1
        upp<- res + qt(1-.05/2, df=quan.df)*se1
      } # end of if (quan.type=="stdnorm") else ...
      
      # browser()  
      if (cimethod=="cloglog") {
        low <- 1-exp(-exp(low))
        upp <- 1-exp(-exp(upp))
        
      } # end of if (cimethod=="cloglog")
      # se and se1 are differengt
      
    } else { # when the hessian is singular
      low <- Inf
      upp <- -Inf
      se1  <- NaN
    } # end of if (class(covmat)[1] != "try-error") else ...
    
    return(list(snp.CIF=snp.CIF, low=low, upp=upp, se=se1))
    
  }, error=function(ee) {print("From function: snp.cif.cb"); ee; browser()})
} # end of snp.cif.cb

#===============================================================================
snp.cif <- function(tt, theta, base.dens, bet, y, P, mod="AFT") {
  # The estimated cif of tt for one competing risk with respect to the base 
  # density specified in base.dens for any of the following regression models: 
  # PH, AFT and PO
  
  #------------------------------
  # Caller:
  
  # - snp.haz.CI
  
  #------------------------------
  # Input:
  
  # - tt          the time
  
  # - theta       the parameter vectors containing
  
  # + phi         the spherical coordinates
  
  # + mu          MU
  
  # + sig         Sigma
  
  # - base.dens   either "exp" or "stdnorm" specifying the base density
  
  # - bet         the vector of regression coefficients (beta) for the snp part
  
  # - y           the design matrix with respect to bet 
  
  # - P           the vector of marginal probabilities based on gamma and x. 

  # - mod         either "PH", "AFT" or "PO" referring to the Proportional Hazards
  #               model, Accelerated Failure Time model and the Porportional Odds 
  #               model repsectively

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
  # Created:       Apr 16 2012
  # Last modified: Jun 04 2015
  tryCatch({

    F_K.aft(tt, theta, base.dens, bet, y) * P # give same result as above

  }, error=function(ee) {print("From function: snp.cif"); ee; browser()})
}
#===============================================================================

logprob <- function(X, Gamma, N=nrow(X)) {
  # Calculate the marginal probabilities for the given Gamma and X
  
  #------------------------------
  # Caller(s):
  # - snp.cropt, em.weights, snp.ecrlogLik.mar, spn.crlogLik
  
  #------------------------------
  # Input:  
  
  # - X           the covariates matrix with each row having the covariates for
  #               each patient
  
  # - Gamma       the coefficents
  
  # - N           the number of rows in X (number of patients)  
  #------------------------------
  # Output:
  
  # A matrix of marginal probabilities, each row for each patient
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  #------------------------------
  # Created:       Mar 13 2012
  # Last modified: Oct 23 2012  
  
  tryCatch({ # To catch unprecedented errors        
    # XGamma <- X %*% t(Gamma) # This is wrong 
    XGamma <- matrix(t(apply(X, 1, function(x) x * Gamma)), nrow=nrow(X))
    XGamma <- cbind(XGamma, 0)
    logP   <- t(apply(XGamma, MARGIN=1, normalize))
  }, error=function(ee) {print("From function: logprob"); print(ee); browser()})
  
} # End of logprob

#===============================================================================

normalize <- function(x) {
  # calculate log(exp(x)/sum(exp(x))) in a numerically stable way
  #------------------------------
  # Created:       Apr 11 2012
  # Last modified: Aug 06 2012
  
  tryCatch({ # To catch unprecedented errors
    m <- max(x)
    tmp <- log(sum(exp(x-m)))
    x-m-tmp
  }, error=function(ee) {print("From function: normalize"); ee; browser()})
  
} # End of normalize

#===============================================================================

simttesnp <- function(n, mu, sig, phi, base.dens, x=NULL, bet=NULL) {
  # simulate snp distributed survival time
  #------------------------------
  # Created:       May 29 2012
  # Last modified: Aug 06 2012
  
  tryCatch({ # To catch unprecedented errors
    if (!is.null(x) & !is.null(bet)) {
      tte <- exp(mu + x %*% bet + sig*rsnp(n, phi, base.dens))
    } else {
      tte <- exp(mu + sig*rsnp(n, phi, base.dens))
    }
    tte
  }, error=function(ee) {print("From function: simttesnp"); print(ee); browser()})
  
} # End of normalize

#===============================================================================

# my.cos <- function(x) {
#   tmp <- cos(x)
#   re  <- ifelse(round(x%%(2*pi), 6)==1.570796 | round(x%%(2*pi), 6)==1.570797 |
#                 round(x%%pi, 6)==1.570796 | round(x%%pi, 6)==1.570797, 0, tmp)
#   re
# }
#===============================================================================

snp.RMuSig <- function(Data, Z, t_m=Inf, Base.dens, method=c("BFGS", "Nelder-Mead", 
                                                             "CG", "L-BFGS-B", "SANN"), 
                       maxstep=5, Ks=NULL, init_params=NULL, maxit=100, optonsig=F,
                       bound.onphi=F, parallelism=T, ncores=NULL, criterion="HQd", 
                       check.hess=T, anal.grad=T) {
  # Find the point (Sig, Mu, phi) that maximizes the log-likelihood with a base 
  # density for only the following regression model AFT
  
  # Note that only time independent covariates are dealt with
  
  #------------------------------
  # Caller: 
  
  #------------------------------
  # Input:
  
  # - Data        a data frame or matrix of the following form:
  
  #   + For data subejct to right-censoring and left-truncation
  #   V   |   Del   |  LT | D
  #   + V is the time points. Delta takes value from 0,...,J; where J is the 
  #     number of competing risks. D is a dummy columns to indicate that we're
  #     not dealing with interval-censoring
  
  #   + For data subject to only interval-censoring
  #   L   |   R   |   Delta
  #   + L, R are the left and right time points. Delta takes value from 0,...,J; 
  #     where J is the number of competing risks. Also when R = infinity or Delta 
  #     equals 0, the case is considered right-censored.
  
  # - Z           a list having:
  #               + X       - a matrix of n x p_x dimension where n is the total
  #                           number of observations and p_x is the number of 
  #                           columns of Gamma described later.
  
  #               + Y       - a sub-list having J matrices. Each matrix is of size  
  #                           N x L corresponding to Bet described later.
  
  # - t_m         must be larger than any time points recorded in a study i.e. 
  #               failure times, right-censoring times, study entry time and 
  #               interval endpoints for interval censoring. When t_m equals
  #               infinity, we need to adjust how we model the probabilities of
  #               having events at t_m. 
  
  # - Base.dens   a vector of length J whose elements can be either "stdnorm" or
  #               "exp".
  
  # - method      the methods used by R's optim if "BFGS", "L-BFGS-B",
  #               if "genoud" then function genoud is used
  
  # - maxstep     the maximum number of iterations allowed
  
  # - Ks          (Optional) a vector containing the degrees of polynomials for 
  #               each competing risk. If this is set to NULL by default, the 
  #               forward procedure will be applied. Otherwise, the function will 
  #               use exactly the polynomial degrees provided and will still use 
  #               the forward procedure. If we don't want to use the forward 
  #               procedure, we simply set maxstep to 1
  
  # - init_params (Optional) a list having user selected initial parameters for 
  #               Gammas, Mu_Sigs, Phis and Bets. Note that if this is set to NOT
  #               NULL, Ks above won't be used.
  
  # - maxit       maximum number of iterations used by optim
  
  # - optonsig    if T and only if method is not "BFGS" and not "nlm", 
  #               optimization will be done on the original scale of sigma with 
  #               the constrain that sigma >= 1e-8
  
  # - bound.onphi if T and only if method is not "BFGS"and not "nlm", 
  #               optimization will be done on phis with the box constrains
  #               -pi/2<= phi <= pi/2
  
  # - parallelism a flag indicating whether or not to use parallelism. If NULL, all
  #               available threads will be utilized
  
  # - ncores      the number of cores to be used by parallelism if any
  
  # - criterion   model selection criterion for stepwise selection: AIC, BICn,
  #               HQn, BICd or HQd
  
  # - check.hess  if T estimated parameters whose hessian matrix has negative 
  #               eigenvalues will be removed
  
  # - anal.grad   if T (default), analytic gradient of -loglik is used

  #------------------------------
  # Output:
  
  # - A list with M sub-list, where M is the number of iterations reached and each
  # sub-list has the following fields: 
  
  #               + Gamma   - a matrix of regression coefficients involving the
  #                           multilogit-model whose number of rows is J-1
  
  #               + Mu_Sig  - a matrix of J x 2 dimension with the 1st column
  #                           storing parameters mu and the 2nd colum storing 
  #                           parameters sigma for each competing risk.
  
  #               + Phi     - a sub-list with J elements corresponding to J 
  #                           vectors phi of each competing risk.
  
  #               + Bet     - a sub-list with J elements corresponding to J 
  #                           vectors beta of each competing risk.
  
  #               + P       - a N x (J+1) matrix having the estimated P_ij(z_i)
  
  #               + mloglik - negative log-likelihood
  
  #               + AIC     - corresponding AIC 
  
  #               + hessian, counts, convergence  (outputs given by optim)
  
  #               + p       - the number of parameters
  
  #               + params  - a vector having all the parameters  
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  # - Unknown Anh's publication
  
  #------------------------------
  # Created:       Mar 21 2013
  # Last modified: Jun 04 2015
  
  # This is inside snp.RMuSig
  
  # Note that this function won't work if any of the CR has no observed events
  
  tryCatch({ # To catch unprecedented errors
    
    # Check if we're dealing with interval-censoring only  
    is.IC  <- ncol(Data) < 4
    
    # Check if there's left-truncation
    is.LT  <- any(na.omit(Data[,3])!=0) & !is.IC
    
    # Get the number of CRs
    J      <- length(Base.dens)
    
    # Only methods other than BFGS and nlm support box constrains on phi or sigmas
    if (method=="BFGS" | method=="nlm") {
      optonsig <- bound.onphi <- F
    }
    
    #------------------------------------------------------------------------------
    # Setting up initial values for mus, sigmas, betas, gammas
    # first use cuminc to get P(D=j,T>t), j=1...J
    if (!is.IC) { # RC and LT case
      ftime     <- Data[,1]
      fstatus   <- Data[,2]
      
    } else { # only IC case
      ftime     <- ( Data[,1] + Data[,2] ) / 2
      fstatus   <- Data[,3]    
      ftime[Data[,2]==Inf] <- Data[Data[,2]==Inf,1]
      
    } # End of ifelse      
    
    # results of this will be used to get initial values for mus, sigmas and gammas
    tmp.cum       <- cuminc(ftime, fstatus, cencode=0)
    
    if (is.null(init_params)) { # When the user doesn't pass in initial pars
      # Get the starting values for mus and sigmas when Ks = 0s
      init.Mu_Sig <- matrix(NA, nrow=J, ncol=2)
      
      if (!is.IC) { # If RC and LT case
        
        for (j in 1:J) {
          Id                 <- Data[,2]==j | Data[,2]==0
          Dat                <- Data[Id,] 
          Dat[Dat[,2]==j, 2] <- 1
          
          # create the weight vectors
          weights            <- ifelse(Dat[,2]==j, 1, 
                                       1-timepoints(tmp.cum, Dat[Dat[,2]==0,1])$est[j,])
          weights            <- weights/sum(weights)
          
          colnames(Dat)      <- c("ttdeath", "evdeath", "starttime", "")
          # formula            <- Surv(starttime, ttdeath, evdeath, 
          #                            type="counting")~1
          # for now survreg doesn't support left truncation so we neglect this
          formula            <- Surv(ttdeath, evdeath)~1
          
          if (Base.dens[j]=="stdnorm") { # In case the base density is stdnorm      
            tmp                <- survreg(formula, data.frame(Dat), 
                                          weights=weights, dist="lognormal")                
          } else { # In case the base density is weibull
            
            tmp                <- survreg(formula, data.frame(Dat), 
                                          weights=weights, dist="weibull")               
          } # End of ifelse
          
          mu                 <- coef(tmp)
          attr(mu, "names")  <- NULL
          sig                <- tmp$scale  
          
          init.Mu_Sig[j,]    <- c(mu, sig)
        } # End of forloop 
        
      } else { # If only IC case
        for (j in 1:J) {
          Id                  <- Data[,3]==j | Data[,3]==0
          Dat                 <- Data[Id,] 
          Dat[Dat[,3]==j, 3]  <- 3
          
          # create the weight vectors
          weights             <- ifelse(Id==j, 1, 
                                        1-timepoints(tmp.cum, Dat[Dat[,3]==0,1])$est[j,])
          weights             <- weights/sum(weights)
          
          colnames(Dat)       <- c("l", "r", "evdeath")
          formula             <- Surv(time=l, time2=r, event=evdeath,
                                      type='interval')~1     
          
          # Right now survreg won't work if there is 0 in the left-points
          idl0        <- Dat$l==0
          Dat$l[idl0] <- Dat$l[idl0] + min(Dat$l[Dat$l!=0])/10000
          Dat$r[idl0] <- Dat$r[idl0] + min(Dat$l[Dat$l!=0])/10000
          
          if (Base.dens[j]=="stdnorm") { # In case the base density is stdnorm      
            tmp                <- survreg(formula, Dat, weights=weights,
                                          dist="lognormal")                  
          } else { # In case the base density is weibull
            tmp                <- survreg(formula, Dat, weights=weights, 
                                          dist="weibull")            
          } # End of ifelse
          
          mu                 <- coef(tmp)
          attr(mu, "names")  <- NULL
          sig                <- tmp$scale
          
          init.Mu_Sig[j,]    <- c(mu, sig)
        } # End of forloop 
        
      } # End of ifelse   
      
      #--------------------------------------------------------------------------      
      # In this variant of getting initial values for the intercepts of Gammas,
      # for each set of initial values for Mus and Sigmas, we optimize the 
      # loglikelihood with only the intercepts of Gammas as variable. This
      # will be done later in the first iteration of the algorithm
      init.Gamma  <- matrix(0, nrow=ifelse(t_m<Inf, J, J-1), ncol=ncol(Z$X))
      
      # However we still need starting values for optim, and this will be taken
      # from the current method of getting initial values for gammas at the first
      # step
      init.P    <- rep(NA, J)
      
      # Potential problems will occur if one or more CR doesn't have any
      # corresponding observed event time (not yet addressed, actually user are not 
      # allowed to put in such data)
      for (j in 1:length(init.P)) {
        init.P[j]    <- max(tmp.cum[[j]]$est)
      }           
      
      if (t_m == Inf) {    
        init.P         <- init.P/sum(init.P)
      } else {
        # for now when t_m < Inf, we still need to consider a new cuminc model
        # with censoring considered as the 3rd cause. Otherwise the P(D=0) at
        # t_m will always be estimated as 0 and this leads to estimates of the
        # gammas which are Inf
        new.fstatus <- ifelse(fstatus==0, 3, fstatus)
        new.cum <- cuminc(ftime, new.fstatus, cencode=0)
        for (j in 1:length(init.P)) {
          init.P[j] <- max(new.cum[[j]]$est)
        }
        init.P[J+1] <- max(new.cum[[J+1]]$est)
        #init.P[J+1]    <- 1-sum(init.P)
      }
      
      init.Gamma[,1] <- log(init.P[-length(init.P)] / init.P[length(init.P)]) # old use my.log
      
    } else { # If the user does put in initial parameters
      
      init.Gamma <- init_params$Gamma
      init.Mu_Sig<- init_params$Mu_Sig
      init.Bet   <- init_params$Bet
      init.Phi   <- init_params$Phi    
      
      if (nrow(init.Gamma) != J & t_m<Inf) {
        stop("The number of rows of init_params$Gamma must equal the # of CRs!") 
      } 
      
      if (nrow(init.Gamma) != (J-1) & t_m==Inf) {
        stop("When t_m = Inf, the number of rows of init_params$Gamma must equal the # of CRs minus 1!") 
      }    
      
      if (ncol(init.Gamma) != ncol(Z$X)) {
        stop("The # of columns of init_params$Gamma must equal ncol(Z$X)")
      }
      
      if (nrow(init.Mu_Sig) != J) {
        stop("The # of rows of init_params$Mu_Sig must equal the # of CRs!")
      }
      
      if (ncol(init.Mu_Sig) != 2) {
        stop("The # of columns of init_params$Mu_Sig must equal 2!")
      }
      
      # if (length(init_params$Bet) != J) {
      #   stop("The length of init_params$Bet must equal the # of CRs!")
      # }
      
      if (length(init_params$Phi) != J) {
        stop("The length of init_params$Phi must equal the # of CRs!")
      }
      
      # for (j in 1:J) {
      #   if (length(init_params$Bet[[j]]) != ncol(Z$Y[[j]])) {
      #     stop("The length of each init_params$Bet[[j]] must equal the # of columns in Z$Y[[j]]")
      #   }
      # }
      
    } # End of ifelse
    
    init.Mu_Sig
    
  }, error=function(ee) {print("From function: snp.RMuSig"); ee; browser()})
  
} # End of snp.cropt
#===============================================================================

slog <- function(x) {
  tryCatch({ # To catch unprecedented errors
    t2    <- 1e-7
    t1    <- 0
    a     <- -299.999999999999886
    b     <- 5667638086.9808321
    c     <- - 28288190434904165    
    
    id1 <- x<=t1 & !is.na(x)
    id2 <- x>t1 & x<t2 & !is.na(x)
    id3 <- x>=t2 & !is.na(x)
    
    sl <- rep(NA, length(x))
    sl[id1] <- a + b*x[id1]
    sl[id2] <- a + b*x[id2] + c*x[id2]^2
    sl[id3] <- log(x[id3])
    
    sl
    
  }, error=function(ee) {print("From function: slog"); ee; browser()})
} # end of slog
#===============================================================================

dslog <- function(x) {
  tryCatch({ # To catch unprecedented errors
    t2    <- 1e-7
    t1    <- 0
    b     <- 5667638086.9808321
    c     <- - 28288190434904165
    
    id1 <- x<=t1
    id2 <- x>t1 & x<t2
    id3 <- x>=t2
    
    sl <- rep(NA, length(x))
    sl[id1] <- b
    sl[id2] <- b + 2*c*x[id2]
    sl[id3] <- 1/x[id3]
    
    sl
  }, error=function(ee) {print("From function: slog"); ee; browser()})
} # end of dslog
#===============================================================================

iwd.snp.onecif <- function(parlist1, parlist2, j=1, wfun=function(tt) {1}, wfix=1,
                           lower=0, upper, ..., subdivisions=500,
                           rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
                           stop.on.error = TRUE) {
  
  # test the same cif coming from 2 different groups
  # this uses integrate and snp-based cifs
  
  # parlist1      a list whose elements are
  #               $base.dens, $Gamma, $Mu_Sig, $Phi and 
  #               $Bet, $Y, $X, $hessian, $K
  #               When a is NULL, $Phi will be used.
  
  # parlist2      the same as parlist1
  
  # wfun          the weight function. By default w is the identity function

  # wfix          the weight constant
  
  # j             the CR of interest
  
  # other parameters are from integrate in package stats  
  
  # created       08 Jan 2014
  # Last modified 30 Dec 2015
  
  
  tryCatch({ # To catch unprecedented errors     
    # browser()
    # get the number of CRs
    parinfo1 <- parinfo2 <- list()
    parinfo1$nc <- parinfo2$nc <- nrow(parlist1$gammas) + 1
    
    parinfo1$gammas <- ncol(parlist1$gammas)
    parinfo1$phis   <- unlist(lapply(parlist1$phis, length))
    parinfo1$bets   <- unlist(lapply(parlist1$bets, length))

    parinfo2$gammas <- ncol(parlist2$gammas)    
    parinfo2$phis   <- unlist(lapply(parlist2$phis, length))
    parinfo2$bets   <- unlist(lapply(parlist2$bets, length))

    nc <- parinfo1$nc
    l1 <- parinfo1$gammas * (nc-1) + 2*nc + sum(parinfo1$phis) + sum(parinfo1$bets)
    l2 <- parinfo2$gammas * (nc-1) + 2*nc + sum(parinfo2$phis) + sum(parinfo2$bets)
        
    n1 <- parlist1$n; n2 <- parlist2$n
    
    Z1   <- Z2 <- list()
    Z1$X <- matrix(rep(1,n1), nrow=n1); Z2$X <- matrix(rep(1,n2), nrow=n2)
    Z1$Y <- Z2$Y <- list(numeric(0), numeric(0))
    
    iwd <- function(p) {# this version is for fixed weight function w
      # p has the parameters from 2 groups      
      p1 <- p[1:l1]
      p2 <- p[(l1+1):(l1+l2)]      
      
      tmp.pl1 <- parvect2parlist(params=p1, parinfo=parinfo1, exp.sig=F)
      tmp.pl2 <- parvect2parlist(params=p2, parinfo=parinfo2, exp.sig=F)
      
      re  <- integrate(f=function(tt) {
        (snp.cif(tt, theta=c(tmp.pl1$Phi[[j]], tmp.pl1$Mu_Sig[j,]), 
                 parlist1$base.dens, tmp.pl1$Bet[[j]], Z1$Y, 
                 exp(logprob(X=Z1$X, tmp.pl1$Gamma))[1,j], 
                 mod="AFT") -
           snp.cif(tt, theta=c(tmp.pl2$Phi[[j]], tmp.pl2$Mu_Sig[j,]), 
                   parlist2$base.dens, tmp.pl2$Bet[[j]], Z2$Y, 
                   exp(logprob(X=Z2$X, tmp.pl2$Gamma))[1,j], 
                   mod="AFT")) * wfun(tt) 
      },  
      lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol, 
      stop.on.error=stop.on.error)
      
      re$value * wfix
    } # end of iwd
    
    p1 <- c(as.vector(t(parlist1$gammas)), as.vector(t(parlist1$mu.sigs)), 
            unlist(parlist1$phis), unlist(parlist1$bets))
    
    p2 <- c(as.vector(t(parlist2$gammas)), as.vector(t(parlist2$mu.sigs)), 
            unlist(parlist2$phis), unlist(parlist2$bets))            
    
    p  <- c(p1, p2)
    
    iwd.snp <- iwd(p)
    
    var <- p.val <- NULL
    
    # get the covariance matrix from the hessian matrices or the sandwich version 
    # if Ks are not NULL, the hessians should correspond to tne NEGATIVE logliks
    if (length(parlist1$K)==0) { 
      covmat1  <- try(solve(parlist1$hessian), silent=T)
    } else {
      hes1inv  <- try(solve(parlist1$hessian), silent=T)
      covmat1  <- try(hes1inv %*% parlist1$K %*% hes1inv, silent=T)
    } # end of if (is.null(parlist1$K)) else...    
    
    if (length(parlist2$K)==0) { 
      covmat2  <- try(solve(parlist2$hessian), silent=T)
    } else {
      hes2inv  <- try(solve(parlist2$hessian), silent=T)
      covmat2  <- try(hes2inv %*% parlist2$K %*% hes2inv, silent=T)
    } # end of if (is.null(parlist2$K)) else...   
    
    if (class(covmat1)[1] != "try-error" & class(covmat2)[1] != "try-error") {
      covmat <- bdiag(covmat1, covmat2)  
      jcb    <- jacobian(iwd, p) # delta method, tested against grad, same result but on avg faster              
      var    <- diag(jcb %*% covmat %*% t(jcb)) # this is infact just 1x1 matrix, so diag is the whole matrix
      
      # below codes are for a ROUGHE second order delta method but is not stable e.g. can give negative variance
      # foo    <- function(p) {sum(grad(iwd, p)^2) - 0}
      # p0     <- optim(p, foo)$par
      # p0     <- optim(p, iwd)$par
      # H      <- hessian(iwd, p)
      # browser()
      # var    <- diag(jcb %*% covmat %*% t(jcb)) + .5 * sum(diag( H%*%covmat%*%H%*%covmat )) - 
      #           2 * diag(jcb %*% covmat %*% H %*% t(jcb)) - diag(p %*% t(jcb)) * sum(diag(t(H)%*%covmat)) #-
                
      # jcb2   <- hessian(iwd, p)
      # var2   <- covmat %*% jcb2
      p.val  <- pnorm(q=-abs(iwd.snp), mean=0, sd=sqrt(var)) * 2
    } # end of if (class(covmat1)[1] != "try-error" & class(covmat2)[1] != "try-error") 
    
    list(iwd=iwd.snp, var=var, p.val=p.val)  
    
  }, error=function(ee) {print("From function: iwd.snp.onecif"); ee; browser()})
} # end of iwd.snp.onecif
#===============================================================================

iwd.snp.onecif.old <- function(parlist1, parlist2, parinfo1, parinfo2, j=1,
                           wfun=function(tt) {1}, parlistw1=NULL, parlistw2=NULL, 
                           wfix=1,
                           lower=0, upper, ..., subdivisions=500,
                           rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
                           stop.on.error = TRUE) {
  
  # test the same cif coming from 2 different groups
  # this uses integrate and snp-based cifs
  
  # parlist1      a list whose elements are
  #               $base.dens, $Gamma, $Mu_Sig, $Phi and 
  #               $Bet, $Y, $X, $hessian, $K
  #               When a is NULL, $Phi will be used.
  
  # parlist2      the same as parlist1
  
  # wfun          the weight function. By default w is the identity function
  
  # parinfo1,2  the corresponding info for extraction of the parameters
  
  # parlistw1     alternatively one can use the snp-based est of the censoring surv in each group
  #               as weights. At the moment we don't allow wfun to depend on covariates 
  
  # parlistw2     see parlistw1
  
  # wfix          the weight constant

  # j             the CR of interest
  
  # other parameters are from integrate in package stats  
  
  # created       08 Jan 2014
  # Last modified 04 Jun 2015
  
  
  tryCatch({ # To catch unprecedented errors     
    
    # get the number of CRs
    J  <- parinfo1$J
    l1 <- parinfo1$Gamma * (J-1) + 2*J + sum(parinfo1$Phi) + sum(parinfo1$Bet)
    l2 <- parinfo2$Gamma * (J-1) + 2*J + sum(parinfo2$Phi) + sum(parinfo2$Bet)
    
    lw1 <- lw2 <- 0
    pw1 <- covmatw1 <- pw2 <- covmat2 <- NULL    
    
    if (is.null(parlistw1) | is.null(parlistw2)) { # if w is not snp-based and not estimated from the data
      
      iwd <- function(p) {# this version is for fixed weight function w
        # p has the parameters from 2 groups      
        p1 <- p[1:l1]
        p2 <- p[(l1+1):(l1+l2)]      
        
        tmp.pl1 <- parvect2parlist(params=p1, parinfo=parinfo1, exp.sig=F)
        tmp.pl2 <- parvect2parlist(params=p2, parinfo=parinfo2, exp.sig=F)
        
        re  <- integrate(f=function(tt) {
          (snp.cif(tt, theta=c(tmp.pl1$Phi[[j]], tmp.pl1$Mu_Sig[j,]), 
                   parlist1$base.dens, tmp.pl1$Bet[[j]], parlist1$Y, 
                   exp(logprob(X=parlist1$X, tmp.pl1$Gamma))[1,j], 
                   mod="AFT") -
             snp.cif(tt, theta=c(tmp.pl2$Phi[[j]], tmp.pl2$Mu_Sig[j,]), 
                     parlist2$base.dens, tmp.pl2$Bet[[j]], parlist2$Y, 
                     exp(logprob(X=parlist2$X, tmp.pl2$Gamma))[1,j], 
                     mod="AFT")) * wfun(tt) 
        },  
        lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol, 
        stop.on.error=stop.on.error)
        
        re$value * wfix
      } # end of iwd
      
    } else { # if w is snp-based censoring model
      lw1 <- 2 + length(parlistw1$Phi)
      iwd <- function(p) {# this version is for fixed weight function w
        # p has the parameters from 2 groups      
        p1 <- p[1:l1]
        p2 <- p[(l1+1):(l1+l2)]      
        pw1 <- p[(l1+l2+1):(l1+l2+lw1)]      
        pw2 <- p[(l1+l2+lw1+1):(l1+l2+lw1+lw2)]      
        
        tmp.pl1 <- parvect2parlist(params=p1, parinfo=parinfo1, exp.sig=F)
        tmp.pl2 <- parvect2parlist(params=p2, parinfo=parinfo2, exp.sig=F)
        
        re  <- integrate(f=function(tt) {
          C1 <- 1 - snp.cif(tt, theta=c(parlistw1$Phi, parlistw1$Mu_Sig), parlistw1$base.dens, 
                            NULL, NULL, P=1)
          C2 <- 1 - snp.cif(tt, theta=c(parlistw2$Phi, parlistw2$Mu_Sig), parlistw2$base.dens, 
                            NULL, NULL, P=1)
          
          (snp.cif(tt, theta=c(tmp.pl1$Phi[[j]], tmp.pl1$Mu_Sig[j,]), 
                   parlist1$base.dens, tmp.pl1$Bet[[j]], parlist1$Y, 
                   exp(logprob(X=parlist1$X, tmp.pl1$Gamma))[1,j], 
                   mod="AFT") -
             snp.cif(tt, theta=c(tmp.pl2$Phi[[j]], tmp.pl2$Mu_Sig[j,]), 
                     parlist2$base.dens, tmp.pl2$Bet[[j]], parlist2$Y, 
                     exp(logprob(X=parlist2$X, tmp.pl2$Gamma))[1,j], 
                     mod="AFT")) * 
            (C1*C2) / (p1*C1 + p2*C2)
        },  
        lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol, 
        stop.on.error=stop.on.error)
        
        re$value * wfix
      } # end of iwd
      
    } # end of if (is.null(parlistw1) | is.null(parlistw2)) else ...
    
    p1 <- c(as.vector(t(parlist1$Gamma)), as.vector(t(parlist1$Mu_Sig)), 
            unlist(parlist1$Phi), unlist(parlist1$Bet))
    
    p2 <- c(as.vector(t(parlist2$Gamma)), as.vector(t(parlist2$Mu_Sig)), 
            unlist(parlist2$Phi), unlist(parlist2$Bet))            
    
    p  <- c(p1, p2, pw1, pw2)
    
    iwd.snp <- iwd(p)
    
    var <- p.val <- NULL
    
    # get the covariance matrix from the hessian matrices or the sandwich version 
    # if Ks are not NULL, the hessians should correspond to tne NEGATIVE logliks
    if (length(parlist1$K)==0) { 
      covmat1  <- try(solve(parlist1$hessian), silent=T)
    } else {
      covmat1  <- try(solve(parlist1$hessian) %*% parlist1$K %*% solve(parlist1$hessian), silent=T)
    } # end of if (is.null(parlist1$K)) else...    
    
    if (length(parlist2$K)==0) { 
      covmat2  <- try(solve(parlist2$hessian), silent=T)
    } else {
      covmat2  <- try(solve(parlist2$hessian) %*% parlist2$K %*% solve(parlist2$hessian), silent=T)
    } # end of if (is.null(parlist2$K)) else...    
    
    if (is.null(parlistw1) | is.null(parlistw2)) {
      
      if (class(covmat1)[1] != "try-error" & class(covmat2)[1] != "try-error") {
        covmat <- bdiag(covmat1, covmat2)  
        jcb    <- jacobian(iwd, p) # delta method              
        var    <- diag(jcb %*% covmat %*% t(jcb)) # this is infact just 1x1 matrix, so diag is the whole matrix
        p.val  <- pnorm(q=-abs(iwd.snp), mean=0, sd=sqrt(var)) * 2
      } # end of if (class(covmat1)[1] != "try-error" & class(covmat2)[1] != "try-error")     
      
    } else {
      
      if (length(parlistw1$K)==0) { 
        covmatw1  <- try(solve(parlistw1$hessian), silent=T)
      } else {
        covmatw1  <- try(solve(parlistw1$hessian) %*% parlistw1$K %*% solve(parlistw1$hessian), silent=T)
      } # end of if (is.null(parlistw1$K)) else...   
      
      if (length(parlistw2$K)==0) { 
        covmatw2  <- try(solve(parlistw2$hessian), silent=T)
      } else {
        covmatw2  <- try(solve(parlistw2$hessian) %*% parlistw2$K %*% solve(parlistw2$hessian), silent=T)
      } # end of if (is.null(parlistw1$K)) else...   
      
      if (class(covmat1)[1] != "try-error" & class(covmat2)[1] != "try-error" & 
           class(covmatw1)[1] != "try-error" & class(covmatw2)[1] != "try-error") {
        covmat <- bdiag(covmat1, covmat2, covmatw1, covmatw2)  
        jcb    <- jacobian(iwd, p) # delta method              
        var    <- diag(jcb %*% covmat %*% t(jcb))
        p.val  <- pnorm(q=-abs(iwd.snp), mean=0, sd=sqrt(var)) * 2
      } # end of if (...)
    } # end of if (is.null(parlistw1) | is.null(parlistw2)) else ...
    
    
    list(iwd=iwd.snp, var=var, p.val=p.val)  
    
  }, error=function(ee) {print("From function: iwd.snp.onecif.old"); ee; browser()})
} # end of iwd.snp.onecif.old
#===============================================================================

iwd.npmle.onecif <- function(etmCIF.re1, etmCIF.re2, j=1, wfun=1, dat1=NULL, dat2=NULL, 
                             useC=F, wfix=1, seed=69, B=100,
                             lower=0, upper, ..., subdivisions=500,
                             rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
                             stop.on.error = TRUE) { # no alternative could be found from survival task view
  # test the same cif coming from 2 different groups
  # this uses integrate and etmCIF results
  
  # etmCIF.re1    result from etmCIF for the first group
  
  # etmCIF.re2    like etmCIF.re1 but for the 2nd group
  
  # j             the CR of interest
  
  # wfun          the weight function. By default
  #               w is the identity function  
  
  # dat1, dat2    each is a matrix having 2 columns corresponding to event times and types
  #               to get km est of the censoring for each group when useC is T, see below
  
  # useC          whether or not one should add to the weight function the censoring from each group, eq(10) Doehler.
  #               If this is F then wkm1 and wkm2 will be disregarded in the calculation of the IWD. But
  #               They are still needed for extracting the censoring times
  
  # wfix          constant weight
  
  # seed          for simulation based p value
  
  # B             the number of resampling times
  
  # other parameters are from function integrate in package stats. Actually integrate is only
  # used when wfun is used. But in any case, upper must be specified. 
  
  # Note that upper is taken to be the minimum of all the maximum times where each 
  # parametric curve has a jump
  
  # Actually here we don't need to calculate the variance. According to PKlein 2007, we simulate the
  # test statistics many times using normal deviates and the est of some elements and get a 
  # Monte Carlo type p-value and maybe 95% CI  
  
  # created       09 Jan 2014
  # Last modified 20 Feb 2014
  
  tryCatch({ # To catch unprecedented errors
    
    set.seed(seed)
    
    J  <- length(summary(etmCIF.re1)[[1]])
    n1 <- nrow(dat1)
    n2 <- nrow(dat2)
    
    p1 <- n1/(n1+n2)
    p2 <- n2/(n1+n2)
    
    # first pool and order the event time b/c the integrand is just a discrete function
    t.ev1 <- summary(etmCIF.re1)[[1]][[j]]$time
    t.ev2 <- summary(etmCIF.re2)[[1]][[j]]$time
    
    max.tev1 <- max(t.ev1)
    max.tev2 <- max(t.ev2)
    
    # not that whether or not the weight is based on censoring, the censoring times are needed anyway
    t.evw1 <- t.evw2 <- max.tevw1 <- max.tevw2 <- NULL
    
    tmpdat1 <- data.frame(dat1)
    colnames(tmpdat1) <- c("ttdead", "evdead")
    tmpdat1[tmpdat1[,2]!=0, 2] <- 1
    tmpdat1[,2] <- (tmpdat1[,2] - 1) * (-1)
    wkm1 <- survfit(Surv(ttdead, evdead, type="right")~1, data=tmpdat1)
    
    t.evw1    <- summary(wkm1)$time  
    max.tevw1 <- max(t.evw1)
    
    
    tmpdat2 <- data.frame(dat2)
    colnames(tmpdat2) <- c("ttdead", "evdead")
    tmpdat2[tmpdat2[,2]!=0, 2] <- 1
    tmpdat2[,2] <- (tmpdat2[,2] - 1) * (-1)  
    wkm2 <- survfit(Surv(ttdead, evdead, type="right")~1, data=tmpdat2)
    
    t.evw2    <- summary(wkm2)$time  
    max.tevw2 <- max(t.evw2)  
    
    t.ev  <- c(0, sort(unique(c(t.ev1, t.ev2, t.evw1, t.evw2))))
    
    c1s <-rep(1, length(t.ev))
    c2s <-rep(p1/(1-p2), length(t.ev)) # note when useC is F, this way will give W=1 in the iwd calculation
    
    if (useC) {
      c1s <- summary(wkm1, times=t.ev)$surv
      c2s <- summary(wkm2, times=t.ev)$surv
    }
    
    upper <- min(upper, max.tev1, max.tev2, max.tevw1, max.tevw2)
    
    iwd.npmle <- 0
    
    if (!is.function(wfun)) {
      # Unlike the snp case, here if wfun is a constant function then
      # one should explicitly specify it as a constant and we don't have to use integrate
      # which can save time
      for (i in 1:(length(t.ev)-1)) {
        # c1  <- summary(wkm1, times=t.ev[i])$surv
        # c2  <- summary(wkm2, times=t.ev[i])$surv
        c1  <- c1s[i]
        c2  <- c2s[i]
        tmp <- (all.cif(t.ev[i], etmCIF.re1, p=NULL, distname="nonparametric", j=j) - 
                  all.cif(t.ev[i], etmCIF.re2, p=NULL, distname="nonparametric", j=j)) *
          (c1*c2) / (c1*p1 + c2*p2) * (t.ev[i+1] - t.ev[i])       
        
        iwd.npmle <- iwd.npmle + tmp
        
      } # end of for (i in 1:(length(t.ev)-1))
      
      iwd.npmle <- iwd.npmle * wfun
      
    } else {
      
      for (i in 1:(length(t.ev)-1)) {
        # c1  <- summary(wkm1, times=t.ev[i])$surv
        # c2  <- summary(wkm2, times=t.ev[i])$surv
        c1  <- c1s[i]
        c2  <- c2s[i]
        tmp <- (all.cif(t.ev[i], etmCIF.re1, p=NULL, distname="nonparametric", j=j) - 
                  all.cif(t.ev[i], etmCIF.re2, p=NULL, distname="nonparametric", j=j)) *
          (c1*c2) / (c1*p1 + c2*p2) * 
          integrate(f=wfun, lower=t.ev[i], upper=t.ev[i+1], ..., subdivisions=subdivisions/(length(t.ev)-1), 
                    rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=stop.on.error)$value    
        
        iwd.npmle <- iwd.npmle + tmp
        
      } # end of for (i in 1:(length(t.ev)-1))    
      
    } # end of if (!is.function(wfun)) else ...
    
    iwd.npmle <- iwd.npmle * wfix # * sqrt(n1*n2/(n1+n2))
    
    # Actually here we don't need to calculate the variance. According to PKlein 2007, we simulate the
    # test statistics many times using normal deviates and the est of some elements and get a 
    # Monte Carlo type p-value and maybe 95% CI
    
    # t.ev <- c(t.ev[t.ev < upper]) # should not include upper => why??
    
    # ns # should be a vector having the numbers of sbjs in each group (maybe n.risk from summary.etmCIF)
    ns <- c(n1, n2)
    
    # next is to extract for each t in t.ev, the terms in C.0.3 and C.0.4: I, Y, N
    I <- array(NA, dim=c(J, length(t.ev), 2)) # CIF
    Y <- matrix(NA, nrow=2, ncol=length(t.ev)) # At risk process
    N <- list()
    N[[1]] <- array(0, dim=c(length(t.ev), n1, J)) # Counting processes
    N[[2]] <- array(0, dim=c(length(t.ev), n2, J))
    
    I[,1,] <- 0
    for (jj in 1:J) {
      I[jj, -1, 1] <- c(0, summary(etmCIF.re1)[[1]][[jj]]$P)[findInterval(t.ev[-1], c(0, summary(etmCIF.re1)[[1]][[jj]]$time))]
      I[jj, -1, 2] <- c(0, summary(etmCIF.re2)[[1]][[jj]]$P)[findInterval(t.ev[-1], c(0, summary(etmCIF.re2)[[1]][[jj]]$time))]
    } # end of for (jj in 1:J)
    
    Y[1,1] <- summary(etmCIF.re1)[[1]][[1]]$n.risk[1]
    Y[2,1] <- summary(etmCIF.re2)[[1]][[1]]$n.risk[1]
    
    Y[1,-1]<- c(Y[1,1], summary(etmCIF.re1)[[1]][[1]]$n.risk)[findInterval(t.ev[-1], c(0, summary(etmCIF.re1)[[1]][[jj]]$time))]
    Y[2,-1]<- c(Y[2,1], summary(etmCIF.re2)[[1]][[1]]$n.risk)[findInterval(t.ev[-1], c(0, summary(etmCIF.re2)[[1]][[jj]]$time))]
    
    dats <- list(dat1, dat2)
    
    for (g in 1:2) { # right now this is not so efficient any sbj has at most 1 event at only one point so this has a lot of 0s
      # take mmemory and time => should mod
      for (jj in 1:J) {
        
        for (i in 1:ns[g]) {          
          
          for (l in 1:length(t.ev)) {
            
            if (dats[[g]][i,1]==t.ev[l] & dats[[g]][i,2]==jj) {
              N[[g]][l, i, jj] <- 1
              break
            } # end of if (dats[[g]][1,]==t.ev[l] & dat[[g]][,2]==g) 
            
          } # end of for (l in 1:length(t.ev))
          
        } # end of for (i in 1:ns[g])
        
      } # end of for (jj in 1:J)
      
    } # end of for (g in 1:2)
    
    #   # browser()
    #   Dk <- function(t, k) { 
    #     # k is the group index, so far we consider only 2 groups
    #   
    #     # G is a list with 2 elements wrt 2 groups, each is a matrix (nrow=1,J; ncol=1,n_k). Each entry is an stdnorm rv    
    #     re <- 0
    #     
    #     II <- I[j, findInterval(t, t.ev), k] 
    #     
    #     for (l in 2:length(t.ev)) { # note that 
    #       
    #       tmp.t <- t.ev[l]
    #       
    #       if (tmp.t <= t) {
    #         
    #         tmp0 <- (1 - sum(I[(1:J)[-j], l-1, k]) ) 
    #         GN <- rep(NA, J)
    #         for (jj in 1:J) GN[jj] <- G[[k]][jj,] %*% N[[k]][l-1, , jj]
    # 
    #         re <- re + (tmp0 * GN[j] + I[j, l-1, k] * sum(GN[(1:J)[-j]]) - II * sum(GN)) / Y[k, l-1] 
    # 
    #       } else {
    #         break
    #       } # end of if (tmp.t <= t) else ...
    #       
    #     } # end of for (l in 2:lenth(t.ev))
    #     
    #     re 
    #     
    #   } # end of Dk    
    
    
    Dk.piece1 <- function(t, k, G) { 
      # k is the group index, so far we consider only 2 groups
      
      # G is a list with 2 elements wrt 2 groups, each is a matrix (nrow=1,J; ncol=1,n_k). Each entry is an stdnorm rv    
      re1 <- re2 <- 0
      
      l <- findInterval(t, t.ev)
      
      tmp.t <- t.ev[l]
      
      if (tmp.t <= t) {
        
        tmp0 <- (1 - sum(I[(1:J)[-j], l-1, k]) ) 
        GN <- rep(0, J)
        for (jj in 1:J) GN[jj] <- G[[k]][jj,] %*% N[[k]][l-1, , jj]
        
        re1 <- re1 + tmp0 * GN[j] + I[j, l-1, k] * sum(GN[(1:J)[-j]]) 
        re1 <- re1 / Y[k, l-1] 
        re2 <- sum(GN) / Y[k, l-1] 
      } else {
        break
      } # end of if (tmp.t <= t) else ...
      
      list(piece1=re1, GN=re2) 
      
    } # end of Dk.piece1   
    
    # browser()
    sim.iwds <- rep(NA, B)
    
#     cl    <- makeCluster(detectCores(), 
#                          methods=F) #, outfile="parallelerror.txt")
#     
#     sim.iwds <-  parLapply(cl, 1:B, function(b) { # set.seed outside does not really affect
#       G <- list()
#       G[[1]] <- matrix(rnorm(J*n1), nrow=J, ncol=n1)
#       G[[2]] <- matrix(rnorm(J*n2), nrow=J, ncol=n2)
#       
#       #D1.pieces <- unlist(lapply(2:length(t.ev), function(i) Dk.piece(t.ev[i], 1)))
#       #D2.pieces <- unlist(lapply(2:length(t.ev), function(i) Dk.piece(t.ev[i], 2)))
#       
#       #D1s <- cumsum(D1.pieces)
#       #D2s <- cumsum(D2.pieces)
#       tmp.pieces1 <- matrix(unlist(lapply(2:length(t.ev), function(i) Dk.piece1(t.ev[i], 1, G))), ncol=length(t.ev)-1, byrow=F)
#       tmp.pieces2 <- matrix(unlist(lapply(2:length(t.ev), function(i) Dk.piece1(t.ev[i], 2, G))), ncol=length(t.ev)-1, byrow=F)
#       
#       pieces1 <- cumsum(tmp.pieces1[1,])
#       GNs1    <- cumsum(tmp.pieces1[2,])
#       pieces2 <- cumsum(tmp.pieces2[1,])
#       GNs2    <- cumsum(tmp.pieces2[2,])
#       
#       # browser()
#       tmp.iwd<- 0
#       
#       if (!is.function(wfun)) {
#         for (i in 2:length(t.ev)) {
#           # c1  <- summary(wkm1, times=t.ev[i-1])$surv # this can be done outside the loop for faster performance
#           # c2  <- summary(wkm2, times=t.ev[i-1])$surv
#           c1  <- c1s[i-1]
#           c2  <- c2s[i-1]
#           # tmp <- (Dk(t.ev[i], 1) - Dk(t.ev[i], 2)) * (c1*c2) / (c1*p1 + c2*p2) * (t.ev[i] - t.ev[i-1])       
#           # tmp <- (D1s[i-1] - D2s[i-1]) * (c1*c2) / (c1*p1 + c2*p2) * (t.ev[i] - t.ev[i-1])       
#           tmp <- (pieces1[i-1] - I[j, i, 1]*GNs1[i-1] - pieces2[i-1] + I[j, i, 2]*GNs2[i-1]) * 
#             (c1*c2) / (c1*p1 + c2*p2) * (t.ev[i] - t.ev[i-1])       
#           
#           tmp.iwd <- tmp.iwd + tmp   
#         } # end of for (i in 2:length(t.ev))
#         tmp.iwd <- tmp.iwd * wfun
#         
#       } else {
#         for (i in 2:length(t.ev)) {
#           # c1  <- summary(wkm1, times=t.ev[i-1])$surv
#           # c2  <- summary(wkm2, times=t.ev[i-1])$surv
#           c1  <- c1s[i-1]
#           c2  <- c2s[i-1]
#           # tmp <- (Dk(t.ev[i], 1) - Dk(t.ev[i], 2)) * (c1*c2) / (c1*p1 + c2*p2) * 
#           #        integrate(f=wfun, lower=t.ev[i-1], upper=t.ev[i], ..., subdivisions=subdivisions/(length(t.ev)-1), 
#           #                  rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=stop.on.error)$value       
#           tmp <- (pieces1[i-1] - I[j, i, 1]*GNs1[i-1] - pieces2[i-1] + I[j, i, 2]*GNs2[i-1]) * (c1*c2) / (c1*p1 + c2*p2) * 
#             integrate(f=wfun, lower=t.ev[i-1], upper=t.ev[i], ..., subdivisions=subdivisions/(length(t.ev)-1), 
#                       rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=stop.on.error)$value       
#           
#           tmp.iwd <- tmp.iwd + tmp   
#         } # end of for (i in 2:length(t.ev))
#       } # end of if (is.null(wfun)) else ...
#       
#       tmp.iwd <- tmp.iwd * wfix # * sqrt(n1*n2/(n1+n2))
#       
#       tmp.iwd
#     })
#     
#     sim.iwds <- unlist(sim.iwds)
#     stopCluster(cl)
    
    for (b in 1:B) {
        # print(b)
        G <- list()
        G[[1]] <- matrix(rnorm(J*n1), nrow=J, ncol=n1)
        G[[2]] <- matrix(rnorm(J*n2), nrow=J, ncol=n2)
        
        #D1.pieces <- unlist(lapply(2:length(t.ev), function(i) Dk.piece(t.ev[i], 1)))
        #D2.pieces <- unlist(lapply(2:length(t.ev), function(i) Dk.piece(t.ev[i], 2)))
        
        #D1s <- cumsum(D1.pieces)
        #D2s <- cumsum(D2.pieces)
        tmp.pieces1 <- matrix(unlist(lapply(2:length(t.ev), function(i) Dk.piece1(t.ev[i], 1, G))), ncol=length(t.ev)-1, byrow=F)
        tmp.pieces2 <- matrix(unlist(lapply(2:length(t.ev), function(i) Dk.piece1(t.ev[i], 2, G))), ncol=length(t.ev)-1, byrow=F)
        
        pieces1 <- cumsum(tmp.pieces1[1,])
        GNs1    <- cumsum(tmp.pieces1[2,])
        pieces2 <- cumsum(tmp.pieces2[1,])
        GNs2    <- cumsum(tmp.pieces2[2,])
        
        # browser()
        tmp.iwd<- 0
        
        if (!is.function(wfun)) {
          for (i in 2:length(t.ev)) {
            # c1  <- summary(wkm1, times=t.ev[i-1])$surv # this can be done outside the loop for faster performance
            # c2  <- summary(wkm2, times=t.ev[i-1])$surv
            c1  <- c1s[i-1]
            c2  <- c2s[i-1]
            # tmp <- (Dk(t.ev[i], 1) - Dk(t.ev[i], 2)) * (c1*c2) / (c1*p1 + c2*p2) * (t.ev[i] - t.ev[i-1])       
            # tmp <- (D1s[i-1] - D2s[i-1]) * (c1*c2) / (c1*p1 + c2*p2) * (t.ev[i] - t.ev[i-1])       
            tmp <- (pieces1[i-1] - I[j, i, 1]*GNs1[i-1] - pieces2[i-1] + I[j, i, 2]*GNs2[i-1]) * 
                   (c1*c2) / (c1*p1 + c2*p2) * (t.ev[i] - t.ev[i-1])       
            
            tmp.iwd <- tmp.iwd + tmp   
          } # end of for (i in 2:length(t.ev))
          tmp.iwd <- tmp.iwd * wfun
          
        } else {
          for (i in 2:length(t.ev)) {
            # c1  <- summary(wkm1, times=t.ev[i-1])$surv
            # c2  <- summary(wkm2, times=t.ev[i-1])$surv
            c1  <- c1s[i-1]
            c2  <- c2s[i-1]
            # tmp <- (Dk(t.ev[i], 1) - Dk(t.ev[i], 2)) * (c1*c2) / (c1*p1 + c2*p2) * 
            #        integrate(f=wfun, lower=t.ev[i-1], upper=t.ev[i], ..., subdivisions=subdivisions/(length(t.ev)-1), 
            #                  rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=stop.on.error)$value       
            tmp <- (pieces1[i-1] - I[j, i, 1]*GNs1[i-1] - pieces2[i-1] + I[j, i, 2]*GNs2[i-1]) * (c1*c2) / (c1*p1 + c2*p2) * 
                   integrate(f=wfun, lower=t.ev[i-1], upper=t.ev[i], ..., subdivisions=subdivisions/(length(t.ev)-1), 
                             rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=stop.on.error)$value       
            
            tmp.iwd <- tmp.iwd + tmp   
          } # end of for (i in 2:length(t.ev))
        } # end of if (is.null(wfun)) else ...
        
        tmp.iwd <- tmp.iwd * wfix 
        
        sim.iwds[b] <- tmp.iwd
        
      } # end of for (b in 1:B)

    if (iwd.npmle < 0)  p.val <- length(sim.iwds[sim.iwds<=iwd.npmle]) / B * 2
    if (iwd.npmle >= 0) p.val <- length(sim.iwds[sim.iwds>=iwd.npmle]) / B * 2
    # browser()
    list(iwd=iwd.npmle, tmax=upper, p.val=p.val)
  }, error=function(ee) {print("From function: iwd.npmle.onecif"); ee; browser()})  
} # end of iwd.npmle.onecif
#===============================================================================

my.timepoints <- function (w, times) {
  if (!is.null(w$Tests)) 
    w <- w[names(w) != "Tests"]
  ng <- length(w)
  # times <- sort(unique(times))
  nt <- length(times)
  storage.mode(times) <- "double"
  storage.mode(nt) <- "integer"
  ind <- matrix(0, ncol = nt, nrow = ng)
  oute <- matrix(NA, ncol = nt, nrow = ng)
  outv <- oute
  storage.mode(ind) <- "integer"
  slct <- rep(TRUE, ng)
  for (i in 1:ng) {
    if (is.null((w[[i]])$est)) {
      slct[i] <- FALSE
    }
    else {
      z <- .Fortran("tpoi", as.double(w[[i]][[1]]), as.integer(length(w[[i]][[1]])), 
                    ind[i, ], times, nt, PACKAGE = "cmprsk")
      ind[i, ] <- z[[3]]
      oute[i, ind[i, ] > 0] <- w[[i]][[2]][z[[3]]]
      if (length(w[[i]]) > 2) 
        outv[i, ind[i, ] > 0] <- w[[i]][[3]][z[[3]]]
    }
  }
  dimnames(oute) <- list(names(w)[1:ng], as.character(times))
  dimnames(outv) <- dimnames(oute)
  list(est = oute[slct, , drop = FALSE], var = outv[slct, , 
                                                    drop = FALSE])
}
#===============================================================================
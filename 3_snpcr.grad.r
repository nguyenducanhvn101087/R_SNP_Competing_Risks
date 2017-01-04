# VERSION 2.3 SNPCR.GRAD - No SLOG

# ANH NGUYEN DUC - OXFORD UNIVERSITY CLINICAL RESEARCH UNIT

# SUPERVISED BY:  MARCEL WOLBERS
# CREATED:        FEB    22  2013
# LAST MODIFIED:  JUN    04  2015
#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

snpcr.grad.loglik <- function(params, Data, parinfo, Z, Base.dens, obs.res=T) { # All t_m related removed
  tryCatch({ # To catch unprecedented errors
    # final output
    grad.loglik <- rep(NA, length(params))
    
    # extract the parameters
    parList <- parvect2parlist(params=params, parinfo=parinfo, exp.sig=T)  
    Gamma   <- parList$Gamma
    Mu_Sig  <- parList$Mu_Sig
    Phi     <- parList$Phi
    Bet     <- parList$Bet
    
    # Get the number of CRs
    J          <- parinfo$nc
    
    # Get the number of observations
    N          <- nrow(Data)
    
    # Check if we're dealing with interval-censoring only  
    is.IC      <- !ncol(Data) < 5
    
    x <- Z$X
    
    # get the gradient of gamma
    logP <- logprob(x, Gamma)
    P    <- exp(logP)
    Q    <- 1-P
    PQ   <- P*Q  
    N    <- nrow(x)
    
    m <- nrow(Gamma)
    n <- ncol(Gamma)

    # calculate grad gam_l logP_j when l<>j
    # we know that grad gam_l logP_j = grad gam_j logP_l
    if (J > 2) {
      com          <- combn(J, J-1)
      grad.gam_l.Pj<- lapply(1:ncol(com), 
                             FUN=function(i) {
                               
                               -matrix(
                                 unlist(
                                   lapply(1:N, function(r) 
                                     prod(P[r, com[,i]])*x[r,])
                                 ), 
                                 nrow=N, byrow=T) 
                             }
      ) 
      
      grad.gam_l.Pj<- matrix(unlist(grad.gam_l.Pj), nrow=N)                        
      
      # now clone the result 
      tmat <- matrix(1:(J*n), nrow=J, byrow=T)
      tid  <- as.vector(t(tmat[as.vector(com),] )) 
      tid  <- tid[1:(m*n*(J-1))] # if m = J-1 i.e. FSNPD model, we don't consider grad gamma_J....
      grad.gam_l.Pj <- grad.gam_l.Pj[,tid] # OK
      
    } else {
      
      grad.gam_l.Pj <- -matrix(unlist(lapply(1:N, function(r) prod(P[r, 1:2])*x[r,])), nrow=N, byrow=T)
      
      if (m==2) { # if m = J-1 i.e. FSNPD model, we don't consider grad gamma_J....
        grad.gam_l.Pj <- cbind(grad.gam_l.Pj, grad.gam_l.Pj)  
      } 
      
    } # not yet cover the case when there is only one gamma (ISNPND, J=2)
    
    # now it is important to get the correct P in the denominator
    # In case there are  gamma1 gamma2 gamma3 and each gamma_i has 2 element, this
    # function calculates  grad gam1 P2, grad gam 1 P3, 
    #                      grad gam2 P1, grad gam 2 P3,
    #                      grad gam3 P1, grad gam 3 P2
    idp <- rep(unlist(lapply(1:m, FUN=function(i){(1:J)[-i]})), each=n) 
    grad.gam_l.lopj <- grad.gam_l.Pj / P[,idp] 
    
    # calculate grad gam_j logP_j
    grad.gam_j.Pj <- unlist(lapply(1:m, FUN=function(i){
      x * PQ[,i]
    }))
    
    grad.gam_j.Pj <- matrix(grad.gam_j.Pj, nrow=N)    
    grad.gam_j.lopj <- grad.gam_j.Pj / P[, rep(1:m, each=n)]
    
    # final gradient, in correct order of, if m=2
    # function calculates  grad gam1 P1, grad gam 1 P2, 
    #                      grad gam2 P1, grad gam 2 P2
    grad.gam.logP <- matrix(NA, ncol=m*J*n, nrow=N)
    
    mn <- m*n
    
    # first insert grad gam_j logP_j
    idgjpj <- unlist(lapply(1:m, function(i){
      ((i-1)*(J+1)*n+1):((i-1)*(J+1)*n+n)
    })) # OK
    
    grad.gam.logP[,idgjpj] <- grad.gam_j.lopj
    
    grad.gam.logP[,setdiff(1:(J*mn), idgjpj)] <- grad.gam_l.lopj 
    
    Ks <- rep(NA, J)  
    

    if (!is.IC) { # if RC and LT
      TT <- CC <- As <- Bs <- YB <- emYB <- TTL <- CCL <- TTm <- CCm<- DBase <- DBasem <- DBaseL <- 
      Jacob.Phi.A      <- BsICC      <- BsICC1    <- BsCC  <- CCExtra<- list()           
      # each element of As is a set of polycoef and polycoef of P_K^2 for Bs
      # each element of TT is tt * exp(-y%*%bet), and for CC is 
      # (log(TT - mu)/sig) for stdnorm or (TT/exp(mu))^(1/sig) for exp
      # TTL and CCL is the same as CC but for LT time
      # TTm and CCm is also the same as TT and CC but for only ONE t_m, but still
      # each individual has different parameter => TTm and CCm for each CRj is still
      # a vector of length N
      
      # each element of Jacob.phi.a is the jacobian matrix of phi w.r.t a(phi) where
      # a is the vector of polynomial coefficients derived from phi
      
      grad.phi.logf_K  <- grad.phi.F_K  <- grad.phi.F_K.L <- 
      matrix(0, nrow=N, ncol=length(unlist(Phi))) 
      
      grad.mu.logf_K <- grad.sig.logf_K <- grad.mYB.logf_K <- matrix(0, nrow=N, ncol=J) 
      
      
      V    <- Data[,1] # put drop=F here will brings about errors
      Del  <- Data[,2]    
      LTT  <- Data[,3] # in the process of calculating grad.gamma... already used LT
      
      # Likewise, if there's any right-censoring, we precalculate F_k(t_i|z_i)
      F_K.t_i <- matrix(NA, nrow=N, ncol=J)
      
      # Similarly, if there's any left-truncation, we precalculate F_k(l_i|z_i)
      F_K.l_i <- matrix(NA, nrow=N, ncol=J)
      
      # We precalculate f_k(t_i|z_i)
      f_K.t_i <- matrix(NA, nrow=N, ncol=J)
      
      Del0    <- Del==0
      LT0     <- LTT!=0
      
      for (j in 1:J) { # For regression sth is wrong inside this loop
        logf_0V   <- NULL
        logf_0Vz  <- NULL
        F_0V      <- NULL
        F_0Vz     <- NULL
        F_0tm     <- NULL      
        F_0tmz    <- NULL
        
        # This is where a potential problem occrus as Phi is just an empty list !
        # Therefore extracting Phi[[j]] will cause subscript out of bound error 
        # and theta becomes undelcaered !!! To avoid this, we use try      
        theta     <- try(c(Phi[[j]], Mu_Sig[j,]), TRUE)
        if (class(theta)=="try-error") {
          theta   <- Mu_Sig[j,]
        }  
        
        # In case the betas of all causes are NULL we expect an "subscript out of
        # bounds" error. The same happens to y
        bet       <- try(Bet[[j]], TRUE)
        if (class(bet)=="try-error") {
          bet     <- NULL
        } else if (any(is.na(bet))) {
          bet     <- NULL
        }
        
        y         <- try(Z$Y[[j]], TRUE)
        if (class(y)=="try-error") {
          y       <- NULL
        } else if (any(is.na(y))) {
          y     <- NULL
        }
        
        base.dens <- Base.dens[j]            
        
        # These are needed for calculation of dev mu, sigma
        Ks[j]     <- parinfo$phis[j]
        
        if (length(y) > 0 & length(bet) > 0) {
          YB[[j]]   <- y %*% bet
          emYB[[j]] <- exp(-YB[[j]])
          TT[[j]]   <- V * emYB[[j]] 
          TTL[[j]]  <- LTT * emYB[[j]] 
         
        } else {
          TT[[j]]   <- V
          TTL[[j]]  <- LTT
          YB[[j]]   <- emYB[[j]] <- numeric(0) 
        } # end of if (length(y) > 0 & length(bet) > 0) else ...

        if (Base.dens[j]=="stdnorm") {
          CCL[[j]] <- rep(-Inf, N)
          CCm[[j]] <- -CCL[[j]]
          
          if (length(emYB[[j]]) > 0) {
            CC[[j]] <- (log(V)   - YB[[j]] - Mu_Sig[j,1]) / Mu_Sig[j,2] 
            
            if (any(CCL[[j]]>0)) CCL[[j]][CCL[[j]]>0, 1]<- (log(LTT[CCL[[j]]>0]) - YB[[j]][YB[[j]], 1] - Mu_Sig[j,1]) / Mu_Sig[j,2]

          } else {
            CC[[j]] <- (log(V)   - Mu_Sig[j,1]) / Mu_Sig[j,2]   
            
            if (any(CCL[[j]]>0)) CCL[[j]][CCL[[j]]>0, 1]<- (log(LTT[CCL[[j]]>0]) - Mu_Sig[j,1]) / Mu_Sig[j,2]

          } # end of if (length(emYB[[j]]) > 0) else... 
          
          CCExtra[[j]] <- CC[[j]]
          
          grad.mu.logf_K[,j] <- 1/Mu_Sig[j,2] * CC[[j]]
          
          grad.sig.logf_K[,j]<- 1/Mu_Sig[j,2] * (CC[[j]]^2 - 1)
          
          grad.mYB.logf_K[,j]<- -CC[[j]]/Mu_Sig[j,2]
          
        } else {
          CCL[[j]] <- rep(0, N)
          CCm[[j]] <- rep(Inf, N)
          
          CC[[j]] <- ( TT[[j]] / exp(Mu_Sig[j,1])) ^ (1 / Mu_Sig[j,2])

          if (any(CCL[[j]]>0)) CCL[[j]][CCL[[j]]>0, 1]<- ( TTL[[j]][CCL[[j]]>0] / exp(Mu_Sig[j,1])) ^ (1 / Mu_Sig[j,2])

          if (length(emYB[[j]]) > 0) {
            CCExtra[[j]] <- (log(V)   - YB[[j]] - Mu_Sig[j,1]) / Mu_Sig[j,2]                         
          } else {
            CCExtra[[j]] <- (log(V)   - Mu_Sig[j,1]) / Mu_Sig[j,2]               
          } # end of if (length(emYB[[j]]) > 0) else... 
          
          # grad.mu.logf_K
          grad.mu.logf_K[,j] <- 1/Mu_Sig[j,2] * (CC[[j]] - 1)
          
          # grad.sig.logf_K
          grad.sig.logf_K[,j]<- 1/Mu_Sig[j,2] * ((CC[[j]]-1) * CCExtra[[j]] - 1)
          
          # grad.mYB.logf_K
          grad.mYB.logf_K[,j]<- (1/Mu_Sig[j,2] ) - CC[[j]]/Mu_Sig[j,2]
          
        } # end of if (Base.dens[j]=="stdnorm") else..
        
        DBase[[j]]  <- dbase(Base.dens[j], CC[[j]])
        DBasem[[j]] <- dbase(Base.dens[j], CCm[[j]])
        DBaseL[[j]] <- dbase(Base.dens[j], CCL[[j]])
        
        if (Ks[j] > 0) {        
          As[[j]] <- polycofsnp(base.dens=Base.dens[j], phi=Phi[[j]])
          Bs[[j]] <- coef(polynomial(coef=As[[j]])^2) 
          
          # this is a matrix, each row i is C_i^0,...C_i^K
          CCI     <- t(sapply(1:length(CC[[j]]), function(i) CC[[j]][i]^(0:Ks[j])))
          CCI1    <- CCI[, 1:(ncol(CCI)-1), drop=F]
          
          # below are matrices
          BsCC[[j]]  <- CCI %*% As[[j]]
          BsICC[[j]] <- CCI[, -1, drop=F] %*% (As[[j]][-1]* (1:Ks[j]))
          BsICC1[[j]]<- CCI1 %*% (As[[j]][-1] * (1:Ks[j]))
          
          # grad.mu.logf_K, grad.sig.logf_K, grad.mYB.logf_K
          if (Base.dens[j]=="stdnorm") {
            grad.mu.logf_K[,j] <- grad.mu.logf_K[,j] - 2/Mu_Sig[j,2] * 
                                  BsICC1[[j]] / BsCC[[j]]       
          
            grad.sig.logf_K[,j]<- grad.sig.logf_K[,j] - 2/Mu_Sig[j,2] * 
                                  BsICC[[j]] / BsCC[[j]] 
             
            grad.mYB.logf_K[,j]<- grad.mYB.logf_K[,j] + 2/Mu_Sig[j,2] * BsICC1[[j]] / BsCC[[j]]
              
          } else {
            grad.mu.logf_K[,j] <- grad.mu.logf_K[,j] - 2/Mu_Sig[j,2] * 
                                  BsICC[[j]] / BsCC[[j]] 
            
            grad.sig.logf_K[,j]<- grad.sig.logf_K[,j] - 2/Mu_Sig[j,2] * 
                                  BsICC[[j]] / BsCC[[j]] * CCExtra[[j]]
            
            grad.mYB.logf_K[,j]<- grad.mYB.logf_K[,j] + 2/Mu_Sig[j,2] * 
                                  BsICC[[j]] / BsCC[[j]]
            
          } # end of if (Base.dens[j]=="stdnorm") else...
          
          # phi
          tmp.Seq    <- seq(from = 0, to = 2*Ks[j], by = 1)
          
          moment.Seq <- moment(Base.dens[j], tmp.Seq)
          
          A          <- outer(1:(Ks[j]+1), 1:(Ks[j]+1), function(x,y) moment.Seq[x+y-1])
          
          B          <- chol(A)
 
          B.inv      <- chol2inv(B) %*% t(B) 
          
          cos.phi    <- cos(Phi[[j]]) 
          sin.phi    <- sin(Phi[[j]])
          
          co         <- cumprod(cos.phi)

          if (Ks[j] > 1) {    
            
            D_phi_c    <- matrix(0, nrow=Ks[j]+1, ncol=Ks[j])
            diag(D_phi_c[1:Ks[j],]) <- co[1:Ks[j]]
            
            foo1 <- function(c, r) {
              - co[r-1] / cos.phi[c] * sin.phi[c] * sin.phi[r]
            } # end of foo1
            
            for (r in 2:Ks[j]) {
              D_phi_c[r, 1:(r-1)] <- apply(X=matrix(1:(r-1), nrow=1), MARGIN=2 ,FUN=foo1, r=r)
            } # end of for (j in 2:K)
            
            foo2 <- function(c, r) {
              - co[r] / cos.phi[c] * sin.phi[c]
            } # end of foo2
            D_phi_c[Ks[j]+1, ] <- apply(X=matrix(1:Ks[j], nrow=1), MARGIN=2 ,FUN=foo2, r=Ks[j])

          } else {
            
            D_phi_c <- matrix(c(cos.phi[1], -sin.phi[1]), ncol=1)
            
          } # end of if (K > 1) else ...
          
          Jacob.Phi.A[[j]] <- B.inv %*% D_phi_c 
          
          if (j==1) idphi <- 1:parinfo$phis[j]       
          if (j>1)  idphi <- (sum(parinfo$phis[1:(j-1)])+1):(sum(parinfo$phis[1:(j)])) 
          
          # grad phi f_K
              
          grad.a.logf_K <- 2* CCI / as.vector(CCI %*% As[[j]])
          
          grad.phi.logf_K[, idphi] <-  grad.a.logf_K %*% Jacob.Phi.A[[j]]
          
          # grad phi F_K
          d         <- iter.ints1(x=CC[[j]], base.dens=Base.dens[j], k=Ks[j])
          
          grad.a.F_K<- -2 * matrix(unlist(lapply(1:(Ks[j]+1), FUN=function(k){d[,k:(k+Ks[j])] %*% As[[j]]})), nrow=N)
          
          grad.phi.F_K[, idphi] <- grad.a.F_K %*% Jacob.Phi.A[[j]]
          
          # grad phi F_K.L
          if (any(LT0)) {
            dl          <- iter.ints1(x=CCL[[j]][LT0], base.dens=Base.dens[j], k=Ks[j])
            
            grad.a.F_K.L<- -2 * matrix(unlist(lapply(1:(Ks[j]+1), FUN=function(k){dl[,k:(k+Ks[j])] %*% As[[j]]})), 
                                       nrow=length(LTT[LT0]))
            
            grad.phi.F_K.L[LT0, idphi] <- grad.a.F_K.L %*% Jacob.Phi.A[[j]]
          }
          
        } # end of if (Ks[j] > 0)
        
        # Precalculate terms that will be used by both logf_K and S_K
        if (length(bet) == 0 | length(y) == 0) {
          # logf_0V  <- logf_0(V, theta, base.dens) # Using an ifelse here as well 
          # as for the other case turns
          # out to be more time consuming      
          
          # As both beta and y are NULL, S_0tm will be only one number, this will
          # cause trouble when calculating S1. To avoid this, we replicate S_0tm
          # S_0tm    <- S_0(t_m, theta, base.dens)
          # S_0tm    <- rep(S_0tm, N)
          
          F_0V     <- NULL
          F_0Vz    <- NULL
          
          e        <- NULL
          by       <- NULL      
          y        <- NULL
          
        } else { # This is w.r.t. AFT case only
          by       <- YB[[j]]
          e        <- 1/emYB[[j]]             
          
        } # End of if-else                   
        
        F_K.t_i[,j] <- ifelse(Del0,
                              F_K.aft(V, theta, base.dens, bet, y, by, e, 
                                      F_0V, F_0Vz, CC[[j]], As[[j]]),
                              NA)     
        
        F_K.l_i[,j] <- ifelse(LT0,
                              F_K.aft(LTT, theta, base.dens, bet, y, by, e,
                                      NULL, NULL, CCL[[j]], As[[j]]),
                              NA)           
        
        f_K.t_i[,j] <- ifelse(Del==j, exp(logf_K.aft(V, theta, base.dens, bet, y, by, e,
                                                     logf_0V, logf_0Vz, cc=CC[[j]], a=As[[j]])
        ),
                              NA)
      } # End of for (j in 1:J)
      
      D <- R <- LT <- matrix(0, ncol=ncol(grad.gam.logP), nrow=N)
      
      for (j in 1:J) { # run through CR (not m as in FSNPDN # gamma=J-1)
        
        for (i in 1:m) {
          tmpid1 <- ( (i-1)*J*n+(j-1)*n+1 ) : ( (i-1)*J*n+(j-1)*n+n ) # positions of grad gamma_i logP_j
          
          D[Del==j, tmpid1] <- grad.gam.logP[Del==j, tmpid1, drop=F]# indexing of D is correct
          
        } # end of for (i in 1:m) 
        
      } # end of for (j in 1:J)
      
      if (any(Del0)) { 
        R[Del0, setdiff(1:(J*mn), idgjpj)] <- -grad.gam_l.Pj[Del0,, drop=F]
        R[Del0, idgjpj] <- -grad.gam_j.Pj[Del0, , drop=F]
        R[Del0, ] <- R[Del0, , drop=F] / 
          (1-apply(F_K.t_i[Del0,,drop=F] * P[Del0, 1:J, drop=F], 1, sum))
      }
      
      if (any(LT0)) {
        LT[LT0, setdiff(1:(J*mn), idgjpj)] <- grad.gam_l.Pj[LT0,, drop=F]
        LT[LT0, idgjpj] <- grad.gam_j.Pj[LT0, , drop=F]
        LT[LT0, ] <- LT[LT0, , drop=F] / 
          (1-apply(F_K.l_i[LT0,,drop=F] * P[LT0, 1:J, drop=F], 1, sum))
      }
      
      
      for (j in 1:m) { # here j runs through gamma not CRs as above
        
        for (i in (1:J)) {
          tmpid3 <- ((j-1)*J*n+(i-1)*n+1):((j-1)*J*n+(i-1)*n+n) # here it must drop
          R[Del0, tmpid3] <- R[Del0, tmpid3, drop=F] * F_K.t_i[Del0, i]
          
          LT[LT0, tmpid3] <- LT[LT0, tmpid3, drop=F] * F_K.l_i[LT0, i]
        }
      } # end of for (j in 1:m)
      
      DRLT <- D + R + LT 
      
      grad.gamma.loglk <- matrix(NA, nrow=N, ncol=mn)
      for (j in 1:m) {
        tmp <- DRLT[, ((j-1)*J*n+1):((j-1)*J*n+J*n)]
        for (i in 1:n)
          grad.gamma.loglk[, ((j-1)*n+i)] <- apply(tmp[, seq(i, (J-1)*n+i, by=n)], MARGIN=1, sum)
      } 
      
      grad.gamma.loglk.vec <- apply(grad.gamma.loglk, MARGIN=2, sum)
      
      
      # Calculate grad.mu.loglk, grad.sig.loglk (or grad.logsig.loglk)
      
      # assume J=2, the columns are dev mu_1 f_K1, dev mu_2 f_K2, 
      grad.mu.F_K   <- grad.mu.F_K.L   <- 
      grad.sig.F_K  <- grad.sig.F_K.L  <- 
      grad.mYB.F_K  <- grad.mYB.F_K.L  <- matrix(0, nrow=N, ncol=J) 
      grad.mu.loglk <- grad.sig.loglk  <- grad.mYB.loglk <- matrix(0, nrow=N, ncol=J)
      grad.bet.loglk<- matrix(0, nrow=N, ncol=length(unlist(Bet)))
      grad.phi.loglk<- matrix(0, nrow=N, ncol=length(unlist(Phi)))
      
      for (j in 1:J) {       
        
        # dev mu_j F_Kj.L dev sig_j F_Kj.L and dev emXBj F_Kj.L only for LT cases
        dev.muj.Fj  <- dev.muj.Fj.L  <- 
        dev.sigj.Fj <- dev.sigj.Fj.L <- 
        dev.mYBj.Fj <- dev.mYBj.Fj.L <- 1
        
        if (Ks[j] > 0) {      
          dev.muj.Fj  <- dev.muj.Fj.L  <- 
          dev.sigj.Fj <- dev.sigj.Fj.L <- 
          dev.mYBj.Fj <- dev.mYBj.Fj.L<- 0
          
          for (i in 0:(2*Ks[j])){ # can vectorize this for loop, but doesnt gain much => use loop for ease of comprehension
            dev.muj.Fj   <- dev.muj.Fj + (CC[[j]])^i * Bs[[j]][i+1]
            
            if (any(LT0)) {
              dev.muj.Fj.L   <- dev.muj.Fj.L + (CCL[[j]][LT0])^i * Bs[[j]][i+1]
            } # end of if (any(LT0))
            
          } # end of for (i in 0:(2*Ks[j]))
          
          # first only need to calculate for mu then assign to sig and emXB
          dev.mYBj.Fj  <- dev.sigj.Fj   <- dev.muj.Fj
          dev.mYBj.Fj.L<- dev.sigj.Fj.L <- dev.muj.Fj.L
          
        }         
        
        dev.muj.Fj   <- (-1/Mu_Sig[j,2]) * dev.muj.Fj 
        dev.muj.Fj.L <- (-1/Mu_Sig[j,2]) * dev.muj.Fj.L # avoid ifelse, if no LT, this calculation is harmless
        
        if (Base.dens[j]=="stdnorm") {
          dev.muj.Fj   <- dev.muj.Fj * DBase[[j]] 
          dev.sigj.Fj  <- dev.sigj.Fj * ( (-CC[[j]]) / Mu_Sig[j,2] * DBase[[j]] ) 
          # wo this (...) we plunge into -Inf * 0
          
          if (length(Bet[[j]]) > 0) 
            dev.mYBj.Fj <- dev.mYBj.Fj / Mu_Sig[j,2] * DBase[[j]]
          
          dev.muj.Fj.L   <- dev.muj.Fj.L * dnorm(CCL[[j]][LT0]) 
          dev.sigj.Fj.L  <- dev.sigj.Fj.L * ( (-CCL[[j]][LT0]) / Mu_Sig[j,2] * DBaseL[[j]][LT0] )
          
          if (length(Bet[[j]]) > 0) 
            dev.mYBj.Fj.L <- dev.mYBj.Fj.L / Mu_Sig[j,2] * DBaseL[[j]][LT0]
          
        } else { # exp
          dev.muj.Fj   <- dev.muj.Fj * ( CC[[j]] * DBase[[j]] ) 
          # wo this (...) we plunge into -Inf * 0
          
          dev.sigj.Fj  <- dev.sigj.Fj * ( CC[[j]] / Mu_Sig[j,2]^2 * 
                                          (Mu_Sig[j,1]-log(TT[[j]])) * DBase[[j]] )
          # wo this (...) we plunge into -Inf * 0
          
          if (length(Bet[[j]]) > 0) 
            dev.mYBj.Fj <- dev.mYBj.Fj * ( CC[[j]] / Mu_Sig[j,2] * DBase[[j]] )
          
          dev.muj.Fj.L   <- dev.muj.Fj.L * ( CCL[[j]][LT0] * dexp(CCL[[j]][LT0]) )
          dev.sigj.Fj.L  <- dev.sigj.Fj.L * ( CCL[[j]][LT0] / Mu_Sig[j,2]^2 * 
                                              (Mu_Sig[j,1]-log(TTL[[j]][LT0])) * DBaseL[[j]][LT0] )
          if (length(Bet[[j]]) > 0) 
            dev.mYBj.Fj.L <- dev.mYBj.Fj.L * ( CCL[[j]][LT0] / Mu_Sig[j,2] * DBaseL[[j]][LT0] )
          
        } # end of if (Base.dens[j]=="stdnorm") else ...
        
        
        grad.mu.F_K[,j]       <- dev.muj.Fj
        grad.sig.F_K[,j]      <- dev.sigj.Fj
        if (length(Bet[[j]]) > 0) grad.mYB.F_K[,j]     <- dev.mYBj.Fj
        
        grad.mu.F_K.L[LT0, j]   <- dev.muj.Fj.L
        grad.sig.F_K.L[LT0, j]  <- dev.sigj.Fj.L
        if (length(Bet[[j]]) > 0) grad.mYB.F_K.L[LT0, j] <- dev.mYBj.Fj.L
             
        
        
        grad.mu.loglk[Del==j, j] <- grad.mu.logf_K[Del==j, j, drop=F]
        
        grad.mu.loglk[Del0, j] <- -P[Del0, j, drop=F] * 
                                  grad.mu.F_K[Del0, j, drop=F] /
                                  (1 - apply(F_K.t_i[Del0, 1:J, drop=F] *
                                               P[Del0, 1:J, drop=F], 
                                             1, sum)
                                  )
        
        grad.mu.loglk[LT0, j] <- P[LT0, j, drop=F] *                                
                                grad.mu.F_K.L[LT0, j, drop=F] /
                                (1 - apply(F_K.l_i[LT0, 1:J, drop=F] *
                                             P[LT0, 1:J, drop=F], 
                                           1, sum)
                                ) + grad.mu.loglk[LT0, j, drop=F]
        
        # sig
        grad.sig.loglk[Del==j, j] <- grad.sig.logf_K[Del==j, j, drop=F]
        
        grad.sig.loglk[Del0, j] <- -P[Del0, j, drop=F] * 
                                  grad.sig.F_K[Del0, j, drop=F] /
                                  (1 - apply(F_K.t_i[Del0, 1:J, drop=F] *
                                               P[Del0, 1:J, drop=F], 
                                             1, sum))
        
        grad.sig.loglk[LT0, j] <- P[LT0, j, drop=F] * 
                                  grad.sig.F_K.L[LT0, j, drop=F] /
                                  (1 - apply(F_K.l_i[LT0, 1:J, drop=F] *
                                               P[LT0, 1:J, drop=F], 
                                             1, sum)
                                  ) + grad.sig.loglk[LT0, j, drop=F]
        
        # mYB
        if (parinfo$bets[j] > 0) {
          
          grad.mYB.loglk[Del==j, j] <- grad.mYB.logf_K[Del==j, j, drop=F]
          
          grad.mYB.loglk[Del0, j] <- -P[Del0, j, drop=F] * 
                                    grad.mYB.F_K[Del0, j, drop=F] /
                                    (1 - apply(F_K.t_i[Del0, 1:J, drop=F] *
                                                 P[Del0, 1:J, drop=F], 
                                               1, sum))
          
          grad.mYB.loglk[LT0, j] <- P[LT0, j, drop=F] * 
                                    grad.mYB.F_K.L[LT0, j, drop=F] /
                                    (1 - apply(F_K.l_i[LT0, 1:J, drop=F] *
                                                 P[LT0, 1:J, drop=F], 
                                               1, sum)
                                    ) + grad.mYB.loglk[LT0, j, drop=F]
          
          if (j==1) idbet <- 1:parinfo$bets[j]       
          if (j>1)  idbet <- (sum(parinfo$bets[1:(j-1)])+1):(sum(parinfo$bets[1:(j)])) 
          grad.bet.loglk[,idbet] <- grad.mYB.loglk[,j] * -Z$Y[[j]]
          
        } # end of if (parinfo$bets[j] > 0)
        
        if (parinfo$phis[j] > 0) {
          # phi
          if (j==1) iphi <- 1:parinfo$phis[j]       
          if (j>1)  iphi <- (sum(parinfo$phis[1:(j-1)])+1):(sum(parinfo$phis[1:(j)])) 
          
          grad.phi.loglk[Del==j, iphi] <- grad.phi.logf_K[Del==j, iphi, drop=F]
          
          grad.phi.loglk[Del0, iphi] <- -as.vector(P[Del0, j, drop=F]) *
                                        grad.phi.F_K[Del0, iphi, drop=F] /
                                        (1 - apply(F_K.t_i[Del0, 1:J, drop=F] *
                                                     P[Del0, 1:J, drop=F], 
                                                   1, sum)) 
          
          grad.phi.loglk[LT0, iphi] <- as.vector(P[LT0, j, drop=F]) * 
                                      grad.phi.F_K.L[LT0, iphi, drop=F] /
                                      (1 - apply(F_K.l_i[LT0, 1:J, drop=F] *
                                                   P[LT0, 1:J, drop=F], 
                                                 1, sum) 
                                      ) + grad.phi.loglk[LT0, iphi, drop=F]
          
        } # end of if (parinfo$phis[j] > 0)
        
        
      } # end of for (j in 1:J)
      
      grad.mu.loglk.vec <- apply(grad.mu.loglk, MARGIN=2, sum) 
      
      #if (!optonsig) 
      grad.sig.loglk <- t(apply(grad.sig.loglk ,MARGIN=1, function(r) r * Mu_Sig[,2]))
      
      grad.sig.loglk.vec  <- apply(grad.sig.loglk, MARGIN=2, sum)
      
      grad.bet.loglk.vec  <- apply(grad.bet.loglk, MARGIN=2, sum)
      
      grad.phi.loglk.vec  <- apply(grad.phi.loglk, MARGIN=2, sum) 
      
    } else { # if only IC 
      LL <- RR <- CCL <- CCR <- CCLT <- As <- Bs <- YB <- emYB <- TTm <- CCm <- TTLT <-
      DBaseL <- DBaseR <- DBasem <- DBaseLT <-
      grad.B.emYB     <- Jacob.Phi.A      <- list() 
      # grad betj loglik <- grad.B.emYB[[j]] * dev exp(-YjBj) loglik(exp(-YjBj))
      # it is expected that for each cRj grad.B.emYB[[j]] is a N x length(Bet[[j]]) matrix
      
      # each element of As is a set of polycoef and polycoef of P_K^2 for Bs
      # each element of TT is tt * exp(-y%*%bet), and for CC is 
      # (log(TT - mu)/sig) for stdnorm or (TT/exp(mu))^(1/sig) for exp
      # TTL and CCL is the same as CC but for LT time
      # TTm and CCm is also the same as TT and CC but for only ONE t_m, but still
      # each individual has different parameter => TTm and CCm for each CRj is still
      # a vector of length N
      
      # each element of Jacob.phi.a is the jacobian matrix of phi w.r.t a(phi) where
      # a is the vector of polynomial coefficients derived from phi
      
      grad.phi.F_K.L  <- grad.phi.F_K.R  <- grad.phi.F_K.LT <- 
      matrix(0, nrow=N, ncol=length(unlist(Phi))) 
      
      L        <- Data[,1]
      R        <- Data[,2]
      Del      <- Data[,3]
      LTT      <- Data[,4]
      
      # For each case, if there's any  RC, we precalculate F_k(t_m|z_i) and 
      # F_k(t_i|z_i)
      F_K.L_i <- F_K.R_i <- F_K.LT_i <- matrix(NA, nrow=N, ncol=J)  
      
      # Either Del==0 or R==Inf can mean RC
      IsRC     <- Del==0 | R==Inf    
      LT0      <- LTT!=0
      
      for (j in 1:J) { 
        F_0L      <- NULL
        F_0Lz     <- NULL
        
        F_0R      <- NULL
        F_0Rz     <- NULL

        # This is where a potential problem occrus as Phi is just an empty list !
        # Therefore extracting Phi[[j]] will cause subscript out of bound error 
        # and theta becomes undelcaered !!! To avoid this, we use try      
        theta     <- try(c(Phi[[j]], Mu_Sig[j,]), TRUE)
        if (class(theta)=="try-error") {
          theta   <- Mu_Sig[j,]
        }  
        
        # In case the betas of all causes are NULL we expect an "subscript out of
        # bounds" error. The same happens to y
        bet       <- try(Bet[[j]], TRUE)
        if (class(bet)=="try-error") {
          bet     <- NULL
        } else if (any(is.na(bet))) {
          bet     <- NULL
        }
        
        y         <- try(Z$Y[[j]], TRUE)
        if (class(y)=="try-error") {
          y       <- NULL
        } else if (any(is.na(y))) {
          y       <- NULL
        }
        
        base.dens <- Base.dens[j]            
        
        # These are needed for calculation of dev mu, sigma
        Ks[j]     <- parinfo$phis[j]
        
        if (length(y) > 0 & length(bet) > 0) {
          YB[[j]]   <- y %*% bet
          emYB[[j]] <- exp(-YB[[j]])
          LL[[j]]   <- L * emYB[[j]] 
          RR[[j]]   <- R * emYB[[j]] 
          TTLT[[j]]  <- LTT * emYB[[j]] 
          grad.B.emYB[[j]] <- -as.vector(emYB[[j]]) * y
        } else {
          LL[[j]]   <- L
          RR[[j]]   <- R
          TTLT[[j]] <- LTT
          YB[[j]]   <- emYB[[j]] <- grad.B.emYB[[j]] <- numeric(0)
        } # end of if (length(y) > 0 & length(bet) > 0) else ...
        
        if (Base.dens[j]=="stdnorm") {
          CCL[[j]] <- (log(LL[[j]])  - Mu_Sig[j,1]) / Mu_Sig[j,2]   
          CCR[[j]] <- (log(RR[[j]])  - Mu_Sig[j,1]) / Mu_Sig[j,2]
          CCLT[[j]]<- (log(TTLT[[j]]) - Mu_Sig[j,1]) / Mu_Sig[j,2]
        } else {
          CCL[[j]] <- ( LL[[j]] / exp(Mu_Sig[j,1])) ^ (1 / Mu_Sig[j,2])
          CCR[[j]] <- ( RR[[j]] / exp(Mu_Sig[j,1])) ^ (1 / Mu_Sig[j,2])
          CCLT[[j]]<- ( TTLT[[j]] / exp(Mu_Sig[j,1])) ^ (1 / Mu_Sig[j,2])
        } # end of if (Base.dens[j]=="stdnorm") else..
        
        DBaseL[[j]] <- dbase(Base.dens[j], CCL[[j]])
        DBaseR[[j]] <- dbase(Base.dens[j], CCR[[j]])
        DBaseLT[[j]]<- dbase(Base.dens[j], CCLT[[j]])
        
        if (Ks[j] > 0) {        
          As[[j]] <- polycofsnp(base.dens=Base.dens[j], phi=Phi[[j]])
          Bs[[j]] <- coef(polynomial(coef=As[[j]])^2) 
          
          tmp.Seq    <- seq(from = 0, to = 2*Ks[j], by = 1)
          
          moment.Seq <- moment(Base.dens[j], tmp.Seq)
          
          A          <- outer(1:(Ks[j]+1), 1:(Ks[j]+1), function(x,y) moment.Seq[x+y-1])
          
          B          <- chol(A)
          
          B.inv      <- chol2inv(B) %*% t(B) 
          
          cos.phi    <- cos(Phi[[j]]) # my.cos(phi)
          sin.phi    <- sin(Phi[[j]])
          
          co         <- cumprod(cos.phi)
          
          if (Ks[j] > 1) {    
            
            D_phi_c    <- matrix(0, nrow=Ks[j]+1, ncol=Ks[j])
            diag(D_phi_c[1:Ks[j],]) <- co[1:Ks[j]]
            
            foo1 <- function(c, r) {
              - co[r-1] / cos.phi[c] * sin.phi[c] * sin.phi[r]
            } # end of foo1
            
            for (r in 2:Ks[j]) {
              D_phi_c[r, 1:(r-1)] <- apply(X=matrix(1:(r-1), nrow=1), MARGIN=2 ,FUN=foo1, r=r)
            } # end of for (j in 2:K)
            
            foo2 <- function(c, r) {
              - co[r] / cos.phi[c] * sin.phi[c]
            } # end of foo2
            D_phi_c[Ks[j]+1, ] <- apply(X=matrix(1:Ks[j], nrow=1), MARGIN=2 ,FUN=foo2, r=Ks[j])
            
          } else {
            
            D_phi_c <- matrix(c(cos.phi[1], -sin.phi[1]), ncol=1)
            
          } # end of if (K > 1) else ...
          
          Jacob.Phi.A[[j]] <- B.inv %*% D_phi_c 
          
          if (j==1) idphi <- 1:parinfo$phis[j]       
          if (j>1)  idphi <- (sum(parinfo$phis[1:(j-1)])+1):(sum(parinfo$phis[1:(j)])) 
          
          
          # grad phi F_K.L
          d.L         <- iter.ints1(x=CCL[[j]], base.dens=Base.dens[j], k=Ks[j])
          
          grad.a.F_K.L<- -2 * matrix(unlist(lapply(1:(Ks[j]+1), FUN=function(k){d.L[,k:(k+Ks[j])] %*% As[[j]]})), nrow=N)
          
          grad.phi.F_K.L[, idphi] <- grad.a.F_K.L %*% Jacob.Phi.A[[j]]
          
          # grad phi F_K.R
          d.R         <- iter.ints1(x=CCR[[j]][!IsRC], base.dens=Base.dens[j], k=Ks[j]) # There everyone is RC then this brings
          # error as d.R will be NULL
          
          grad.a.F_K.R<- -2 * matrix(unlist(lapply(1:(Ks[j]+1), FUN=function(k){d.R[,k:(k+Ks[j])] %*% As[[j]]})), 
                                     nrow=length(CCR[[j]][!IsRC]))
          
          grad.phi.F_K.R[!IsRC, idphi] <- grad.a.F_K.R %*% Jacob.Phi.A[[j]]
          
          # grad phi F_K.L
          if (any(LT0)) {
            dl          <- iter.ints1(x=CCLT[[j]][LT0], base.dens=Base.dens[j], k=Ks[j])
            
            grad.a.F_K.LT<- -2 * matrix(unlist(lapply(1:(Ks[j]+1), FUN=function(k){dl[,k:(k+Ks[j])] %*% As[[j]]})), 
                                        nrow=length(LTT[LT0]))
            
            grad.phi.F_K.LT[LT0, idphi] <- grad.a.F_K.LT %*% Jacob.Phi.A[[j]]
          }
          
        } # end of if (Ks[j] > 0)
        
        # Precalculate terms that will be used by both logf_K and S_K
        if (length(bet) == 0 | length(y) == 0) {
          
          e        <- NULL
          by       <- NULL      
          y        <- NULL
          
        } else { # This is w.r.t. AFT case only
          by       <- YB[[j]]
          e        <- 1/emYB[[j]]             
          
        } # End of if-else                   
        
        F_K.L_i[,j] <- F_K.aft(L, theta, base.dens, bet, y, by, e,
                               F_0L, F_0Lz, CCL[[j]], As[[j]])                            
        
        F_K.R_i[,j] <- ifelse(!IsRC,
                              F_K.aft(R, theta, base.dens, bet, y, by, e,
                                      F_0R, F_0Rz, CCR[[j]], As[[j]]),
                              NA)
        
        F_K.LT_i[,j] <- ifelse(LT0,
                               F_K.aft(LTT, theta, base.dens, bet, y, by, e,
                                       NULL, NULL, CCL[[j]], As[[j]]),
                               NA)   

      } # End of for (j in 1:J)
      
      D <- PR <- LT <- matrix(0, ncol=ncol(grad.gam.logP), nrow=N)
      
      for (j in 1:J) { # run through CR (not m as in FSNPDN # gamma=J-1)
        
        for (i in 1:m) {
          tmpid1 <- ( (i-1)*J*n+(j-1)*n+1 ) : ( (i-1)*J*n+(j-1)*n+n ) # positions of grad gamma_i logP_j
          
          D[Del==j, tmpid1] <- grad.gam.logP[Del==j, tmpid1, drop=F]# indexing of D is correct
          
        } # end of for (i in 1:m) 
        
      } # end of for (j in 1:J) 
      
      if (any(IsRC)) { 
        PR[IsRC, setdiff(1:(J*mn), idgjpj)] <- -grad.gam_l.Pj[IsRC,, drop=F]
        PR[IsRC, idgjpj] <- -grad.gam_j.Pj[IsRC, , drop=F]
        PR[IsRC, ] <- PR[IsRC, , drop=F] / 
          (1-apply(F_K.L_i[IsRC,,drop=F] * P[IsRC, 1:J, drop=F], 1, sum))
      }    
      
      if (any(LT0)) {
        LT[LT0, setdiff(1:(J*mn), idgjpj)] <- grad.gam_l.Pj[LT0,, drop=F]
        LT[LT0, idgjpj] <- grad.gam_j.Pj[LT0, , drop=F]
        LT[LT0, ] <- LT[LT0, , drop=F] / 
          (1-apply(F_K.LT_i[LT0,,drop=F] * P[LT0, 1:J, drop=F], 1, sum))
      }
      
      for (j in 1:m) { # here j runs through gamma not CRs as above
        
        for (i in (1:J)) {
          tmpid3 <- ((j-1)*J*n+(i-1)*n+1):((j-1)*J*n+(i-1)*n+n) # here it must drop
          PR[IsRC, tmpid3] <- PR[IsRC, tmpid3, drop=F] * F_K.L_i[IsRC, i]
          
          LT[LT0, tmpid3] <- LT[LT0, tmpid3, drop=F] * F_K.LT_i[LT0, i]
        } # end of for (i in (1:J))
        
      } # end of for (j in 1:m)
      
      DR <- D + PR + LT# 
      
      grad.gamma.loglk <- matrix(NA, nrow=N, ncol=mn)
      for (j in 1:m) {
        tmp <- DR[, ((j-1)*J*n+1):((j-1)*J*n+J*n)]
        for (i in 1:n)
          grad.gamma.loglk[, ((j-1)*n+i)] <- apply(tmp[, seq(i, (J-1)*n+i, by=n)], MARGIN=1, sum)
      } 
      
      grad.gamma.loglk.vec <- apply(grad.gamma.loglk, MARGIN=2, sum) 
      
      # Calculate grad.mu.loglk, grad.sig.loglk (or grad.logsig.loglk)
      
      # assume J=2, the columns are dev mu_1 f_K1, dev mu_2 f_K2, 
      grad.mu.F_K.L   <- grad.mu.F_K.R   <- grad.mu.F_K.LT <-
      grad.sig.F_K.L  <- grad.sig.F_K.R  <- grad.sig.F_K.LT <- 
      grad.emYB.F_K.L <- grad.emYB.F_K.R <- grad.emYB.F_K.LT <- matrix(0, nrow=N, ncol=J) 
      grad.mu.loglk <- grad.sig.loglk<- grad.emYB.loglk <- matrix(0, nrow=N, ncol=J)
      grad.bet.loglk<- matrix(0, nrow=N, ncol=length(unlist(Bet)))
      grad.phi.loglk<- matrix(0, nrow=N, ncol=length(unlist(Phi)))
      
      for (j in 1:J) {       
        
        # dev mu_j F_Kj.L dev sig_j F_Kj.L and dev emXBj F_Kj.L 
        # dev mu_j F_Kj.R dev sig_j F_Kj.R and dev emXBj F_Kj.R 
        dev.muj.Fj.L  <- dev.muj.Fj.R  <- dev.muj.Fj.LT  <- 
          dev.sigj.Fj.L <- dev.sigj.Fj.R <- dev.sigj.Fj.LT <-
          dev.emYBj.Fj.L<- dev.emYBj.Fj.R<- dev.emYBj.Fj.LT<- 1
        
        if (Ks[j] > 0) {      
          dev.muj.Fj.L  <- dev.muj.Fj.R  <- dev.muj.Fj.LT  <-
            dev.sigj.Fj.L <- dev.sigj.Fj.R <- dev.sigj.Fj.LT <-
            dev.emYBj.Fj.L<- dev.emYBj.Fj.R<- dev.emYBj.Fj.LT<- 0
          
          for (i in 0:(2*Ks[j])){ # can vectorize this for loop, but doesnt gain much => use loop for ease of comprehension
            dev.muj.Fj.L   <- dev.muj.Fj.L + (CCL[[j]])^i * Bs[[j]][i+1]
            
            if (any(!IsRC)) {
              dev.muj.Fj.R   <- dev.muj.Fj.R + (CCR[[j]][!IsRC])^i * Bs[[j]][i+1]
            } # end of if (any(LT0))
            
            if (any(LT0)) {
              dev.muj.Fj.LT   <- dev.muj.Fj.LT + (CCLT[[j]][LT0])^i * Bs[[j]][i+1]
            } # end of if (any(LT0))
            
          } # end of for (i in 0:(2*Ks[j]))
          
          # first only need to calculate for mu then assign to sig and emXB
          dev.emYBj.Fj.L<- dev.sigj.Fj.L <- dev.muj.Fj.L
          dev.emYBj.Fj.R<- dev.sigj.Fj.R <- dev.muj.Fj.R
          #dev.mYBj.Fj.LT<- dev.sigj.Fj.LT <- dev.muj.Fj.LT # suspectedly wrong
          dev.emYBj.Fj.LT<- dev.sigj.Fj.LT <- dev.muj.Fj.LT

          # is no Bet, dev.emYBj... 
          # won't be used, but we 
          # catch this later
          
        } # end of if (Ks[j] > 0)       
        
        dev.muj.Fj.L <- (-1/Mu_Sig[j,2]) * dev.muj.Fj.L 
        dev.muj.Fj.LT<- (-1/Mu_Sig[j,2]) * dev.muj.Fj.LT
        if (any(!IsRC)) dev.muj.Fj.R <- (-1/Mu_Sig[j,2]) * dev.muj.Fj.R 
        
        if (Base.dens[j]=="stdnorm") {
          dev.muj.Fj.L <- dev.muj.Fj.L * DBaseL[[j]]          
          dev.sigj.Fj.L<- dev.sigj.Fj.L * ( (-CCL[[j]]) / Mu_Sig[j,2] * DBaseL[[j]] )
          
          if (length(Bet[[j]]) > 0) {
            dev.emYBj.Fj.L <- dev.emYBj.Fj.L / emYB[[j]] / Mu_Sig[j,2] * DBaseL[[j]]
            dev.emYBj.Fj.LT <- dev.emYBj.Fj.LT / emYB[[j]][LT0] / Mu_Sig[j,2] * DBaseLT[[j]][LT0] 
          }
          
          if (any(!IsRC)) {
            dev.muj.Fj.R   <- dev.muj.Fj.R * dnorm(CCR[[j]][!IsRC]) 
            
            dev.sigj.Fj.R  <- dev.sigj.Fj.R * ( (-CCR[[j]][!IsRC]) / Mu_Sig[j,2] * 
                                                DBaseR[[j]][!IsRC] )
            
            if (length(Bet[[j]]) > 0) 
              dev.emYBj.Fj.R <- dev.emYBj.Fj.R / emYB[[j]][!IsRC] / Mu_Sig[j,2] * DBaseR[[j]][!IsRC]
          } # end of if (any(!IsRC))
          
          dev.muj.Fj.LT   <- dev.muj.Fj.LT * dnorm(CCLT[[j]][LT0]) 
          
          dev.sigj.Fj.LT  <- dev.sigj.Fj.LT * ( (-CCLT[[j]][LT0]) / Mu_Sig[j,2] * DBaseLT[[j]][LT0] )
          
        } else { # exp 
          dev.muj.Fj.L <- dev.muj.Fj.L * ( CCL[[j]] * DBaseL[[j]] )
          
          dev.sigj.Fj.L<- dev.sigj.Fj.L * ( CCL[[j]] / Mu_Sig[j,2]^2 * 
                                            (Mu_Sig[j,1]-log(LL[[j]])) * DBaseL[[j]] )
          
          if (length(Bet[[j]]) > 0) {
            dev.emYBj.Fj.L <- dev.emYBj.Fj.L * ( CCL[[j]] / emYB[[j]] / Mu_Sig[j,2] * 
                                                 DBaseL[[j]] )

            dev.emYBj.Fj.LT <- dev.emYBj.Fj.LT * ( CCLT[[j]][LT0] / emYB[[j]][LT0] / Mu_Sig[j,2] * DBaseLT[[j]][LT0] )
          }
          
          if (any(!IsRC)) {
            dev.muj.Fj.R   <- dev.muj.Fj.R * ( CCR[[j]][!IsRC] * dexp(CCR[[j]][!IsRC]) )
            dev.sigj.Fj.R  <- dev.sigj.Fj.R * ( CCR[[j]][!IsRC] / Mu_Sig[j,2]^2 * 
              (Mu_Sig[j,1]-log(RR[[j]][!IsRC])) * DBaseR[[j]][!IsRC] )
            if (length(Bet[[j]]) > 0) 
              dev.emYBj.Fj.R <- dev.emYBj.Fj.R * ( CCR[[j]][!IsRC] / 
                                                   emYB[[j]][!IsRC] / Mu_Sig[j,2] * 
                                                   DBaseR[[j]][!IsRC] )
          } # end of if (any(!IsRC))
          
          dev.muj.Fj.LT   <- dev.muj.Fj.LT * ( CCLT[[j]][LT0] * dexp(CCLT[[j]][LT0]) )
          dev.sigj.Fj.LT  <- dev.sigj.Fj.LT * ( CCLT[[j]][LT0] / Mu_Sig[j,2]^2 * 
                                                (Mu_Sig[j,1]-log(TTLT[[j]][LT0])) * DBaseLT[[j]][LT0] )
          
        } # end of if (Base.dens[j]=="stdnorm") else ...
        
        
        grad.mu.F_K.L[,j]     <- dev.muj.Fj.L
        grad.sig.F_K.L[,j]    <- dev.sigj.Fj.L
        if (length(Bet[[j]]) > 0) grad.emYB.F_K.L[,j]   <- dev.emYBj.Fj.L
        
        if (any(!IsRC)) {
          grad.mu.F_K.R[!IsRC, j]   <- dev.muj.Fj.R
          grad.sig.F_K.R[!IsRC, j]  <- dev.sigj.Fj.R
          if (length(Bet[[j]]) > 0) grad.emYB.F_K.R[!IsRC, j] <- dev.emYBj.Fj.R                          
        } # end of if (any(!IsRC))
        
        grad.mu.F_K.LT[LT0, j]   <- dev.muj.Fj.LT
        grad.sig.F_K.LT[LT0, j]  <- dev.sigj.Fj.LT

        if (length(Bet[[j]]) > 0) grad.emYB.F_K.LT[LT0, j] <- dev.emYBj.Fj.LT
        
        # mu    
        grad.mu.loglk[Del==j, j] <- (grad.mu.F_K.R[Del==j, j, drop=F] - grad.mu.F_K.L[Del==j, j, drop=F]) / 
          (F_K.R_i[Del==j, j, drop=F] - F_K.L_i[Del==j, j, drop=F])  
        
        grad.mu.loglk[IsRC, j] <- -P[IsRC, j, drop=F] * 
          grad.mu.F_K.L[IsRC, j, drop=F] /
          (1 - apply(F_K.L_i[IsRC, 1:J, drop=F] *
                       P[IsRC, 1:J, drop=F], 
                     1, sum)
          )

        grad.mu.loglk[LT0, j] <- P[LT0, j, drop=F] *                                
          grad.mu.F_K.LT[LT0, j, drop=F] /
          (1 - apply(F_K.LT_i[LT0, 1:J, drop=F] *
                       P[LT0, 1:J, drop=F], 
                     1, sum)
          ) + grad.mu.loglk[LT0, j, drop=F]
        

        # sig
        grad.sig.loglk[Del==j, j] <- (grad.sig.F_K.R[Del==j, j, drop=F] - grad.sig.F_K.L[Del==j, j, drop=F]) /
          (F_K.R_i[Del==j, j, drop=F] - F_K.L_i[Del==j, j, drop=F])
        
        grad.sig.loglk[IsRC, j] <- -P[IsRC, j, drop=F] * 
          grad.sig.F_K.L[IsRC, j, drop=F] /
          (1 - apply(F_K.L_i[IsRC, 1:J, drop=F] *
                       P[IsRC, 1:J, drop=F], 
                     1, sum)
          )
        
        grad.sig.loglk[LT0, j] <- P[LT0, j, drop=F] * 
          grad.sig.F_K.LT[LT0, j, drop=F] /
          (1 - apply(F_K.LT_i[LT0, 1:J, drop=F] *
                       P[LT0, 1:J, drop=F], 
                     1, sum)
          ) + grad.sig.loglk[LT0, j, drop=F]
        
        
        # emYB
        if (parinfo$bets[j] > 0) {
          
          grad.emYB.loglk[Del==j, j] <- (grad.emYB.F_K.R[Del==j, j, drop=F] - grad.emYB.F_K.L[Del==j, j, drop=F]) /
            (F_K.R_i[Del==j, j, drop=F] - F_K.L_i[Del==j, j, drop=F])  
          
          grad.emYB.loglk[IsRC, j] <- -P[IsRC, j, drop=F] * 
            grad.emYB.F_K.L[IsRC, j, drop=F] /
            (1 - apply(F_K.L_i[IsRC, 1:J, drop=F] *
                         P[IsRC, 1:J, drop=F], 
                       1, sum))
          
          # In the old version emYB was mYB and is proven wrong
          grad.emYB.loglk[LT0, j] <- P[LT0, j, drop=F] * 
            grad.emYB.F_K.LT[LT0, j, drop=F] /
            (1 - apply(F_K.LT_i[LT0, 1:J, drop=F] *
                         P[LT0, 1:J, drop=F], 
                       1, sum)
            ) + grad.emYB.loglk[LT0, j, drop=F]
          
          if (j==1) idbet <- 1:parinfo$bets[j]       
          if (j>1)  idbet <- (sum(parinfo$bets[1:(j-1)])+1):(sum(parinfo$bets[1:(j)])) 
          grad.bet.loglk[,idbet] <- grad.emYB.loglk[,j] * grad.B.emYB[[j]]
          
        } # end of if (parinfo$bets[j] > 0)
        
        if (parinfo$phis[j] > 0) {
          # phi
          if (j==1) iphi <- 1:parinfo$phis[j]       
          if (j>1)  iphi <- (sum(parinfo$phis[1:(j-1)])+1):(sum(parinfo$phis[1:(j)])) 
          
          grad.phi.loglk[Del==j, iphi] <- (grad.phi.F_K.R[Del==j, iphi, drop=F] - grad.phi.F_K.L[Del==j, iphi, drop=F]) / 
            as.vector(F_K.R_i[Del==j, j, drop=F] - F_K.L_i[Del==j, j, drop=F])
          
          grad.phi.loglk[IsRC, iphi] <- -as.vector(P[IsRC, j, drop=F]) * 
            grad.phi.F_K.L[IsRC, iphi, drop=F] /
            (1 - apply(F_K.L_i[IsRC, 1:J, drop=F] *
                         P[IsRC, 1:J, drop=F], 
                       1, sum)) 
          
          grad.phi.loglk[LT0, iphi] <- as.vector(P[LT0, j, drop=F]) * 
            grad.phi.F_K.LT[LT0, iphi, drop=F] /
            (1 - apply(F_K.LT_i[LT0, 1:J, drop=F] *
                         P[LT0, 1:J, drop=F], 
                       1, sum) 
            ) + grad.phi.loglk[LT0, iphi, drop=F]
          
        } # end of if (parinfo$phis[j] > 0)
           
      } # end of for (j in 1:J)
      
      grad.mu.loglk.vec <- apply(grad.mu.loglk, MARGIN=2, sum) 
      
      #if (!optonsig) 
      grad.sig.loglk <- t(apply(grad.sig.loglk ,MARGIN=1, function(r) r * Mu_Sig[,2]))
      
      grad.sig.loglk.vec  <- apply(grad.sig.loglk, MARGIN=2, sum)
      
      grad.bet.loglk.vec  <- apply(grad.bet.loglk, MARGIN=2, sum)
      
      grad.phi.loglk.vec  <- apply(grad.phi.loglk, MARGIN=2, sum)
      
    } # end of if (!is.IC) else ...
    
    if (!obs.res) {
      re <- c(grad.gamma.loglk.vec, as.vector(rbind(grad.mu.loglk.vec, grad.sig.loglk.vec)), 
               grad.phi.loglk.vec, grad.bet.loglk.vec) 
    } else {
      grad.musig.loglk <- NULL
      for (j in 1:J) grad.musig.loglk <- cbind(grad.musig.loglk, grad.mu.loglk[,j], grad.sig.loglk[,j])
      
      re <- cbind(grad.gamma.loglk, grad.musig.loglk, grad.phi.loglk, grad.bet.loglk)
    } # end of if (retype="vec") else ...

    re
    
  }, error=function(ee) {print("From function: snpcr.grad.loglik"); ee; browser()})  
  
} # end of snpcr.grad.loglik
#==========================================================================================

snp.crlogLik <- function(params, Data, parinfo, Z, Base.dens, obs.res=T) {
  # The log-likelihood with a fixed base density for each competing risk under 
  # ONLY the AFT model. The case of left-truncation and right-censoring and the 
  # case of only interval-censoring are handled separately. This is based on
  # snp.logLik
  
  # Note that only time independent covariates are dealt with
  
  #------------------------------
  # Caller: snp.sub.cropt
  
  #------------------------------
  # Input:
  
  # - Data        a data frame or matrix of the following form:
  
  #   + For data subject to right-censoring and left-truncation
  #   V   |   Del   |  LT | D
  #   + V is the time points. Delta takes value from 0,...,J; where J is the 
  #     number of competing risks. D is a dummy columns to indicate that we're
  #     not dealing with interval-censoring
  
  #   + For data subject to only interval-censoring
  #   L   |   R   |   Delta
  #   + L, R are the left and right time points. Delta takes value from 0,...,J; 
  #     where J is the number of competing risks. Also when R = infinity or Delta 
  #     equals 0, the case is considered right-censored.
  
  # - params      a vectorized version of Params from snp.cropt
  
  # - parinfo     information about the specific length of each element in params
  
  # - Z           a list having:
  #               + X       - a matrix of n x p_x dimension where n is the total
  #                           number of observations and p_x is the number of 
  #                           columns of Gamma above.
  
  #               + Y       - a sub-list having J matrices. Each matrix is of size  
  #                           N x L corresponding to Bet above.
  
  # - Base.dens   a vector of length J whose elements can be either "stdnorm" or
  #               "exp".

  # - obs.res     logical, default T, indicates if the results should be vector of 
  #               loglikelihood contributions for each obs or a single number
  
  #------------------------------
  # Output:
  
  # - The value of the log-likelihood
  
  #------------------------------
  # References:
  
  # - ZhangDavidian_SmoothRegression_Biometrics2008
  
  # - ZhangDavidian_SmoothRegression_Appendix
  
  # - 'Smooth' inference for survival functions with arbitrarily censored data
  
  # - Unknown Anh's publication
  
  #------------------------------
  # Created:       May 10 2011
  # Last modified: Jun 04 2015
  
  # Note that this function won't work if any of the CR has no observed events
  
  tryCatch({ # To catch unprecedented errors
    
    # Get the number of CRs
    J          <- parinfo$nc
    
    # Get the number of observations
    N          <- nrow(Data)
    
    # Check if we're dealing with interval-censoring only  
    is.IC      <- !ncol(Data) < 5
    
    # Extract parameters from params. Note as we assumed Gamma_0 = 0, the number
    # of rows is only J. This is a J x parinfo$Gamma matrix. When t_m is Inf, we
    # also remove Gamma_0 and instead let Gamma_J be 0. Thus this is only a 
    # (J-1) x parinfo$Gamma matrix.
    parList    <- parvect2parlist(params, parinfo, exp.sig=T) # If
    # we don't optimize on sigmas but their logs => must exponentiate the pars
    # which are log sigmas. Before exp.sig=T
    Gamma      <- parList$Gamma
    Mu_Sig     <- parList$Mu_Sig
    Phi        <- parList$Phi
    Bet        <- parList$Bet
    
    # First calculate P_j(z_i) for all j and i as they'll almost be needed
    logP       <- logprob(Z$X, Gamma, N)  
    P          <- exp(logP)
    
    loglk      <- NA
    
    if(!is.IC) { # In case it's not interval-censoring i.e RC and LT
      V    <- Data[,1] # put drop=F here will brings about errors
      Del  <- Data[,2]    
      LT   <- Data[,3]
      
      # Likewise, if there's any right-censoring, we precalculate F_k(t_i|z_i)
      F_K.t_i <- matrix(NA, nrow=N, ncol=J)
      
      # Similarly, if there's any left-truncation, we precalculate F_k(l_i|z_i)
      F_K.l_i <- matrix(NA, nrow=N, ncol=J)
      
      Del0    <- Del==0
      LT0     <- LT!=0
      
      # S1, S2, and S3 are N vectors having respectively the individual
      # contributions from: exact death, RC and LT    
      S1      <- S2 <- S3 <- rep(0, N)
      
      # Take care of the terms for failure-cases and conveniently get the needed
      # elements in F_K.t_i and F_K.l_i             
      for(j in 1:J) { 
        logf_0V   <- NULL
        logf_0Vz  <- NULL
        F_0V      <- NULL
        F_0Vz     <- NULL
        a         <- NULL
        # This is where a potential problem occrus as Phi is just an empty list !
        # Therefore extracting Phi[[j]] will cause subscript out of bound error 
        # and theta becomes undelcaered !!! To avoid this, we use try      
        theta     <- try(c(Phi[[j]], Mu_Sig[j,]), TRUE)
        if(class(theta)=="try-error") {
          theta   <- Mu_Sig[j,]
          a       <- NULL
        } else {
          if(length(Phi[[j]]) > 0)
            a       <- polycofsnp(Base.dens[j], Phi[[j]])
        }     
        
        # In case the betas of all causes are NULL we expect an "subscript out of
        # bounds" error. The same happens to y
        bet       <- try(Bet[[j]], TRUE)
        if(class(bet)=="try-error") {
          bet     <- NULL
        } else if(any(is.na(bet)) | length(bet)==0) {
          bet     <- NULL
        } 
        
        
        y         <- try(Z$Y[[j]], TRUE)
        if(class(y)=="try-error") {
          y       <- NULL
        } else if(any(is.na(y)) | length(y)==0) {
          y       <- NULL
        } 
        
        base.dens <- Base.dens[j]            
        
        # Precalculate terms that will be used by both logf_K and S_K
        if(length(bet) == 0 | length(y) == 0) {
          # logf_0V  <- logf_0(V, theta, base.dens) # Using an ifelse here as well 
          # as for the other case turns
          # out to be more time consuming      
          
          # As both beta and y are NULL, S_0tm will be only one number, this will
          # cause trouble when calculating S1. To avoid this, we replicate S_0tm
          
          S_0V     <- NULL
          S_0Vz    <- NULL
          
          e        <- NULL
          by       <- NULL      
          y        <- NULL
          
          mu  <- Mu_Sig[j,1]
          sig <- Mu_Sig[j,2]
          
          if(base.dens=="stdnorm") {          
            cc   <- (log(V)-mu)/sig
            ccL  <- (log(LT)-mu)/sig
            
          } else {
            cc   <- (V/exp(mu))^(1/sig)
            ccL  <- (LT/exp(mu))^(1/sig)
            
          } # end of if(base.dens=="stdnorm") else...  
          
        } else { # This is w.r.t. AFT case only, only modifiy here 1200 31 Jan 2013
          by       <- y %*% bet
          e        <- exp(by)
          
          mu  <- Mu_Sig[j,1]
          sig <- Mu_Sig[j,2]
          
          if(base.dens=="stdnorm") {          
            cc   <- (log(V)-by-mu)/sig
            ccL  <- (log(LT)-by-mu)/sig
            
          } else {
            cc   <- (V/exp(mu+by))^(1/sig)
            ccL  <- (LT/exp(mu+by))^(1/sig)
            
          } # end of if(base.dens=="stdnorm") else...                       
          
        } # End of if-else                   
        IDel        <- Del==j
        
        S1[IDel]    <- logf_K.aft(tt=V[IDel], theta=theta, base.dens=base.dens, bet=bet, x=y[IDel,], 
                                  bx=by[IDel,], e=e[IDel,], logf_0tt=logf_0V[IDel], logf_0z=logf_0Vz[IDel], 
                                  S_0tt=NULL, cc[IDel], a=a) + logP[IDel,j]          
        
        F_K.t_i[,j] <- ifelse(Del0,
                              F_K.aft(tt=V, theta=theta, base.dens=base.dens, bet=bet, x=y, 
                                      bx=by, e=e, 
                                      F_0tt=F_0V, F_0z=F_0Vz, cc=cc, a=a),
                              NA)     
        
        F_K.l_i[,j] <- ifelse(LT0,
                              F_K.aft(tt=LT, theta=theta, base.dens=base.dens, bet=bet, x=y, 
                                      bx=by, e=e,
                                      F_0tt=NULL, F_0z=NULL, cc=ccL, a=a),
                              NA)                                                                       
      } # End of for    
      
      S2[Del0] <- log(1-
                        apply(F_K.t_i[Del0,,drop=F] * 
                                P[Del0, 1:J], 
                              1, sum)
      )
      
      if(any(LT0)) {
        
        S3[LT0]  <- log(1-
                          apply(F_K.l_i[LT0,,drop=F] * 
                                  P[LT0, 1:J], 
                                1, sum)
        )
      }
      
      # Note that eventhough sum(S1) + sum(S2) + sum(S3) is faster and 
      # theoretically should give the same result, but many times one of the sums
      # will be too large and the final result ignores the smaller ones
      loglk           <- S1 + S2 - S3
      if(!obs.res) loglk  <- sum(loglk)
      
    } else { # The IC case
      # Actually here our formula hasn't incorporated the case of exact death yet.
      # This requries an extra term like f()^{I(L==R)}. However, when the user
      # decides to use IC model, it's almost impossible to have L==R
      L        <- Data[,1]
      R        <- Data[,2]
      Del      <- Data[,3]
      LT       <- Data[,4]
      
      # For each case, if there's any  RC, we precalculate 
      # F_k(t_i|z_i)
      F_K.L_i  <- matrix(NA, nrow=N, ncol=J)  
      
      # Similarly, if there's any left-truncation, we precalculate F_k(lt_i|z_i)
      F_K.lt_i <- matrix(NA, nrow=N, ncol=J)
      
      # Either Del==0 or R==Inf can mean RC
      IsRC     <- Del==0 | R==Inf    
      LT0      <- LT!=0
      
      # S1 and S2 are N vectors having respectively the individual
      # contributions from: real interval-censored data and RC data
      S1       <- S2 <- S3 <- rep(0, N)    
      
      # Take care of the terms for really interval-censored data and conveniently 
      # get the needed elements in F_K.t_i     
      for(j in 1:J) {
        
        F_0L      <- NULL
        F_0Lz     <- NULL
        
        F_0R      <- NULL
        F_0Rz     <- NULL
        
        F_0tm     <- NULL      
        F_0tmz    <- NULL
        
        a         <- NULL
        
        theta     <- try(c(Phi[[j]], Mu_Sig[j,]), TRUE)
        if(class(theta)=="try-error") {
          theta   <- Mu_Sig[j,]
          a       <- NULL
        } else {
          if(length(Phi[[j]]) > 0)
            a       <- polycofsnp(Base.dens[j], Phi[[j]])
        }
        
        # In case the betas of all causes are NULL we expect an "subscript out of
        # bounds" error. The same happens to y  
        bet       <- try(Bet[[j]], TRUE)
        if(class(bet)=="try-error") {
          bet     <- NULL
        } else if(any(is.na(bet)) | length(bet)==0) {
          bet     <- NULL
        } 
        
        
        y         <- try(Z$Y[[j]], TRUE)
        if(class(y)=="try-error") {
          y       <- NULL
        } else if(any(is.na(y)) | length(y)==0) {
          y       <- NULL
        } 
        
        base.dens <- Base.dens[j] 
        
        if(is.null(bet) | is.null(y)) {    
          y        <- e <- by <- NULL
          
          mu  <- Mu_Sig[j,1]
          sig <- Mu_Sig[j,2]   
          
          if(base.dens=="stdnorm") {                      
            ccL  <- (log(L)-mu)/sig
            ccR  <- (log(R)-mu)/sig
            ccLT <- (log(LT)-mu)/sig
            
          } else {            
            ccL  <- (L/exp(mu))^(1/sig)
            ccR  <- (R/exp(mu))^(1/sig)
            ccLT <- (LT/exp(mu))^(1/sig)
            
          } # end of if(base.dens=="stdnorm") else... 
          
        } else { # This is w.r.t. AFT case only
          by       <- y %*% bet
          e        <- exp(by)
          
          mu  <- Mu_Sig[j,1]
          sig <- Mu_Sig[j,2]             
          
          if(base.dens=="stdnorm") {                      
            ccL  <- (log(L)-by-mu)/sig
            ccR  <- (log(R)-by-mu)/sig
            ccLT <- (log(LT)-by-mu)/sig
            
          } else {            
            
            ccL  <- (L/exp(mu+by))^(1/sig)
            ccR  <- (R/exp(mu+by))^(1/sig)
            ccLT <- (LT/exp(mu+by))^(1/sig)
            
          } # end of if(base.dens=="stdnorm") else... 
          
        } # End of if-else            
        
        Ind         <- Del==j & !IsRC
        
        
        S1[Ind]     <- log(P[Ind,j]) + # used to use my.log
          logS_Kdif.aft(L=L[Ind], R=R[Ind], theta=theta, base.dens=base.dens, bet=bet, 
                        y[Ind,,drop=F], by[Ind,,drop=F], e[Ind,,drop=F], F_0L[Ind], 
                        F_0Lz[Ind], ccL[Ind], ccR[Ind], a)  # used to use my.log
        
        F_K.L_i[,j] <- ifelse(IsRC,
                              F_K.aft(L, theta, base.dens, bet, y, by, e,
                                      F_0L, F_0Lz, ccL, a), 
                              NA)    
        
        F_K.lt_i[,j]<- ifelse(LT0,
                              F_K.aft(tt=LT, theta=theta, base.dens=base.dens, bet=bet, x=y, 
                                      bx=by, e=e,
                                      F_0tt=NULL, F_0z=NULL, cc=ccL, a=a),
                              NA) 
        
      } # End of forloop
      S2[IsRC] <- log(1-
                        apply(F_K.L_i[IsRC,,drop=F] * 
                                P[IsRC, 1:J], 
                              1, sum)
      )  
      
      if(any(LT0)) {
        
        S3[LT0]  <- log(1-
                          apply(F_K.lt_i[LT0,,drop=F] * 
                                  P[LT0, 1:J], 
                                1, sum)
        )
      }
      
      # Note that eventhough sum(S1) + sum(S2) is faster and theoretically should 
      # give the same result, but many times one of the sums will be too large and 
      # the final result ignores the smaller ones
      loglk          <- S1 + S2 - S3
      if(!obs.res) loglk <- sum(loglk)                 
      
    } # End of ifelse
    
    # if(is.na(loglk)) loglk <- -Inf
    
    c(loglk)
    
  }, error=function(ee) {cat("From function: snp.crlogLik"); ee; browser()})   
  
} # End of snp.crlogLik
OurMethodOnePairEta0 <- function(Y, X, Gamma_beta, Gamma_e,
                      eps=1e-7, max.iter=5000, crit.loglik=1e-4){
  ## Initialization
  n = nrow(Y)
  ## Preprocessing
  Y.vec = as.vector(Y)
  # Eigenvalue decomposition for X%*%t(X)
  K = X%*%t(X)
  K = (K + t(K))/2
  K.eig = eigen(K)
  U = K.eig$vectors
  d = K.eig$values
  d = ifelse(d<eps, eps, d)
  UTY = t(U) %*% Y
  UTY1 = UTY[,1]
  UTY2 = UTY[,2]
  
  loglikeli_old = 0
  loglikelis = NULL
  for(iter in 1: max.iter){
    # Generalized eigen-value decomposition
    edsi = 1 / sqrt( diag(Gamma_e) ) # Gamma_e.diag.sqrt.inverse
    Gamma = matrix(c( (edsi[1]^2)*Gamma_beta[1,1],
                 edsi[1]*edsi[2]*Gamma_beta[2,1],
                 edsi[1]*edsi[2]*Gamma_beta[1,2],
                 (edsi[2]^2)*Gamma_beta[2,2] ), 2, 2)
    Gamma.eigen = eigen(Gamma)
    lambda = Gamma.eigen$values
    lambda = ifelse(lambda<eps, eps, lambda)
    Phi = edsi * Gamma.eigen$vectors
    Phi.inv = t(Gamma.eigen$vectors / edsi)
    
    # Calculation of log-likelihood
    temp1 = 1/(lambda[1] * d + 1)
    temp2 = 1/(lambda[2] * d + 1)     
    log.Omega.det = drop(-sum(log(c(temp1, temp2))) + n * sum(log(diag(Gamma_e))))
    PhiUTY.vec = matrix(0, (2*n), 1)
    PhiUTY.vec[1:n, 1] = Phi[1,1]*UTY1 + Phi[2,1]*UTY2
    PhiUTY.vec[(n+1):(2*n), 1] = Phi[1,2]*UTY1 + Phi[2,2]*UTY2
    loglikeli = drop(-1/2 * log.Omega.det - 1/2 * t(PhiUTY.vec) %*% (c(temp1, temp2)*PhiUTY.vec) )
    loglikelis = c(loglikelis, loglikeli)
    
    if(abs((loglikeli-loglikeli_old)/loglikeli_old) < crit.loglik) break
    if(iter>1 & loglikeli < loglikeli_old) message("Likelihood decreasing")
    loglikeli_old = loglikeli
    
    # evaluation of M1 and M2
    M1 = Phi %*% diag(c(sum(d*temp1), sum(d*temp2))) %*% t(Phi)
    M1 = (M1 + t(M1))/2
    
    tr1 = sum(temp1)
    tr2 = sum(temp2)
    M2.diag = c( tr1*(Phi[1,1]^2)+tr2*(Phi[1,2]^2), tr1*(Phi[2,1]^2)+tr2*(Phi[2,2]^2) )

    # evaluation of N1 and N2    
    temp3 = t(U) %*% Y %*% Phi * matrix(c(temp1, temp2), ncol=2, nrow = n, byrow=F)
    N1 = diag(sqrt(d)) %*% temp3 %*% diag(lambda) %*% Phi.inv
    N2 = temp3 %*% Phi.inv
    N11 = t(N1) %*% N1
    N22 = t(N2) %*% N2
    N22.diag = diag(N22)
    
    # update function for Gamma_beta and Gamma_e
    L1 = t(chol(M1))
    L1.inv = solve(L1)
    mid_1 = t(L1) %*% N11 %*% L1
    mid_1 = (mid_1 + t(mid_1))/2
    eig1 = eigen(mid_1)
    eig1.values = ifelse(eig1$values<eps, eps, eig1$values)
    Gamma_beta = t(L1.inv) %*% (eig1$vectors %*% diag(sqrt(eig1.values)) %*% t(eig1$vectors)) %*% L1.inv
    Gamma_beta = (Gamma_beta + t(Gamma_beta))/2
    
    L2.diag = sqrt(M2.diag)
    Gamma_e = diag( sqrt(N22.diag)/L2.diag )
  }
#   print(iter)
#   print(loglikelis)
#   print(Gamma_beta)
#   print(Gamma_e)
  
  Ge.eigen = eigen(Gamma_e)
  Ge.values = Ge.eigen$values
  Ge.values = ifelse(Ge.values<eps, eps, Ge.values)
  Ge.inv = Ge.eigen$vectors %*% diag(1/Ge.values) %*% t(Ge.eigen$vectors)
  list(
    loglikelis = loglikelis,
    loglikeli = loglikeli,
    Gamma_beta = Gamma_beta,
    Gamma_e = Gamma_e,
    precision.pair = Ge.inv
  )
}


OurMethodOnePairRankTwoEta0 <- function(Y, Zij, U, d, XXT, Gamma_beta, Gamma_e,
                           eps=1e-7, max.iter=5000, crit.loglik=1e-4){
  ## Initialization
  n = nrow(Y)
  ## Preprocessing
  Y.vec = as.vector(Y)
  UTZij = t(U) %*% Zij
  UTZijT = t(UTZij)
  UTY = t(U) %*% Y
  UTY1 = UTY[,1]
  UTY2 = UTY[,2]

  d.sqrt = sqrt(d)
  V = UTZij / d.sqrt
  VTV = t(V)%*%V
  I_VTV = diag(2) - VTV # I_2-V^TV
  I_VTV.eigen = eigen(I_VTV)
  I_VTV.values = I_VTV.eigen$values
  I_VTV.values[I_VTV.values<eps] = eps
  I_VTV.sqrt = I_VTV.eigen$vectors %*% (sqrt(I_VTV.values)*t(I_VTV.eigen$vectors))
  VH = V %*% solve(VTV) %*% ( I_VTV.sqrt - diag(2) )

  loglikeli_old = 0
  loglikelis = NULL
  for(iter in 1: max.iter){   
    # Generalized eigen-value decomposition
    edsi = 1 / sqrt( diag(Gamma_e) ) # Gamma_e.diag.sqrt.inverse
    Gamma = matrix(c( (edsi[1]^2)*Gamma_beta[1,1], # Gamma_e^{-1/2} %*% Gamma_beta %*% Gamma_e^{-1/2}
                 edsi[1]*edsi[2]*Gamma_beta[2,1],
                 edsi[1]*edsi[2]*Gamma_beta[1,2],
                 (edsi[2]^2)*Gamma_beta[2,2] ), 2, 2)
    Gamma.eigen = eigen(Gamma)
    lambda = Gamma.eigen$values
    lambda = ifelse(lambda<eps, eps, lambda)
    Phi = edsi * Gamma.eigen$vectors
    Phi.inv = t(Gamma.eigen$vectors / edsi)
    
    # Preprocessing
    l1DI = (lambda[1] * d + 1) # diagonal of (\lambda_1 D+I_n)
    l1DI.inv = 1 / l1DI # diagonal of (\lambda_1 D+I_n)^{-1}
    l1DI.inv.UTZ = l1DI.inv * UTZij # (\lambda_1 D+I_n)^{-1} U^T Z_{ij}
    UTZT.l1DI.inv = t(l1DI.inv.UTZ) # (U^T Z_{ij})^T (\lambda_1 D+I_n)^{-1} 
    UTZT.l1DI.inv.UTZ = UTZijT %*% l1DI.inv.UTZ # (U^T Z_ij)^T (\lambda_1 D+I_n)^{-1} U^T Z_{ij}
    
    l2DI = (lambda[2] * d + 1) # diagonal of (\lambda_2 D+I_n)
    l2DI.inv = 1 / l2DI # diagonal of (\lambda_2 D+I_n)^{-1}
    l2DI.inv.UTZ = l2DI.inv * UTZij # (\lambda_2 D+I_n)^{-1} U^T Z_{ij}
    UTZT.l2DI.inv = t(l2DI.inv.UTZ) # (U^T Z_{ij})^T (\lambda_2 D+I_n)^{-1}
    UTZT.l2DI.inv.UTZ = UTZijT %*% l2DI.inv.UTZ # (U^T Z_ij)^T (\lambda_2 D+I_n)^{-1} U^T Z_{ij}

    ## Get \Omega^{-1}
    # center part: (A+BCD)^{-1} = A^{-1} - A^{-1}B (C^{-1}+DA^{-1}B)^{-1} DA^{-1}
    # [ 1/lambda_1 I_2 - (U^TZ)^T (lambda_1 D+I_n)^{-1} U^TZ ]^{-1}, inverse of 2*2 matrix
    l1.inv = 1/lambda[1]
    l1I_T_1_1 = solve( diag(c(l1.inv, l1.inv)) - UTZT.l1DI.inv.UTZ ) 
#     center.11.inv = diag( l1DI.inv ) + l1DI.inv.UTZ %*% l1I_T_1_1 %*% UTZT.l1DI.inv
    center.11.inv = l1DI.inv.UTZ %*% l1I_T_1_1 %*% UTZT.l1DI.inv

    l2.inv = 1/lambda[2]
    l2I_T_1_1 = solve( diag(c(l2.inv, l2.inv)) - UTZT.l2DI.inv.UTZ )
#     center.22.inv = diag( l2DI.inv ) + l2DI.inv.UTZ %*% l2I_T_1_1 %*% UTZT.l2DI.inv
    center.22.inv = l2DI.inv.UTZ %*% l2I_T_1_1 %*% UTZT.l2DI.inv
    
    UTYPhi.vec1 = Phi[1,1]*UTY1 + Phi[2,1]*UTY2
    UTYPhi.vec2 = Phi[1,2]*UTY1 + Phi[2,2]*UTY2
    
    # Calculation of log-likelihood
    det.1 = det( diag(2)-lambda[1]*UTZT.l1DI.inv.UTZ )
    det.2 = det( diag(2)-lambda[2]*UTZT.l2DI.inv.UTZ )
    Omega.log.det = sum(log(l1DI)) + log(det.1) + sum(log(l2DI)) + log(det.2) + n*sum(log(diag(Gamma_e)))
    loglikeli = -1/2 * Omega.log.det -1/2 *
        drop( t(UTYPhi.vec1)%*%center.11.inv%*%UTYPhi.vec1 +
              t(UTYPhi.vec1)%*%(l1DI.inv*UTYPhi.vec1) +
              t(UTYPhi.vec2)%*%center.22.inv%*%UTYPhi.vec2 +
              t(UTYPhi.vec2)%*%(l2DI.inv*UTYPhi.vec2)
            )
    loglikelis = c(loglikelis, loglikeli)
    
    if(abs((loglikeli-loglikeli_old)/loglikeli_old) < crit.loglik) break
    if(iter>1 & loglikeli < loglikeli_old) message("Likelihood decreasing")
    loglikeli_old = loglikeli

    ## Evaluation of M_\beta and M_\epsilon
    Mb.diag.1 = sum(d*l1DI.inv) + sum(diag( l1I_T_1_1 %*% UTZT.l1DI.inv %*% (d*l1DI.inv.UTZ) )) -
                sum(diag( UTZT.l1DI.inv.UTZ )) - sum(diag( UTZT.l1DI.inv.UTZ %*% l1I_T_1_1 %*% UTZT.l1DI.inv.UTZ ))
    Me.diag.1 = sum(l1DI.inv) + sum(diag( l1I_T_1_1 %*% UTZT.l1DI.inv %*% l1DI.inv.UTZ ))
              
    Mb.diag.2 = sum(d*l2DI.inv) + sum(diag( l2I_T_1_1 %*% UTZT.l2DI.inv %*% (d*l2DI.inv.UTZ) )) -
                sum(diag( UTZT.l2DI.inv.UTZ )) - sum(diag( UTZT.l2DI.inv.UTZ %*% l2I_T_1_1 %*% UTZT.l2DI.inv.UTZ ))
    Me.diag.2 = sum(l2DI.inv) + sum(diag( l2I_T_1_1 %*% UTZT.l2DI.inv %*% (l2DI.inv.UTZ) ))
    
    M.beta = Phi %*% (c(Mb.diag.1, Mb.diag.2)*t(Phi))
    M.beta = (M.beta+t(M.beta))/2

    M.epsilon.diag = c( Me.diag.1*(Phi[1,1]^2) + Me.diag.2*(Phi[1,2]^2), Me.diag.1*(Phi[2,1]^2) + Me.diag.2*(Phi[2,2]^2) )

    # Evaluation of N_\beta and N_\epsilon
    G1 = center.11.inv %*% UTYPhi.vec1 + (l1DI.inv*UTYPhi.vec1)
    G2 = center.22.inv %*% UTYPhi.vec2 + (l2DI.inv*UTYPhi.vec2)
    G = matrix(c(G1, G2), nrow=n, ncol=2)
    
    DGLP = (d.sqrt*G) %*% (lambda*Phi.inv) # D^{1/2} %*% G %*% lambda %*% Phi^{-1}
    N.beta = DGLP + VH %*% (t(V) %*% DGLP)
    N.epsilon = G %*% Phi.inv

    NTN.beta = t(N.beta) %*% N.beta
    NTN.epsilon = t(N.epsilon) %*% N.epsilon
    NTN.epsilon.diag = diag(NTN.epsilon)

    # update function for Gamma_beta and Gamma_e
    L1 = t(chol(M.beta))
    L1.inv = solve(L1)
    mid_1 = t(L1) %*% NTN.beta %*% L1
    mid_1 = (mid_1 + t(mid_1))/2
    eig1 = eigen(mid_1)
    eig1.values = ifelse(eig1$values<eps, eps, eig1$values)
    Gamma_beta = t(L1.inv) %*% (eig1$vectors %*% diag(sqrt(eig1.values)) %*% t(eig1$vectors)) %*% L1.inv
    Gamma_beta = (Gamma_beta + t(Gamma_beta))/2
    
    L2.diag = sqrt(M.epsilon.diag)
    Gamma_e = diag( sqrt(NTN.epsilon.diag)/L2.diag )
  }
#   print(iter)
#   print(loglikelis)
#   print(Gamma_beta)
#   print(Gamma_e)
  
  Ge.eigen = eigen(Gamma_e)
  Ge.values = Ge.eigen$values
  Ge.values = ifelse(Ge.values<eps, eps, Ge.values)
  Ge.inv = Ge.eigen$vectors %*% diag(1/Ge.values) %*% t(Ge.eigen$vectors)
  list(
    loglikelis = loglikelis,
    loglikeli = loglikeli,
    Gamma_beta = Gamma_beta,
    Gamma_e = Gamma_e,
    precision.pair = Ge.inv
  )
}



get_se_rho <- function(Gamma_e,Gamma_e_cov){
  d = c(-1/2 * Gamma_e[1,1]^{-3/2} * Gamma_e[2,2]^{-1/2} * Gamma_e[1,2], 
      -1/2 * Gamma_e[1,1]^{-1/2} * Gamma_e[2,2]^{-3/2} * Gamma_e[1,2],
      Gamma_e[1,1]^{-1/2} * Gamma_e[2,2]^{-1/2})
  d = matrix(d, nrow=1, ncol=3)
  return(sqrt(d %*% Gamma_e_cov %*% t(d)))
}
InferWald <- function(Omega.inv, n, K, Y.vec, Gamma_e){
#   infer.start=Sys.time()###########
  Omega.inv.11 = Omega.inv[1:n, 1:n]
  Omega.inv.12 = Omega.inv[1:n, (n+1):(2*n)]
  Omega.inv.22 = Omega.inv[(n+1):(2*n), (n+1):(2*n)]
  Omega.inv.11K = Omega.inv.11 %*% K
  Omega.inv.12K = Omega.inv.12 %*% K
  Omega.inv.22K = Omega.inv.22 %*% K
  
  term0 = 1/2*diag(2*n) - Omega.inv %*% Y.vec %*% t(Y.vec)
  terms = replicate(6, matrix(0, 2*n, 2*n), simplify=F)
  terms[[1]][1:n, 1:n] = Omega.inv.11K
  terms[[1]][(n+1):(2*n), 1:n] = Omega.inv.12K
  terms[[2]][1:n, 1:n] = Omega.inv.11 
  terms[[2]][(n+1):(2*n), 1:n] = Omega.inv.12
  terms[[3]][1:n, (n+1):(2*n)] = Omega.inv.12K 
  terms[[3]][(n+1):(2*n), (n+1):(2*n)] = Omega.inv.22K
  terms[[4]][1:n, (n+1):(2*n)] = Omega.inv.12  
  terms[[4]][(n+1):(2*n), (n+1):(2*n)] = Omega.inv.22 
  terms[[5]][1:n, 1:n] = Omega.inv.12K   
  terms[[5]][1:n, (n+1):(2*n)] = Omega.inv.11K
  terms[[5]][(n+1):(2*n), 1:n] = Omega.inv.22K
  terms[[5]][(n+1):(2*n), (n+1):(2*n)] = Omega.inv.12K
  terms[[6]][1:n, 1:n] = Omega.inv.12  
  terms[[6]][1:n, (n+1):(2*n)] = Omega.inv.11 
  terms[[6]][(n+1):(2*n), 1:n] = Omega.inv.22 
  terms[[6]][(n+1):(2*n), (n+1):(2*n)] = Omega.inv.12
  
  # Fisher information matrix
  FI = matrix(0, 6, 6)
  for(i in 1:6){
    for(j in 1:i){
      FI[i,j] = FI[j,i] = -sum(diag(term0 %*% terms[[i]] %*% terms[[j]]))  
    }
  }
  eta = drop(Gamma_e[1,2])
  rho = drop(Gamma_e[1,2] / sqrt(Gamma_e[1,1]*Gamma_e[2,2]))
  COV = tryCatch({
      solve(FI)  
    }, error=function(e){
      message("error:",e)
      return(NA)
    }, warning=function(w){
      message("warning:",w)
      return(NA)
    }
  )
  
  if(any(is.na(COV)) | COV[6,6] < 0){
    se.eta = NA
    pval.eta = 1
    se.rho = NA
    pval.rho = 1
  }
  else{
    se.eta = drop(sqrt(COV[6,6]))
    pval.eta = pchisq((eta^2)/(se.eta^2), 1, lower.tail=F) # Wald test
    se.rho = get_se_rho(Gamma_e, COV[c(2,4,6), c(2,4,6)])
    pval.rho = pchisq((rho^2)/(se.rho^2), 1, lower.tail=F)
  }
#   infer.end=Sys.time()###########
#   cat('Inference: ', infer.end-infer.start, '\n')###########
  
  list(
      COV = COV,
      eta = eta,
      se.eta = se.eta,
      pval.eta = pval.eta,
      rho = rho,
      se.rho = se.rho,
      pval.rho = pval.rho
  )
}




OurMethodOnePair <- function(Y, X, Gamma_beta=NULL, Gamma_e=NULL,
                    infer='llr', fix.eta=FALSE, eps=1e-7, max.iter=5000, crit.loglik=1e-4){
  ## Initialization
  n = nrow(Y)
  if(is.null(Gamma_beta)) Gamma_beta = cov(Y)/2
  if(is.null(Gamma_e)) Gamma_e = cov(Y)/2
  
  ## Preprocessing
  Y.vec = as.vector(Y)
  # Eigenvalue decomposition for X%*%t(X)
  K = X%*%t(X)
  K = (K + t(K))/2
  K.eig = eigen(K)
  U = K.eig$vectors
  d = K.eig$values
  d = ifelse(d<eps, eps, d)
  UTY = t(U) %*% Y
  UTY1 = UTY[,1]
  UTY2 = UTY[,2]
  
  loglikeli_old = 0
  loglikelis = NULL
  for(iter in 1: max.iter){
    # Generalized eigen-value decomposition
    Gamma_e.eig = eigen(Gamma_e)
    Gamma_e.eig$values = ifelse(Gamma_e.eig$values<eps, eps, Gamma_e.eig$values)

    Gamma_e_sqrt = Gamma_e.eig$vectors %*% (sqrt(Gamma_e.eig$values) * t(Gamma_e.eig$vectors) )
    Gamma_e_sqrt_inv = Gamma_e.eig$vectors %*% (sqrt(1/Gamma_e.eig$values) * t(Gamma_e.eig$vectors) )
    G = Gamma_e_sqrt_inv %*% Gamma_beta %*% Gamma_e_sqrt_inv
    
    eig_Gamma = eigen(G)
    lambda = eig_Gamma$values
    lambda = ifelse(lambda<eps, eps, lambda)
    Phi = Gamma_e_sqrt_inv %*% eig_Gamma$vectors
    Phi.inv = t(eig_Gamma$vectors) %*% Gamma_e_sqrt
    
    # Calculation of log-likelihood
    temp1 = 1/(lambda[1] * d + 1)
    temp2 = 1/(lambda[2] * d + 1)     
    log.Omega.det = drop(-sum(log(c(temp1, temp2))) + n * sum(log(Gamma_e.eig$values)))
    PhiUTY.vec = matrix(0, (2*n), 1)
    PhiUTY.vec[1:n, 1] = Phi[1,1]*UTY1 + Phi[2,1]*UTY2
    PhiUTY.vec[(n+1):(2*n), 1] = Phi[1,2]*UTY1 + Phi[2,2]*UTY2
    loglikeli = drop(-1/2 * log.Omega.det - 1/2 * t(PhiUTY.vec) %*% (c(temp1, temp2)*PhiUTY.vec) )
    loglikelis = c(loglikelis, loglikeli)
    
    if(abs((loglikeli-loglikeli_old)/loglikeli_old) < crit.loglik) break
    if(iter>1 & loglikeli < loglikeli_old) message("Likelihood decreasing")
    loglikeli_old = loglikeli
    
    # evaluation of M1 and M2
    M1 = Phi %*% (c(sum(d*temp1), sum(d*temp2)) * t(Phi) )
    M1 = (M1 + t(M1))/2
    M2 = Phi %*% (c(sum(temp1), sum(temp2)) * t(Phi) )
    M2 = (M2 + t(M2))/2
    

    # evaluation of N1 and N2    
    temp3 = t(U) %*% Y %*% Phi * matrix(c(temp1, temp2), ncol=2, nrow = n, byrow=F)
    N1 = (sqrt(d) * temp3 ) %*% (lambda * Phi.inv )
    N2 = temp3 %*% Phi.inv
    N11 = t(N1) %*% N1
    N22 = t(N2) %*% N2
    
    # update function for Gamma_beta and Gamma_e
    L1 = t(chol(M1))
    L1.inv = solve(L1)
    mid_1 = t(L1) %*% N11 %*% L1
    mid_1 = (mid_1 + t(mid_1))/2
    eig1 = eigen(mid_1)
    eig1.values = ifelse(eig1$values<eps, eps, eig1$values)
    Gamma_beta = t(L1.inv) %*% (eig1$vectors %*% (sqrt(eig1.values) * t(eig1$vectors)) ) %*% L1.inv
    Gamma_beta = (Gamma_beta + t(Gamma_beta))/2

    L2 = t(chol(M2))
    L2.inv = solve(L2)
    mid_2 = t(L2) %*% N22 %*% L2
    mid_2 = (mid_2 + t(mid_2))/2
    eig2 = eigen(mid_2) 
    eig2.values = ifelse(eig2$values<eps, eps, eig2$values)
    Gamma_e = t(L2.inv) %*% (eig2$vectors %*% (sqrt(eig2.values) * t(eig2$vectors)) ) %*% L2.inv
    Gamma_e = (Gamma_e + t(Gamma_e))/2
  }
#   print(iter)
#   print(loglikelis)
#   print(Gamma_beta)
#   print(Gamma_e)
  
  Ge.eigen = eigen(Gamma_e)
  Ge.values = Ge.eigen$values
  Ge.values = ifelse(Ge.values<eps, eps, Ge.values)
  Ge.inv = Ge.eigen$vectors %*% ( (1/Ge.values) * t(Ge.eigen$vectors) )
  if(infer=='llr'){
    return(
      list(
        loglikelis = loglikelis,
        loglikeli = loglikeli,
        Gamma_beta = Gamma_beta,
        Gamma_e = Gamma_e,
        precision.pair = Ge.inv,
        eta = drop(Gamma_e[1,2]),
        rho = drop(Gamma_e[1,2] / sqrt(Gamma_e[1,1]*Gamma_e[2,2]))
      )
    )
  }
  
  
  ## inference
  K.PhiU = kronecker(Phi, U)
  Omega.inv = K.PhiU %*% (c(temp1, temp2)*t(K.PhiU))
  infer.wald = InferWald(Omega.inv, n, K, Y.vec, Gamma_e)
  
  list(
    loglikelis = loglikelis,
    Gamma_beta = Gamma_beta,
    Gamma_e = Gamma_e,
    precision.pair = Ge.inv,
    COV = infer.wald$COV,
    eta = infer.wald$eta,
    se.eta = infer.wald$se.eta,
    pval.eta = infer.wald$pval.eta,
    rho = infer.wald$rho,
    se.rho = infer.wald$se.rho,
    pval.rho = infer.wald$pval.rho
  )
}




OurMethodOnePairRankTwo <- function(Y, Zij, U, d, XXT, Gamma_beta=NULL, Gamma_e=NULL,
                         infer='llr', fix.eta=FALSE, eps=1e-7, max.iter=5000, crit.loglik=1e-4){
  ## Initialization
  n = nrow(Y)
  if(is.null(Gamma_beta)) Gamma_beta = cov(Y)/2
  if(is.null(Gamma_e)) Gamma_e = cov(Y)/2

  ## Preprocessing
  Y.vec = as.vector(Y)
  UTZij = t(U) %*% Zij
  UTZijT = t(UTZij)
  UTY = t(U) %*% Y
  UTY1 = UTY[,1]
  UTY2 = UTY[,2]

  d.sqrt = sqrt(d)
  V = UTZij / d.sqrt
  VTV = t(V)%*%V
  I_VTV = diag(2) - VTV # I_2-V^TV
  I_VTV.eigen = eigen(I_VTV)
  I_VTV.values = I_VTV.eigen$values
  I_VTV.values[I_VTV.values<eps] = eps
  I_VTV.sqrt = I_VTV.eigen$vectors %*% (sqrt(I_VTV.values)*t(I_VTV.eigen$vectors))
  VH = V %*% solve(VTV) %*% ( I_VTV.sqrt - diag(2) )

  loglikeli_old = 0
  loglikelis = NULL
  
  for(iter in 1: max.iter){   
    # Generalized eigen-decomposition
    Gamma_e.eig = eigen(Gamma_e)
    Gamma_e.eig$values = ifelse(Gamma_e.eig$values<eps, eps, Gamma_e.eig$values)
    Gamma_e = Gamma_e.eig$vectors %*% diag(Gamma_e.eig$values) %*% t(Gamma_e.eig$vectors)

    Gamma_e_sqrt = Gamma_e.eig$vectors %*% diag(sqrt(Gamma_e.eig$values)) %*% t(Gamma_e.eig$vectors)
    Gamma_e_sqrt_inv = Gamma_e.eig$vectors %*% diag(sqrt(1/Gamma_e.eig$values)) %*% t(Gamma_e.eig$vectors)
    Gamma_e_sqrt_inv = (Gamma_e_sqrt_inv + t(Gamma_e_sqrt_inv))/2
    G = Gamma_e_sqrt_inv %*% Gamma_beta %*% Gamma_e_sqrt_inv
    
    eig_Gamma = eigen((G+t(G))/2)
    lambda = eig_Gamma$values
    lambda = ifelse(lambda<eps, eps, lambda)
    Phi = Gamma_e_sqrt_inv %*% eig_Gamma$vectors
    Phi.inv = t(eig_Gamma$vectors) %*% Gamma_e_sqrt
    
    # Preprocessing
    l1DI = (lambda[1] * d + 1) # diagonal of (\lambda_1 D+I_n)
    l1DI.inv = 1 / l1DI # diagonal of (\lambda_1 D+I_n)^{-1}
    l1DI.inv.UTZ = l1DI.inv * UTZij # (\lambda_1 D+I_n)^{-1} U^T Z_{ij}
    UTZT.l1DI.inv = t(l1DI.inv.UTZ) # (U^T Z_{ij})^T (\lambda_1 D+I_n)^{-1} 
    UTZT.l1DI.inv.UTZ = UTZijT %*% l1DI.inv.UTZ # (U^T Z_ij)^T (\lambda_1 D+I_n)^{-1} U^T Z_{ij}
    
    l2DI = (lambda[2] * d + 1) # diagonal of (\lambda_2 D+I_n)
    l2DI.inv = 1 / l2DI # diagonal of (\lambda_2 D+I_n)^{-1}
    l2DI.inv.UTZ = l2DI.inv * UTZij # (\lambda_2 D+I_n)^{-1} U^T Z_{ij}
    UTZT.l2DI.inv = t(l2DI.inv.UTZ) # (U^T Z_{ij})^T (\lambda_2 D+I_n)^{-1}
    UTZT.l2DI.inv.UTZ = UTZijT %*% l2DI.inv.UTZ # (U^T Z_ij)^T (\lambda_2 D+I_n)^{-1} U^T Z_{ij}

    ## Get \Omega^{-1}
    # center part: (A+BCD)^{-1} = A^{-1} - A^{-1}B (C^{-1}+DA^{-1}B)^{-1} DA^{-1}
    # [ 1/lambda_1 I_2 - (U^TZ)^T (lambda_1 D+I_n)^{-1} U^TZ ]^{-1}, inverse of 2*2 matrix
    l1.inv = 1/lambda[1]
    l1I_T_1_1 = solve( diag(c(l1.inv, l1.inv)) - UTZT.l1DI.inv.UTZ ) 
#     center.11.inv = diag( l1DI.inv ) + l1DI.inv.UTZ %*% l1I_T_1_1 %*% UTZT.l1DI.inv
    center.11.inv = l1DI.inv.UTZ %*% l1I_T_1_1 %*% UTZT.l1DI.inv

    l2.inv = 1/lambda[2]
    l2I_T_1_1 = solve( diag(c(l2.inv, l2.inv)) - UTZT.l2DI.inv.UTZ )
#     center.22.inv = diag( l2DI.inv ) + l2DI.inv.UTZ %*% l2I_T_1_1 %*% UTZT.l2DI.inv
    center.22.inv = l2DI.inv.UTZ %*% l2I_T_1_1 %*% UTZT.l2DI.inv
    
    UTYPhi.vec1 = Phi[1,1]*UTY1 + Phi[2,1]*UTY2
    UTYPhi.vec2 = Phi[1,2]*UTY1 + Phi[2,2]*UTY2
    
    # Calculation of log-likelihood
    det.1 = det( diag(2)-lambda[1]*UTZT.l1DI.inv.UTZ )
    det.2 = det( diag(2)-lambda[2]*UTZT.l2DI.inv.UTZ )
    Omega.log.det = sum(log(l1DI)) + log(det.1) + sum(log(l2DI)) + log(det.2) + n*sum(log(Gamma_e.eig$values))
    loglikeli = -1/2 * Omega.log.det -1/2 *
        drop( t(UTYPhi.vec1)%*%center.11.inv%*%UTYPhi.vec1 +
              t(UTYPhi.vec1)%*%(l1DI.inv*UTYPhi.vec1) +
              t(UTYPhi.vec2)%*%center.22.inv%*%UTYPhi.vec2 +
              t(UTYPhi.vec2)%*%(l2DI.inv*UTYPhi.vec2)
            )
    loglikelis = c(loglikelis, loglikeli)
    
    if(abs((loglikeli-loglikeli_old)/loglikeli_old) < crit.loglik) break
    if(iter>1 & loglikeli < loglikeli_old) message("Likelihood decreasing")
    loglikeli_old = loglikeli

    ## Evaluation of M_\beta and M_\epsilon
    Mb.diag.1 = sum(d*l1DI.inv) + sum(diag( l1I_T_1_1 %*% UTZT.l1DI.inv %*% (d*l1DI.inv.UTZ) )) -
                sum(diag( UTZT.l1DI.inv.UTZ )) - sum(diag( UTZT.l1DI.inv.UTZ %*% l1I_T_1_1 %*% UTZT.l1DI.inv.UTZ ))
    Me.diag.1 = sum(l1DI.inv) + sum(diag( l1I_T_1_1 %*% UTZT.l1DI.inv %*% l1DI.inv.UTZ ))
              
    Mb.diag.2 = sum(d*l2DI.inv) + sum(diag( l2I_T_1_1 %*% UTZT.l2DI.inv %*% (d*l2DI.inv.UTZ) )) -
                sum(diag( UTZT.l2DI.inv.UTZ )) - sum(diag( UTZT.l2DI.inv.UTZ %*% l2I_T_1_1 %*% UTZT.l2DI.inv.UTZ ))
    Me.diag.2 = sum(l2DI.inv) + sum(diag( l2I_T_1_1 %*% UTZT.l2DI.inv %*% (l2DI.inv.UTZ) ))
    
    M.beta = Phi %*% (c(Mb.diag.1, Mb.diag.2)*t(Phi))
    M.beta = (M.beta+t(M.beta))/2

    M.epsilon = Phi %*% (c(Me.diag.1, Me.diag.2)*t(Phi))
    M.epsilon = (M.epsilon+t(M.epsilon))/2  

    # Evaluation of N_\beta and N_\epsilon
    G1 = center.11.inv %*% UTYPhi.vec1 + (l1DI.inv*UTYPhi.vec1)
    G2 = center.22.inv %*% UTYPhi.vec2 + (l2DI.inv*UTYPhi.vec2)
    G = matrix(c(G1, G2), nrow=n, ncol=2)
    
    DGLP = (d.sqrt*G) %*% (lambda*Phi.inv) # D^{1/2} %*% G %*% lambda %*% Phi^{-1}
    N.beta = DGLP + VH %*% (t(V) %*% DGLP)
    N.epsilon = G %*% Phi.inv

    NTN.beta = t(N.beta) %*% N.beta
    NTN.epsilon = t(N.epsilon) %*% N.epsilon

    # update function for Gamma_beta and Gamma_e
    L1 = t(chol(M.beta))
    L1.inv = solve(L1)
    mid_1 = t(L1) %*% NTN.beta %*% L1
    mid_1 = (mid_1 + t(mid_1))/2
    eig1 = eigen(mid_1)
    eig1.values = ifelse(eig1$values<eps, eps, eig1$values)
    Gamma_beta = t(L1.inv) %*% (eig1$vectors %*% diag(sqrt(eig1.values)) %*% t(eig1$vectors)) %*% L1.inv
    Gamma_beta = (Gamma_beta + t(Gamma_beta))/2

    L2 = t(chol(M.epsilon))
    L2.inv = solve(L2)
    mid_2 = t(L2) %*% NTN.epsilon %*% L2
    mid_2 = (mid_2 + t(mid_2))/2
    eig2 = eigen(mid_2) 
    eig2.values = ifelse(eig2$values<eps, eps, eig2$values)
    Gamma_e = t(L2.inv) %*% (eig2$vectors %*% diag(sqrt(eig2.values)) %*% t(eig2$vectors)) %*% L2.inv 
    Gamma_e = (Gamma_e + t(Gamma_e))/2

  }
#   print(iter)
#   print(loglikelis)
#   print(Gamma_beta)
#   print(Gamma_e)
  
  Ge.eigen = eigen(Gamma_e)
  Ge.values = Ge.eigen$values
  Ge.values = ifelse(Ge.values<eps, eps, Ge.values)
  Ge.inv = Ge.eigen$vectors %*% diag(1/Ge.values) %*% t(Ge.eigen$vectors)
  if(infer=='llr'){
    return(
      list(
        loglikelis = loglikelis,
        loglikeli = loglikeli,
        Gamma_beta = Gamma_beta,
        Gamma_e = Gamma_e,
        precision.pair = Ge.inv,
        eta = drop(Gamma_e[1,2]),
        rho = drop(Gamma_e[1,2] / sqrt(Gamma_e[1,1]*Gamma_e[2,2]))
      )
    )
  }

  ## inference by wald test based on fisher information
  K = XXT
  U.center.11.inv.UT = U %*% center.11.inv %*% t(U)
  U.center.22.inv.UT = U %*% center.22.inv %*% t(U)
  Omega.inv.11 = (Phi[1,1]^2) * U.center.11.inv.UT + (Phi[1,2]^2) * U.center.22.inv.UT
  Omega.inv.12 = (Phi[1,1]*Phi[2,1]) * U.center.11.inv.UT + (Phi[1,2]*Phi[2,2]) * U.center.22.inv.UT
  Omega.inv.22 = (Phi[2,1]^2) * U.center.11.inv.UT + (Phi[2,2]^2) * U.center.22.inv.UT
  Omega.inv = matrix(0, (2*n), (2*n))
  Omega.inv[1:n, 1:n] = Omega.inv.11
  Omega.inv[1:n, (n+1):(2*n)] = Omega.inv[(n+1):(2*n), 1:n] = Omega.inv.12
  Omega.inv[(n+1):(2*n), (n+1):(2*n)] = Omega.inv.22
  infer.wald = InferWald(Omega.inv, n, K, Y.vec, Gamma_e)

  list(
      loglikelis = loglikelis,
      Gamma_beta = Gamma_beta,
      Gamma_e = Gamma_e,
      precision.pair = Ge.inv,
      COV = infer.wald$COV,
      eta = infer.wald$eta,
      se.eta = infer.wald$se.eta,
      pval.eta = infer.wald$pval.eta,
      rho = infer.wald$rho,
      se.rho = infer.wald$se.rho,
      pval.rho = infer.wald$pval.rho
  )
}




## data: a n*p data matrix
## infer: option of different tests of inference: 'llr' for likelihood ratio test, 'wald' for wald test based on Fisher Information Matrix
## low.rank: whether to use low rank update to shrink the time of eigendecomposition of XX^T, default to be TRUE when sample size larger than 1000
## scale.var: whether to scale the variance of X to 1/p, default to be T(RUE)
## eps: a small term to avoid numerical problems, default to be 1e-7
OurMethod <- function(data, scale.var=T, low.rank=NULL, infer='llr', eps=1e-7, crit.loglik=1e-4){
  all.start=Sys.time()

  ## preprocessing
  data = as.matrix(data)
  N=dim(data)[1]
  P=dim(data)[2]
  data.c = scale(data, center=T, scale=F) # centering: center each column of data to have zero mean
  if(scale.var){
    Z = scale(data.c, center = T, scale = T)
  }
  else{
    Z = data.c
  }
  Z = Z / sqrt(P-2)
  if(is.null(low.rank)){
    low.rank = ifelse(N>=1000, T, F)
  }
  
  ## Eigen-decomposition of ZZ^T for low rank update
  if(low.rank){
    eigen.start = Sys.time()
    ZZT = Z%*%t(Z)
    ZZT = (ZZT+t(ZZT))/2
    ZZT.eigen = eigen(ZZT)
    U = ZZT.eigen$vectors
    d = ZZT.eigen$values
    d = ifelse(d<eps, eps, d)
    eigen.end = Sys.time()
    exe.time = eigen.end-eigen.start
    cat('Eigen: ', exe.time, '\n')
  }

  ## Initialization
  precision.est = matrix(0,P,P)
  precision.est.diag = matrix(0,P,P) # precision.est.diag[i,j] means estimation of precision[i,i] when paired with j
  precision.pval = matrix(0,P,P) # p-value of edges
  beta.var = rep(0, P)
  epsilon.var = rep(0, P)
  eta.est = matrix(0,P,P)
  eta.se = matrix(0,P,P)
  rho.est = matrix(0,P,P)
  rho.se = matrix(0,P,P)
  rho.pval = matrix(0,P,P)
  sigma.a2 = matrix(0,P,P)
  sigma.b2 = matrix(0,P,P)
  
  for(i in 1:(P-1)){
    for(j in (i+1):P){
      ## consider edge between node i and j each time
#       cat("pair: ", i, " and ", j, "\n")
      
      ## Get data for each pair
      Y = data.c[, c(i,j)]
      if(low.rank == TRUE){
        Zij = Z[, c(i,j)]
        XXT = ZZT - Zij%*%t(Zij)
      }
      else{
        X = Z[, -c(i,j)]
      }
      ## Initialization of Gamma
      if(i==1){
        Gamma_beta = cov(Y)/2
        Gamma_e = cov(Y)/2
      }
      else{
        Gamma_beta = matrix(0,2,2)
        Gamma_beta[1,1] = beta.var[i]
        Gamma_beta[2,2] = beta.var[j]
        Gamma_e = matrix(0,2,2)
        Gamma_e[1,1] = epsilon.var[i]
        Gamma_e[2,2] = epsilon.var[j]
      }
      ## For inference based on likelihood ratio test, first estimate when eta = 0
      if(infer == 'llr'){
        Gamma_e[1,2] = Gamma_e[2,1] = 0
        if(low.rank == TRUE){
          pair.eta0 = OurMethodOnePairRankTwoEta0(Y, Zij, U, d, XXT, Gamma_beta, Gamma_e, eps=eps, crit.loglik=crit.loglik)
        }
        else{
          pair.eta0 = OurMethodOnePairEta0(Y, X, Gamma_beta, Gamma_e, eps=eps, crit.loglik=crit.loglik)
        }
        Gamma_beta = pair.eta0$Gamma_beta
        Gamma_e = pair.eta0$Gamma_e
      }
      ## Estimate
      if(low.rank == TRUE){
        pair.result = OurMethodOnePairRankTwo(Y, Zij, U, d, XXT, Gamma_beta, Gamma_e,
                                 infer=infer, fix.eta=FALSE, eps=eps, crit.loglik=crit.loglik)
      }
      else{
        pair.result = OurMethodOnePair(Y, X, Gamma_beta, Gamma_e, infer=infer, fix.eta=FALSE, eps=eps, crit.loglik=crit.loglik)
      }
      precision.est[i,j] = precision.est[j,i] = pair.result$precision.pair[1,2]
      precision.est.diag[i,j] = pair.result$precision.pair[1,1]
      precision.est.diag[j,i] = pair.result$precision.pair[2,2]
      if(i == 1){
        ## Collect variance of Gammas for warm start in following iterations
        beta.var[i] = pair.result$Gamma_beta[1,1]
        beta.var[j] = pair.result$Gamma_beta[2,2]
        epsilon.var[i] = pair.result$Gamma_e[1,1]
        epsilon.var[j] = pair.result$Gamma_e[2,2]
      }
      
      if(infer == 'llr'){
        LLR = -2 * (pair.eta0$loglikeli - pair.result$loglikeli)
        pval.eta = pchisq(LLR, 1, lower.tail = F)
      }
      else if(infer == 'wald'){
        pval.eta = pair.result$pval.eta
        eta.se[i,j] = eta.se[j,i] = pair.result$se.eta
        rho.se[i,j] = rho.se[j,i] = pair.result$se.rho
        rho.pval[i,j] = rho.pval[j,i] = pair.result$pval.rho
      }
      precision.pval[i,j] = precision.pval[j,i] = pval.eta
      eta.est[i,j] = eta.est[j,i] = pair.result$eta
      rho.est[i,j] = rho.est[j,i] = pair.result$rho
      sigma.a2[i,j] = sigma.a2[j,i] = pair.result$Gamma_e[1,1]
      sigma.b2[i,j] = sigma.b2[j,i] = pair.result$Gamma_e[2,2]
    }
  }

  for(i in 1:P){
    diag.ests=precision.est.diag[i,]
    diag.avg=mean(diag.ests[diag.ests>0])
    precision.est[i,i]=diag.avg
  }
    
  edge.bonferroni=matrix(0, P, P)
  edge.bonferroni[precision.pval<( 0.05/(P*(P-1)/2) )]=1
  diag(edge.bonferroni)=0
  edge.fdr=ifelse(p.adjust(precision.pval, "BH") < 0.1, 1, 0)
  edge.fdr=matrix(edge.fdr, nrow=P, ncol=P)
  diag(edge.fdr)=0
  
  all.end=Sys.time()
  exe.time=all.end-all.start
    
  list(
    precision.est=precision.est,
    precision.pval=precision.pval,
    edge.bonferroni=edge.bonferroni,
    edge.fdr=edge.fdr,
    eta.est=eta.est,
    eta.se=eta.se,
    rho.est=rho.est,
    rho.se=rho.se,
    partialcor.pval = rho.pval,
    sigma.a2 = sigma.a2,
    sigma.b2 = sigma.b2,
    exe.time=exe.time
  )
}




OurMethod_oneEdge <- function(data, i, j, scale.var=T, infer='llr', eps=1e-7){
  all.start=Sys.time()
  ## preprocessing
  data = as.matrix(data)
  N=dim(data)[1]
  P=dim(data)[2]
  data.c = scale(data, center=T, scale=F) # centering: center each column of data to have zero mean
  if(scale.var){
    Z = scale(data.c, center = T, scale = T)
  }
  else{
    Z = data.c
  }
  Z = Z / sqrt(P-2)
  ## Get data for each pair
  Y = data.c[, c(i,j)]
  X = Z[, -c(i,j)]
  ## Initialization of Gamma
  Gamma_beta = cov(Y)/2
  Gamma_e = cov(Y)/2
  ## For inference based on likelihood ratio test, first estimate when eta = 0
  if(infer == 'llr'){
    Gamma_e[1,2] = Gamma_e[2,1] = 0
    pair.eta0 = OurMethodOnePairEta0(Y, X, Gamma_beta, Gamma_e, eps=eps)
    Gamma_beta = pair.eta0$Gamma_beta
    Gamma_e = pair.eta0$Gamma_e
  }
  ## Estimate
  pair.result = OurMethodOnePair(Y, X, Gamma_beta, Gamma_e, infer=infer, fix.eta=FALSE, eps=eps)
  ## Inference
  if(infer == 'llr'){
    LLR = -2 * (pair.eta0$loglikeli - pair.result$loglikeli)
    pval.eta = pchisq(LLR, 1, lower.tail = F)
  }
  else if(infer == 'wald'){
    pval.eta = pair.result$pval.eta
    eta.se = pair.result$se.eta
    rho.se = pair.result$se.rho
    rho.pval = pair.result$pval.rho
  }
  rho.est = pair.result$rho

  all.end=Sys.time()
  exe.time=all.end-all.start
  list(
    loglikelis = pair.result$loglikelis,
    precision.est = pair.result$precision.pair[1,2],
    precision.pval = pval.eta,
    Gamma_beta = pair.result$Gamma_beta,
    Gamma_e = pair.result$Gamma_e,
    eta.est = pair.result$eta,
    eta.se = pair.result$se.eta,
    rho.est=rho.est,
    rho.se = rho.se,
    partialcor.pval = drop(rho.pval),
    exe.time = exe.time
  )
}
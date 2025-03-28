#' FLAG for one pair of random variables, using random effects model.
#' This is a repeated function for FLAG.
#'
#' @param Y Matrix, with size n*2.
#' @param X Matrix, with size n*(p-2).
#' @param Gamma_beta Matrix, with size 2*2.
#' @param Gamma_e Matrix, with size 2*2.
#' @param infer Character, option of different tests of inference where 'llr' for likelihood ratio test and 'wald' for Wald test based on Fisher Information Matrix.
#' @param fix.eta Logical, whether to fix eta, default to be FALSE.
#' @param eps Numeric, a small term to avoid numerical problems, default to be 1e-7.
#' @param max.iter Integer, the maximum number of iterations, default to be 5000.
#' @param crit.loglik Numeric, the criteria of the change ratio of log likelihood to stop.
#'
#' @return List,
#' the list of log likelihood during iterations,
#' the estimated Gamma_beta matrix with size 2*2, in the random effects model,
#' the estimated Gamma_epsilon matrix with size 2*2, in the random effects model,
#' the estimated 2*2 submatrix of the precision matrix,
#' the covariance matrix, which is the inverse of the Fisher information matrix, inferred by the Wald test,
#' the estimated off-diagonal element eta in the matrix Gamma_epsilon,
#' the standard error of eta,
#' the p-value of eta,
#' the estimated partial correlation rho,
#' the standard error of rho,
#' the p-value of rho.
#'
#' @export
#'
#' @importFrom stats cov
#'
#' @examples
#' \donttest{
#' FlagOnePair(Y, X)
#' }
FlagOnePair <- function(Y, X, Gamma_beta=NULL, Gamma_e=NULL, infer='llr', fix.eta=FALSE, eps=1e-7, max.iter=5000, crit.loglik=1e-4){
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


### For each pair of variables, use rank-two update.
FlagOnePairRankTwo <- function(Y, Zij, U, d, XXT, Gamma_beta=NULL, Gamma_e=NULL, infer='llr', fix.eta=FALSE, eps=1e-7, max.iter=5000, crit.loglik=1e-4){
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

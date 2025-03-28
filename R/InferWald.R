#' Get the standard error of rho.
#'
#' @param Gamma_e Matrix, with size 2*2.
#' @param Gamma_e_cov Matrix, with size 2*2.
#'
#' @return Numeric.
#' @export
#'
#' @importFrom stats pchisq
#'
#' @examples
#' \donttest{
#' GetSeRho(Gamma_e,Gamma_e_cov)
#' }
GetSeRho <- function(Gamma_e,Gamma_e_cov){
  d = c(-1/2 * Gamma_e[1,1]^{-3/2} * Gamma_e[2,2]^{-1/2} * Gamma_e[1,2],
      -1/2 * Gamma_e[1,1]^{-1/2} * Gamma_e[2,2]^{-3/2} * Gamma_e[1,2],
      Gamma_e[1,1]^{-1/2} * Gamma_e[2,2]^{-1/2})
  d = matrix(d, nrow=1, ncol=3)
  return(sqrt(d %*% Gamma_e_cov %*% t(d)))
}

#' Infer by the Wald test.
#'
#' @param Omega.inv Matrix, with size (2n)*(2n).
#' @param n Integer.
#' @param K Matrix, with size n*n.
#' @param Y.vec Matrix, with size (2*n)*1, the vectorized Y.
#' @param Gamma_e Matrix, with size 2*2.
#'
#' @return List,
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
#' @examples
#' \donttest{
#' InferWald(Omega.inv, n, K, Y.vec, Gamma_e)
#' }
InferWald <- function(Omega.inv, n, K, Y.vec, Gamma_e){
#   infer.start=Sys.time()
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
    se.rho = GetSeRho(Gamma_e, COV[c(2,4,6), c(2,4,6)])
    pval.rho = pchisq((rho^2)/(se.rho^2), 1, lower.tail=F)
  }
#   infer.end=Sys.time()
#   cat('Inference: ', infer.end-infer.start, '\n')

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

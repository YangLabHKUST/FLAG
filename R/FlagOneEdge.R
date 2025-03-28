#' Use FLAG to infer one edge. Given n*p data matrix, when we only interest in the conditional dependence between i-th and j-th variables.
#'
#' @param data Matrix, with size n*p.
#' @param i integer, the index of one element.
#' @param j integer, the index of another element.
#' @param scale.var Logical, whether to scale the variance of X to 1/p, default to be T(RUE).
#' @param infer Character, option of different tests of inference where 'llr' for likelihood ratio test and 'wald' for Wald test based on Fisher Information Matrix.
#' @param eps Numeric, a small term to avoid numerical problems, default to be 1e-7.
#'
#' @return List,
#' the list of log likelihood during iterations,
#' the estimated precision value,
#' the p-value of precision value estimation,
#' the estimated Gamma_beta matrix with size 2*2, in the random effects model,
#' the estimated Gamma_epsilon matrix with size 2*2, in the random effects model,
#' the estimated off-diagonal element eta in the matrix Gamma_epsilon,
#' the standard error of eta,
#' the estimated partial correlation rho,
#' the standard error of rho,
#' the p-value of rho,
#' the execution time.
#'
#' @export
#'
#' @importFrom stats cov
#'
#' @examples
#' \donttest{
#' FlagOneEdge(matrix, i, j)
#' }
FlagOneEdge <- function(data, i, j, scale.var=T, infer='llr', eps=1e-7){
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
    pair.eta0 = FlagOnePairEta0(Y, X, Gamma_beta, Gamma_e, eps=eps)
    Gamma_beta = pair.eta0$Gamma_beta
    Gamma_e = pair.eta0$Gamma_e
  }
  ## Estimate
  pair.result = FlagOnePair(Y, X, Gamma_beta, Gamma_e, infer=infer, fix.eta=FALSE, eps=eps)
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

#' FLAG is the main function to fulfill the whole process.
#'
#' @param data Matrix, with size n*p.
#' @param scale.var Logical, whether to scale the variance of X to 1/p, default to be T(RUE).
#' @param low.rank Logical, whether to use low rank update to shrink the time of eigen-decomposition of XX^T, default to be TRUE when sample size larger than 1000.
#' @param infer Character, option of different tests of inference where 'llr' for likelihood ratio test and 'wald' for Wald test based on Fisher Information Matrix.
#' @param eps Numeric, a small term to avoid numerical problems, default to be 1e-7.
#' @param crit.loglik Numeric, the criteria of the change ratio of log likelihood to stop.
#'
#' @return List,
#' the estimated precision matrix,
#' the p-value of precision matrix estimation,
#' the edge existence using Bonferroni correction,
#' the edge existence using false discovery rate,
#' the matrix of estimated eta,
#' the standard error or estimated eta,
#' the matrix of estimated partial correlation rho,
#' the standard error or estimated partial correlation rho,
#' the p-value of partial correlation matrix estimation,
#' the matrix of estimated sigma_a^2,
#' the standard error or estimated sigma_b^2,
#' the execution time.
#' @export
#'
#' @importFrom stats p.adjust
#' @importFrom stats cov
#'
#' @examples
#' \donttest{
#' FLAG(matrix)
#' }
#'
FLAG <- function(data, scale.var=T, low.rank=NULL, infer='llr', eps=1e-7, crit.loglik=1e-4){
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
          pair.eta0 = FlagOnePairRankTwoEta0(Y, Zij, U, d, XXT, Gamma_beta, Gamma_e, eps=eps, crit.loglik=crit.loglik)
        }
        else{
          pair.eta0 = FlagOnePairEta0(Y, X, Gamma_beta, Gamma_e, eps=eps, crit.loglik=crit.loglik)
        }
        Gamma_beta = pair.eta0$Gamma_beta
        Gamma_e = pair.eta0$Gamma_e
      }
      ## Estimate
      if(low.rank == TRUE){
        pair.result = FlagOnePairRankTwo(Y, Zij, U, d, XXT, Gamma_beta, Gamma_e,
                                 infer=infer, fix.eta=FALSE, eps=eps, crit.loglik=crit.loglik)
      }
      else{
        pair.result = FlagOnePair(Y, X, Gamma_beta, Gamma_e, infer=infer, fix.eta=FALSE, eps=eps, crit.loglik=crit.loglik)
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

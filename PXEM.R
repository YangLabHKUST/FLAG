PXEM <- function(Y, X, max.iter=500, crit.loglik=1e-7){  
  start=Sys.time()

  ## Initialization
  n = nrow(Y)
  p = ncol(X)
  Gamma_b = cov(Y)/2
  Gamma_e = cov(Y)/2
  # Eigenvalue decomposition for X%*%t(X)
  XXT= X%*%t(X)
  XTX = t(X)%*%X
  XTX.eigen = eigen(XTX)
  V = XTX.eigen$vectors
  q = XTX.eigen$values
  # q = ifelse(q<eps, eps, q)

  ## Preprocessing
  Y1 = Y[,1]
  Y2 = Y[,2]
  XTY1 = t(X)%*%Y1
  XTY2 = t(X)%*%Y2
  VTXTY1 = t(V) %*% XTY1
  VTXTY2 = t(V) %*% XTY2

  ELBO.old = 0
  ELBOs = NULL

  for(iter in 1: max.iter){
    delta = 1

    # E-step:
    Gamma_b.inv = solve(Gamma_b)
    Gamma_e.inv = solve(Gamma_e)
    Gei11 = Gamma_e.inv[1,1]
    Gei12 = Gei21 = Gamma_e.inv[1,2]
    Gei22 = Gamma_e.inv[2,2]

    h11 = Gei11*q + Gamma_b.inv[1,1]
    h12 = h21 = Gei12*q + Gamma_b.inv[1,2]
    h22 = Gei22*q + Gamma_b.inv[2,2]
    ad_cb = h11*h22 - h21*h12
    hi11 = h22/(ad_cb)
    hi12 = hi21 = -h12/(ad_cb)
    hi22 = h11/(ad_cb)

    mu1 = ( V * rep( (Gei11*hi11 + Gei21*hi12), rep(p, p)) ) %*% VTXTY1 +
          ( V * rep( (Gei12*hi11 + Gei22*hi12), rep(p, p)) ) %*% VTXTY2
    mu2 = ( V * rep( (Gei11*hi21 + Gei21*hi22), rep(p, p)) ) %*% VTXTY1 +
          ( V * rep( (Gei12*hi21 + Gei22*hi22), rep(p, p)) ) %*% VTXTY2

    mu1.mu1 = sum(mu1^2)
    mu1.mu2 = t(mu2)%*%mu1
    mu2.mu2 = sum(mu2^2)
    YXmu1.YXmu1 = sum((Y1 - X%*%mu1)^2)
    YXmu1.YXmu2 = t(Y2 - X%*%mu2) %*% (Y1 - X%*%mu1)
    YXmu2.YXmu2 = sum((Y2 - X%*%mu2)^2)

    Qfunc = - n/2*log(det(Gamma_e)) - p/2*log(det(Gamma_b)) -
            1/2* ( Gei11*YXmu1.YXmu1 + 2*Gei12*YXmu1.YXmu2 + Gei22*YXmu2.YXmu2 ) -
            1/2* ( Gamma_b.inv[1,1]*mu1.mu1 + 2*Gamma_b.inv[1,2]*mu1.mu2 + Gamma_b.inv[2,2]*mu2.mu2 )
    ELBO = Qfunc - sum(log(ad_cb))/2
    ELBOs = c(ELBOs, ELBO)
    if(abs((ELBO-ELBO.old)/ELBO.old) < crit.loglik) break
    if(iter>1 & ELBO < ELBO.old) message("Likelihood decreasing")
    ELBO.old = ELBO

    # M-step:
    qhi11 = q*hi11
    qhi12 = q*hi12
    qhi22 = q*hi22
    GeXTY1 = Gei11*XTY1 + Gei12*XTY2
    GeXTY2 = Gei21*XTY1 + Gei22*XTY2
    delta = drop( ( t(GeXTY1)%*%mu1 + t(GeXTY2)%*%mu2 ) /
                  ( Gei11*t(mu1)%*%XTX%*%mu1 + 2*Gei12*t(mu1)%*%XTX%*%mu2 + Gei22*t(mu2)%*%XTX%*%mu2 +
                    sum(Gei11*qhi11 + 2*Gei12*qhi12 + Gei22*qhi22) ) )
    YdXmu1.YdXmu1 = sum((Y1 - delta* X%*%mu1)^2)
    YdXmu1.YdXmu2 = t(Y2 - delta* X%*%mu2) %*% (Y1 - delta* X%*%mu1)
    YdXmu2.YdXmu2 = sum((Y2 - delta* X%*%mu2)^2)
    U.11.tr = YdXmu1.YdXmu1 + (delta^2)* sum( qhi11 )
    U.12.tr = U.21.tr = YdXmu1.YdXmu2 + (delta^2)* sum( qhi12 )
    U.22.tr = YdXmu2.YdXmu2 + (delta^2)* sum( qhi22 )
    W.11.tr = mu1.mu1 + sum( hi11 )
    W.12.tr = W.21.tr = mu1.mu2 + sum( hi12 )
    W.22.tr = mu2.mu2 + sum( hi22 )
    U.tr = matrix(c(U.11.tr, U.12.tr, U.21.tr, U.22.tr), 2)
    W.tr = matrix(c(W.11.tr, W.12.tr, W.21.tr, W.22.tr), 2)

    Gamma_b = (delta^2) * W.tr / p
    Gamma_e = U.tr / n

  }
  end=Sys.time()

  list(
    iterations = iter,
    ELBOs = ELBOs,
    Gamma_beta = Gamma_b,
    Gamma_epsilon = Gamma_e,
    exe.time = end-start
  )
}
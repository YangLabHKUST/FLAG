library(ggplot2)
library(igraph)

############# Visualization #############
get.roc.labels <- function(roc.list, digits=2){
  methods = names(roc.list)
  labels = c()
  for (method in methods){
    roc.tmp = roc.list[method]
    auc.tmp = roc.tmp[[1]]$auc
    labels = c(labels, paste0(method, ", AUC=", round(auc.tmp, digits)))
  }
  labels
}

sp100.significant <- function(method, pvals){
  P = nrow(pvals)
  S = matrix('', P, P)
  if(is.null(method)) return(S)
  fdr = get.fdr(pvals)
  for(i in 1:P){
    for(j in 1:P){
      if(fdr[i,j] == 1) S[i,j] = paste(S[i,j],'*',sep='')
      if(pvals[i,j] < 0.05/(P*(P-1)/2)) S[i,j] = paste(S[i,j],'*',sep='')
    }
  }
  diag(S)=''
  S
}
plot.sp100.edge.sign <- function(dims, prec, title, method=NULL, edge=NULL, pvals=NULL){
  P = length(dims)
  df = data.frame(X=rep(dims, P), Y=rep(dims, each=P))
  df$X = factor(df$X, levels = dims)
  df$Y = factor(df$Y, levels = dims)
  pc = prec2pc(prec)
  diag(pc)=0
  if(is.null(edge)) df$Z = c(pc)
  else df$Z = c(pc*edge)
  df$S = c(sp100.significant(method, pvals))
#   options(repr.plot.width=12, repr.plot.height=12)
  ggplot(data = df, aes(Y, X, fill = Z))+ geom_tile(aes(fill = Z)) + coord_fixed() + 
    scale_fill_gradient2(low = "#008837", high = "#7b3294", mid = "white", midpoint = 0, limits=c(-1,1), name='corr') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 7, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(),axis.ticks = element_blank(),
      legend.justification = c(1, 0), legend.position = c(0.97, 0.99), legend.direction = "horizontal") +
    geom_text(aes(Y, X, label = S), color = "black", size = 3)+ggtitle(title)+
    geom_hline(yintercept=c(9.5, 20.5, 31.5, 34.5, 49.5, 63.5, 75.5, 91.5, 94.5, 96.5)) +
    geom_vline(xintercept=c(9.5, 20.5, 31.5, 34.5, 49.5, 63.5, 75.5, 91.5, 94.5, 96.5))
}

### Brain data
## precision matrix
brain.significant <- function(S, method, pval, edge){
  if(method!='bggm'){
    P = nrow(pval)
    fdr = get.fdr(pval)
  }
  else P = nrow(edge)
  
  for(i in 1:P){
    for(j in 1:P){
      if(method=='bggm'){
        if(edge[i,j]==1) S[i,j] = paste(S[i,j],'*',sep='')
      }
      else{
        if(fdr[i,j] == 1) S[i,j] = paste(S[i,j],'*',sep='')
        if(pval[i,j] < 0.05/(P*(P-1)/2)) S[i,j] = paste(S[i,j],'*',sep='')
      }
    }
  }
  S
}
brain.heatmap <- function(dims, mat, title, method=NULL, pval=NULL, edge=NULL, tri=FALSE){
#   print(range(get.tri(mat)))
  P = length(dims)
  df = data.frame(X=rep(dims, P), Y=rep(dims, each=P))
  df$X = factor(df$X, levels = dims)
  df$Y = factor(df$Y, levels = dims)
  diag(mat)=0
  S = round(mat,2)
  
  if(tri){
    mat[lower.tri(mat)]=0
  }
  df$Z = c(mat)
  
  if(!is.null(method)) S=brain.significant(S, method, pval, edge)
  diag(S)=''
  if(tri){
    S[upper.tri(S)]=''
  }
  df$S = c(S)
  ggplot(df, aes(X, Y, fill= Z)) + geom_tile(aes(fill = Z)) + ggtitle(title) + coord_fixed() +
    scale_fill_gradient2(low = "#008837", high = "#7b3294", mid = "white", midpoint = 0, limits=c(-1,1), name='corr') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(),axis.ticks = element_blank(),
      legend.justification = c(1, 0), legend.position = c(0.97, 0.99), legend.direction = "horizontal") +
    geom_text(aes(Y, X, label = S), color = "black", size = 4)
}


## graph
get.net<-function(pval, mat){
  bc = get.bonferroni(pval)
  fdr = get.fdr(pval)
  edge = fdr*mat
  edge[bc==1] = -edge[bc==1]
  diag(edge)=0
  net <- graph_from_adjacency_matrix(edge, mode = "undirected",weighted = TRUE)
  E(net)$color[E(net)$weight>0] <- 'blue'
#   E(net)$lty[E(net)$weight>0] <- 2
  E(net)$color[E(net)$weight<0] <- 'red'
#   E(net)$lty[E(net)$weight<0] <- 1
  if(!is.null(E(net)$weight)) E(net)$weight = abs(E(net)$weight)
  net
}
brain.graph <- function(dims, mat, title, method=NULL, pval=NULL, path=NULL, filename=NULL){
  set.seed(0)
  mat = abs(mat)
  diag(mat)=0
  if(method=='om' | method=='fastggm') net = get.net(pval, mat)
  else net = graph_from_adjacency_matrix(mat, mode = "undirected",weighted = TRUE)
  if(method=='om.meta') E(net)$color='red'
  options(repr.plot.width=4, repr.plot.height=4)
  if(is.null(path)){
    plot(net, vertex.label = dims, layout = layout_in_circle,
       vertex.size = 1, vertex.color='red', vertex.label.family = "Helvetica", vertex.label.color = "black",
       vertex.label.dist = 2, vertex.label.cex= 1.3, vertex.shape = "circle", edge.width = 6*E(net)$weight, main=title)
  }
  if(!is.null(path)){
    if(is.null(filename)) jpeg(file = paste(c(path, title,'.jpeg'), collapse = ""))
    else jpeg(file = paste(c(path, filename,'.jpeg'), collapse = ""))
    plot(net, vertex.label = dims, layout = layout_in_circle,
       vertex.size = 1, vertex.color='red', vertex.label.family = "Helvetica", vertex.label.color = "black",
       vertex.label.dist = 2, vertex.label.cex= 1.3, vertex.shape = "circle", edge.width = 6*E(net)$weight, main=title)
    dev.off()
  }
}


############# Result Evaluation #############
## matrix operations
prec2pc<-function(prec){
  pc = -cov2cor(prec)
  diag(pc) = 1
  pc
}
get.tri <- function(matrix){
  matrix[ lower.tri( matrix ) ]
}
get.2hub <- function(matrix, p){
  matrix[1:p, 1:p]
}
get.hubEdge <- function(matrix, p, P=-1){
  if(P==-1){
    P = ncol(matrix)
  }
  matrix[1:p, (p+1):P]
}
get.block <- function(matrix, p, P=-1){
  if(P==-1){
    P = ncol(matrix)
  }
  matrix[(p+1):P, (p+1):P]
}
remove.diag <- function(matrix, val=0){
  diag(matrix)=val
  matrix
}

## parameter estimation
cal.err <- function(est, true, relative=FALSE, s1=1, e1=nrow(true), s2=1, e2=ncol(true)) {
  err.f=norm( est[s1:e1, s2:e2]-true[s1:e1, s2:e2], type='F' ) # forbenius norm
  err.f/norm(true[s1:e1, s2:e2], type='F')
#   err.f.re=err.f/norm( true, type='F' ) # relative forbenius norm error
#   list(err.f, err.f.re)
}

## graph recovery
# get edge (adjacency matrix)
get.bonferroni<-function(pval, thr=0.05, diag.val=0){
  P = ncol(pval)
  edge.bonferroni=matrix(0, P, P)
  edge.bonferroni[pval<( thr/(P*(P-1)/2) )]=1
  diag(edge.bonferroni)=diag.val
  edge.bonferroni
}
get.fdr<-function(pval, thr=0.1, prt=FALSE, diag.val=0){
  fdr=ifelse(p.adjust(pval, "BH") < thr, 1, 0)
  P = ncol(pval)
  fdr=matrix(fdr, nrow=P, ncol=P)
  diag(fdr)=diag.val
  if(prt) cat(sum(fdr)/2, '\n')
  fdr
}
get.adj<-function(pre, eps=1e-3, diag.val=0){
    P = ncol(pre)
    adj=matrix(0,P,P)
    adj[abs(pre)>eps]=1
    diag(adj)=diag.val
    adj
}
cal.rec <- function(est.p, true, getadj=TRUE) {
  if(getadj) est = get.adj(est.p)
  else est = est.p
  TP=sum( est * true ) # True Positive
#   TN=sum( (1-est) * (1-true) ) # True Negative
  FP=sum( est * (1-true) ) # False Positive
#   FN=sum( (1-est) * true ) # False Negative
  list(TP=TP, FP=FP)
  
#   MISR=(FN+FP)/P/(P-1) # misspecification rate
#   MCC=(TP*TN-FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) # Matthews correlation coefficient in[-1,1], 1 for perfect classifier
#   F1=2*TP/(2*TP+FP+FN)
#   list(MISR=MISR, F1=F1, MCC=MCC)
}


############# Run Methods #############
run.mle <- function(Z){
  solve(cov(Z))
}
## CLIME
# run.clime <- function(Z, lambda=0.2){
run.clime <- function(Z, lambda=0.1){
  clime1=fastclime::fastclime(Z)
  clime2=fastclime::fastclime.selector(clime1$lambdamtx, clime1$icovlist, lambda)
  clime2
}
## ANT method, by FastGGM algorithm
run.ant <- function(Z){
  FastGGM::FastGGM(Z)
}
## Bayesian GGM
run.bggm <- function(Z, prt = FALSE){
  start=Sys.time()
  bggm.fit = BGGM::estimate(Z, type = "continuous", analytic = FALSE, progress=FALSE)
  bggm.prec = BGGM::precision(bggm.fit)[[1]]
  bggm.edge = BGGM::select(bggm.fit)[[1]]
  end=Sys.time()
  if(prt) print(end-start)
  return(list(precision=bggm.prec, edge = bggm.edge))
}
## GLasso
# run.glasso <- function(Z, K=10, lambda=NULL){
run.glasso <- function(Z, K=5, lambda=NULL){
  if(!is.null(lambda)){
    GLASSO = glasso::glasso(cov(Z), rho=lambda)
    GLASSO$lambda = lambda
    return(GLASSO)
  }
  start=Sys.time()
#   glasso.CV = CVglasso::CVglasso(X=Z, lam = 10^seq(-6, 2, 0.2), K=K, trace='none')
#   glasso.CV = CVglasso::CVglasso(X=Z, lam = 10^seq(-3, 2, 0.2), K=K, trace='none')
  glasso.CV = CVglasso::CVglasso(X=Z, lam = 10^seq(-2, 2, 0.2), K=K, trace='none')
  end=Sys.time()
  GLASSO = glasso::glasso(cov(Z), rho=glasso.CV$Tuning[2])
  GLASSO$lambda = glasso.CV$Tuning[2]
  GLASSO$tune.time = end-start
  GLASSO
}
## Hub GLasso: input data Z
# precision: HGLASSO$Theta; hub nodes: HGLASSO$hubind
run.hglasso <- function(Z, lambda1=0.3, lambda2s=seq(0.05,0.5,by=0.05),lambda3=1, prt=FALSE, BICcriterion=NULL){
  for(lambda2 in lambda2s){
      res1=hglasso::hglasso(cov(Z), lambda1, lambda2, lambda3)
      BICcriterion=c(BICcriterion, hglasso::hglassoBIC(res1, cov(Z))$BIC)
  }
  lambda2=lambda2s[which(BICcriterion==min(BICcriterion))]
  HGLASSO=hglasso::hglasso(cov(Z), lambda1, lambda2[1], lambda3)
  if(prt) cat('HGL: ', lambda2[1], ',', HGLASSO$hubind)
  HGLASSO
}
run.dsglasso <- function(Z){
  SILGGM::SILGGM(Z,method='D-S_GL')
}


run.jgl <- function(Z.list, penalty='group', lam.list=10^seq(-0.6, 1.2, 0.2)){
  cal.aic<- function(N.list, S.list, Theta.list){
    sum = 0
    len = length(N.list)
    for(i in 1:len){
      N = N.list[i]
      S = S.list[[i]]
      Theta = as.matrix(Theta.list[[i]])
      Theta.adj = get.adj(Theta, diag.val=1)
#       sum = sum + N*sum(diag( S%*%Theta )) - N*log(det(Theta)) + 2*sum(Theta.adj)
      sum = tryCatch({
        sum + N*sum(diag( S%*%Theta )) - N*log(det(Theta)) + 2*sum(Theta.adj)
      }, error=function(e){
#         message("error:",e)
#         cat('\t',i,'\n')
#         print(S)
#         print(Theta)
        return(1e9)
      }, warning=function(w){
#         message("warning:",w)
#         cat('\t',i,'\n')
#         print(S)
#         print(Theta)
        return(1e9)
      })
    }
    sum
  }
  tune.start=Sys.time()
  N.list = c() # sample size
  S.list = list() # sample covariance
  Z.len = length(Z.list)
  for(i in 1:Z.len){
    N.list = c(N.list, nrow(Z.list[[i]]))
    S.list[[i]] = cov(Z.list[[i]])
  }
  lam1.list = lam2.list = lam.list
  lam.len = length(lam1.list)
  aic.mat = matrix(0, lam.len, lam.len)
  for(i in 1:lam.len){
    for(j in 1:lam.len){
      lam1 = lam1.list[i]
      lam2 = lam2.list[j]
#       cat(lam1, lam2, '\n')
      fgl = JGL::JGL(Z.list, penalty=penalty, weights=N.list, lambda1=lam1, lambda2=lam2)
      aic.mat[i,j] = cal.aic(N.list, S.list, fgl$theta)
    }
  }
  ind.i = which(aic.mat==min(aic.mat), arr.ind=TRUE)[1]
  ind.j = which(aic.mat==min(aic.mat), arr.ind=TRUE)[2] 
  tune.end=Sys.time()

  exe.start=Sys.time()
  jgl = JGL::JGL(Z.list, penalty=penalty, weights=N.list, lambda1=lam1.list[ind.i], lambda2=lam2.list[ind.j])
  exe.end=Sys.time()
  
  jgl$lam1 = lam1.list[ind.i]
  jgl$lam2 = lam2.list[ind.j]
  jgl$tune.time = tune.end-tune.start
  jgl$exe.time = exe.end-exe.start
  jgl
}

# source('/home/yqianai/graphical_model/precision_matrix_estimation/Algorithm/our_method_v2.R')
# OM=OurMethod(Z)
# cat('OM:', OM$exe.time, '\n')
# ANT=FastGGM::FastGGM(Z)



############# QQ Plot #############
qqplot1 <- function(dat.group, comp1, comp.labels=NULL, max=5){
  ci = 0.95
  n <- nrow(dat.group)
  df <- data.frame(pval = c(sort(dat.group[, comp1])), Method = c(rep(comp1, n)))
  df$pval[-log10(df$pval)>=max] = 10^(-1*max)
  df <- data.frame(Method = df$Method,
               observed = -log10(df$pval),
               expected = rep(-log10(ppoints(n)), length(unique(df$Method))),
               clower   = rep(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))),
               cupper   = rep(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))))
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  if(is.null(comp.labels)){
    df$Method = factor(df$Method, levels = comp1)
  }else{
    df$Method = factor(df$Method, levels = comp1, labels = comp.labels)
  }
  
  plot = ggplot(df) +
    geom_point(aes(expected, observed, color=Method), size = 2, shape = 16) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) + ylab(log10Po) + xlim(0,max) + ylim(0,max) +
    theme(legend.position="none",
          legend.title=element_blank(),
          legend.text = element_text(hjust = 0, size=15),
          axis.title = element_text(size = 15, color = "black"),
          axis.text  = element_text(size = 10,color = "black"),
          plot.title = element_text(face="bold", size = 20))
  return(plot)
}

qqplot2 <- function(dat.group, comp1, comp2, comp3, comp.labels=NULL,max=40, strong=NULL, mark=T){
  ci = 0.95
  n  <- nrow(dat.group)
  df <- data.frame(pval = c(sort(dat.group[, comp1]), sort(dat.group[, comp2])), Method = c(rep(comp1, n), rep(comp2, n)))
  df$pval[-log10(df$pval)>=max] = 10^(-1*max)
  df <- data.frame(Method = df$Method,
               observed = -log10(df$pval),
               expected = rep(-log10(ppoints(n)), length(unique(df$Method))),
               clower   = rep(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))),
               cupper   = rep(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))))
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  if(is.null(comp.labels)){
    df$Method = factor(df$Method, levels = c(comp1, comp2)
    )
  }else{
    df$Method = factor(df$Method, levels = c(comp1, comp2), labels = comp.labels
    )
  }
  
  plot = ggplot(df, aes(x=expected, y=observed)) +
    geom_point(data = df, aes(x=expected, y=observed,color=Method),size = 3, shape = 16) +
#     geom_point(data = df, aes(x=expected, y=observed,color=Method),size = 3, shape = 16, alpha=0.4,colour=c(rep('blue',n),rep('green',n))) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, color="#999999") +
    geom_line(aes(expected, clower), linetype = 2, color="#999999") +
    xlab(log10Pe) + ylab(log10Po) + ylim(0, max)+ xlim(0, max)+
    theme(legend.position = c(1,0.025),
          legend.justification = c(0.95,0),
          legend.text = element_text(hjust = 0, size=20),
          legend.title = element_blank(),
          axis.title = element_text(size = 20, color = "black"),
          axis.text  = element_text(size = 15,color = "black"),
          plot.title = element_text(hjust=-0.1, size = 20, face = "bold"))
    if(!is.null(strong)) plot= plot + geom_point(data=subset(df, Method==strong),aes(color=Method), size=2,shape=16)
  return(plot)
}

qqplot3 <- function(dat.group, comp1, comp2, comp3, comp.labels=NULL,max=40, strong=NULL, mark=T){
  ci = 0.95
  n  <- nrow(dat.group)
  df <- data.frame(pval = c(sort(dat.group[, comp1]), sort(dat.group[, comp2]), sort(dat.group[, comp3])),
              Method = c(rep(comp1, n), rep(comp2, n), rep(comp3, n)))
  df$pval[-log10(df$pval)>=max] = 10^(-1*max)
  df <- data.frame(Method = df$Method,
               observed = -log10(df$pval),
               expected = rep(-log10(ppoints(n)), length(unique(df$Method))),
               clower   = rep(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))),
               cupper   = rep(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = (n):1)),  length(unique(df$Method))))
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  if(is.null(comp.labels)){
    df$Method = factor(df$Method, levels = c(comp1, comp2, comp3)
    )
  }else{
    df$Method = factor(df$Method, levels = c(comp1, comp2, comp3), labels = comp.labels
    )
  }
  
  plot = ggplot(df, aes(x=expected, y=observed)) +
    geom_point(data = df, aes(x=expected, y=observed,color=Method),size = 3, shape = 16) +
#     geom_point(data = df, aes(x=expected, y=observed,color=Method),size = 3, shape = 16, alpha=0.4,colour=c(rep('blue',n),rep('green',n))) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, color="#999999") +
    geom_line(aes(expected, clower), linetype = 2, color="#999999") +
    xlab(log10Pe) + ylab(log10Po) + ylim(0, max)+ xlim(0, max)+
    theme(legend.position = c(1,0.025),
          legend.justification = c(0.95,0),
          legend.text = element_text(hjust = 0, size=20),
          legend.title = element_blank(),
          axis.title = element_text(size = 20, color = "black"),
          axis.text  = element_text(size = 15,color = "black"),
          plot.title = element_text(hjust=-0.1, size = 20, face = "bold"))
    if(!is.null(strong)) plot= plot + geom_point(data=subset(df, Method==strong),aes(color=Method), size=2,shape=16)
  return(plot)
}
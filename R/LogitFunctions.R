#' @title Logistic mixed effect regression with EM algorithm based on Polya-Gamma distributed latent variables.
#'
#' @description
#' Performs logistic mixed effect regression on fixed effects \code{X.fe} and \code{X.re}.
#' Algorithm employs Polya-Gamma latent variables to transform the phenotype into Gaussian framework.
#' Parameter coefficients are estimated using an expectation-maximisation algorithm.
#' 
#' @param y A boolean phenotype.
#' @param X.fe Design matrix of fixed effects.
#' @param X.re Design matrix of random effects.
#' @param c Covariance of fixed effects.
#' @param phi Covariance of random effects (genetic additive variance).
#'  @param inc.m boolean for whether the intercept should be estimated?
#' 
#' @return Logistic regression parameters for the fixed effects and random effects
#' specified in \code{X.re} and \code{X.fe}.
#' @export
logit.PG.EM.mm.gen <- function(y, X.re, X.fe, phi,c,  inc.m=F, n=rep(1, length(y)),
                               m0=rep(0, ncol(X.fe)))
{
  #cov FE = c
  #cov ME = delta/M
  
  M=nrow(X.fe)
  P0=matrix(0, nrow=ncol(X.fe), ncol=ncol(X.fe))
  diag(P0)=c
  y = as.numeric(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  #Block X
  
  ##Block Cov
  
  k = ((y-n/2))
  #inc.m=F
  if(inc.m) {
    X    = cbind(1,X.fe,X.re)
  }else {X    = cbind(X.fe,X.re)}
  P.a   = ncol(X.re)
  P.b   = ncol(X.fe)
  P.ab  = P.a + P.b
  P.all = P.ab + inc.m
  PP=matrix(0, nrow=P.all, ncol=P.all)
  P1=diag(phi/M, nrow=ncol(X.re), ncol=ncol(X.re))
  N = nrow(X)
  a.idc = inc.m+1:P.b
  b.idc = inc.m+1:P.a + P.b
  m.idc = which(inc.m)
  
  #Z = colSums(X * (y-n/2));
  #length(Z)
  m1=rep(0, ncol(X.re))
  
  mu=rep(0,P.all)
  dbm=rnorm(P.all,mu,sd=3)
  #save(init,file='Init.RData')
  #dbm=init
  betaold=rnorm(P.all,mu)
  PP[a.idc,a.idc]=P0
  PP[b.idc,b.idc]=P1
  if(inc.m){
    #Add covariance of intercept
    PP[m.idc, m.idc] = c
  }
  
  while ( abs(sum(betaold-dbm))>1E-5){
    print(abs(sum(betaold-dbm)))
    ## draw w
    #Length of number of individuals
    psi = drop(X %*% dbm)
    #  str(X)
    
    #length(psi)
    #Compute expected value of w
    #E-step on w
    w=drop((n/(2*psi))*tanh(psi/2))
    
    S = t(X) %*% (X * w);
    
    betaold=dbm
    B=t(X)%*%k
    A=S+solve(PP)
    
    dbm= solve(A,B);
    
  }
  dbm
}

#' @title Logistic mixed effect regression with Gibbs sampling algorithm based on Polya-Gamma distributed latent variables.
#'
#' @description
#' Performs logistic mixed effect regression on fixed effects \code{X.fe} and \code{X.re}.
#' Algorithm employs Polya-Gamma latent variables to transform the phenotype into Gaussian framework.
#' Parameter coefficients are estimated using a Gibbs sampler that returns a posterior for 
#' estimated fixed effects. Adapted from mixed-effect model in BayesLogit package. Point estimates for
#' the mixed effect estimates could be retrieved from posterior means.
#' 
#' @param y A boolean phenotype.
#' @param X.fe Design matrix of fixed effects.
#' @param X.re Design matrix of random effects.
#' @param c Covariance of fixed effects.
#' @param phi Covariance of random effects (genetic additive variance).
#' @param inc.m boolean for whether the intercept should be estimated?
#' 
#' @return Logistic regression parameter posteriors for the fixed effect coefficients
#' specified in \code{X.fe}.
#' @export
logit.PG.EM.mm.gen.gibbs <- function(y, X.re, X.fe, phi,c,  inc.m=F, n=rep(1, length(y)),
                                     m0=rep(0, ncol(X.fe)),kappa=0,samp=500,burn=500,verbose=100, P0=matrix(0, nrow=ncol(X.fe), ncol=ncol(X.fe)))
{
  #kappa=Inf
  #phi=0.007
  M=nrow(X.fe)
  diag(P0)=c
  y = as.numeric(y)
  X.re = as.matrix(X.re)
  X.fe = as.matrix(X.fe)
  #Block X
  
  ##Block Cov
  
  k = ((y-n/2))
  #inc.m=F
  if(inc.m) {
    X    = cbind(1,X.fe,X.re)
  }else {X    = cbind(X.fe,X.re)}
  P.a   = ncol(X.re)
  P.b   = ncol(X.fe)
  P.ab  = P.a + P.b
  P.all = P.ab + inc.m
  PP=matrix(0, nrow=P.all, ncol=P.all)
  P1=diag(phi/M, nrow=ncol(X.re), ncol=ncol(X.re))
  N = nrow(X)
  
  #str(X)
  a.idc = inc.m+1:P.b
  b.idc = inc.m+1:P.a + P.b
  m.idc = which(inc.m)
  
  #Z = colSums(X * (y-n/2));
  
  
  #length(Z)
  m1=rep(0, ncol(X.re))
  #Vector of means
  #mu=c(1,m0,m1)
  mu=rep(0,P.all)
  #Parameter vector
  #Initialise
  dbm=rnorm(P.all,mu)
  #save(init,file='Init.RData')
  #dbm=init
  
  start.time = proc.time()
  ## Sample
  PP[a.idc,a.idc]=P0
  PP[b.idc,b.idc]=P1
  if(inc.m){
    #Add covariance of intercept
    PP[m.idc, m.idc] = c
  }
  Z = colSums(X * (y-n/2));
  if(kappa!=0){
    out <- list(w = matrix(nrow=samp, ncol=N),
                dbm = matrix(nrow=samp, ncol=P.ab+1))
  }else{
    out <- list(w = matrix(nrow=samp, ncol=N),
                dbm = matrix(nrow=samp, ncol=P.ab))
    #  a.idc
  }
  #samp=10
  #burn=10
  for ( j in 1:(samp+burn) ){
    ## draw w
    #Length of number of individuals
    psi = drop(X %*% dbm)
    #E-step on w
    #w=drop((n/(2*psi))*tanh(psi/2))
    w = rpg.devroye(N, n, psi);
    
    S = t(X) %*% (X * w);
    Vw=solve(S+solve(PP))
    
    pm = Vw %*% as.vector(Z);
    dbm = pm + Vw %*% rnorm(P.all);
    if (j>burn) {
      out$w[j-burn,]   = w
      out$dbm[j-burn,] = dbm
      
    }
    if (j %% verbose == 0) { print(paste("LogitPG MM 2: Iteration", j)); }
  }
  if (inc.m) out$fe = out$dbm[,c(m.idc, a.idc)]
  else out$fe = cbind(0,out$dbm[,a.idc])
  colnames(out$fe) = c("Intercept", colnames(X.fe))
  out
  
} 
#' @title Logistic regression with EM algorithm based on Polya-Gamma distributed latent variables.
#'
#' @description
#' \code{logit.EM.R} returns logistic regression parameters for the fixed effects specified in \code{X}.
#'
#' 
#' @param y A binomial phenotype.
#' @param X Design matrix of fixed effects.
#' @param c Covariance of fixed effects.
#' 
#' @return Logistic regression parameters for the fixed effects specified in \code{X}.
#' @export
logit.EM.R <- function(y, X, c=1)
{
  #X=X.fe
  n=rep(1, length(y))
  m0=rep(0, ncol(X))
  P0=matrix(0, nrow=ncol(X), ncol=ncol(X))
  diag(P0)=c
  ## X: n by p matrix
  X = as.matrix(X);
  y = as.numeric(y)
  p = ncol(X)
  N = nrow(X)
  w = rep(0,N)
  ## w = w.known;
  mu=rep(0,p)
  betaold=rnorm(length(mu),0)
  beta=rnorm(length(mu),0)
  
  while (     abs(sum(betaold-beta))>1*10E-10){
    
    psi = drop(X%*%beta)
    w=drop((n/(2*psi))*tanh(psi/2))
    k = ((y-n/2))
    S = t(X) %*%diag(w) %*% (X ) ;
    betaold=beta
    beta = solve(S+chol2inv(chol(P0)),t(X)%*%k+chol2inv(chol(P0))%*%mu);
    print(abs(sum(betaold-beta)))
    
  }
  return(beta)  
  
} 
#A function for 
beta_calc_wb=function(w,X,prec_mat,k,mu,K,phi){
  pi_=diag(w)
  X.t=t(X)
  inv.pi=chol2inv(chol(pi_))
  Kpi=K%*%pi_
  M=ncol(X)
  #Woodbury-matrix 
  inv.part=pi_-pi_%*%solve(diag(M/phi,nrow=nrow(K))+Kpi)%*%Kpi
  sum_comp=X.t%*%inv.part%*%X
  A=(sum_comp+prec_mat)
  sum_comp2=X.t%*%inv.part%*%inv.pi%*%(k)
  B=sum_comp2+prec_mat%*%mu
  #Solve a system of linear equations
  beta=solve(A,B)
  beta
}

gamma_calc_wb=function(w,X.re,k,phi,c){
  mu=rep(0,ncol(X.re))
  M=ncol(X.re)
  prec_mat=solve(diag(phi/M,nrow=M))
  pi_=diag(w)
  B_=(diag(c,nrow=nrow(pi_)))
  X.t=t(X.re)
  inv.pi=chol2inv(chol(pi_))
  inv.part=solve(inv.pi+B_)
  sum_comp=(X.t%*%(inv.part))%*%X.re
  A=(sum_comp+prec_mat)
  sum_comp2=X.t%*%inv.part%*%inv.pi%*%(k)
  B=sum_comp2+prec_mat%*%mu
  gamma=solve(A,B)
  gamma
}

#' @title Alternative EM-based logistic mixed effect regression on Polya-Gamma distributed latent variables.
#'
#' @description
#' Performs logistic mixed effect regression on fixed effects \code{X.fe} and \code{X.re}.
#' Algorithm employs Polya-Gamma latent variables to transform the phenotype into Gaussian framework.
#' Parameter coefficients are from iterative adjustments to Polya-Gamma omega (latent variables),
#' beta (fixed effects) and gamma (random effects).
#' Omega expected values come from expecepted values of Polya-Gamma distribution, beta and omega
#' values use maximisation of Gaussian parameters from the maximum likelihood estimate
#' of Bayesian linear regression. Note! It is extremely slow and not included in benchmarking.
#'  
#' 
#' @param y A boolean phenotype.
#' @param X.fe Design matrix of fixed effects.
#' @param X.re Design matrix of random effects.
#' @param c Covariance of fixed effects.
#' @param phi Covariance of random effects (genetic additive variance).
#'
#' @return Logistic regression parameter posteriors for the fixed effect coefficients
#' specified in \code{X.fe}.
#' @export
logit.PG.EM.fastlmm.ksisamp <- function(y, X.re, X.fe,phi,c=1,n=rep(1, length(y)),
                                        m0=rep(0, ncol(X.fe)), P0=matrix(0, nrow=ncol(X.fe), ncol=ncol(X.fe))){
  y = as.numeric(y)
  X.re = as.matrix(X.re)
  #Knorm=norm.kinship(K)
  #norm.kinship=apply(X.re,2,function(x){
  #  #x=mat_f[,1]
  #  pj=mean(x)/2
  #  (x-mean(x))/(sqrt(pj*(1-pj)))
  #})
  #K=as.matrix(norm.kinship)%*%as.matrix(t(norm.kinship))  
  
  K=as.matrix(X.re)%*%as.matrix(t(X.re))  

  X.fe = as.matrix(X.fe)
  k = (y-(n/2))
  diag(P0)=c
  ## X: n by p matrix
  X = as.matrix(cbind(X.fe,X.re));
  y = as.numeric(y)
  p = ncol(X.fe)
  N = nrow(X.fe)
  ## w = w.known;
  mu=rep(0,p)
  betaold=rnorm(p)

  prec_mat=solve(P0)
  beta=rnorm(p)
  #gamma=rnorm(dim(X.re)[2])
  gamma=rep(0,dim(X.re)[2])
  ## Sample
  while ( abs(sum(betaold-beta))>1E-3){
   
    psi=X%*%as.matrix(c(beta,gamma))
    w=drop((n/(2*psi))*tanh(psi/2))
    rem=(abs(sum(betaold-beta)))
    wrong.way.mat=c(rem,wrong.way.mat[1:19])
    betaold=beta

    beta=beta_calc_wb(w,X.fe,prec_mat,k,mu,K,phi)
    gamma=gamma_calc_wb(w,X.re,k,phi,c)

    print(sum(betaold-beta))
    }
  
  #delta
  list(beta=beta,phi=phi)
}

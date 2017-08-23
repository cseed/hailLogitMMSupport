#' @title A simple data generation method for use benchmarking
#'
#' @description
#' Generate a dataset with 3 fixed effects and with an arbitrary number
#' of genotypes and some relatedness.
#' 
#' @param nx Number of variants.
#' @param probset1  Probabilities of genotypes for first half of dataset ie (c(0.1,0.2,0.5)) .
#' @param probset2 Probabilities of genotypes for second half of dataset ie (c(0.1,0.2,0.5)) .
#' @param X.re Design matrix of random effects.
#' @param c Covariance of fixed effects.
#' @param phi Covariance of random effects (genetic additive variance).
#' 
#' @return A list of phenotypes ($y), random effects ($X.re), fixed effects ($X.fe), probabilites
#' for the phenotype ($probs).
#' @export
Data=function(nx,beta=c(0.4,0.6,-0.4),phi=1,probset1=c(0.4,0.3,0.3),probset2=c(0.1,0.2,0.5)){
  
  mat1=data.frame(Z=matrix(replicate(nx,sample(c(0,1,2),nx/2,
                                               prob = probset1,replace=T) ),nrow=nx))
  mat2=data.frame(Z=matrix(replicate(nx,sample(c(0,1,2),nx/2,
                                               prob = probset2,replace=T) ),nrow=nx))
  
  mat_f=as.matrix(cbind(mat1,mat2))
  colnames(mat_f)=paste('Z.',1:ncol(mat_f),sep='')
  
  ##Normalization by Patterson (Population Structure and Eigenanalysis)...
  
  #norm.kinship=apply(mat_f,2,function(x){
    #x=mat_f[,1]
  #  pj=mean(x)/2
  #  (x-mean(x))/(sqrt(pj*(1-pj)))
  #})
  #K=as.matrix(norm.kinship)%*%as.matrix(t(norm.kinship))
  
  K=as.matrix(mat_f)%*%as.matrix(t((mat_f)))
  
  num=rnorm(nrow(mat_f))
  num2=rnorm(nrow(mat_f))+0.01*(1:nrow(mat_f))
  
  ksi=mvrnorm(mu=cbind(1,num,num2)%*%beta,Sigma=phi*K)
  probs=exp(ksi)/((exp(ksi)+1))
  
  y=sapply(probs,function(z){
    rbinom(n=1,size=1,prob=z)
  })
  
  X.fe=as.matrix(model.matrix(y~num+num2))
  colnames(X.fe)[1]='Inter'
  y=y
  X.re=as.matrix(mat_f)
  list(y=data.frame(y),X.fe=X.fe,X.re=X.re,probs=probs)
}

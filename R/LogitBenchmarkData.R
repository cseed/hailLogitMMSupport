#' @title A simple data generation method for use benchmarking
#'
#' @description
#' A general file for simple test and use cases for 
#' all implemented versions of logistic (mixed) models1) GLM - standard classic logistic regression R implementation for fixed effect estimate benchmarks
#' 2) logit.EM.R - EM version of logistic regression with Polya-Gamma latent variables
#' 3) logit.PG.EM.mm.gen.gibbs - Gibbs sampling version of logistic mixed model regression with Polya-Gamma latents
#' 4) logit.PG.EM.mm.gen - EM version of logistic mixed model regression with Polya-Gamma latents
#' 
#' The sample data has been created from Balding-Nichols model
#' and was exported as a file as given in Jupyter Notebook sampleNotebook.ipynb.
#' 3 fixed effects (V1,V2,V3) were extracted from sampleAnnot.tsv.
#' Phenotype column Pheno.
#' 
#' @param data a datset for benchmarking.
#' @param neff number of fixed effects
#' @param c Covariance of fixed effects.
#' @param phi Covariance of random effects (genetic additive variance).
#' 
#' @examples
#' generateBenchmarks(system.file('extdata','sampleData.tsv',package='hailLogitMMSupport'))
#' 
#' Dat=Data(200)
#' generateBenchmarks(Dat)
#' 
#' @return A data frame of fixed effect parameters benched with four aforementioned models.
#' 
#' @export
generateBenchmarks=function(data=system.file('extdata','sampleData.tsv',package='hailLogitMMSupport')
        ,neff=3,phi=0.007,c=1.0){
if(!is.character(data)){
  BenchData=data
  y=BenchData$y
  X.fe=BenchData$X.fe
  X.re=BenchData$X.re
  dv_=paste(colnames(X.fe),collapse='+')
  dataset=data.frame(cbind(X.fe,y))
  nm=colnames(y)
  y=unlist(y)
  #GLM
  #Model 1
  modGlm=glm(as.formula(paste(nm,' ~ ',dv_,' +0',sep='')),data=dataset,family = 'binomial')
  
}else{
BenchData=read.table(data)
vars=paste('V',1:neff,sep='')
dv_=paste(vars,collapse='+')
colnames(BenchData)=c('ID','Variant','GT',vars,'Pheno')
dataset=reshape2::dcast(as.formula(paste('ID + ',dv_,' + Pheno ~ Variant',sep='')),value.var='GT',data=BenchData)
#head(dataset)
y=dataset$Pheno
X.fe=dataset[,2:(neff+1)]
X.re=dataset[,(neff+3):ncol(dataset)]
#GLM
#Model 1
modGlm=glm(as.formula(paste('Pheno~ ',dv_,' +0',sep='')),data=dataset,family = 'binomial')
}

#Logit EM / PG
#Model 2
modEmr=logit.EM.R(y,X.fe)

#Model 3
Gibbs_Base2=logit.PG.EM.mm.gen.gibbs(y, X.re, X.fe,phi=phi,c=c,samp=1000, burn=1000, verbose=100,kappa=0)
#Model 4
modmm.PGEM=logit.PG.EM.mm.gen(y, X.re, X.fe, phi=phi,c=c,inc.m=F)

result_=data.frame(
GLM=coef(modGlm),
GLM_EM=modEmr[,1],
GLMM_GIBBS=colMeans(Gibbs_Base2$fe)[-1],
GLMM_EM=modmm.PGEM[1:neff])
return(result_)
}

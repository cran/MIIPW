## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(MIIPW)
data("srdata1")
head(srdata1)
apply(srdata1,2,anyNA)
mice::md.pattern(srdata1[,-c(1,2)],plot = TRUE)

## -----------------------------------------------------------------------------
formula<-C6kine~ActivinRIB+ActivinRIIA+ActivinRIIAB+Adiponectin+AgRP+ALCAM
pMat<-mice::make.predictorMatrix(srdata1[names(srdata1)%in%all.vars(formula)])
m1<-MeanScore(data=srdata1,
formula<-formula,id='ID',
visit='Visit',family='gaussian',init.beta = NULL,
init.alpha=NULL,init.phi=1,tol=.00001,weights = NULL,
corstr = 'exchangeable',maxit=50,m=2,pMat=pMat)
summary_meanscore(m1)

## ----eval=FALSE---------------------------------------------------------------
#  m2<-SIPW(data=srdata1,formula<-formula,id='ID',
#  visit='Visit',family='gaussian',corstr = 'exchangeable',maxit=5)
#  
#  m3<-AIPW(data=srdata1,
#  formula<-formula,id='ID',
#  visit='Visit',family='gaussian',init.beta = NULL,
#  init.alpha=NULL,init.phi=1,tol=.00001,weights = NULL,
#  corstr = 'exchangeable',maxit=50,m=3,pMat=pMat)
#  
#  m4<-miSIPW(data=srdata1,
#  formula<-formula,id='ID',
#  visit='Visit',family='gaussian',init.beta = NULL,
#  init.alpha=NULL,init.phi=1,tol=0.001,weights = NULL,
#  corstr = 'exchangeable',maxit=50,m=2,pMat=pMat)
#  
#  m1<-miAIPW(data=srdata1,
#  formula<-formula,id='ID',
#   visit='Visit',family='gaussian',init.beta = NULL,
#  init.alpha=NULL,init.phi=1,tol=.00001,weights = NULL,
#  corstr = 'exchangeable',maxit=4,m=2,pMat=pMat)
#  

## -----------------------------------------------------------------------------
m1<-MeanScore(data=srdata1,
             formula<-formula,id='ID',
             visit='Visit',family='gaussian',init.beta = NULL,
             init.alpha=NULL,init.phi=1,tol=.00001,weights = NULL,
             corstr = 'exchangeable',maxit=50,m=2,pMat=pMat)
 m11<-MeanScore(data=srdata1,
             formula<-formula,id='ID',
             visit='Visit',family='gaussian',init.beta = NULL,
             init.alpha=NULL,init.phi=1,tol=.00001,weights = NULL,
            corstr = 'independent',maxit=50,m=2,pMat=pMat)
QICmiipw(model.R=m1,model.indep=m11,family="gaussian")
##


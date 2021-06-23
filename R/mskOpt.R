#' @title
#' Estimates of parameter corresponding to minimum AIC
#' @description provides estimates of parameter of linear regression model of the
#' response variable corresponding to the minimum AIC value using mean score method with different covariance structure
#' @details It calculates the AIC value for the linear regression model
#' \deqn{Y_{ij}=\beta_0+\beta_1x_{1ij}+\beta_2 x_{2ij}+...+\beta_p x_{pij}+e_{ij}}
#' using mean score method with different covariance structures and gives the estimates of
#' parameter for minimum AIC value
#'
#' @param data balanced longitudinal data set where each subject's outcome has
#'  been measured at same time points and number of visits for each patient is same,
#' covariance structure of the outcome variable like "unstuctured","compound","ToE",
#' "AR1","markov","independence"
#' @param Dep column name of dependent variable in the dataset
#' @param Id column name of id of subjects in the dataset
#' @param Time column name of timepoints in the dataset
#' @param m starting column number of covariates
#' @param n ending column number of covariates
#' @return estimated parameter value for multiple linear regression model for that covariance structure with minimum AIC value.
#' @import matlib
#'
#' @examples
#' data(srdata)
#' mskopt(Dep="C6kine",Id="ID",Time="Visit",m=5,n=10,data=srdata)
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma
#'
#'
mskopt<-function(Dep,Id,Time,m,n,data){
  InD<-Dep
  dat1<-data
  cvstr1=c("independence","compound","markov","AR1","ToE","unstructured")
  cat("\n Estimated parameters for minimum AIC value   \n")
  aIC<-c()
  for(i in 1:length(cvstr1)){
    aIC[i]<-mskall(data=dat1,cvstr=cvstr1[i],Dep=InD,Id=Id,Time=Time,m=m,n=n)$aic
  }
  d<-data.frame(cvstr1,aIC)
  params<-mskall(data=dat1,cvstr=d$cvstr1[which.min(d$aIC)],Dep=InD,Id=Id,Time=Time,m=m,n=n)$coff
  rslt<-list()
  rslt$params<-params
  rslt$cvs<-d$cvstr1[which.min(d$aIC)]
  rslt$Aic<-d
  rslt

}

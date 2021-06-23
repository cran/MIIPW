#' @title
#' Bayesian analysis of repeated measurement data using
#' linear regression with an unstructured covariance matrix
#'
#'
#' @description Bayesian analysis is performed using MCMC and uses a
#' linear regression with an unstructured covariance matrix and the prior distributions are
#' uninformative normal distributions for the regression coefficients and an uninformative Wishart for
#' the 3 by 3 precision matrix of the outcome variable measured at 3 timepoints.
#' @details The mean reponse model is \deqn{E(Y_{ij})=\beta_0+\beta_1t_j}  with unstructured
#' covariance \eqn{\Sigma} where i refers to ith individual and j is the timepoint.
#'
#' @param m starting column number
#' @param n ending column number
#' @param n.chains number of MCMC chains
#' @param data  dataset whose first row is the respective time points at which observations are taken
#' @return posterior distribution results of the parameters
#' @import R2jags
#' @examples
#'  ##
#'  data(repeatdata)
#'  byrgu(m=1,n=3,n.chains=1,data=repeatdata)
#'  ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma


byrgu<-function(m,n,n.chains,data){

  data1<-function(d){
    data1<-matrix(ncol=ncol(d),nrow=nrow(d))
    for(i in 1:nrow(d))
    {
      for(j in 1:ncol(d)){
        data1[i,j]=ifelse(is.na(d[i,j])==TRUE,0.0001,d[i,j])
      }
    }
    data1
  }
  #m,n are the number of column from which to which column we want our output
  age<-unname(unlist(data[1,c(m:n)]))
  data<-data1(data[-1,])
  data<-data[,m:n]
  N1 <- nrow(data)
  M1 <- n-m+1
  Y <-data
  # length of age1 variable must be same with number of column in the dataset
  # length of age variable must be same with number of column we are considering



  R <- diag(1,ncol(data))

  mdata<-list("N1","M1","Y",
              "age",
              "R" )
  mreg101<-function(){
    # prior distribution for the regression coefficients.
    beta1 ~ dnorm(0.0, 0.001)
    beta2 ~ dnorm(0.0, 0.001)
    for(i in 1:N1){Y[i,1:M1]~dmnorm(mu1[],Omega[,])}
    for (j in 1:M1){mu1[j]<-beta1+beta2*age[j]}
    # non-informative precision matrix
    Omega[1:M1,1:M1]~dwish(R[,],M1)
    Sigma[1:M1,1:M1]<-inverse(Omega[,])
    for(i in 1:M1){
      for(j in (i+1):M1){
        rho[i,j]<-Sigma[i,j]/sqrt(Sigma[i,i]*Sigma[j,j])
      }}

  }


  jagsfit1 <- jags( model.file=mreg101,
                    data =mdata,
                    n.chains=n.chains,parameters.to.save =c('beta1','beta2'))
  jagsfit1



}

utils::globalVariables(c("beta1","beta2","inverse","Omega","jags"))




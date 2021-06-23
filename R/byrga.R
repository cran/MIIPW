#' @title
#' Bayesian analysis of mean response model with autoregressive  covariance
#' matrix
#'
#' @description Bayesian analysis is performed using MCMC and uses a
#' linear regression with an autoregressive covariance matrix for the response
#' @details The model for the response is \deqn{Y_{ij}=X_{ij}'\beta+e_{ij}}
#' and \deqn{e_{ij}=\rho e_{ij-1}+u_{ij}} ,\eqn{u_{ij}\sim N(0,1/\tau);\rho} is the correlation
#' coefficient where i refers to ith individual and j is the timepoint.
#'
#' @param m starting column number
#' @param n ending column number
#' @param n.chains number of MCMC chains
#' @param data dataset whose first row is the respective time points at which observations are taken
#' @return posterior distribution  results of the parameters
#' @import R2jags
#' @examples
#' ##
#' data(repeatdata)
#' byrga(m=1,n=3,n.chains=1,data=repeatdata)
#' ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma

byrga<-function(m,n,n.chains,data){
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


mdata<-list("N1","M1","Y",
            "age"
            )




mreg102<-function(){
  beta1 ~ dnorm(0.0, 0.001)
  beta2 ~ dnorm(0.0, 0.001)
  for(i in 1:N1){for (j in 2:M1){Y[i,j]~dnorm(mu[i,j],tau)}}
  for(i in 1:N1){for (j in 2:M1){mu[i,j]<-beta1*(1-rho)
  +beta2*(age[j]-rho*age[j-1])+rho*Y[i,j-1]}}
  for(i in 1:N1){Y[i,1]~dnorm(mu[i,1],tau)}
  rho~dbeta(1,1)
  for(i in 1:N1){mu[i,1]<-beta1+beta2*age[1]}
  tau~dgamma(.01,.01)
  sigma<-1/tau
}
jagsfit1 <- jags( model.file=mreg102,
                  data =mdata,
                  n.chains=n.chains,parameters.to.save =c('beta1','beta2') )
jagsfit1
}

utils::globalVariables(c("rho"))


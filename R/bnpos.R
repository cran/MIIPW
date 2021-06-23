#' @title
#' Bayesian analysis of generalised nonlinear mixed model with one random effect
#' @description Bayesian analysis of generalized nonlinear mixed model where response follows poisson distribution using
#' MCMC
#' @details Here the response variable \eqn{Y_{ij}} has poisson distribution
#'  mean and variance given one random effect as \eqn{E(Y_{ij}|b_i)=Var(Y_{ij}|b_i)
#'  } where link function is 
#'  \deqn{log(E(Y_{ij}|b_i))=\beta_1+b_{1i}+\beta_2 exp(-\beta_{3}x_{ij})}
#'  where i is the ith subject and j is the timepoint and \eqn{b_i\sim N(0,\sigma_1^2)} are independent.
#'
#' @param m starting column number
#' @param n ending column number
#' @param n.chains number of MCMC chains
#' @param data dataset with integer entries with missing values. first row represents proportion,size of dataset 
#' should be 11 by 6. Inside the function we are taking first 6 column of the propdata in this package for the example given.
#'
#'
#' @return posterior distribution result of parameters
#' @import R2jags
#' @references Broemeling, Lyle D. Bayesian methods for repeated measures. CRC Press, 2015.
#' @examples
#'  ##
#'  data(propdata)
#'  bnpos(m=1,n=3,n.chains=1,data=propdata)
#'  ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma
#'

bnpos<-function(m,n,n.chains,data)
  {
  cl<-match.call()
  cl1<-as.list(cl[-1])
  
  if(cl1$data==deparse(substitute(propdata))){
    data<-propdata[,c(1:6)]
  }else{
    data<-data
  }
  
  
  
s<-unname(unlist(data[1,][m:n]))
data<-data[-1,m:n]
M <- ncol(data)
N <- nrow(data)
Y <-data

mdata<-list("N","M","Y",
            "s")
mreg109<-function(){
  beta[1] ~ dnorm(0,.001)
  beta[2] ~ dnorm(0,.001)
  beta[3] ~ dnorm(0,.001)
  for(i in 1:N){for (j in 1:M){Y[i,j]~dpois(mu[i,j])}}
  for(i in 1:N){for (j in 1:M){log(mu[i,j])<-b[i]+beta[1]
  +beta[2]*exp(-beta[3]*s[j])}}
  for(i in 1:N){b[i]~dnorm(0,tau1)}
  tau1~dgamma(.001,.001)
  # sigma1 is the variance of the random effects
  sigma1<-1/tau1
}
jagsfit1 <- jags( model.file=mreg109,
                  data =mdata,
                  n.chains=n.chains,parameters.to.save =c('beta','sigma1') )
jagsfit1



}

utils::globalVariables(c("b","propdata","log<-"))
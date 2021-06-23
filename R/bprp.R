#' @title
#' Bayesian analysis of generalised nonlinear mixed model with one random effect
#' using MCMC
#'
#' @description provides Bayesian analysis of generalized nonlinear mixed model where response follows Poisson distribution using
#' MCMC
#' @details The response variable \eqn{Y_{ij}} follows poisson distribution
#' with conditional mean and variance \eqn{E(Y_{ij}|b_i)=Var(Y_{ij}|b_i)}
#' \eqn{b_i} are independent \eqn{N(0,\sigma^2)} and link function is 
#' \deqn{logit(\theta_{ij})=\beta_1+\beta_2(1-x_{ij}^{\beta_3+b_i})}
#' where \eqn{\theta_{ij}=P(Y_{ij}=1/b_i)},
#' \eqn{x_{ij}} is the proportion 
#'
#' @param m starting column number
#' @param n ending column number
#' @param n.chains  number of MCMC chains
#' @param data dataset with entries 0,1. first row represents proportion;size of dataset 
#' should be 11 by 6. Inside the function we are taking first 6 column of the propdata in this package for the given example 
#'
#' @return posterior distribution result of parameters.
#' @import R2jags
#' @references Broemeling, Lyle D. Bayesian methods for repeated measures. CRC Press, 2015.
#' @examples
#' ##
#' data(propdata)
#' bprp(m=1,n=4,n.chains=1,data=propdata)
#' ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma



bprp<-function(m,n,n.chains,data)
  {
  cl<-match.call()
  cl1<-as.list(cl[-1])
  
  if(cl1$data==deparse(substitute(propdata))){
    data<-propdata[,c(1:6)]
  }else{
    data<-data
  }
  
  
  
  
  
  #m,n are the number of column from which to which column we want our output


  x<-unname(unlist(data[1,][m:n]))
  data<-data[-1,m:n]
  M <- ncol(data)
  N <- nrow(data)
  Y <-data

  mdata<-list("N","M","Y",
              "x")

mreg108<-function(){
  # prior distribution for the regression coefficients.
  beta[1]~dnorm(11,.1)
  beta[2]~dnorm(-23,.1)
  beta[3]~dunif(.7,1)
  for(i in 1:N){for(j in 1:M){Y[i,j]~dpois(theta[i,j])}}
  for(i in 1:N){for (j in 1:M)
  {logit(theta[i,j])<-beta[1]+beta[2]*(1-
                                         pow(x[j],beta[3]+b1[i]))}}
  for (i in 1:N){b1[i]~dnorm(0,tau)}
  # prior distribution for the precision of the random
  #effect
  tau~dgamma(3,2)
  sigma<-1/tau
}
jagsfit1 <- jags( model.file=mreg108,
                  data =mdata,
                  n.chains=n.chains,parameters.to.save =c('beta','sigma') )
jagsfit1


}


utils::globalVariables(c("logit<-","pow","b1","tau","propdata"))

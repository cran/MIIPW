#' @title
#' Bayesian analysis of generalised linear mixed model for Poisson outcome variable with two random effect
#' @description provides bayesian analysis of generalised linear mixed model with log link function 
#' for categorical response using MCMC 
#' @details The response variable \eqn{Y_{ij}} follows poisson distribution
#' ,mean and variance given random effects \eqn{E(Y_{ij}|b_{1i},b_{2i})=Var(Y_{ij}|b_{1i},b_{2i})}
#' with link function \eqn{log(\mu_{ij})=\beta_1+b_{1i}+(\beta_2+b_{2i})(X_{ij}-\beta_3)}
#' where i is the ith subject and j is the timepoint.
#'
#' @param m  starting column number
#' @param n ending column number
#' @param n.chains number of MCMC chains
#' @param data dataset with integer entries(including NA values), first row represents time points.
#' In this function we are using observations at four timepoints, function takes first four columns of the countdata
#' for the example given in the package
#'
#' @return posterior distribution result of parameters
#' @import R2jags
#' @references Broemeling, Lyle D. Bayesian methods for repeated measures. CRC Press, 2015.
#' @examples
#' ##
#' data(countdata)
#' bgnml(m=1,n=2,n.chains = 1,data=countdata)
#' ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma

bgnml<-function(m,n,n.chains,data)
{
  cl<-match.call()
  cl1<-as.list(cl[-1])
  
  if(cl1$data==deparse(substitute(countdata))){
    data<-countdata[,c(1:4)]
  }else{
    data<-data
  }
  
  x<-unname(unlist(data[1,][m:n]))
  data<-data[-1,m:n]
  M <- ncol(data)
  N <- nrow(data)
  Y <-data
  mdata<-list("N","M","Y",
              "x")
  mreg1010<-function(){
    beta[1] ~ dnorm(0,.01)
    beta[2] ~ dnorm(0,.01)
    beta[3] ~ dnorm(0,.01)
    for(i in 2:N){for (j in 1:M){Y[i,j]~dpois(mu[i,j])}}
    for(i in 2:N){for (j in 1:M){Z[i,j]~dpois(mu[i,j])}}
    for(i in 2:N){for (j in 1:M){log(mu[i,j])<-(beta[1]+b1[i])+
      (beta[2]+b2[i])*(x[j]-beta[3])}}
    for(i in 1:N){b1[i]~dnorm(0,tau1)}
    for(i in 1:N){b2[i]~dnorm(0,tau2)}
    tau1~dgamma(3,2)
    sigma1<-1/tau1
    tau2~dgamma(3,2)
    sigma2<-1/tau2
  }
  
  jagsfit1 <- jags( model.file=mreg1010,
                    data =mdata,
                    n.chains=n.chains,parameters.to.save =c('beta','b1','b2','sigma1','sigma2'))
  jagsfit1
  
}
utils::globalVariables(c("countdata"))
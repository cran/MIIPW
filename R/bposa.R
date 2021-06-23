#' @title
#'Bayesian analysis of generalized mixed linear model using MCMC
#'
#' @description  provides Bayesian analysis of generalized mixed linear model where the repeated measure(with missing value) has poisson distrbution
#'using MCMC
#' @details The model for this function is
#'\deqn{Y_{ij}\sim Poisson(\mu_{ij})}
#'with link function
#'\deqn{log(\mu_{ij}|b_{1i},b_{2i})=\beta_1+\beta_2t_j+b_{1i}+b_{2i}t_j}
#'where the \eqn{b_{1i},b_{2i}}, i = 1,2,…,N are independent and have a two-dimensional
#'normal distribution with a 2 by 1 mean vector 0 and unknown 2 by 2 covariance matrix \eqn{\Sigma}
#'
#' @param m starting column number
#' @param n ending column number
#' @param n.chains number of MCMC chains
#' @param data dataset whose first row is the respective time points at which observations(integers) are taken
#' where timepoints are the respective column names,dimension have to be 26 by 3.
#'
#' @return posterior distribution result of the parameters
#' @import R2jags
#' @references Broemeling, Lyle D. Bayesian methods for repeated measures. CRC Press, 2015.
#' @examples
#' ##
#' data(catdata)
#' bposa(m=1,n=2,n.chains=1,data=catdata)
#' ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma

bposa<-function(m,n,n.chains,data)
  {
  cl<-match.call()
  cl1<-as.list(cl[-1])
  
  if(cl1$data==deparse(substitute(catdata))){
    data<-catdata[c(1:26),c(9:11)]
  }else{
    data<-data
  }
  
  data1=function(d){
    data1=matrix(ncol=ncol(d),nrow=nrow(d))
    for(i in 1:nrow(d))
    {
      for(j in 1:ncol(d)){
        data1[i,j]=ifelse(is.na(d[i,j])==TRUE,1,d[i,j])
      }
    }
    data1
  }
  #putting 1 in place of NA, if we put less than equal to 0 then it willnot work
  #m,n are the number of column from which to which column we want our output
  tim<-unname(unlist(data[1,c(m:n)]))
  data<-data1(data)
  data<-data[,m:n]
  N <- nrow(data)
  M <- n-m+1
  Y <-data

  mdata<-list("N","M","Y",
              "tim")
 mreg105A<-function(){
  for (i in 1:N){for (j in 1:M){Y[i,j]~dpois(mu[i,j])
    log(mu[i,j])<-beta[1]+beta[2]*tim[j]+b1[i]+b2[i]*tim[j]}}
  for(i in 1:2){beta[i]~dnorm(0,.0001)}
  for (i in 1:N){b1[i]~dnorm(0,tau1)}
  for (i in 1:N){b2[i]~dnorm(0,tau2)}
  tau1~dgamma(0.001,.001)
  tau2~dgamma(0.001,.001)
  sigma1<-1/tau1
  sigma2<-1/tau2
}
 jagsfit1 <- jags( model.file=mreg105A,
                   data =mdata,
                   n.chains=n.chains,parameters.to.save =c('beta','sigma1') )
 jagsfit1


}

utils::globalVariables(c("b2","tau1","tau2","catdata","log<-"))

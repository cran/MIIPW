#' @title
#' Bayesian analysis of non linear mean response model with random effect compound covariance
#' matrix
#'
#' @description Provides bayesian analysis of random effect model for repeated measurement data with missing values
#' using MCMC for compound symmetry covariance structure,where age follows normal distribution
#' @details The model for the response is \deqn{Y_{ij}=\beta_1+b_i+\beta_{2}exp(-\beta_{3}x_{ij})+e_{ij}}
#',where \eqn{e_{ij}} are independent \eqn{ N(0,\sigma^2)} and independent of n random effects \eqn{b_i\sim N(0,\sigma_b^2)}
#' ,where i refers to ith individual and j is the timepoint.
#'
#'
#'
#' @param m starting column number
#' @param n ending column number
#' @param n.chains number of MCMC chains
#' @param data dataset with first row represents proportion,dimension have to be 9 by 11.
#' @return posterior distribution result of parameters
#' @import R2jags
#' @references Broemeling, Lyle D. Bayesian methods for repeated measures. CRC Press, 2015.
#' @examples
#' ##
#' data(propdata)
#' byran(m=1,n=4,n.chains=1,data=propdata)
#' ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma

byran<-function(m,n,n.chains,data)
  {
  
  #if(c("age")%in%names(data)){data<-data[,-which(colnames(data)=="age")]
  #colnames(data)<-c("age","1","2","3","4")
  #}
  #else{data<-data}
  cl<-match.call()
  cl1<-as.list(cl[-1])
  
  if(cl1$data==deparse(substitute(propdata))){
    data<-propdata[c(1:9),c(7:17)]
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

mreg104=function()
  {
  beta[1] ~ dnorm(0,.001)
  beta[2] ~ dnorm(0,.001)
  beta[3] ~ dnorm(0,.001)
  for(i in 1:N){for (j in 1:M){Y[i,j]~dnorm(mu[i,j],tau)}}
  for(i in 1:N){for (j in 1:M)
  {mu[i,j]<-b[i]+beta[1]+beta[2]*exp(-beta[3]*s[j])}}
  for(i in 1:N){b[i]~dnorm(0,tau1)}
  tau~dgamma(.001,.001)
  vr<-sigma+sigma1
  sigma<-1/tau
  tau1~dgamma(.001,.001)
  sigma1<-1/tau1
  rho<-sigma1/(sigma+sigma1)

}
jagsfit1 <- jags( model.file=mreg104,
                  data =mdata,
                  n.chains=n.chains,parameters.to.save =c('beta','sigma1','mu') )
jagsfit1


}
utils::globalVariables(c("b","propdata"))

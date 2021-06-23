#' @title
#'Bayesian analysis of generalized mixed linear model using MCMC
#'
#' @description  provides Bayesian analysis of generalized mixed linear model where the repeated measure(with missing value) has poisson distrbution
#'using MCMC
#' @details The model for this function is
#'\deqn{Y_{ij}\sim Poisson(\mu_{ij})}
#'with link function
#'\deqn{log(\mu_{ij}|b_{1i},b_{2i})=\beta_1+\beta_2t_j+\beta_3X_1+\beta_4 X_2+\beta_5 X_3+\beta_6 X_4
#'+b_{1i}+b_{2i}t_j}
#'where the \eqn{b_{1i},b_{2i}}, i = 1,2,…,N are independent and have a two-dimensional
#'normal distribution with a 2 by 1 mean vector 0 and unknown 2 by 2 covariance matrix \eqn{\Sigma}
#' @param m starting column number
#' @param n ending column number
#' @param n.chains number of MCMC chains
#' @param data dataset whose first row is the respective time points at which observations(integers) are taken
#' where timepoints are the respective column names,dimension should be 26 by 3 and design matrix with column names X1,X2,X3,X4
#'
#' @return posterior distribution result of the parameter
#' @import R2jags
#' @references Broemeling, Lyle D. Bayesian methods for repeated measures. CRC Press, 2015.
#' @examples
#' ##
#' data(catdata)
#' bpois(m=1,n=3,n.chains=1,data=catdata)
#' ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma

bpois<-function(m,n,n.chains,data)
  {
  cl<-match.call()
  cl1<-as.list(cl[-1])
  
  if(cl1$data==deparse(substitute(catdata))){
    data<-catdata[c(1:26),c(9:11)]
    X<-catdata[c(1:25),c(15:18)]
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
  #m,n are the number of column from which to which column we want our output
  tim<-unname(unlist(data[1,c(m:n)]))
  data<-data1(data[-1,])
  data<-data[,m:n]
  N <- nrow(data)
  M <- n-m+1
  Y <-data
  X1<-X[,colnames(X)=="X1"]
  X2<-X[,colnames(X)=="X2"]
  X3<-X[,colnames(X)=="X3"]
  X4<-X[,colnames(X)=="X4"]


  mdata<-list("N","M","Y",
              "tim","X1","X2","X3","X4")

mreg105<-function(){
  for (i in 1:N){for (j in 1:M){Y[i,j]~dpois(mu[i,j])
    log(mu[i,j])<-beta[1]+beta[2]*tim[j]+beta[3]*X1[i]+beta[4]*X2[i]+beta[5]*X3[i]+beta[6]*X4[i]+b1[i]+b2[i]*tim[j]}}
  for(i in 1:6){beta[i]~dnorm(0,.0001)}
  for (i in 1:N){b1[i]~dnorm(0,tau1)}
  for (i in 1:N){b2[i]~dnorm(0,tau2)}
  tau1~dgamma(0.001,.001)
  tau2~dgamma(0.001,.001)
  sigma1<-1/tau1
  sigma2<-1/tau2
  d34<-beta[3]-beta[4]
  d35<-beta[3]-beta[5]
  d36<-beta[3]-beta[6]
  d45<-beta[4]-beta[5]
  d46<-beta[4]-beta[6]
  d56<-beta[5]-beta[6]
}

jagsfit1 <- jags( model.file=mreg105,
                  data =mdata,
                  n.chains=n.chains,parameters.to.save =c('beta','sigma1','sigma2'))
jagsfit1

}
utils::globalVariables(c("b2","tau1","tau2","catdata","log<-"))

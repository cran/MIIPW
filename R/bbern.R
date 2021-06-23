#' @title
#' Bayesian analysis of generalized mixed linear model using MCMC
#'
#' @description Provides bayesian analysis of generalized mixed linear model where the repeated measure(with missing value) follows bernoulli distrbution
#'using MCMC
#' @details The model for the response variable is given by
#' \deqn{Y_{ij}\sim Bernoulli(\mu_{ij})}
#' where link function is
#' \deqn{logit(\mu_{ij})=\beta_1+\beta_2 t_j+\beta_3 t_{j}^2+\beta_4 age_{i}+\beta_{5}gen_{i}
#' \beta_{6}grp_{i}+b_{1i}}
#' where i is the ith individual and j is the timepoint.
#' @param m starting column number
#' @param n ending column number
#' @param nc number of MCMC chains
#' @param data dataset with entries 0,1. column names are age,grp(group),gen(gender),0,1,2,3,4 (time points)
#'
#'
#' @return posterior distribution results of parameters.
#' @import R2jags
#' @references  Broemeling, Lyle D. Bayesian methods for repeated measures. CRC Press, 2015.
#' @examples
#' ##
#' data(catadata)
#' bbern(m=1,n=3,nc=1,data=catadata)
#' ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma


bbern<-function(m,n,nc,data){
  n.chains<-nc
  grp<-data[,which(colnames(data)=="grp")]
  z1<-which(colnames(data)=="age")
  z2<-which(colnames(data)=="grp")
  z3<-which(colnames(data)=="gen")
  t<-as.numeric(colnames(data[,-c(z1,z2,z3)]))[m:n]

  #m,n are the number of column from which to which column we want our output
  age<-data[,which(colnames(data)=="age")]
  gen<-data[,which(colnames(data)=="gen")]
  data<-data[,-c(z1,z2,z3)]
  data<-data[,m:n]
  N <- nrow(data)
  M <- n-m+1
  Y <-data


  mdata<-list("N","M","Y",
              "t","gen","grp","age")

 mreg106<-function(){
  for(i in 1:N){for (j in 1:M){Y[i,j]~dbern(mu[i,j])}}
  for(i in 1:N){for(j in 1:M){logit(mu[i,j])<-beta[1]+beta[2]*t[j]+beta[3]*t[j]*t[j]+beta[4]*age[i]+beta[5]*gen[i]+
    beta[6]*grp[i]+b1[i]}}
  for(i in 1:N){for (j in 1:M){p[i,j]<-(exp(beta[1]+beta[2]*t[j]+beta[3]*t[j]*t[j]+beta[4]*age[i]+beta[5]*gen[i]+
                                              beta[6]*grp[i]+b1[i]))/(1+exp(beta[1]+beta[2]*t[j]+beta[3]*t[j]*t[j]+beta[4]*age[i]+beta[5]*gen[i]+beta[6]*grp[i]
                                                                            +b1[i]))}}
  for(i in 1:6){beta[i]~dnorm(0,.0001)}
  for(i in 1:N){b1[i]~dnorm(0,tau1)}
  tau1~dgamma(.0001,.0001)
  sig1<-1/tau1
}

 jagsfit1 <- jags( model.file=mreg106,
                   data =mdata,
                   n.chains=n.chains,parameters.to.save =c('beta','sig1') )
 jagsfit1


}


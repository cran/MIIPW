#' @title
#' Bayesian analysis of mean response model with autoregressive  covariance
#' matrix
#'
#' @description  Bayesian analysis of mean response over time and age (quadratic trends)
#' using MCMC for AR1 covariance structure,where age follows normal distribution
#'
#' @details The model for the response is \deqn{Y_{ij}=\beta_1+\beta_{2}t_{ij}
#' +\beta_{3}t_{ij}^2+\beta_{4} age_i+e_{ij}} \deqn{e_{ij}=\rho e_{ij-1}+u_{ij}} ,\eqn{u_{ij}\sim N(0,1/\tau);\rho} is the correlation
#' coefficient where i refers to ith individual and j is the timepoint. Missing values of covariate age is
#' imputed assuming age follows normal distribution.
#'
#'
#' @param m starting column number
#' @param n ending column number
#' @param n.chains number of MCMC chains
#' @param data dataset with first column is age( there are missing values in age),and columns other than age are observation at
#' four different timepoints, where timepoints are the respective column names.
#'
#' @references Broemeling, Lyle D. Bayesian methods for repeated measures. CRC Press, 2015.
#' @return posterior distribution result of parameters
#' @import R2jags
#'
#' @examples
#' ##
#' data(agedata)
#' bygan(m=1,n=3,n.chains=1,data=agedata)
#' ##
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma
#' @seealso Brgan


bygan<-function(m,n,n.chains,data)
  {
  #if(c("age")%in%names(data)){data<-data[,-which(colnames(data)=="age")]
  #colnames(data)<-c("age","1","2","3","4")
  #}
  #else{data<-data}
  cl<-match.call()
  cl1<-as.list(cl[-1])
  
  if(cl1$data==deparse(substitute(agedata))){
    data<-agedata[,-1]
    colnames(data)<-c("age","1","2","3","4")
  }else{
    data<-data
  }
  
  data1<-function(d){
    data1=matrix(ncol=ncol(d),nrow=nrow(d))
    for(i in 1:nrow(d)){
      for(j in 1:ncol(d)){
        data1[i,j]=ifelse(is.na(d[i,j])==TRUE,0.0001,d[i,j])
      }
    }
    data1
  }
  #m,n are the number of column from which to which column we want our output
  age<-data$age#here length is same as number of row of dataset
  time<-as.numeric(colnames(data[,-which(colnames(data)=="age")]))[m:n]
  data<-data1(data[,-which(colnames(data)=="age")])
  data<-data[,m:n]
  N <- nrow(data)
  M <- n-m+1
  Y <-data

  mdata<-list("N","M","Y",
              "age","time")
mreg103B<-function(){

  for (i in 1:4){beta[i] ~ dnorm(0.0, 0.001)}
  for (i in 1:N){age[i]~dnorm(vu, taua)}
  vu~dnorm(.001,.001)
  taua~dgamma(.001,.001)
  sigmaa<-1/taua
  for(i in 1:N){for (j in 2:M){Y[i,j]~dnorm(mu[i,j],tau)}}
  # future values for hematocrit denoted by Z
  for(i in 1:N){for (j in 1:M){Z[i,j]~dnorm(mu[i,j],tau)}}
  for(i in 1:N){for (j in 2:M){mu[i,j]<-beta[1]*(1-rho)+
    beta[2]*(time[j]-rho*time[j-1])+beta[3]*(time[j]*time[j]-rho
                                             *time[j-1]*time[j-1])+rho*Y[i,j-1]+beta[4]*age[i]*(1-rho)}}
  for(i in 1:N){Y[i,1]~dnorm(mu[i,1],tau)}
  for(i in 1:N){mu[i,1]<-beta[1]+beta[2]*time[1]+beta[3]*time[1]*time[1]+beta[4]*age[i]}
  rho~dbeta(1,1)
  tau~dgamma(.001,.001)
  sigma<-1/tau
}

jagsfit1 <- jags( model.file=mreg103B,data =mdata,
                  n.chains=n.chains,parameters.to.save =c('beta','sigma') )
jagsfit1

}

utils::globalVariables(c("taua","rho","log-","agedata"))

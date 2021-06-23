#' @title
#' Mean score method for missing covariate value in linear regression model for repeated measurement data
#'
#' @description provides estimates of parameter from linear regression model using meanscore method for repeated measurement data.
#' @details  Mean score method is used for getting the missing score function value in the estimating
#' equation.
#'
#' @param data balanced longitudinal data set where each subject's outcome has been measured at same time points and
#' number of visits for each patient is samiliar covariance structure
#' @param cvstr "unstuctured","compound","ToE","AR1","markov","independence"
#' @param Dep column name of dependent variable in the dataset
#' @param Id column name of id of subjects in the dataset
#' @param Time column name of timepoints in the dataset
#' @param m starting column number of covariates
#' @param n ending column number of covariates
#' @return estimated parameter value for multiple linear regression model
#' @import matlib
#' @examples
#' data(srdata)
#' mskall(cvstr="ToE",Dep="C6kine",Id="ID",Time="Visit",m=5,n=10,data=srdata)
#' @export
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma
#'
mskall<-function(cvstr="unstructured",Dep,Id,Time,m,n,data)
{
  InD<-Dep
  dat1<-data
  dat2<-dat1[,-c(which(colnames(dat1)==c(Id)),which(colnames(dat1)==c(InD)),which(colnames(dat1)==c(Time)))][c(m:n)]
  dat1<-data.frame(id=dat1[,which(colnames(dat1)==c(Id))],t=dat1[,which(colnames(dat1)==c(Time))]
                   ,y=dat1[,which(colnames(dat1)==c(InD))],dat2)
  dat1$y[is.na(dat1$y)]<-mean(dat1$y,na.rm=TRUE)
  dat1<-as.data.frame(dat1)

  id1<-dat1$id
  id2<-dat1$id[1]
  tm<-dat1$t[dat1$id==dat1$id[1]]
  id3<-length(id1[id1==id2])
  f<-unique(dat1$id)
  ##unstuctured covarince matrix

  if(is.null(cvstr)==T){
    ustr<-matrix(dat1$y,byrow = F,nrow=length(id1)/id3,ncol=id3)
    sig1<-cov(ustr)
  }else if(cvstr=="unstructured"){
    ustr<-matrix(dat1$y,byrow = F,nrow=length(id1)/id3,ncol=id3)
    sig1<-cov(ustr)
  }else if(cvstr=="compound"){
    ustr<-matrix(dat1$y,byrow = F,nrow=length(id1)/id3,ncol=id3)
    vr<-sum(unlist(lapply(as.data.frame(ustr),function(x){var(x)}))*(length(f)-1))
    sig<-vr/((nrow(ustr)-1)*ncol(ustr))
    ro<-cor(ustr[,1],ustr[,2])
    ##correlation coefficient is taken as the correlation of first two time points
    sig1<-(ro*matrix(1,nrow=id3,ncol=id3)+(1-ro)*diag(1,nrow=id3,ncol=id3))*sig
  }else if(cvstr=="AR1"){
    ##autoregressive covarince matrix
    ustr<-matrix(dat1$y,byrow = F,nrow=length(id1)/id3,ncol=id3)
    vr<-sum(unlist(lapply(as.data.frame(ustr),function(x){var(x)}))*(length(f)-1))
    sig<-vr/((nrow(ustr)-1)*ncol(ustr))
    ro<-cor(ustr[,1],ustr[,2])
    sig1<-matrix(nrow=id3,ncol=id3)
    for(i in 1:id3){
      for(j in 1:id3 ){
        if(i!=j){
          sig1[i,j]=sig*(ro^(max(i,j)-min(i,j)))}
        else{sig1[i,j]=sig}
      }}
  }else if(cvstr=="ToE"){
    ##Torplitz or Banded Covarince Matrix
    ustr<-matrix(dat1$y,byrow = F,nrow=length(id1)/id3,ncol=id3)
    vr<-sum(unlist(lapply(as.data.frame(ustr),function(x){var(x)}))*(length(f)-1))
    sig<-vr/((nrow(ustr)-1)*ncol(ustr))
    sig1<-matrix(nrow=id3,ncol=id3)
    for(i in 1:id3){
      for(j in 1:id3 ){
        if(i!=j){
          sig1[i,j]=sig*(cor(ustr[,i],ustr[,j]))}
        else{sig1[i,j]=sig}
      }
    }
  }else if(cvstr=="independence"){
    ustr<-matrix(dat1$y,byrow = F,nrow=length(id1)/id3,ncol=id3)
    sig1<-matrix(nrow=id3,ncol=id3)
    for(i in 1:id3){
      for(j in 1:id3 ){
        if(i==j){
          sig1[i,j]=cov(ustr[,i],ustr[,j])}
        else{sig1[i,j]=0}
      }
    }
  }else if(cvstr=="markov"){
    ##autoregressive covarince matrix
    ustr<-matrix(dat1$y,byrow = F,nrow=length(id1)/id3,ncol=id3)
    vr<-sum(unlist(lapply(as.data.frame(ustr),function(x){var(x)}))*(length(f)-1))
    sig<-vr/((nrow(ustr)-1)*ncol(ustr))
    ro<-cor(ustr[,1],ustr[,2])
    sig1<-matrix(nrow=id3,ncol=id3)
    for(i in 1:id3){
      for(j in 1:id3 ){
        if(i!=j){
          sig1[i,j]=sig*(ro^(abs(tm[i]-tm[j])))}
        else{sig1[i,j]=sig}
      }}

  }
  ms<-function(dat){
    d<-function(x){
      dlt1<-c()
      for(j in 1:ncol(x)){
        if(anyNA.data.frame(x[,j])==T){
          dlt1[j]<-0
        }
        else{
          dlt1[j]<-1
        }
      }
      dlt1
    }


    dlt<-d(t(dat))
    dlt1<-d(dat)
    if(sum(dlt)!=0){
      cc<-as.matrix(t(dat)[,which(dlt==1)])
      if(0%in%dlt){
        nc<-as.matrix(t(dat)[,-which(dlt==1)])
        ncc<-as.matrix(nc[which(dlt1==1),])
        ccc<-as.matrix(cc[which(dlt1==1),])

        phik<-function(u,w){
          ##u is for missing (y,z,w)
          ##w is for complete (y,z,w)
          v=matrix(nrow=ncol(w),ncol=ncol(u))
          for(i in 1:ncol(w)){
            for(j in 1:ncol(u)){
              if(identical(u[,j],w[,i])==T){
                v[i,j]=1
              }else{
                v[i,j]=0
              }
            }

          }
          k=list()
          for(j in 1:ncol(v)){
            k[[j]]=which(v[,j]==1)

          }
          k}
        k<-phik(u=ncc,w=ccc)
        y1<- cc[1,]
        dsg=rbind(1,as.matrix(cc[-1,]))
        k2=list()
        for(j in 1:length(k)){
          if(length(k[[j]])!=0){
            k2[[j]]=as.matrix(dsg[,k[[j]]])%*%Ginv(as.matrix(sig1[k[[j]],k[[j]]]))%*%as.matrix(y1[k[[j]]])/length(k[[j]])}
          else{k2[[j]]=0}}

        k3=list()
        for(j in 1:length(k)){
          if(length(k[[j]])!=0){
            k3[[j]]=as.matrix(dsg[,k[[j]]])%*%Ginv(as.matrix(sig1[k[[j]],k[[j]]]))%*%t(as.matrix(dsg[,k[[j]]]))/length(k[[j]])
          }else{k3[[j]]=0}
        }
        if(length(k2)!=0){
          h2=Reduce("+",k2)
        }else{h2=0}
        if(length(k3)!=0){h3=Reduce("+",k3)}else{ h3=0}
        xy1<-as.matrix(dsg)%*%Ginv(as.matrix(sig1[which(dlt==1),which(dlt==1)]))%*%as.matrix(cc[1,])+h2
        xx1<-as.matrix(dsg)%*%Ginv(as.matrix(sig1[which(dlt==1),which(dlt==1)]))%*%t(dsg)+h3

      }else{
        dsg=rbind(1,as.matrix(cc[-1,]))
        xy1<-as.matrix(dsg)%*%Ginv(as.matrix(sig1[which(dlt==1),which(dlt==1)]))%*%as.matrix(cc[1,])
        xx1<-as.matrix(dsg)%*%Ginv(as.matrix(sig1[which(dlt==1),which(dlt==1)]))%*%t(dsg)
      }}else{xx1<-NULL;xy1<-NULL}
    list(xx1,xy1)
  }

  l1<-list()
  l2<-list()
  dat1<-data.frame(dat1)
  for(i in 1:length(unique(dat1$id))){

    l1[i]<-ms(subset(dat1,dat1$id==unique(dat1$id)[i])[,-c(1,2)])[1]
    l2[i]<-ms(subset(dat1,dat1$id==unique(dat1$id)[i])[,-c(1,2)])[2]
  }
  l1[sapply(l1,is.null)]<-NULL
  l2[sapply(l2,is.null)]<-NULL
  xx<-Reduce("+",l1)
  xy<-Reduce("+",l2)
  n<-nrow(xx)
  a<-matrix(unlist(xx),nrow=n,ncol=n,byrow = T)
  b<-matrix(unlist(xy),nrow=n,ncol=1,byrow = T)
  parm<-Ginv(a)%*%b
  rslt<-list()
  #cat("\n Estimated parameters using Mean Square method :\n")
  name1<-which(colnames(dat1)=="id")
  name2<-which(colnames(dat1)=="t")
  name3<-which(colnames(dat1)=="y")
  coefficients<-c("intercept",names(dat1)[-c(name1,name2,name3)])
  estimate<-parm

  rslt$coff<-data.frame(coefficients,estimate)
  ##AIC
  ##
  #print("AIC :")
  fit<-cbind(1,dat1[,-c(name1,name2,name3)])
  resid<-na.omit(dat1$y-as.matrix(fit)%*%parm)
  Aic<- nrow(na.omit(dat1))*(log(2*pi)+1+log((sum(resid^2)/nrow(na.omit(dat1)))))+((n+1)*2)
  rslt$aic<-Aic
  ##BIC
  #BIC = n*ln(MSE) + log(n) * k
  #print("BIC :")
  Bic<-nrow(na.omit(dat1))*log((sum(resid^2)/nrow(na.omit(dat1))))+n*log(nrow(na.omit(dat1)))
  rslt$Bic<-Bic
  rslt

}
utils::globalVariables(c("cov","var","cor","na.omit"))

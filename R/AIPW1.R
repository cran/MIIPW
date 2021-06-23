#' @title
#' Estimate of linear regression parameter from AIPW1
#' @description provides augmented inverse probability weighted estimates of parameters for linear regression model
#' of response variable using different covariance structure
#' @details It uses the inverse probability weighted method to reduce the bias
#' due to missing covariate in linear regression model. The estimating equation is
#' \deqn{\sum_{i=1}^{k}\sum_{j=1}^{n}(\frac{\delta_{ij}}{\pi_{ij}}S(Y_{ij},\mathbf{X}_{ij},\mathbf{X}'_{ij})+(1-\frac{\delta_{ij}}{\pi_{ij}})\phi(\mathbf{V}=\mathbf{v}))=0}
#' where \eqn{\delta_{ij}=1} if there is missing value in covariates and 0 otherwise,
#' \eqn{\mathbf{X}} is fully observed all subjects and \eqn{\mathbf{X}'} is partially missing,
#'  where \eqn{\mathbf{V}=(Y,\mathbf{X})}
#'
#'
#' @param data balanced longitudinal data set where each subject's outcome has been measured at same time points and number
#'  of visits for each patient is similar
#' covariance structure of the outcome variable like "unstuctured","compound","ToE","AR1"
#' ,"markov","independence"
#' @param cvstr "unstuctured","compound","ToE","AR1","markov",
#' "independence"
#' @param Dep column name of dependent variable in the dataset
#' @param Id column name of id of subjects in the dataset
#' @param Time column name of timepoints in the dataset
#' @param prob vector of probabilities not having missing value in covariate, which must be known by the
#' user from previous studies. In the example we consider 4 observations for each subject,so we create a vector of 4 and applied in the function.
#' @param m starting column number of covariates
#' @param n ending column number of covariates
#' @return estimated parameter value for multiple linear regression model,AIC,BIC
#' @import matlib
#' @examples
#'
#' data(srdata)
#' aipw1(cvstr="ToE",Dep="C6kine",Id="ID",Time="Visit",m=5,n=10,prob=rep(0.1,4),data=srdata)
#' @author Atanu Bhattacharjee,Bhrigu Kumar Rajbongshi and Gajendra K Vishwakarma
#' @export
#'

aipw1<-function(cvstr="unstructured",Dep,Id,Time,m,n,prob,data)
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
  ipw2<-function(dat){
    d<-function(x){
      dlt<-c()
      for(j in 1:ncol(x)){
        if(anyNA.data.frame(x[,j])==T){
          dlt[j]=0
        }
        else{
          dlt[j]=1
        }
      }
      dlt
    }
    dlt1=d(dat)
    dlt=d(t(dat))
    if(sum(dlt)!=0){
      cc=as.matrix(t(dat)[,which(dlt==1)])
      if(0%in%dlt){
        nc=as.matrix(t(dat)[,-which(dlt==1)])
        ncc=matrix(nc[which(dlt1==1),],nrow=sum(dlt1))
        ccc=matrix(cc[which(dlt1==1),],nrow=sum(dlt1))
        ph<-function(u,w){
          ##u is for including missing (y,z,w)
          ##w is for without missing (y,z,w)
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
        k<-ph(u=ncc,w=ccc)

        dsg=rbind(1,as.matrix(cc[-1,]))
        dsg<-as.matrix(dsg)
        fcc=data.frame(ncc,ccc)

        fl=ph(w=as.matrix(data.frame(ncc,ccc)),u=as.matrix(ccc))
        fl1=ph(w=as.matrix(ccc),u=as.matrix(ccc))
        invp<-c()
        for(i in 1:ncol(dsg)){
          invp[i]=prob[i]
          #length(fl1[[i]])/length(fl[[i]])

        }

        ##
        inv1<-unlist(lapply(invp,function(x){1/x}))
        inv2<-unlist(lapply(invp,function(x){1-1/x}))
        y1<-cc[1,]
        k2=list()
        for(j in 1:length(k)){
          if(length(k[[j]])!=0){
            k2[[j]]=as.matrix(dsg[,k[[j]]])%*%diag(inv2[k[[j]]],nrow=length(k[[j]]),ncol=length(k[[j]]))%*%Ginv(as.matrix(sig1[k[[j]],k[[j]]]))%*%as.matrix(y1[k[[j]]])/length(k[[j]])
          }else{
            k2[[j]]=0
          }
        }

        k3=list()
        for(j in 1:length(k)){
          if(length(k[[j]])!=0){
            k3[[j]]=as.matrix(dsg[,k[[j]]])%*%diag(inv2[k[[j]]],nrow=length(k[[j]]),ncol=length(k[[j]]))%*%Ginv(as.matrix(sig1[k[[j]],k[[j]]]))%*%t(as.matrix(dsg[,k[[j]]]))/length(k[[j]])
          }else{k3[[j]]=0}
        }
        if(length(k2)!=0){
          h2=Reduce("+",k2)
        }else{h2=0}
        if(length(k3)!=0){h3=Reduce("+",k3)}else{ h3=0}
        xy1<-as.matrix(dsg)%*%diag(inv1,ncol=ncol(as.matrix(dsg)),nrow=ncol(as.matrix(dsg)))%*%Ginv(as.matrix(sig1[which(dlt==1),which(dlt==1)]))%*%as.matrix(cc[1,])+h2
        #k<-list()
        #for(i in 1:ncol(dsg)){
        # k[[i]]<-dsg[,i]%*% t(dsg[,i])*inv1[i]
        #}
        xx1<-as.matrix(dsg)%*%diag(inv1,ncol=ncol(as.matrix(dsg)),nrow=ncol(as.matrix(dsg)))%*%Ginv(as.matrix(sig1[which(dlt==1),which(dlt==1)]))%*%t(as.matrix(dsg))+h3

      } else {
        dsg=rbind(1,as.matrix(cc[-1,]))
        xy1<-as.matrix(dsg)%*%Ginv(as.matrix(sig1[which(dlt==1),which(dlt==1)]))%*%as.matrix(cc[1,])
        xx1<-as.matrix(dsg)%*%Ginv(as.matrix(sig1[which(dlt==1),which(dlt==1)]))%*%t(as.matrix(dsg))
      }}else{xx1<-NULL;xy1<-NULL}
    ## beta value for identity link function
    list(xx1,xy1)
  }


  l1<-list()
  l2<-list()
  dat1<-data.frame(dat1)
  for(i in 1:length(unique(dat1$id))){

    l1[i]<-ipw2(subset(dat1,dat1$id==unique(dat1$id)[i])[,-c(1,2)])[1]
    l2[i]<-ipw2(subset(dat1,dat1$id==unique(dat1$id)[i])[,-c(1,2)])[2]
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
  #cat("\nEstimated parameters using IPW2 :\n")
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


source("my_methods/my_method.R")
source("R/scores.R")
source("R/amputation.R")
source("R/amputation_mar.R")


library(mice)
library(missForest)


######################amputation################################
ratio<-0.2
datasets<-c("pumadyn32nm", "scm1d", "scm20d", "parkinsons")
n<-2000

for (data in datasets){
  
  X0 <- readRDS(paste0("data/datasets/complete/", data, ".RDS"))
  N<-nrow(X0)
  
  ###Only choose numerical variables ####
  X0<-X0[, apply(X0,2, function(x) length(unique(x)))/N >0.7 ]
  
  
  set.seed(1)
  indextrain<-sample(1:N, size=n, replace=F)
  indextest<-sample(setdiff(1:N, indextrain), size=n, replace=F)
  
  Xtrain<-X0[indextrain,]
  Xtest<-X0[indextest,]
    

  d<-ncol(X0)
  npattern<-sample(3:max(round(n/100),3),size=1)
  
  patterns<-matrix(sample(c(0,1), size=npattern*(d-1), replace=T ), nrow=npattern, ncol=d-1, byrow = T)
  patterns<-cbind(patterns, rep(1, nrow(patterns)))
  
  ##add fully observed pattern
  if (all(rowSums(patterns)< d)){
    
    patterns<-rbind(patterns, rep(1,d))
    
  }
  
  
  X.NA <-mar(Xtrain,ratio=ratio, by.patterns=T, patterns=patterns)
  M<- is.na(X.NA)*1
  unique(M)
  nrow(unique(M))
  mean(is.na(X.NA))
  
  colnames(Xtest)<-paste0("X", 1:ncol(X0))
  colnames(X.NA)<-paste0("X", 1:ncol(X0))
  
  saveRDS(X.NA, file=paste0("results/amputedsplit/", "mar.", ratio, ".1.", data, ".RDS"))
  saveRDS(Xtest, file=paste0("data/datasets/split/", "test.", ratio, ".1.", data, ".RDS"))
}
#############################################################

datasets<-c("parkinsons")

ratio<-0.2
S<-1

methods<-c("missForest","mice_rf", "mice_cart") #"missForest","mice_rf"


imputations<-list()
ediffres<-matrix(NaN, nrow=length(datasets), ncol=length(methods)+1)

colnames(ediffres)<-c(methods, "wgf")
rownames(ediffres)<-datasets

###After amputation

for (data in datasets){
  
  set.seed(1)
  
  #Here we use Xtest
  X <- readRDS(paste0("data/datasets/split/", "test.", ratio, ".1.", data, ".RDS"))
  X.NA <- readRDS(paste0("results/amputedsplit/", "mar.", ratio, ".1.", data, ".RDS"))
  
  colnames(X)<-paste0("X", 1:ncol(X))
  colnames(X.NA)<-paste0("X", 1:ncol(X))
  
  
  if ("knn" %in% methods){  imputations[["knn"]]<-impute.knn(as.matrix(X.NA))$data}
  if ("missForest" %in% methods){imputations[["missForest"]]<-missForest(X.NA)$ximp}
  if ("mice_cart" %in% methods){  blub <- mice(X.NA, method = "cart", m = 1,remove.collinear = FALSE)
  imputations[["mice_cart"]]<-mice::complete(blub, action="all")[[1]]}
  if ("mice_rf" %in% methods){  blub <- mice(X.NA, method = "rf", m = 1, remove.collinear = FALSE)
  imputations[["mice_rf"]]<-mice::complete(blub, action="all")[[1]]}
  
  
  imputations[["wgf"]]<-readRDS(paste0("results/imputedsplit/", "wgf.mar.", ratio, ".1.", data, ".RDS"))
  
  n<-nrow(X)
  
  
  ##Continue here!!!
  ediff<-matrix(NaN, nrow=S, ncol=length(methods)+1)
  colnames(ediff)<-c(methods, "wgf")
  for (s in 1:S){
    
    set.seed(s+1)
    
    Xboot<-X #X[sample(1:n, size=n, replace = T),]
    colnames(Xboot)<-paste0("X",1:ncol(X))
    
    ediff[s,"wgf"]<- - energy_std(imputations[["wgf"]],Xboot) #-eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
    for (method in c(methods)){
      
      
      Ximp<-imputations[[method]]
      
      colnames(Ximp)<-paste0("X",1:ncol(X))
      ediff[s,method]<- -energy_std(Ximp,Xboot) #-eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
    }
    
  }
  
  ediffres[data,]<-colMeans(ediff)
}


# par(mfrow=c(3,1))
# 
# hist(Xwgfimpute[,1])
# hist(Xboot[,1])
# hist(Xmiceimp[,1])

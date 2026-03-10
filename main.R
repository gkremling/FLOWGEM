
source("my_methods/my_method.R")
source("R/scores.R")
source("R/amputation.R")
source("R/amputation_mar.R")


library(mice)
library(missForest)

######################amputation################################
ratio<-0.2
datasets<-c("fat", "windspeed", "airfoil_self_noise", "enb", "scm1d", "yeast", "stock", "energy", "wine", "gas", "solar", "forest", "housing", "pumadyn32nm", "parkinsons")


for (data in datasets){
  
X <- readRDS(paste0("data/datasets/complete/", data, ".RDS"))
  

set.seed(1)

#X.NA <-mar(X,ratio=ratio)
 
n<-nrow(X)
d<-ncol(X)
npattern<-sample(3:max(round(n/100),3),size=1)

patterns<-matrix(sample(c(0,1), size=npattern*(d-1), replace=T ), nrow=npattern, ncol=d-1, byrow = T)
patterns<-cbind(patterns, rep(1, nrow(patterns)))

##add fully observed pattern
if (all(rowSums(patterns)< d)){
  
  patterns<-rbind(patterns, rep(1,d))
  
}


X.NA <-mar(X,ratio=ratio, by.patterns=T, patterns=patterns)
M<- is.na(X.NA)*1
unique(M)
cat(nrow(unique(M)))
mean(is.na(X.NA))

  
  
saveRDS(X.NA, file=paste0("results/amputed/", "mar.", ratio, ".1.", data, ".RDS"))
}
#############################################################

datasets<-c("fat", "enb", "windspeed")

ratio<-0.2
S<-10

methods<-c("missForest", "mice_cart", "mice_rf")


imputations<-list()
ediffres<-matrix(NaN, nrow=length(datasets), ncol=length(methods)+1)

colnames(ediffres)<-c(methods, "wgf")
rownames(ediffres)<-datasets

###After amputation

for (data in datasets){

set.seed(1)

X <- readRDS(paste0("data/datasets/complete/", data, ".RDS"))
X.NA <- readRDS(paste0("results/amputed/", "mar.", ratio, ".1.", data, ".RDS"))



if ("knn" %in% methods){  imputations[["knn"]]<-impute.knn(as.matrix(X.NA))$data}
if ("missForest" %in% methods){imputations[["missForest"]]<-missForest(X.NA)$ximp}
if ("mice_cart" %in% methods){  blub <- mice(X.NA, method = "cart", m = 1, remove.collinear = FALSE)
imputations[["mice_cart"]]<-mice::complete(blub, action="all")[[1]]}
if ("mice_rf" %in% methods){  blub <- mice(X.NA, method = "rf", m = 1, eps = 0, remove.collinear = FALSE)
imputations[["mice_rf"]]<-mice::complete(blub, action="all")[[1]]}


imputations[["wgf"]]<-readRDS(paste0("results/imputed/", "wgf.mar.", ratio, ".1.", data, ".RDS"))

n<-nrow(X)


##Continue here!!!
ediff<-matrix(NaN, nrow=S, ncol=length(methods)+1)
colnames(ediff)<-c(methods, "wgf")
for (s in 1:S){

set.seed(s+1)
  
Xboot<-X[sample(1:n, size=n, replace = T),]
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

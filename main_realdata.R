

source("R/scores.R")
source("R/amputation.R")
source("R/amputation_mar.R")


library(mice)
library(missForest)


######################amputation################################
ratio <- 0.6
#datasets<-c("pumadyn32nm", "scm1d", "scm20d", "parkinsons")


datasets <- c(
  "pumadyn32nm",
  "scm1d",
  "scm20d",
  "parkinsons",
  "allergens",
  "concrete",
  "windspeed",
  "forest",
  "housing",
  "stock"
)

table_dataXNA <- data.frame(dataset = character(),
                            n = integer(),
                            d = integer())
table_dataXtest <- data.frame(dataset = character(),
                              n = integer(),
                              d = integer())



for (data in datasets) {
  X0 <- readRDS(paste0("data/datasets/complete/", data, ".RDS"))
  N <- nrow(X0)
  
  
  ###Only choose numerical variables ####
  X0 <- X0[, apply(X0, 2, function(x)
    length(unique(x))) / N > 0.1 , drop = F]
  
  if (ncol(X0) <= 5) {
    next
  }
  
  
  set.seed(1)
  n <- round(N / 2) - 1
  indextrain <- sample(1:N, size = n , replace = F)
  indextest <- sample(setdiff(1:N, indextrain),
                      size = n,
                      replace = F)
  
  Xtrain <- X0[indextrain, ]
  Xtest <- X0[indextest, ]
  
  
  d <- ncol(X0)
  npattern <- sample(3:max(round(n / 100), 3), size = 1)
  
  patterns <- matrix(
    sample(c(0, 1), size = npattern * (d - 1), replace = T),
    nrow = npattern,
    ncol = d - 1,
    byrow = T
  )
  patterns <- cbind(patterns, rep(1, nrow(patterns)))
  
  ##add fully observed pattern
  if (all(rowSums(patterns) < d)) {
    patterns <- rbind(patterns, rep(1, d))
    
  }
  
  
  # Separate the all-ones (complete) pattern
  all_ones <- apply(patterns, 1, function(r)
    all(r == 1))
  incomplete_patterns <- patterns[!all_ones, , drop = FALSE]
  n_pats <- nrow(incomplete_patterns)
  
  n_total <- nrow(Xtrain)
  n_complete <- round(n_total / (n_pats + 1))  # ~1/21 of rows stay complete
  
  # Randomly assign rows to the "complete" pattern
  complete_idx <- sample(n_total, n_complete)
  Xtrain_to_ampute <- Xtrain[-complete_idx, , drop = FALSE]
  
  # Ampute the remaining rows with equal frequency across incomplete patterns
  tmp <- ampute(
    Xtrain_to_ampute,
    patterns  = incomplete_patterns,
    freq      = rep(1 / n_pats, n_pats),
    prop      = 0.99,
    # ampute all but 1/(n_pats+1) fraction
    mech      = "MAR",
    bycases   = TRUE
  )
  
  # Reassemble, preserving original row order
  X.NA <- rbind(Xtrain[complete_idx, , drop = FALSE], tmp$amp)
  X.NA <- X.NA[order(c(complete_idx, seq_len(nrow(Xtrain))[-complete_idx])), , drop = FALSE]
  
  M <- is.na(X.NA) * 1
  unique(M)
  nrow(unique(M))
  mean(is.na(X.NA))
  
  nrow(X.NA[complete.cases(X.NA), ])
  nrow(X.NA[!complete.cases(X.NA), ])
  
  colnames(Xtest) <- paste0("X", 1:ncol(X0))
  colnames(X.NA) <- paste0("X", 1:ncol(X0))
  
  saveRDS(X.NA,
          file = paste0("results/amputedsplit/", "mar.", ratio, ".1.", data, ".RDS"))
  saveRDS(Xtest,
          file = paste0("data/datasets/split/", "test.", ratio, ".1.", data, ".RDS"))
  
  
  table_dataXNA <- rbind(table_dataXNA, data.frame(
    dataset = data,
    n = nrow(X.NA),
    d = ncol(X.NA)
  ))
  table_dataXtest <- rbind(table_dataXtest, data.frame(
    dataset = data,
    n = nrow(Xtest),
    d = ncol(Xtest)
  ))
  
}



# Print in R
print(table_dataXNA)

# LaTeX output
library(xtable)
xtable(table_dataXNA,
       caption = "Dataset dimensions after filtering",
       label = "tab:datasets",
       digits = 0)


# Print in R
print(table_dataXtest)

# LaTeX output
library(xtable)
xtable(
  table_dataXtest,
  caption = "Dataset dimensions after filtering",
  label = "tab:datasets",
  digits = 0
)


#############################################################
# 
# #######################Gaussian EM X0#######################################
# 
# 
# library(EMgaussian)
# 
# datasets <- c("pumadyn32nm", "scm20d", "parkinsons", "scm1d")
# 
# 
# ratio <- 0.6
# S <- 1
# 
# 
# imputations <- list()
# ediffres <- matrix(NaN, nrow = length(datasets), ncol = length(methods))
# 
# rownames(ediffres) <- datasets
# 
# 
# for (data in datasets) {
#   set.seed(1)
#   
#   X <- readRDS(paste0("data/datasets/split/", "test.", ratio, ".1.", data, ".RDS"))
#   X.NA <- readRDS(paste0("results/amputedsplit/", "mar.", ratio, ".1.", data, ".RDS"))
#   
#   tmp <- em.cov(X.NA)
#   X0 <- rmvnorm(n = nrow(X.NA),
#                 mu = tmp$mu,
#                 sigma = tmp$S)
#   
#   saveRDS(X0,
#           file = paste0("results/imputedsplit/", "X0.", ratio, ".1.", data, ".RDS"))
#   
#   ediffres[data, ] <- -energy_std(X0, X)
#   
#   
# }
# 
# 
# #######################Gaussian EM X0#######################################







######################Run mice methods#######################################

#datasets<-c("pumadyn32nm", "scm20d") #"gas", "scm1d","parkinsons"

datasets <- c("parkinsons") #"gas", "scm1d","parkinsons"

ratio <- 0.6
S <- 1

methods <- c("mice_norm") #"missForest","mice_rf", "mice_cart",


imputations <- list()



for (data in datasets) {
  set.seed(1)
  
  #Here we use Xtest
  X <- readRDS(paste0("data/datasets/split/", "test.", ratio, ".1.", data, ".RDS"))
  X.NA <- readRDS(paste0("results/amputedsplit/", "mar.", ratio, ".1.", data, ".RDS"))
  
  colnames(X) <- paste0("X", 1:ncol(X))
  colnames(X.NA) <- paste0("X", 1:ncol(X))
  
  
  if ("missForest" %in% methods) {
    imputations[["missForest"]] <- missForest(X.NA)$ximp
    saveRDS(
      object = imputations[["missForest"]],
      paste0(
        "results/",
        "imputations",
        "missForest",
        ratio,
        ".",
        data,
        ".RDS"
      )
    )
  }
  if ("mice_cart" %in% methods) {
    blub <- mice(
      X.NA,
      method = "cart",
      m = 1,
      remove.collinear = FALSE,
      eps = 0
    )
    imputations[["mice_cart"]] <- mice::complete(blub, action = "all")[[1]]
    saveRDS(
      object = imputations[["mice_cart"]],
      paste0(
        "results/",
        "imputations",
        "mice_cart",
        ratio,
        ".",
        data,
        ".RDS"
      )
    )
  }
  if ("mice_rf" %in% methods) {
    blub <- mice(
      X.NA,
      method = "rf",
      m = 1,
      remove.collinear = FALSE,
      eps = 0
    )
    imputations[["mice_rf"]] <- mice::complete(blub, action = "all")[[1]]
    
    saveRDS(object = imputations[["mice_rf"]],
            paste0(
              "results/",
              "imputations",
              "mice_rf",
              ratio,
              ".",
              data,
              ".RDS"
            ))
  }
  
  if ("mice_norm" %in% methods) {
    blub <- mice(
      X.NA,
      method = "norm.nob",
      m = 1,
      remove.collinear = FALSE,
      eps = 0
    )
    imputations[["mice_norm"]] <- mice::complete(blub, action = "all")[[1]]
    
    saveRDS(
      object = imputations[["mice_norm"]],
      paste0(
        "results/",
        "imputations",
        "mice_norm",
        ratio,
        ".",
        data,
        ".RDS"
      )
    )
  }
  
  # saveRDS(imputations, file = paste0("results/",  "imputations.mar.", ratio,".", data, paste0(methods, collapse="_"), ".RDS"))
  #
  #
  # n<-nrow(X)
  #
  # ediff<-matrix(NaN, nrow=S, ncol=length(methods)+1)
  # colnames(ediff)<-c(methods, "wgf")
  # for (s in 1:S){
  #
  #   set.seed(s+1)
  #
  #   Xboot<-X #X[sample(1:n, size=n, replace = T),]
  #   colnames(Xboot)<-paste0("X",1:ncol(X))
  #
  #   ediff[s,"wgf"]<- - energy_std(imputations[["wgf"]],Xboot) #-eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
  #   for (method in c(methods)){
  #
  #
  #     Ximp<-imputations[[method]]
  #
  #     colnames(Ximp)<-paste0("X",1:ncol(X))
  #     ediff[s,method]<- -energy_std(Ximp,Xboot) #-eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
  #   }
  #
  # }
  #
  # ediffres[data,]<-colMeans(ediff)
}


######################Run mice methods#######################################

######################Evaluate imputations#######################################

######################################
###After saving imputations ####
###################################



#datasets<-c("pumadyn32nm", "scm20d", "parkinsons", "scm1d", "allergens", "concrete", "windspeed", "forest", "housing", "stock")

datasets <- c(
  "pumadyn32nm",
  "scm20d",
   "scm1d",
  "parkinsons",
  "allergens",
  "concrete",
  "windspeed",
  "forest",
  "housing",
  "stock"
)


datasets <- c(
  "pumadyn32nm",
  "scm20d",
  "parkinsons",
  "allergens",
  "concrete",
  "windspeed",
   "housing",
  "forest",
  "stock"
)




table_dataXtest <- data.frame(dataset = character(),
                              n = integer(),
                              d = integer())



ratio <- 0.6
S <- 1

## Ours: "wgf_new"
## Python methods: "gain", "hyperimpute" 'miri', "knn", "mice"
## MICE R methods: 'mice_cart', "mice_rf", "mice_norm"

methods <- c('flowgem', "gain", "hyperimpute", "mice", "miri")#

imputations <- list()
ediffres <- matrix(NaN, nrow = length(datasets), ncol = length(methods))

colnames(ediffres) <- c(methods)
rownames(ediffres) <- datasets


for (data in datasets) {
  X <- readRDS(paste0("data/datasets/split/", "test.", ratio, ".1.", data, ".RDS"))
  
  sig <- cor(X)
  print(max(abs(sig[upper.tri(sig)])))
  
  n <- nrow(X)
  d<- ncol(X)
  
  for (method in methods) {
    imputations[[method]] <- readRDS(paste0("results/", "imputations", method, ratio, ".", data, ".RDS"))
    
  }
  
  ediff <- matrix(NaN, nrow = S, ncol = length(methods))
  colnames(ediff) <- c(methods)
  for (s in 1:S) {
    set.seed(s + 1)
    
    Xboot <- as.matrix(X) #X[sample(1:n, size=n, replace = T),]
    colnames(Xboot) <- paste0("X", 1:ncol(X))
    
    for (method in c(methods)) {
      Ximp <- imputations[[method]]
      Ximp <- as.data.frame(apply(Ximp, 2, as.numeric))
      
      colnames(Ximp) <- paste0("X", 1:ncol(X))
      if (any(is.nan(Ximp))==F){
      ediff[s, method] <- energy_std(Xboot, Ximp) #eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
      }else{
        ediff[s, method]<-NA
      }
      
    }
    
  }
  
  ediffres[data, ] <- colMeans(ediff)
  
  table_dataXtest <- rbind(table_dataXtest, 
                           data.frame(
                             dataset = data,
                             n = n,
                             d = d
                           ))
  
  
}

colnames(ediffres)[colnames(ediffres) == "wgf"] <- "\\name"

tablefinal<-cbind(table_dataXtest, ediffres)
tablefinal <- tablefinal[order(tablefinal$n, decreasing = TRUE), ]


rownames(tablefinal)<-NULL

# LaTeX output
library(xtable)
xtable(
  tablefinal,
  caption = "Datasets and results",
  label = "tab:datasets",
  digits = 2,
  include.rownames       = FALSE
)



######################Evaluate imputations#######################################


######################Evaluate specific data sets#######################################

######################################
###After saving imputations ####
###################################



datasets <- c(
  "pumadyn32nm",
  "scm20d",
  "scm1d",
  "parkinsons",
  "allergens",
  "concrete",
  "windspeed",
  "forest",
  "housing",
  "stock"
)



dataset <- "parkinsons"

ratio <- 0.6
S <- 1

## Ours: "wgf_new"
## Python methods: "gain", "hyperimpute" 'miri', "knn", "mice"
## MICE R methods: 'mice_cart', "mice_rf", "mice_norm"

methods <- c('wgf_new', "gain", "hyperimpute", "mice")

imputations <- list()
ediffres <- matrix(NaN, nrow = length(datasets), ncol = length(methods))

colnames(ediffres) <- c(methods)
rownames(ediffres) <- datasets


X0 <- readRDS(paste0("data/datasets/complete/", data, ".RDS"))
N <- nrow(X0)


###Only choose numerical variables ####
X0 <- X0[, apply(X0, 2, function(x)
  length(unique(x))) / N > 0.1 , drop = F]
set.seed(1)
n <- round(N / 2) - 1
indextrain <- sample(1:N, size = n , replace = F)
indextest <- sample(setdiff(1:N, indextrain), size = n, replace = F)



Xtrain <- X0[indextrain, ]
Xtest <- X0[indextest, ]

#Xtest2 <- readRDS(paste0("data/datasets/split/", "test.", ratio, ".1.", data, ".RDS"))
X.NA <- readRDS(paste0("results/amputedsplit/", "mar.", ratio, ".1.", data, ".RDS"))

sig <- cor(X)
print(max(abs(sig[upper.tri(sig)])))

n <- nrow(X)

for (method in methods) {
  imputations[[method]] <- readRDS(paste0("results/", "imputations", method, ratio, ".", data, ".RDS"))
  
}

head(Xtrain)
head(imputations[["wgf_new"]])


edifftest <- matrix(NaN, nrow = S, ncol = length(methods))
colnames(edifftest) <- c(methods)
for (s in 1:S) {
  set.seed(s + 1)
  
  Xboot <- as.matrix(Xtest) #X[sample(1:n, size=n, replace = T),]
  colnames(Xboot) <- paste0("X", 1:ncol(X))
  
  for (method in c(methods)) {
    Ximp <- imputations[[method]]
    Ximp <- as.data.frame(apply(Ximp, 2, as.numeric))
    
    colnames(Ximp) <- paste0("X", 1:ncol(X))
    edifftest[s, method] <- -energy_std(Ximp, Xboot) #-eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
  }
  
}

edifftrain <- matrix(NaN, nrow = S, ncol = length(methods))
colnames(edifftrain) <- c(methods)
for (s in 1:S) {
  set.seed(s + 1)
  
  Xboot <- as.matrix(Xtrain) #X[sample(1:n, size=n, replace = T),]
  colnames(Xboot) <- paste0("X", 1:ncol(X))
  
  for (method in c(methods)) {
    Ximp <- imputations[[method]]
    Ximp <- as.data.frame(apply(Ximp, 2, as.numeric))
    
    colnames(Ximp) <- paste0("X", 1:ncol(X))
    edifftrain[s, method] <- -energy_std(Ximp, Xboot) #-eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
  }
  
}



######################Evaluate specific data sets#######################################

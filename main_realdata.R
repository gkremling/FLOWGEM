###################### AMPUTATION ######################

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
  ## Read in the dataset
  X0 <- readRDS(paste0("datasets/complete/", data, ".RDS"))
  N <- nrow(X0)
  
  
  ## Filter for continuous variables
  X0 <- X0[, apply(X0, 2, function(x)
    length(unique(x))) / N > 0.1 , drop = F]
  
  if (ncol(X0) <= 5) {
    next
  }
  
  
  set.seed(1)
  
  ## Split dataset into train (amputed) and test (fully observed)
  n <- round(N / 2) - 1
  indextrain <- sample(1:N, size = n , replace = F)
  indextest <- sample(setdiff(1:N, indextrain),
                      size = n,
                      replace = F)
  
  Xtrain <- X0[indextrain, ]
  Xtest <- X0[indextest, ]
  
  ## Sample random number of missingess patterns
  d <- ncol(X0)
  npattern <- sample(3:max(round(n / 100), 3), size = 1)
  
  ## Sample missigness patterns (excl. all zeros)
  patterns <- matrix(
    sample(c(0, 1), size = npattern * (d - 1), replace = T),
    nrow = npattern,
    ncol = d - 1,
    byrow = T
  )
  patterns <- cbind(patterns, rep(1, nrow(patterns)))
  
  ## Add fully observed pattern if not included yet
  if (all(rowSums(patterns) < d)) {
    patterns <- rbind(patterns, rep(1, d))
    
  }
  
  
  ## Separate the all-ones (complete) pattern
  all_ones <- apply(patterns, 1, function(r)
    all(r == 1))
  incomplete_patterns <- patterns[!all_ones, , drop = FALSE]
  n_pats <- nrow(incomplete_patterns)
  
  ## Randomly assign rows to the "complete" pattern
  n_total <- nrow(Xtrain)
  n_complete <- round(n_total / (n_pats + 1))
  complete_idx <- sample(n_total, n_complete)
  Xtrain_to_ampute <- Xtrain[-complete_idx, , drop = FALSE]
  
  ## Ampute the remaining rows with equal frequency across incomplete patterns
  tmp <- ampute(
    Xtrain_to_ampute,
    patterns  = incomplete_patterns,
    freq      = rep(1 / n_pats, n_pats),
    prop      = 0.99,
    # ampute all but 1/(n_pats+1) fraction
    mech      = "MAR",
    bycases   = TRUE
  )
  
  ## Reassemble, preserving original row order
  X.NA <- rbind(Xtrain[complete_idx, , drop = FALSE], tmp$amp)
  X.NA <- X.NA[order(c(complete_idx, seq_len(nrow(Xtrain))[-complete_idx])), , drop = FALSE]
  
  ## Print some statistics about the missingness
  M <- is.na(X.NA) * 1
  unique(M)
  nrow(unique(M))
  mean(is.na(X.NA))
  
  nrow(X.NA[complete.cases(X.NA), ])
  nrow(X.NA[!complete.cases(X.NA), ])
  
  ## Save the resulting datatset with missing values
  colnames(Xtest) <- paste0("X", 1:ncol(X0))
  colnames(X.NA) <- paste0("X", 1:ncol(X0))
  
  saveRDS(X.NA,
          file = paste0("datasets/split_amputed/", "mar.", data, ".RDS"))
  saveRDS(Xtest,
          file = paste0("datasets/split_test/", "test.", data, ".RDS"))
  
  ## Create a table with n and d for each dataset
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






###################### EVALUATE IMPUTATIONS ######################

#### only reasonable to run AFTER running main_realdata.ipynb ####

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

table_dataXtest <- data.frame(dataset = character(),
                              n = integer(),
                              d = integer())

## Ours: "flowgem"
## Python methods: "gain", "hyperimpute" "miri", "mice", "knewimp", "missdiff", "bayes"

methods <- c('flowgem', "gain", "hyperimpute", "mice", "miri")#

imputations <- list()
ediffres <- matrix(NaN, nrow = length(datasets), ncol = length(methods))

colnames(ediffres) <- c(methods)
rownames(ediffres) <- datasets


energy_dist<-function (X, X_imp){ 
  energy::eqdist.e(rbind(X, X_imp), c(nrow(X), nrow(X_imp)))
}


energy_std <- function(original_data, imputed_data) {
  
  scaled_original <- scale(original_data)
  
  scaled_imputed <- sapply(1:ncol(imputed_data), function(i) {
    (imputed_data[, i] - attr(scaled_original, "scaled:center")[i])/ 
      attr(scaled_original, "scaled:scale")[i]
  })
  
  return(energy_dist(X = scaled_original, X_imp = scaled_imputed))
}


for (data in datasets) {
  X <- readRDS(paste0("datasets/split_test/", "test.", data, ".RDS"))
  
  sig <- cor(X)
  print(max(abs(sig[upper.tri(sig)])))
  
  n <- nrow(X)
  d <- ncol(X)
  
  for (method in methods) {
    imputations[[method]] <- readRDS(paste0("results/", "imputations.", method, ".", data, ".RDS"))
  }
  
  set.seed(2)
  
  Xboot <- as.matrix(X) #X[sample(1:n, size=n, replace = T),]
  colnames(Xboot) <- paste0("X", 1:ncol(X))
  
  for (method in c(methods)) {
    Ximp <- imputations[[method]]
    Ximp <- as.data.frame(apply(Ximp, 2, as.numeric))
    
    colnames(Ximp) <- paste0("X", 1:ncol(X))
    if (any(is.nan(Ximp))==F){
      ediffres[data, method] <- energy_std(Xboot, Ximp) #eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
    } else{
      ediffres[data, method] <- NA
    }
  }
  
  table_dataXtest <- rbind(table_dataXtest, 
                           data.frame(
                             dataset = data,
                             n = n,
                             d = d
                           ))
  
  
}

tablefinal <- cbind(table_dataXtest, ediffres)
tablefinal <- tablefinal[order(tablefinal$n, decreasing = TRUE), ]


rownames(tablefinal) <- NULL

# LaTeX output
library(xtable)
xtable(
  tablefinal,
  caption = "Datasets and results",
  label = "tab:datasets",
  digits = 2,
  include.rownames       = FALSE
)

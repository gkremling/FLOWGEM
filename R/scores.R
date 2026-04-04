
# ------------------------------------------------------------- IMPLEMENTATIONS



energy_dist<-function (X, X_imp){ 
  energy::eqdist.e(rbind(X, X_imp), c(nrow(X), nrow(X_imp)))
}


energy_std <- function(original_data, imputed_data, ...) {
  
  scaled_original <- scale(original_data)
  
  scaled_imputed <- sapply(1:ncol(imputed_data), function(i) {
    (imputed_data[, i] - attr(scaled_original, "scaled:center")[i])/ 
      attr(scaled_original, "scaled:scale")[i]
  })
  
  return(energy_dist(X = scaled_original, X_imp = scaled_imputed))
}

mae <- function(original_data, imputed_data, amputed_data) {
  observed <- original_data[is.na(amputed_data)]
  imputed <- imputed_data[is.na(amputed_data)]
  mean(abs(observed - imputed))
}

rmse <- function(original_data, imputed_data, amputed_data) {
  observed <- original_data[is.na(amputed_data)]
  imputed <- imputed_data[is.na(amputed_data)]
  sqrt(mean((observed - imputed)^2))
}

nrmse <- function(original_data, imputed_data, amputed_data) {
  observed <- original_data[is.na(amputed_data)]
  imputed <- imputed_data[is.na(amputed_data)]
  sqrt(mean((observed - imputed)^2) / var(observed))
}

energy <- function(original_data, imputed_data, ...) {
  safe_score(
    energy_dist(X = original_data, X_imp = imputed_data),
    original_data,
    imputed_data
  )  
}




feature_wise_wasserstein <- function(original_data, imputed_data, ...) {
  safe_score({
    mean(sapply(1:ncol(original_data), function(ith_col) 
      transport::wasserstein1d(as.matrix(original_data)[, ith_col], 
                               as.matrix(imputed_data)[, ith_col])), na.rm = TRUE)
  }, original_data, imputed_data)
}

KLD <- function(original_data, imputed_data, ...) {
  safe_score({ 
    median(FNN::KL.dist(as.matrix(original_data), 
                        as.matrix(imputed_data), 
                        floor(sqrt(nrow(original_data)))), na.rm = TRUE)
  }, original_data, imputed_data)
}

entropic_wasserstein <- function(original_data, imputed_data, ...) {
  safe_score({ 
    T4transport::sinkhorn(as.matrix(original_data), 
                          as.matrix(imputed_data), 
                          lambda = 1)$distance
  }, original_data, imputed_data)
}


sliced_wasserstein <- function(original_data, imputed_data, ...) {
  safe_score({ 
    T4transport::swdist(as.matrix(original_data), 
                        as.matrix(imputed_data))$distance
  }, original_data, imputed_data)
}

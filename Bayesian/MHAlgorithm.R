

# -----------------------------------------------------------------------------
# 2.  Log-likelihood: x_i | Sigma ~ N(mu, Sigma),  mu fixed
#     S  = pre-computed scatter matrix = sum_i (x_i - mu)(x_i - mu)^T
# -----------------------------------------------------------------------------

# Log-prior density
log_prior_Sigma <- function(Sigma, Psi0, nu0) {
  d <- nrow(Sigma)
  
  # # Inverse-Wishart density 
  # log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  # 
  # log_dens <- -(nu0 + d + 1)/2 * log_det_Sigma
  # log_dens <- log_dens - 0.5 * sum(diag(Psi0 %*% solve(Sigma)))
  # ## Normalizing constants:
  # log_dens<- log_dens  + nu0 / 2 * determinant(Psi0, logarithm = TRUE)$modulus
  
  log_dens<-log(diwish(Sigma,  v = nu0, S = Psi0))
  
  return(log_dens)
}



log_lik_Sigma <- function(Sigma, z, mu, data) {
  #log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  #-n / 2 * log_det_Sigma - 0.5 * sum(diag(solve(Sigma) %*% S))
  #sum(dmvnorm(x=X, mean=mu, sigma=Sigma, log=T))
  index<-unique(z)
  loglik<-0
  for (k in index){
    loglik<-loglik+ loglikGauss(data[z==k,, drop=F], Sigma, mu[k,])
  }  
    
  return(loglik)
}



# # -----------------------------------------------------------------------------
# # 3.  Log unnormalised target:  log pi(Sigma) + log f(x | Sigma)
# # -----------------------------------------------------------------------------
# log_target <- function(Sigma, mu, data, Psi0, nu0) {
#   # Guard against non-positive-definite proposals
#   eigs <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
#   if (any(eigs <= 0)) return(-Inf)
#   log_prior_Sigma(Sigma, Psi0, nu0) + log_lik_Sigma(Sigma, mu, data)
# }



# -----------------------------------------------------------------------------
# 2.  Log-likelihood: x_i | Sigma ~ N(mu, Sigma),  mu fixed
#     S  = pre-computed scatter matrix = sum_i (x_i - mu)(x_i - mu)^T
# -----------------------------------------------------------------------------




# Sigma -> unconstrained vector (log of diagonal, raw off-diagonal)
Sigma_to_theta <- function(Sigma) {
  L <- t(chol(Sigma))                      # lower Cholesky
  d <- nrow(L)
  theta <- numeric(d * (d + 1) / 2)
  idx <- 1
  for (j in 1:d) {
    for (i in j:d) {
      theta[idx] <- if (i == j) log(L[i, j]) else L[i, j]
      idx <- idx + 1
    }
  }
  theta
}

# unconstrained vector -> Sigma
theta_to_Sigma <- function(theta, d) {
  L <- matrix(0, d, d)
  idx <- 1
  for (j in 1:d) {
    for (i in j:d) {
      L[i, j] <- if (i == j) exp(theta[idx]) else theta[idx]
      idx <- idx + 1
    }
  }
  L %*% t(L)
}

# Log Jacobian of the transformation: sum of (d - j + 1) * log(L[j,j])
# needed to correct the prior from Sigma-space to theta-space
log_jacobian <- function(theta, d) {
  # diagonal entries are at positions 1, 2, ..., d in the packed vector
  # (i.e. indices 1, 1+d, 1+d+(d-1), ... but easier to re-extract)
  log_diag <- numeric(d)
  idx <- 1
  for (j in 1:d) {
    log_diag[j] <- theta[idx]   # log(L[j,j])
    idx <- idx + (d - j + 1)
  }
  sum((d - seq_len(d) + 2) * log_diag)
}



# Log target on the unconstrained scale (prior x likelihood x Jacobian)
log_target_theta <- function(theta, z, mu, X, Psi0, nu0) {
  d <- ncol(X)
  Sigma <- theta_to_Sigma(theta, d)
  eigs  <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigs <= 0)) return(-Inf)
  log_prior_Sigma(Sigma, Psi0, nu0) +
    log_lik_Sigma(Sigma,z, mu, X) +
    log_jacobian(theta, d)
}

# --- MH sampler on unconstrained scale ---
MHsamplerSigma <- function(data, z, mu, Psi0, nu0,
                           n_iter    = 10000,
                           burnin    = 2000,
                           Sigma_init = NULL) {
  
  
  # mu is a matrix of dim (max(z), d)
  
  n <- nrow(data); d <- ncol(data)
  k <- d * (d + 1) / 2          # dimension of unconstrained vector
  
  step_sd <- min(2.38 / sqrt(d * (d + 1) / 2), 0.001)
  

  
  # Initialise at MLE
  if (is.null(Sigma_init)) {
    result <- optimize_mu_Sigma_constrained(X.NA)
    #muinit<- matrix(result$mu_hat, nrow=max(z), ncol=d, byrow = T)
    Sigma_init<-result$Sigma_hat
  }
  
  theta_cur <- Sigma_to_theta(Sigma_init)
  ltgt_cur  <- log_target_theta(theta_cur, z, mu, data, Psi0, nu0)
  #(theta, z, mu, X, Psi0, nu0)
  
  samples  <- array(NA_real_, dim = c(d, d, n_iter))
  log_tgt  <- numeric(n_iter)
  accepted <- logical(n_iter)
  
  for (t in seq_len(n_iter)) {
    
    # Symmetric Gaussian random walk — correction term = 0
    theta_prop <- theta_cur + rnorm(k, sd = step_sd)
    ltgt_prop  <- log_target_theta(theta_prop,z, mu, data, Psi0, nu0)
    
    log_alpha <- ltgt_prop - ltgt_cur     # no correction term needed!
    
    if (log(runif(1)) < log_alpha) {
      theta_cur <- theta_prop
      ltgt_cur  <- ltgt_prop
      accepted[t] <- TRUE
    }
    
    samples[, , t] <- theta_to_Sigma(theta_cur, d)
    log_tgt[t]     <- ltgt_cur
  }
  
  ar <- mean(accepted[(burnin + 1):n_iter])
  #cat(sprintf("Post-burn-in acceptance rate: %.1f%%\n", 100 * ar))
  
  list(samples     = samples,
       log_target  = log_tgt,
       accepted    = accepted,
       accept_rate = ar,
       burnin      = burnin)
}


# ##Usage
# n_iter    = 200
# burnin    = 50
# fit<-MHsamplerSigma(X.NA,z=sample(c(1,2), size=n, replace=T),mu= matrix(mustar, nrow=2, ncol=d), Psi0=diag(ncol(data)), nu0=ncol(data)+1, n_iter=n_iter, burnin=burnin)
# post_samples <- fit$samples[, , (burnin + 1):n_iter]
# Sigma_hat    <- apply(post_samples, c(1, 2), mean)



# =============================================================================
# MH Sampler for mu | X, Sigma
# =============================================================================

log_prior_mu <- function(mu, mu0, tau0_sq) {
  d <- length(mu)
  # N(mu0, tau0^2 * I)  — keeping normalising constant for correctness
  Rfast::dmvnorm(x=mu, mu=mu0, sigma=tau0_sq*diag(d), log=T)
}

log_lik_mu <- function(mu, Sigma, data) {
  loglikGauss(data, Sigma, mu)
}

log_target_mu <- function(mu, Sigma, data, mu0, tau0_sq) {
  log_prior_mu(mu, mu0, tau0_sq) + log_lik_mu(mu, Sigma, data)
}

MHsamplerMu <- function(data, Sigma, mu0, tau0_sq,
                        step_sd   = NULL,
                        n_iter    = 10000,
                        burnin    = 2000,
                        mu_init   = NULL) {
  
  n <- nrow(data)
  d <- ncol(data)
  stopifnot(burnin < n_iter)
  
  # Default step size: Roberts & Rosenthal optimal scaling
  if (is.null(step_sd)) step_sd <- min(2.38 / sqrt(d), 0.001)
  
  # Initialise at MLE = column means
  if (is.null(mu_init)){
    result <- optimize_mu_Sigma_constrained(X.NA)
    mu_init<- result$mu_hat
  } 
  
  # Storage
  samples  <- matrix(NA_real_, nrow = n_iter, ncol = d)
  log_tgt  <- numeric(n_iter)
  accepted <- logical(n_iter)
  
  mu_cur   <- mu_init
  ltgt_cur <- log_target_mu(mu_cur, Sigma, data, mu0, tau0_sq)
  
  for (t in seq_len(n_iter)) {
    
    # Symmetric Gaussian random walk — no correction term needed
    mu_prop   <- mu_cur + rnorm(d, sd = step_sd)
    ltgt_prop <- log_target_mu(mu_prop, Sigma, data, mu0, tau0_sq)
    
    log_alpha <- ltgt_prop - ltgt_cur
    
    if (log(runif(1)) < log_alpha) {
      mu_cur   <- mu_prop
      ltgt_cur <- ltgt_prop
      accepted[t] <- TRUE
    }
    
    samples[t, ] <- mu_cur
    log_tgt[t]   <- ltgt_cur
  }
  
  ar <- mean(accepted[(burnin + 1):n_iter])
  #cat(sprintf("Post-burn-in acceptance rate: %.1f%%\n", 100 * ar))
  
  list(samples     = samples,
       log_target  = log_tgt,
       accepted    = accepted,
       accept_rate = ar,
       burnin      = burnin)
}


# ##Usage
# n_iter    = 200
# burnin    = 50
# fit<-MHsamplerMu(X.NA, Sigmastar, mu0=rep(0,ncol(X.NA)), tau0_sq=2, n_iter=n_iter, burnin=burnin)
# post_samples <- fit$samples
# mu_hat    <- colMeans(fit$samples)


# fixed p setting
source("utils.R")


if (!require(combinat)) {
  install.packages("combinat")
  library(combinat)
}

# two step estimation for fixed p
two_step_regression <- function(X,y,n,p,s,q){
  # estimate Psi space
  Sigma_hat = t(X)%*%X/n
  svd_Sigma = svd(Sigma_hat)
  # sigma_square_hat = mean(svd_Sigma$d[(q+1):p])
  # Psi_hat = t(svd_Sigma$u[1:p,1:q]%*%diag((svd_Sigma$d[1:q] - sigma_square_hat)**0.5))
  
  # estimate SIV
  B_Psi = svd_Sigma$u[1:p,(q+1):p]
  SIV = X%*%B_Psi
  # least square
  X_hat = SIV%*%solve(t(SIV)%*%SIV,t(SIV)%*%X)
  # the second step regression
  beta_hat = sparse_lm(X_hat,y,p,s)
  return(beta_hat)
}

# l_0 sparsity regression
sparse_lm <- function(X, y, p, s) {
  best_model <- NULL
  best_rss <- Inf
  best_idx <- NULL
  
  # search for all model with sparsity smaller than s
  for (k in 1:s) {
    subsets <- combn(p, k, simplify = FALSE)
    for (idx in subsets) {
      X_sub <- X[, idx, drop = FALSE]
      fit <- lm(y ~ X_sub - 1)
      rss <- sum(residuals(fit)^2)
      if (rss < best_rss) {
        best_rss <- rss
        best_model <- list(fit = fit, vars = idx, rss = rss)
        best_idx <- idx
      }
    }
  }
  
  coef = best_model$fit$coefficients
  beta_hat = rep(0,p)
  beta_hat[best_idx] = coef
  return(beta_hat)
}

# our method
transfer_estimation_fixed_p <- function(X_s,y_s,X_t,y_t,s,q){
  n_s = dim(X_s)[1]
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  # estimate beta_s
  beta_s <- two_step_regression(X_s,y_s,n_s,p,s,q)
  
  # bias correction step
  # lasso_initial = Lasso_func(X_s,y_s)
  beta_initial = beta_s
  Z_s = y_s - X_s%*%beta_initial
  
  cov_z = X_t%*%t(X_t)/n_t**2
  svd_cov_z = svd(cov_z)
  cov_z_pesudo_inverse = svd_cov_z$u[1:n_t,1:p]%*%diag((svd_cov_z$d[1:p])**(-1))%*%t(svd_cov_z$v[1:n_t,1:p])
  Z_t =  cov_z_pesudo_inverse%*%(X_t/n_t)%*%(t(X_s)%*%Z_s/n_s)
  #estimate eta by y_t-Z_t-X_t*beta_initial~X_t
  eta = Lasso_func(X_t, y_t-Z_t-X_t%*%beta_initial)$beta

  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta)
  return(ret_list)
}

# baseline 1 for fixed p
baseline1_estimation_fixed_p <- function(X_s,y_s,X_t,y_t,s,s_t,q){
  n_s = dim(X_s)[1]
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  # estimate beta_s
  beta_s <- two_step_regression(X_s,y_s,n_s,p,s,q)
  
  # estimate beta_t
  beta_t <- two_step_regression(X_t,y_t,n_t,p,s_t,q)
  
  ret_list = list('eta'=beta_t-beta_s,'beta_s'=beta_s,'beta_t'=beta_t)
  return(ret_list)
}

# baseline 2 for fixed p
baseline2_estimation_fixed_p <- function(X_s,y_s,X_t,y_t,s,q){
  n_s = dim(X_s)[1]
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  # estimate beta_s
  beta_s <- two_step_regression(X_s,y_s,n_s,p,s,q)
  
  # estimate eta by y_t-X_t*beta_s~X_t
  # eta = two_step_regression(X_t,y_t-X_t%*%beta_s,n_t,p,s_t,q)
  eta = Lasso_func(X_t, y_t-X_t%*%beta_s)$beta

  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta)
  return(ret_list)
}


# our method
multi_transfer_estimation_fixed_p <- function(X_s,y_s,X_t,y_t,s,q){
  n_s = dim(X_s)[1]/4
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  # estimate beta_s
  beta_s1 <- two_step_regression(X_s[1:n_s,1:p],y_s[1:n_s],n_s,p,s,q)
  beta_s2 <- two_step_regression(X_s[(n_s+1):(2*n_s),1:p],y_s[(n_s+1):(2*n_s)],n_s,p,s,q)
  beta_s3 <- two_step_regression(X_s[(2*n_s+1):(3*n_s),1:p],y_s[(2*n_s+1):(3*n_s)],n_s,p,s,q)
  beta_s4 <- two_step_regression(X_s[(3*n_s+1):(4*n_s),1:p],y_s[(3*n_s+1):(4*n_s)],n_s,p,s,q)
  beta_s <- (beta_s1 + beta_s2 + beta_s3 + beta_s4) / 4
  
  # bias correction step
  # lasso_initial = Lasso_func(X_s,y_s)
  beta_initial = beta_s
  Z_s = y_s - X_s%*%beta_initial
  
  cov_z = X_t%*%t(X_t)/n_t**2
  svd_cov_z = svd(cov_z)
  cov_z_pesudo_inverse = svd_cov_z$u[1:n_t,1:p]%*%diag((svd_cov_z$d[1:p])**(-1))%*%t(svd_cov_z$v[1:n_t,1:p])
  Z_t =  cov_z_pesudo_inverse%*%(X_t/n_t)%*%(t(X_s)%*%Z_s/(n_s*4))
  #estimate eta by y_t-Z_t-X_t*beta_initial~X_t
  eta = Lasso_func(X_t, y_t-Z_t-X_t%*%beta_initial)$beta
  
  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta)
  return(ret_list)
}

# baseline 1 for fixed p
multi_baseline1_estimation_fixed_p <- function(X_s,y_s,X_t,y_t,s,s_t,q){
  n_s = dim(X_s)[1]/4
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  # estimate beta_s
  beta_s1 <- two_step_regression(X_s[1:n_s,1:p],y_s[1:n_s],n_s,p,s,q)
  beta_s2 <- two_step_regression(X_s[(n_s+1):(2*n_s),1:p],y_s[(n_s+1):(2*n_s)],n_s,p,s,q)
  beta_s3 <- two_step_regression(X_s[(2*n_s+1):(3*n_s),1:p],y_s[(2*n_s+1):(3*n_s)],n_s,p,s,q)
  beta_s4 <- two_step_regression(X_s[(3*n_s+1):(4*n_s),1:p],y_s[(3*n_s+1):(4*n_s)],n_s,p,s,q)
  beta_s <- (beta_s1 + beta_s2 + beta_s3 + beta_s4) / 4
  
  # estimate beta_t
  beta_t <- two_step_regression(X_t,y_t,n_t,p,s_t,q)
  
  ret_list = list('eta'=beta_t-beta_s,'beta_s'=beta_s,'beta_t'=beta_t)
  return(ret_list)
}

# baseline 2 for fixed p
multi_baseline2_estimation_fixed_p <- function(X_s,y_s,X_t,y_t,s,q){
  n_s = dim(X_s)[1]/4
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  # estimate beta_s
  beta_s1 <- two_step_regression(X_s[1:n_s,1:p],y_s[1:n_s],n_s,p,s,q)
  beta_s2 <- two_step_regression(X_s[(n_s+1):(2*n_s),1:p],y_s[(n_s+1):(2*n_s)],n_s,p,s,q)
  beta_s3 <- two_step_regression(X_s[(2*n_s+1):(3*n_s),1:p],y_s[(2*n_s+1):(3*n_s)],n_s,p,s,q)
  beta_s4 <- two_step_regression(X_s[(3*n_s+1):(4*n_s),1:p],y_s[(3*n_s+1):(4*n_s)],n_s,p,s,q)
  beta_s <- (beta_s1 + beta_s2 + beta_s3 + beta_s4) / 4
  
  # estimate eta by y_t-X_t*beta_s~X_t
  # eta = two_step_regression(X_t,y_t-X_t%*%beta_s,n_t,p,s_t,q)
  eta = Lasso_func(X_t, y_t-X_t%*%beta_s)$beta
  
  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta)
  return(ret_list)
}


source("utils_nonlinear.R")
library(withr)
# our method
transfer_estimation <- function(X_s,y_s,X_t,y_t,
                                true_beta_s = NULL,
                                rho=0.5,trim_form = 'standard',linear=TRUE){
  n_s = dim(X_s)[1]
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  X = rbind(X_s,X_t)
  y = c(y_s,y_t)
  
  # trim operator
  Q_s = Q_trim(X_s,rho=rho,trim_form = trim_form)
  Q_t = Q_trim(X_t,rho=rho,trim_form = trim_form)
  QX_s = Q_s%*%X_s
  Qy_s = Q_s%*%y_s
  QX_t = Q_t%*%X_t
  Qy_t = Q_t%*%y_t

  # estimate beta_s by Qy_s~QX_s
  beta_s = Lasso_func(QX_s,Qy_s)$beta
  # beta_initial <- Lasso_func(QX_s,Qy_s)$beta
  # lasso_initial = Lasso_func(X_s,y_s)
  # beta_initial = 0.1*beta_initial + 0.9*lasso_initial$beta
  
  # bias correction step
  
  # Add noise
  beta_initial <- withr::with_preserve_seed({
    get_beta_init(
      method = 'noise',
      X_s = X_s, y_s = y_s,
      X_t = X_t, y_t = y_t,
      k = NULL,
      screen_domain = 'source',
      noise_sigma = 0.05
    )
  })

  
  # # Ridge
  # beta_initial <- withr::with_preserve_seed({
  #   # ridge via glmnet: alpha = 0
  #   ridge_fit <- cv.glmnet(
  #     x         = X_s,
  #     y         = y_s,
  #     alpha     = 0,
  #     intercept = FALSE
  #   )
  #   as.numeric(coef(ridge_fit, s = ridge_fit$lambda.min)[-1])
  # })
  
  # screening
  # beta_initial <- withr::with_preserve_seed({
  #   get_beta_init(
  #     method        = "ols_screen_nonlinear",
  #     X_s           = X_s, y_s = y_s,
  #     X_t           = X_t, y_t = y_t,
  #     k             = 160,
  #     screen_domain = "source"
  #   )
  # })
  
  ## zero
  #beta_initial <- rep(0,p)
  
  Z_s <- as.numeric(y_s - X_s %*% beta_initial)
  # Z_s_lasso <- as.numeric(y_s - X_s %*% lasso_initial$beta)
  # delta_beta <- as.numeric(beta_initial - lasso_initial$beta)
  # cat("beta_initial = ", beta_initial[1:10], "\n")
  # cat("non-zero lasso_beta = ", sum(lasso_initial$beta!=0), "\n")
  # cat("non-zero initial_beta = ", sum(beta_initial!=0), "\n")
  # cat("||delta_beta||2 = ", sqrt(sum(delta_beta^2)), "\n")
  # cat("||Z_s (oracle)||_2  =",
  #     sqrt(sum(Z_s^2)), "\n")
  # cat("||Z_s (lasso) ||_2  =",
  #     sqrt(sum(Z_s_lasso^2)), "\n")
  
  Z_t = find_Zt(t(X_t)/n_t, t(X_s)%*%Z_s/n_s)$Z_t
  #estimate eta by y_t-Z_t-X_t*beta_initial~X_t
  eta = Lasso_func(X_t, y_t-Z_t-X_t%*%beta_initial)$beta
  
  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta)
  
}

# baseline 1: no transfer
baseline1_estimation <- function(X_s,y_s,X_t,y_t,rho=0.5,trim_form = 'standard'){
  
  # trim operator
  Q_s = Q_trim(X_s,rho=rho,trim_form = trim_form)
  Q_t = Q_trim(X_t,rho=rho,trim_form = trim_form)
  # QX,Qy
  QX_s = Q_s%*%X_s
  Qy_s = Q_s%*%y_s
  QX_t = Q_t%*%X_t
  Qy_t = Q_t%*%y_t
  
  # estimate beta_s by Qy_s~QX_s
  beta_s = Lasso_func(QX_s,Qy_s)$beta
  
  # estimate beta_t by Qy_t~QX_t
  beta_t = Lasso_func(QX_t,Qy_t)$beta
  
  ret_list = list('eta'=beta_t-beta_s,'beta_s'=beta_s,'beta_t'=beta_t)
}

# baseline 2
# transfer the coefficient beta 
# Transfer learning for high-dimensional linear regression: Prediction, estimation and minimax optimality. Sai Li, T. Tony Cai & Hongzhe Li. JRSSB 2021
baseline2_estimation <- function(X_s,y_s,X_t,y_t,rho=0.5,trim_form = 'standard'){
  # trim operator
  Q_s = Q_trim(X_s,rho=rho,trim_form = trim_form)
  Q_t = Q_trim(X_t,rho=rho,trim_form = trim_form)
  # QX,Qy
  QX_s = Q_s%*%X_s
  Qy_s = Q_s%*%y_s
  QX_t = Q_t%*%X_t
  Qy_t = Q_t%*%y_t
  
  # estimate beta_s by Qy_s~QX_s
  beta_s = Lasso_func(QX_s,Qy_s)$beta
  
  # estimate eta by Qy_t-QX_t*beta_s~QX_t
  eta = Lasso_func(QX_t,Qy_t-QX_t%*%beta_s)$beta
  
  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta)
}

# baseline 3: transfer the structure of X 
baseline3_estimation <- function(X_s,y_s,X_t,y_t,rho=0.5,trim_form = 'standard'){
  n_s = dim(X_s)[1]
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  X = rbind(X_s,X_t)
  y = c(y_s,y_t)
  
  # trim operator
  Q = Q_trim(X,rho=rho,trim_form = trim_form)
  # do QX, QY together
  QX = Q%*%X
  Qy = Q%*%y
  
  # split to source and target 
  QX_s = QX[1:n_s,1:p]
  Qy_s = Qy[1:n_s]
  QX_t = QX[(n_s+1):(n_s+n_t),1:p]
  Qy_t = Qy[(n_s+1):(n_s+n_t)]
  
  # estimate beta_s by Qy_s~QX_s
  beta_s = Lasso_func(QX_s,Qy_s)$beta
  
  # estimate beta_t by Qy_t~QX_t
  beta_t = Lasso_func(QX_t,Qy_t)$beta
  
  ret_list = list('eta'=beta_t-beta_s,'beta_s'=beta_s,'beta_t'=beta_t)
}

# baseline 4: transfer both the structure of X and the coefficient beta
baseline4_estimation <- function(X_s,y_s,X_t,y_t,rho=0.5,trim_form = 'standard'){
  n_s = dim(X_s)[1]
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  X = rbind(X_s,X_t)
  y = c(y_s,y_t)
  
  # trim operator
  Q = Q_trim(X,rho=rho,trim_form = trim_form)
  # do QX, QY together
  QX = Q%*%X
  Qy = Q%*%y
  
  # split to source and target 
  QX_s = QX[1:n_s,1:p]
  Qy_s = Qy[1:n_s]
  QX_t = QX[(n_s+1):(n_s+n_t),1:p]
  Qy_t = Qy[(n_s+1):(n_s+n_t)]
  
  # estimate beta_s by Qy_s~QX_s
  beta_s = Lasso_func(QX_s,Qy_s)$beta
  
  # estimate beta_t by Qy_t-QX_t*beta_s~QX_t
  eta = Lasso_func(QX_t,Qy_t-QX_t%*%beta_s)$beta
  
  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta)
}


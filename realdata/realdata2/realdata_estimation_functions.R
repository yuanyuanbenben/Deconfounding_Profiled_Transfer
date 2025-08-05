
if (!require(glmnet)) {
  install.packages("glmnet")
  library(glmnet)
}

if (!require(lpSolve)) {
  install.packages("lpSolve")
  library(lpSolve)
}

# find trim operator Q
Q_trim <- function(X,rho=0.5,trim_form='standard',quant=2){
  UDV_list = svd(X)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)
  if (trim_form == 'standard'){
    tau = quantile(D, rho)
    Dtilde = pmin(D, tau)
    Q = diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*%  t(U)
    return(Q)
  }
  if (trim_form == 'fixed'){
    tau =D[quant]
    Dtilde = pmin(D, tau)
    Q = diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*%  t(U)
    return(Q)
  }
  if (trim_form == 'zero'){
    gap_list = (D[1:(round(length(D)/2)-1)] - D[2:round(length(D)/2)])/D[2:round(length(D)/2)]
    q_hat = which.max(gap_list)
    Q = diag(nrow(X)) - U[,1:q_hat] %*%t(U[,1:q_hat])
    return(Q)
  }
  if (trim_form =='no_trim'){
    Q = diag(nrow(X))
    return(Q)
  }
}

# lasso for y~X
Lasso_func <- function(X,y){
  fit = cv.glmnet(x=X, y=y,intercept = FALSE)
  beta = as.matrix((coef(fit,S =fit$lambda.min)[-1]))
  result = list("betas" = beta,"res" = y-X %*% beta)
  return(result)
}

# find Zt = argmin\|X Z_t-y\|_{\inf}
find_Zt <- function(X,y){
  p <- dim(X)[1]
  n <- dim(X)[2]
  
  
  # Objective: Minimize t
  # Decision variables: [beta_1, beta_2, ..., beta_n,beta_1_, beta_2_, ..., beta_n_, t]
  f.obj <- c(rep(0, 2*n), 1)  # Coefficients for the objective function
  
  # Constraints:
  # 1. X * beta - t <= y
  # 2. -X * beta - t <= -y
  
  # Construct A_ub matrix
  A_upper <- cbind(X,-X,-1)    # X * beta+ - X * beta- - t
  A_lower <- cbind(-X,X,-1)   # -X * beta+ + X * beta- - t
  A_ub <- rbind(A_upper, A_lower)
  
  # Construct b_ub vector
  b_ub <- c(y, -y)
  
  # Solve the linear program
  lp_result <- lp(
    direction = "min",
    objective.in = f.obj,
    const.mat = A_ub,
    const.dir = rep("<=", 2 * p),
    const.rhs = b_ub,
    compute.sens = FALSE,
    dense.const = TRUE
  )
  
  # Check if the optimization was successful
  if (lp_result$status == 0) {
    # Extract the solution
    solution <- lp_result$solution
    beta_opt <- solution[1:n] - solution[(n+1):(2*n)]
    t_opt <- solution[2*n + 1]
    
    cat("Optimal Zt:\n")
    #print(beta_opt)
    cat("\nMinimum infinity norm:", t_opt, "\n")
  } else {
    cat("Optimization failed. Status code:", lp_result$status, "\n")
    cat("Message:", lp_result$message, "\n")
  }
  return(list('Z_t'=beta_opt,'loss'=t_opt))
}

# our method
transfer_estimation <- function(X_s,y_s,X_t,y_t,rho=0.5,trim_form = 'standard'){
  n_s = dim(X_s)[1]
  n_t = dim(X_t)[1]
  # p = dim(X_s)[2]
  # X = rbind(X_s,X_t)
  # y = c(y_s,y_t)
  
  # trim operator
  
  # Q = Q_trim(X,rho=rho,trim_form = trim_form)
  # # do QX, QY together
  # QX = Q%*%X
  # Qy = Q%*%y
  #
  # # split to source and target
  # QX_s = QX[1:n_s,1:p]
  # Qy_s = Qy[1:n_s]
  # QX_t = QX[(n_s+1):(n_s+n_t),1:p]
  # Qy_t = Qy[(n_s+1):(n_s+n_t)]
  
  # trim operator
  Q_s = Q_trim(X_s,rho=rho,trim_form = trim_form)
  Q_t = Q_trim(X_t,rho=rho,trim_form = trim_form)
  # Q_s_fix = Q_trim(X_s,rho=rho,trim_form = 'fixed',quant=round(n_t/2))
  # Q_s_fix = Q_alter_trim(X_t,rho=rho)
  # QX,Qy
  QX_s = Q_s%*%X_s
  Qy_s = Q_s%*%y_s
  QX_t = Q_t%*%X_t
  Qy_t = Q_t%*%y_t
  # QX_s_fixed = Q_s_fix%*%X_s
  # Qy_s_fixed = Q_s_fix%*%y_s
  
  
  # estimate beta_s by Qy_s~QX_s
  beta_s = Lasso_func(QX_s,Qy_s)$beta
  
  # bias correction step
  lasso_initial = Lasso_func(X_s,y_s)
  # beta_initial = beta_s
  beta_initial = lasso_initial$beta
  Z_s = y_s - X_s%*%beta_initial
  Z_t = find_Zt(t(X_t)/n_t, t(X_s)%*%Z_s/n_s)$Z_t
  #estimate eta by y_t-Z_t-X_t*beta_initial~X_t
  # print(y_t-Z_t-X_t%*%beta_initial)
  lasso_model = Lasso_func(X_t, y_t-Z_t-X_t%*%beta_initial)
  eta = lasso_model$beta
  res1 = mean((y_t - X_t%*%(beta_s+eta))**2)
  res2 = mean(abs(y_t - X_t%*%(beta_s+eta)))
  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta,'res1'=res1,'res2'=res2)
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
  res1 = mean((y_t - X_t%*%(beta_t))**2)
  res2 = mean(abs(y_t - X_t%*%(beta_t)))
  ret_list = list('eta'=beta_t-beta_s,'beta_s'=beta_s,'beta_t'=beta_t,'res1'=res1,'res2'=res2)
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
  res1 = mean((y_t - X_t%*%(beta_s+eta))**2)
  res2 = mean(abs(y_t - X_t%*%(beta_s+eta)))
  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta,'res1'=res1,'res2'=res2)
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
  res1 = mean((y_t - X_t%*%(beta_t))**2)
  res2 = mean(abs(y_t - X_t%*%(beta_t)))
  ret_list = list('eta'=beta_t-beta_s,'beta_s'=beta_s,'beta_t'=beta_t,'res1'=res1,'res2'=res2)
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
  res1 = mean((y_t - X_t%*%(beta_s+eta))**2)
  res2 = mean(abs(y_t - X_t%*%(beta_s+eta)))
  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta,'res1'=res1,'res2'=res2)
}









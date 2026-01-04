source("utils.R")

# our method
transfer_estimation <- function(X_s,y_s,X_t,y_t,rho=0.5,trim_form = 'standard',linear=TRUE){
  n_s = dim(X_s)[1]
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  X = rbind(X_s,X_t)
  y = c(y_s,y_t)
  
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
  beta_initial = lasso_initial$beta
  Z_s = lasso_initial$res
  Z_t = find_Zt(t(X_t)/n_t, t(X_s)%*%Z_s/n_s)$Z_t
  #estimate eta by y_t-Z_t-X_t*beta_initial~X_t
  eta = Lasso_func(X_t, y_t-Z_t-X_t%*%beta_initial)$beta
  
  # v^{(k)} 
  eta1 = Lasso_func(QX_t,Qy_t-QX_t%*%beta_initial-Q_t%*%Z_t)$beta
  v_value = max(abs(t(X_t)%*%(y_t-Z_t-X_t%*%beta_initial-X_t%*%eta1)))/n_t/log(p)
  print(v_value)
  ret_list = list('eta'=eta,'beta_s'=beta_s,'beta_t'=beta_s+eta,'v_value'=v_value)
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
  # print(sum(y_t-X_t%*%(eta+beta_s))**2)
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

# baseline 5: oracle baseline, transfer both the structure of X and the coefficient beta with known confounding structure 
baseline5_estimation <- function(X_s,y_s,X_t,y_t,Psi,rho=0.5,trim_form = 'standard'){
  n_s = dim(X_s)[1]
  n_t = dim(X_t)[1]
  p = dim(X_s)[2]
  X = rbind(X_s,X_t)
  y = c(y_s,y_t)
  
  # projection operator
  H_ = find_HE(X,Psi,n=n_s+n_t,q=3)$H
  # Q = diag(1,n_s+n_t,n_s+n_t) - H_%*%solve(t(H_)%*%H_,t(H_))
  X_H = cbind(X,H_)
  X_H_s = X_H[1:n_s,]
  X_H_t = X_H[(n_s+1):(n_s+n_t),]
  # trim operator
  Q_s = Q_trim(X_H_s,rho=rho,trim_form = trim_form)
  Q_t = Q_trim(X_H_t,rho=rho,trim_form = trim_form)
  # QX,Qy
  QX_s = Q_s%*%X_s
  Qy_s = Q_s%*%y_s
  QX_t = Q_t%*%X_t
  Qy_t = Q_t%*%y_t
  
  # estimate beta_s by Qy_s~QX_s
  beta_s = Lasso_func(QX_s,Qy_s)$beta
  # estimate beta_t by Qy_t-QX_t*beta_s~QX_t
  eta = Lasso_func(QX_t,Qy_t-QX_t%*%beta_s)$beta
  
  ret_list = list('eta'=eta[1:p,1],'beta_s'=beta_s[1:p,1],'beta_t'=beta_s[1:p,1]+eta[1:p,1])
}




















setwd("/home/yuanyuanbenben/project_transfer_learning/transfer_learning_confounder/multi_source")

source("utils.R")
source("estimate_function.R")
library(foreach)
library(doParallel)
# setting
source_num = 4
n_s = 150
n_t = 100
p = 2000
q = 3
s = 3
s_delta = 1
sigma = 1
# control the confounder structure
alpha_Psi = 1
# control the effect of confounder to model
alpha_phi = 1
# trim operator
trim_form = 'standard' # 'zero', 'no_trim'

# parameters
beta = matrix(rep(c(2,0),times = c(s,p-s)),p,1,byrow = TRUE)
eta = matrix(rep(c(-4,0,2),times = c(s_delta,p-2*s_delta,s_delta)),p,1,byrow = TRUE)
beta_s = matrix(rep(beta,source_num),p,source_num)
# for (i in 1:source_num) {
#   beta_s[i,i] = 3
# }
beta_t = rowMeans(beta_s) + eta
alpha_index = c(1,0.75,0.5,0.25)


cl.cores = detectCores(logical = F)
cl <- makeCluster(102)
registerDoParallel(cl)

rep=100
results_out <- foreach(repeat_index = 1:rep,.verbose=TRUE,.combine = rbind,.packages=c("glmnet","lpSolve"))%dopar%{
  set.seed(repeat_index)
  
  # confounder structure 
  # we can relax the spike condition for confounder structure
  Psi = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE) / (n_s*source_num)**0.5 * log(n_s*source_num) /2 * alpha_Psi
  phi_ = matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE) * alpha_phi
  
  ## source data
  X_s_ = matrix(0,source_num*n_s,p)
  y_s_ = rep(0,source_num*n_s)
  for (i in 1:source_num){
    alpha=alpha_index[i]
    phi = phi_*alpha + matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE) * alpha_phi * (1-alpha)
    source = generate_dataset(Psi,phi,beta_s[1:p,i],n_s,p,s,q,sigmaE=2,sigma=sigma,pert=1)
    X_s_[((i-1)*n_s+1):(i*n_s),1:p] = source$X
    y_s_[((i-1)*n_s+1):(i*n_s)] = as.vector(source$Y)
  }
 
  ## target data
  target = generate_dataset(Psi,phi_,beta_t,n_t,p,s,q,sigmaE=2,sigma=sigma,pert=1)
  X_t = target$X
  y_t = as.vector(target$Y)
  
  rho_index = c(1,1.1,1.2,1.3,1.4,1.6,1.8,2,2.5,3)
  
  # our method (no selection)
  result = transfer_estimation(X_s_, y_s_, X_t, y_t,trim_form = trim_form)
  eta_hat <- as.vector(result$eta)
  beta_t_hat <- as.vector(result$beta_t)
  beta_s_hat <- as.vector(result$beta_s)
  cat("\nLoss of our method (no selection)\n")
  # cat("eta_l1_loss: ",sum(abs(eta_hat_ - eta)),"\n")
  # cat("eta_l2_loss: ",(sum((eta_hat_ - eta)**2))**0.5,"\n")
  cat("beta_t_l1_loss: ",sum(abs(beta_t_hat - beta_t)),"\n")
  cat("beta_t_l2_loss: ",(sum((beta_t_hat - beta_t)**2))**0.5,"\n")
  # cat("beta_s_l1_loss: ",sum(abs(beta_s_hat_ - rowMeans(beta_s))),"\n")
  # cat("beta_s_l2_loss: ",(sum((beta_s_hat_ - rowMeans(beta_s))**2))**0.5,"\n")
  
  ret = rep(0,2*length(rho_index)+2)
  ret[length(rho_index)*2+1] = sum(abs(beta_t_hat - beta_t))
  ret[length(rho_index)*2+2] = (sum((beta_t_hat - beta_t)**2))**0.5
  
  v_value_list = rep(0,source_num)
  for (i in 1:source_num){
    result = transfer_estimation_selection(X_s_[((i-1)*n_s+1):(i*n_s),1:p], y_s_[((i-1)*n_s+1):(i*n_s)], X_t, y_t,trim_form = trim_form)$v_value
    v_value_list[i] = result
  }
  print(v_value_list)
  v_min = min(v_value_list)
  loss_l1 = rep(0,length(rho_index))
  loss_l2 = rep(0,length(rho_index))
  residual1 = rep(0,length(rho_index))
  residual2 = rep(0,length(rho_index))
  for (rho_i in 1:length(rho_index)){
    rho = rho_index[rho_i]
    v_index = v_value_list <= rho * v_min
    # selected source
    select_num = sum(v_index)
    X_s = matrix(0,select_num*n_s,p)
    y_s = matrix(0,select_num*n_s)
    j = 0
    for (i in 1:source_num){
      if (v_index[i]){
        j = j + 1
        X_s[((j-1)*n_s+1):(j*n_s),1:p] = X_s_[((i-1)*n_s+1):(i*n_s),1:p]
        y_s[((j-1)*n_s+1):(j*n_s)] = y_s_[((i-1)*n_s+1):(i*n_s)]
      }
    }
    
    # our method 
    result = transfer_estimation(X_s, y_s, X_t, y_t,trim_form = trim_form)
    eta_hat <- as.vector(result$eta)
    beta_t_hat <- as.vector(result$beta_t)
    beta_s_hat <- as.vector(result$beta_s)
    res1 = result$res1
    res2 = result$res2
    # cat("\nLoss of our method","rho=",rho,"\n")
    # cat("eta_l1_loss: ",sum(abs(eta_hat - eta)),"\n")
    # cat("eta_l2_loss: ",(sum((eta_hat - eta)**2))**0.5,"\n")
    # cat("beta_t_l1_loss: ",sum(abs(beta_t_hat - beta_t)),"\n")
    # cat("beta_t_l2_loss: ",(sum((beta_t_hat - beta_t)**2))**0.5,"\n")
    # cat("beta_s_l1_loss: ",sum(abs(beta_s_hat - rowMeans(beta_s))),"\n")
    # cat("beta_s_l2_loss: ",(sum((beta_s_hat - rowMeans(beta_s))**2))**0.5,"\n")
    # cat("residual",res,"\n")
    loss_l1[rho_i] = sum(abs(beta_t_hat - beta_t))
    loss_l2[rho_i] = (sum((beta_t_hat - beta_t)**2))**0.5
    residual1[rho_i] = res1
    residual2[rho_i] = res2
  }
  # selected_rho_index1 = which.min(residual1)
  # selected_rho_index2 = which.min(residual2)
  ret[1:length(rho_index)] = loss_l1
  ret[(1+length(rho_index)):(length(rho_index)*2)] = loss_l2
  # ret[(length(rho_index)*2 +1):(length(rho_index)*3)] = residual1
  # ret[(length(rho_index)*3 +1):(length(rho_index)*4)] = residual2
  # ret[length(rho_index)*4+1] = rho_index[selected_rho_index1]
  # ret[length(rho_index)*4+2] = loss_l1[selected_rho_index1]
  # ret[length(rho_index)*4+3] = loss_l2[selected_rho_index1]
  # ret[length(rho_index)*4+4] = rho_index[selected_rho_index2]
  # ret[length(rho_index)*4+5] = loss_l1[selected_rho_index2]
  # ret[length(rho_index)*4+6] = loss_l2[selected_rho_index2]
  ret
}
write.csv(results_out,file=paste('output/selection_loss_',p,'_',alpha_phi,'_linear.csv',sep = ''))
stopCluster(cl) 


p = 2000
alpha_phi=1
rets <- read.csv(file=paste('output/selection_loss_',p,'_',alpha_phi,'_linear.csv',sep = ''))[2:23]
rets_means <- colMeans(rets)
print(rets_means[c(11,13,15,17,18,19,20,22)])
print(rets_means[c(1,3,5,7,8,9,10,21)])

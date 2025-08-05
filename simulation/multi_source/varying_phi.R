
source("utils.R")
source("estimate_function.R")

library(ggplot2)
library(patchwork)

phi_list <- seq(1,4,by=1)
df_results <- data.frame()

repeat_times <- 100

for (phi_cur in phi_list) {
  # We set p = p_cur
  alpha_phi <- phi_cur
  
  l1_loss_our_eta = 0
  l2_loss_our_eta = 0
  l1_loss_baseline1_eta = 0
  l2_loss_baseline1_eta = 0
  l1_loss_baseline2_eta = 0
  l2_loss_baseline2_eta = 0
  l1_loss_baseline3_eta = 0
  l2_loss_baseline3_eta = 0
  l1_loss_baseline4_eta = 0
  l2_loss_baseline4_eta = 0
  
  l1_loss_our_beta = 0
  l2_loss_our_beta = 0
  l1_loss_baseline1_beta = 0
  l2_loss_baseline1_beta = 0
  l1_loss_baseline2_beta = 0
  l2_loss_baseline2_beta = 0
  l1_loss_baseline3_beta = 0
  l2_loss_baseline3_beta = 0
  l1_loss_baseline4_beta = 0
  l2_loss_baseline4_beta = 0
  
  
  # Repeat 100 times
  for (repeat_index in 1:repeat_times) {
    set.seed(repeat_index)
    
    # setting
    source_num = 4
    n_s = 150
    n_t = 100
    p = 1500
    q = 3
    s = 5
    s_delta = 1
    sigma = 1
    # control the confounder structure
    alpha_Psi = 1
    # control the effect of confounder to model
    alpha_phi = phi_cur
    # trim operator
    trim_form = 'standard' # 'zero', 'no_trim'
    
    # parameters
    beta = matrix(rep(c(2,0),times = c(s,p-s)),p,1,byrow = TRUE)
    eta = matrix(rep(c(-4,0,2),times = c(s_delta,p-2*s_delta,s_delta)),p,1,byrow = TRUE)
    beta_s = matrix(rep(beta,source_num),p,source_num)
    for (i in 1:source_num) {
      beta_s[i,i] = 1.5
    }
    beta_t = rowMeans(beta_s) + eta
    
    # confounder structure 
    # we can relax the spike condition for confounder structure
    Psi = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE) / n_s**0.5 * log(n_s) /2 * alpha_Psi
    phi = matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE) * alpha_phi
    
    ## source data
    X_s = matrix(0,source_num*n_s,p)
    y_s = rep(0,source_num*n_s)
    for (i in 1:source_num){
      source = generate_dataset(Psi,phi,beta_s[1:p,i],n_s,p,s,q,sigmaE=2,sigma=sigma,pert=1)
      X_s[((i-1)*n_s+1):(i*n_s),1:p] = source$X
      y_s[((i-1)*n_s+1):(i*n_s)] = as.vector(source$Y)
    }
    
    ## target data
    target = generate_dataset(Psi,phi,beta_t,n_t,p,s,q,sigmaE=2,sigma=sigma,pert=1)
    X_t = target$X
    y_t = as.vector(target$Y)
    
    # our method
    result = transfer_estimation(X_s, y_s, X_t, y_t,trim_form = trim_form)
    eta_hat <- as.vector(result$eta)
    beta_t_hat <- as.vector(result$beta_t)
    beta_s_hat <- as.vector(result$beta_s)
    
    # baseline 1: no transfer
    baseline1_result = baseline1_estimation(X_s, y_s, X_t, y_t,trim_form = trim_form)
    eta_hat_baseline1 <- as.vector(baseline1_result$eta)
    beta_t_hat_baseline1 <- as.vector(baseline1_result$beta_t)
    beta_s_hat_baseline1 <- as.vector(baseline1_result$beta_s)
    
    # baseline 2
    # transfer the coefficient beta 
    # Transfer learning for high-dimensional linear regression: Prediction, estimation and minimax optimality. Sai Li, T. Tony Cai & Hongzhe Li. JRSSB 2021
    baseline2_result = baseline2_estimation(X_s, y_s, X_t, y_t,trim_form = trim_form)
    eta_hat_baseline2 <- as.vector(baseline2_result$eta)
    beta_t_hat_baseline2 <- as.vector(baseline2_result$beta_t)
    beta_s_hat_baseline2 <- as.vector(baseline2_result$beta_s)
    
    # baseline 3: transfer the structure of X 
    baseline3_result = baseline3_estimation(X_s, y_s, X_t, y_t,trim_form = trim_form)
    eta_hat_baseline3 <- as.vector(baseline3_result$eta)
    beta_t_hat_baseline3 <- as.vector(baseline3_result$beta_t)
    beta_s_hat_baseline3 <- as.vector(baseline3_result$beta_s)
    
    # baseline 4: transfer both the structure of X and the coefficient beta
    baseline4_result = baseline4_estimation(X_s, y_s, X_t, y_t,trim_form = trim_form)
    eta_hat_baseline4 <- as.vector(baseline4_result$eta)
    beta_t_hat_baseline4 <- as.vector(baseline4_result$beta_t)
    beta_s_hat_baseline4 <- as.vector(baseline4_result$beta_s)
    
    l1_loss_our_eta = l1_loss_our_eta+sum(abs(eta_hat - eta))
    l2_loss_our_eta = l2_loss_our_eta+(sum((eta_hat - eta)**2))
    l1_loss_baseline1_eta = l1_loss_baseline1_eta+sum(abs(eta_hat_baseline1 - eta))
    l2_loss_baseline1_eta = l2_loss_baseline1_eta+(sum((eta_hat_baseline1 - eta)**2))
    l1_loss_baseline2_eta = l1_loss_baseline2_eta+sum(abs(eta_hat_baseline2 - eta))
    l2_loss_baseline2_eta = l2_loss_baseline2_eta+(sum((eta_hat_baseline2 - eta)**2))
    l1_loss_baseline3_eta = l1_loss_baseline3_eta+sum(abs(eta_hat_baseline3 - eta))
    l2_loss_baseline3_eta = l2_loss_baseline3_eta+(sum((eta_hat_baseline3 - eta)**2))
    l1_loss_baseline4_eta = l1_loss_baseline4_eta+sum(abs(eta_hat_baseline4 - eta))
    l2_loss_baseline4_eta = l2_loss_baseline4_eta+(sum((eta_hat_baseline4 - eta)**2))
    
    l1_loss_our_beta = l1_loss_our_beta+sum(abs(beta_t_hat - beta_t))
    l2_loss_our_beta = l2_loss_our_beta+(sum((beta_t_hat - beta_t)**2))
    l1_loss_baseline1_beta = l1_loss_baseline1_beta+sum(abs(beta_t_hat_baseline1 - beta_t))
    l2_loss_baseline1_beta = l2_loss_baseline1_beta+(sum((beta_t_hat_baseline1 - beta_t)**2))
    l1_loss_baseline2_beta = l1_loss_baseline2_beta+sum(abs(beta_t_hat_baseline2 - beta_t))
    l2_loss_baseline2_beta = l2_loss_baseline2_beta+(sum((beta_t_hat_baseline2 - beta_t)**2))
    l1_loss_baseline3_beta = l1_loss_baseline3_beta+sum(abs(beta_t_hat_baseline3 - beta_t))
    l2_loss_baseline3_beta = l2_loss_baseline3_beta+(sum((beta_t_hat_baseline3 - beta_t)**2))
    l1_loss_baseline4_beta = l1_loss_baseline4_beta+sum(abs(beta_t_hat_baseline4 - beta_t))
    l2_loss_baseline4_beta = l2_loss_baseline4_beta+(sum((beta_t_hat_baseline4 - beta_t)**2))
  }
  
  # After 100 repeats, compute averages
  l1_loss_our_eta_mean  <- l1_loss_our_eta / repeat_times
  l2_loss_our_eta_mean  <- l2_loss_our_eta / repeat_times
  l1_loss_baseline1_eta_mean <- l1_loss_baseline1_eta / repeat_times
  l2_loss_baseline1_eta_mean <- l2_loss_baseline1_eta / repeat_times
  l1_loss_baseline2_eta_mean <- l1_loss_baseline2_eta / repeat_times
  l2_loss_baseline2_eta_mean <- l2_loss_baseline2_eta / repeat_times
  l1_loss_baseline3_eta_mean <- l1_loss_baseline3_eta / repeat_times
  l2_loss_baseline3_eta_mean <- l2_loss_baseline3_eta / repeat_times
  l1_loss_baseline4_eta_mean <- l1_loss_baseline4_eta / repeat_times
  l2_loss_baseline4_eta_mean <- l2_loss_baseline4_eta / repeat_times
  
  
  
  l1_loss_our_beta_mean  <- l1_loss_our_beta / repeat_times
  l2_loss_our_beta_mean  <- l2_loss_our_beta / repeat_times
  l1_loss_baseline1_beta_mean <- l1_loss_baseline1_beta / repeat_times
  l2_loss_baseline1_beta_mean <- l2_loss_baseline1_beta / repeat_times
  l1_loss_baseline2_beta_mean <- l1_loss_baseline2_beta / repeat_times
  l2_loss_baseline2_beta_mean <- l2_loss_baseline2_beta / repeat_times
  l1_loss_baseline3_beta_mean <- l1_loss_baseline3_beta / repeat_times
  l2_loss_baseline3_beta_mean <- l2_loss_baseline3_beta / repeat_times
  l1_loss_baseline4_beta_mean <- l1_loss_baseline4_beta / repeat_times
  l2_loss_baseline4_beta_mean <- l2_loss_baseline4_beta / repeat_times
  
  # Add one row per method to df_results
  
  # our method
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "our_method",
    metric  = "eta_L1",
    value   = l1_loss_our_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "our_method",
    metric  = "eta_L2",
    value   = l2_loss_our_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "our_method",
    metric  = "beta_L1",
    value   = (l1_loss_our_beta_mean)
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "our_method",
    metric  = "beta_L2",
    value   = (l2_loss_our_beta_mean)
  ))
  
  
  ## baseline1
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_1",
    metric  = "eta_L1",
    value   = l1_loss_baseline1_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_1",
    metric  = "eta_L2",
    value   = l2_loss_baseline1_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_1",
    metric  = "beta_L1",
    value   = (l1_loss_baseline1_beta_mean)
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_1",
    metric  = "beta_L2",
    value   = (l2_loss_baseline1_beta_mean)
  ))
  
  
  ## baseline 2
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_2",
    metric  = "eta_L1",
    value   = l1_loss_baseline2_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_2",
    metric  = "eta_L2",
    value   = l2_loss_baseline2_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi      = phi_cur,
    method  = "baseline_2",
    metric  = "beta_L1",
    value   = (l1_loss_baseline2_beta_mean)
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_2",
    metric  = "beta_L2",
    value   = (l2_loss_baseline2_beta_mean)
  ))
  
  ## baseline 3
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_3",
    metric  = "eta_L1",
    value   = l1_loss_baseline3_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_3",
    metric  = "eta_L2",
    value   = l2_loss_baseline3_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_3",
    metric  = "beta_L1",
    value   = (l1_loss_baseline3_beta_mean)
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_3",
    metric  = "beta_L2",
    value   = (l2_loss_baseline3_beta_mean)
  ))
  
  
  ## baseline 4
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_4",
    metric  = "eta_L1",
    value   = l1_loss_baseline4_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_4",
    metric  = "eta_L2",
    value   = l2_loss_baseline4_eta_mean
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_4",
    metric  = "beta_L1",
    value   = (l1_loss_baseline4_beta_mean)
  ))
  df_results <- rbind(df_results, data.frame(
    phi       = phi_cur,
    method  = "baseline_4",
    metric  = "beta_L2",
    value   = (l2_loss_baseline4_beta_mean)
  ))
}

save(df_results, file = "varying_phi_multi_source_linear.Rdata")




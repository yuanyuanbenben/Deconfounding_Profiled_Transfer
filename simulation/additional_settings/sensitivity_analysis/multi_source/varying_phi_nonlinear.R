rm(list=ls())

source("utils_nonlinear.R")
source("estimate_function_nonlinear_ini.R")

library(ggplot2)
library(patchwork)

phi_list <- seq(1,4,by=1)
df_results <- data.frame()

repeat_times <- 100

for (phi_cur in phi_list) {
  # We set p = phi_cur
  metric_names <- c(
    "our_eta_L1","our_eta_L2","our_beta_L1","our_beta_L2",
    "b1_eta_L1","b1_eta_L2","b1_beta_L1","b1_beta_L2",
    "b2_eta_L1","b2_eta_L2","b2_beta_L1","b2_beta_L2",
    "b3_eta_L1","b3_eta_L2","b3_beta_L1","b3_beta_L2",
    "b4_eta_L1","b4_eta_L2","b4_beta_L1","b4_beta_L2"
  )
  perrep <- setNames(rep(list(numeric(repeat_times)), length(metric_names)), metric_names)
  
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
      # source <- generate_dataset_1(Psi, phi, beta_s[1:p,i], n_s, p, s, q, sigmaE = 2,
      #                                  sigma = sigma, pert = 1,
      #                                  nonlin_type="tanh",
      #                                  nonlin_strength=4)
      X_s[((i-1)*n_s+1):(i*n_s),1:p] = source$X
      y_s[((i-1)*n_s+1):(i*n_s)] = as.vector(source$Y)
    }
    
    ## target data
    target = generate_dataset(Psi,phi,beta_t,n_t,p,s,q,sigmaE=2,sigma=sigma,pert=1)
    # target <- generate_dataset_1(Psi, phi, beta_t, n_t, p, s, q, sigmaE = 2,
    #                                  sigma = sigma, pert = 1,
    #                                  nonlin_type="tanh",
    #                                  nonlin_strength=4)
    X_t = target$X
    y_t = as.vector(target$Y)
    
    # our method
    # result = transfer_estimation(X_s, y_s, X_t, y_t,trim_form = trim_form)
    # eta_hat <- as.vector(result$eta)
    # beta_t_hat <- as.vector(result$beta_t)
    # beta_s_hat <- as.vector(result$beta_s)
    
    ## keep seed
    
    result = with_preserve_seed(
      transfer_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    )
    eta_hat    <- as.vector(result$eta)
    beta_t_hat <- as.vector(result$beta_t)
    beta_s_hat <- as.vector(result$beta_s)
    burn_rng_like_original_ours(X_s, y_s, X_t, y_t, trim_form = trim_form)
    
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
    
    perrep[["our_eta_L1"]][repeat_index]  <- sum(abs(eta_hat - as.vector(eta)))
    perrep[["our_eta_L2"]][repeat_index]  <- sum((eta_hat - as.vector(eta))^2)
    perrep[["our_beta_L1"]][repeat_index] <- sum(abs(beta_t_hat - as.vector(beta_t)))
    perrep[["our_beta_L2"]][repeat_index] <- sum((beta_t_hat - as.vector(beta_t))^2)
    
    ## baseline 1
    perrep[["b1_eta_L1"]][repeat_index]   <- sum(abs(eta_hat_baseline1 - as.vector(eta)))
    perrep[["b1_eta_L2"]][repeat_index]   <- sum((eta_hat_baseline1 - as.vector(eta))^2)
    perrep[["b1_beta_L1"]][repeat_index]  <- sum(abs(beta_t_hat_baseline1 - as.vector(beta_t)))
    perrep[["b1_beta_L2"]][repeat_index]  <- sum((beta_t_hat_baseline1 - as.vector(beta_t))^2)
    
    ## baseline 2
    perrep[["b2_eta_L1"]][repeat_index]   <- sum(abs(eta_hat_baseline2 - as.vector(eta)))
    perrep[["b2_eta_L2"]][repeat_index]   <- sum((eta_hat_baseline2 - as.vector(eta))^2)
    perrep[["b2_beta_L1"]][repeat_index]  <- sum(abs(beta_t_hat_baseline2 - as.vector(beta_t)))
    perrep[["b2_beta_L2"]][repeat_index]  <- sum((beta_t_hat_baseline2 - as.vector(beta_t))^2)
    
    ## baseline 3
    perrep[["b3_eta_L1"]][repeat_index]   <- sum(abs(eta_hat_baseline3 - as.vector(eta)))
    perrep[["b3_eta_L2"]][repeat_index]   <- sum((eta_hat_baseline3 - as.vector(eta))^2)
    perrep[["b3_beta_L1"]][repeat_index]  <- sum(abs(beta_t_hat_baseline3 - as.vector(beta_t)))
    perrep[["b3_beta_L2"]][repeat_index]  <- sum((beta_t_hat_baseline3 - as.vector(beta_t))^2)
    
    ## baseline 4
    perrep[["b4_eta_L1"]][repeat_index]   <- sum(abs(eta_hat_baseline4 - as.vector(eta)))
    perrep[["b4_eta_L2"]][repeat_index]   <- sum((eta_hat_baseline4 - as.vector(eta))^2)
    perrep[["b4_beta_L1"]][repeat_index]  <- sum(abs(beta_t_hat_baseline4 - as.vector(beta_t)))
    perrep[["b4_beta_L2"]][repeat_index]  <- sum((beta_t_hat_baseline4 - as.vector(beta_t))^2)
  }
  
  # After 100 repeats, compute averages
  to_row <- function(method, metric_label, vec, repeat_times) {
    m  <- mean(vec)
    s  <- stats::sd(vec)
    se <- s / sqrt(repeat_times)
    data.frame(
      phi    = phi_cur,
      method = method,
      metric = metric_label,              # "eta_L1","eta_L2","beta_L1","beta_L2"
      mean   = m,
      sd     = se,
      ci_lo  = m - 1.96 * se,
      ci_hi  = m + 1.96 * se
    )
  }
  
  map <- rbind(
    c("our_method","eta_L1","our_eta_L1"),
    c("our_method","eta_L2","our_eta_L2"),
    c("our_method","beta_L1","our_beta_L1"),
    c("our_method","beta_L2","our_beta_L2"),
    c("baseline_1","eta_L1","b1_eta_L1"),
    c("baseline_1","eta_L2","b1_eta_L2"),
    c("baseline_1","beta_L1","b1_beta_L1"),
    c("baseline_1","beta_L2","b1_beta_L2"),
    c("baseline_2","eta_L1","b2_eta_L1"),
    c("baseline_2","eta_L2","b2_eta_L2"),
    c("baseline_2","beta_L1","b2_beta_L1"),
    c("baseline_2","beta_L2","b2_beta_L2"),
    c("baseline_3","eta_L1","b3_eta_L1"),
    c("baseline_3","eta_L2","b3_eta_L2"),
    c("baseline_3","beta_L1","b3_beta_L1"),
    c("baseline_3","beta_L2","b3_beta_L2"),
    c("baseline_4","eta_L1","b4_eta_L1"),
    c("baseline_4","eta_L2","b4_eta_L2"),
    c("baseline_4","beta_L1","b4_beta_L1"),
    c("baseline_4","beta_L2","b4_beta_L2")
  )
  
  for (k in seq_len(nrow(map))) {
    method_k <- map[k, 1]
    metric_k <- map[k, 2]
    key_k    <- map[k, 3]
    df_results <- rbind(
      df_results,
      to_row(method_k, metric_k, perrep[[key_k]], repeat_times)
    )
  }
}

compare_methods_summary <- function(df_results) {
  stopifnot(all(c("phi","method","metric","mean","sd","ci_lo","ci_hi") %in% names(df_results)))
  library(dplyr)
  library(tidyr)
  
  df_results |>
    dplyr::mutate(
      method = factor(method,
                      levels = c("our_method","baseline_1","baseline_2","baseline_3","baseline_4"),
                      labels = c("Our","B1: Target-only","B2: β-transfer",
                                 "B3: X-structure","B4: Both"))
    ) |>
    dplyr::group_by(phi, metric, method) |>
    dplyr::summarise(
      mean = mean(mean),
      sd   = mean(sd),
      ci_lo = mean(ci_lo),
      ci_hi = mean(ci_hi),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols = c(phi, metric),
      names_from = method,
      values_from = mean,
      names_glue = "{method}_mean"
    ) |>
    dplyr::arrange(metric, p)
}

#df_summary <- compare_methods_summary(df_results)
#View(df_summary)

save(df_results, file = "varying_phi_multi_source_nonlinear_with_sd.Rdata")




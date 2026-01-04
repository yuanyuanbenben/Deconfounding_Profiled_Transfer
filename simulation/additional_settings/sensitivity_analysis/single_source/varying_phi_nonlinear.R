rm(list = ls())

source("utils_nonlinear.R")
source("estimate_function_nonlinear_ini.R")

library(ggplot2)
library(patchwork)

phi_list <- seq(1, 4, by = 1)
repeat_times <- 100

## Tidy results container (will hold mean/sd/ci for phi × method × metric)
df_results <- data.frame()

for (phi_cur in phi_list) {
  
  ## per-rep storage for all method × metric (20 vectors)
  metric_names <- c(
    "our_eta_L1","our_eta_L2","our_beta_L1","our_beta_L2",
    "b1_eta_L1","b1_eta_L2","b1_beta_L1","b1_beta_L2",
    "b2_eta_L1","b2_eta_L2","b2_beta_L1","b2_beta_L2",
    "b3_eta_L1","b3_eta_L2","b3_beta_L1","b3_beta_L2",
    "b4_eta_L1","b4_eta_L2","b4_beta_L1","b4_beta_L2"
  )
  perrep <- setNames(rep(list(numeric(repeat_times)), length(metric_names)), metric_names)
  
  for (repeat_index in 1:repeat_times) {
    set.seed(repeat_index)
    
    ## ----- settings -----
    n_s   <- 600
    n_t   <- 100
    p     <- 1500
    q     <- 3
    s_beta <- 5   # sparsity of beta_s (renamed from s)
    s_eta  <- 1   # sparsity of eta     (renamed from s_delta)
    sigma <- 1
    
    # confounder structure & effect
    alpha_Psi <- 1
    alpha_phi <- phi_cur
    trim_form <- "standard"   # 'standard' | 'zero' | 'no_trim'
    
    ## ----- parameters (sparse beta_s and eta) -----
    beta_s <- matrix(rep(c(1, 0), times = c(s_beta, p - s_beta)), p, 1, byrow = TRUE)
    eta    <- matrix(rep(c(-2, 0, 1), times = c(s_eta, p - 2 * s_eta, s_eta)), p, 1, byrow = TRUE)
    beta_t <- beta_s + eta
    
    ## ----- confounder structure (Psi, phi) -----
    Psi <- matrix(rnorm(q * p, 0, 1), q, p, byrow = TRUE) / sqrt(n_s) * (log(n_s) / 2) * alpha_Psi
    phi <- matrix(rnorm(q, 0, 1), q, 1, byrow = TRUE) * alpha_phi
    
    ## ----- source data -----
    # source_dat <- generate_dataset(Psi, phi, beta_s, n_s, p, s_beta, q, sigmaE = 2, sigma = sigma, pert = 1)
    source_dat <- generate_dataset_1(Psi, phi, beta_s, n_s, p, s, q, sigmaE = 2,
                                     sigma = sigma, pert = 1,
                                     nonlin_type="tanh",
                                     nonlin_strength=4)
    X_s <- source_dat$X
    y_s <- as.vector(source_dat$Y)
    
    ## ----- target data -----
    # target_dat <- generate_dataset(Psi, phi, beta_t, n_t, p, s_beta, q, sigmaE = 2, sigma = sigma, pert = 1)
    target_dat <- generate_dataset_1(Psi, phi, beta_t, n_t, p, s, q, sigmaE = 2,
                                     sigma = sigma, pert = 1,
                                     nonlin_type="tanh",
                                     nonlin_strength=4)
    X_t <- target_dat$X
    y_t <- as.vector(target_dat$Y)
    
    ## ----- Fit all methods -----
    # our method
    #result <- transfer_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    #eta_hat    <- as.vector(result$eta)
    #beta_t_hat <- as.vector(result$beta_t)
    
    ## keep seed
    
    result = with_preserve_seed(
      transfer_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    )
    eta_hat    <- as.vector(result$eta)
    beta_t_hat <- as.vector(result$beta_t)
    beta_s_hat <- as.vector(result$beta_s)
    burn_rng_like_original_ours(X_s, y_s, X_t, y_t, trim_form = trim_form)
    
    
    # baseline 1
    baseline1_result <- baseline1_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    eta_hat_baseline1    <- as.vector(baseline1_result$eta)
    beta_t_hat_baseline1 <- as.vector(baseline1_result$beta_t)
    
    # baseline 2
    baseline2_result <- baseline2_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    eta_hat_baseline2    <- as.vector(baseline2_result$eta)
    beta_t_hat_baseline2 <- as.vector(baseline2_result$beta_t)
    
    # baseline 3
    baseline3_result <- baseline3_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    eta_hat_baseline3    <- as.vector(baseline3_result$eta)
    beta_t_hat_baseline3 <- as.vector(baseline3_result$beta_t)
    
    # baseline 4
    baseline4_result <- baseline4_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    eta_hat_baseline4    <- as.vector(baseline4_result$eta)
    beta_t_hat_baseline4 <- as.vector(baseline4_result$beta_t)
    
    ## ----- record per-rep metrics (L1 sums & L2 SSEs) -----
    ## our method
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
  } # end repeats
  
  ## ---- summarize mean, SD, 95% CI for this phi_cur ----
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
} # end phi_list loop

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

#df_summary1 <- compare_methods_summary(df_results)
#View(df_summary1)

save(df_results, file = "varying_phi_single_source_nonlinear_with_sd_ini2_tan.Rdata")

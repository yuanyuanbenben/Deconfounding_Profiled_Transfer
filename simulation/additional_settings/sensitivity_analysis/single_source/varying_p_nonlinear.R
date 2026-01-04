rm(list = ls())

source("utils_nonlinear.R")
source("estimate_function_nonlinear_ini.R")

library(ggplot2)
library(patchwork)

## Grid & repetitions
p_list <- c(500, 1000, 1500, 2000)
repeat_times <- 100

# p_list <- c(1000)
# repeat_times <- 1

## Results container (tidy; will hold mean / sd / CI per p × method × metric)
df_results <- data.frame()

## ===========================
## Main loop over dimensions p
## ===========================
for (p_cur in p_list) {
  
  metric_names <- c(
    "our_eta_L1","our_eta_L2","our_beta_L1","our_beta_L2",
    "b1_eta_L1","b1_eta_L2","b1_beta_L1","b1_beta_L2",
    "b2_eta_L1","b2_eta_L2","b2_beta_L1","b2_beta_L2",
    "b3_eta_L1","b3_eta_L2","b3_beta_L1","b3_beta_L2",
    "b4_eta_L1","b4_eta_L2","b4_beta_L1","b4_beta_L2"
  )
  perrep <- setNames(rep(list(numeric(repeat_times)), length(metric_names)), metric_names)
  
  ## =========================
  ## Repeat the data-generating
  ## =========================
  for (repeat_index in 1:repeat_times) {
    ## If you want identical RNG across p for each repeat_index, keep this seed here.
    ## If instead you want unique randomness across p, place a single set.seed() before the p-loop.
    set.seed(repeat_index)
    
    ## ----- settings -----
    n_s <- 600
    n_t <- 100
    p    <- p_cur
    q    <- 3
    s    <- 5
    #s_beta <- s
    s_delta <- 1
    sigma <- 1
    
    ## Control structure and confounding
    alpha_Psi <- 1   # confounder structure strength
    alpha_phi <- 8   # confounder effect on Y
    #alpha_phi <- 4
    trim_form <- "standard"  # 'standard' | 'zero' | 'no_trim'
    
    ## ----- parameters (sparse beta_s and eta) -----
    beta_s <- matrix(rep(c(2, 0), times = c(s, p - s)), p, 1, byrow = TRUE)
    eta    <- matrix(rep(c(-4, 0, 2), times = c(s_delta, p - 2 * s_delta, s_delta)), p, 1, byrow = TRUE)
    beta_t <- beta_s + eta
    beta_s_true <- beta_s
    
    ## ----- confounder structure (Psi, phi) -----
    Psi <- matrix(rnorm(q * p, 0, 1), q, p, byrow = TRUE) / sqrt(n_s) * (log(n_s) / 2) * alpha_Psi
    phi <- matrix(rnorm(q, 0, 1), q, 1, byrow = TRUE) * alpha_phi
    
    ## ----- source data -----
    # source_dat <- generate_dataset(Psi, phi, beta_s, n_s, p, s, q, sigmaE = 2, sigma = sigma, pert = 1)
    source_dat <- generate_dataset_1(Psi, phi, beta_s, n_s, p, s, q, sigmaE = 2,
                                     sigma = sigma, pert = 1,
                                     nonlin_type="tanh",
                                     nonlin_strength=4)
    X_s <- source_dat$X
    y_s <- as.vector(source_dat$Y)
    
    ## ----- target data -----
    # target_dat <- generate_dataset(Psi, phi, beta_t, n_t, p, s, q, sigmaE = 2, sigma = sigma, pert = 1)
    target_dat <- generate_dataset_1(Psi, phi, beta_t, n_t, p, s, q, sigmaE = 2,
                                     sigma = sigma, pert = 1,
                                     nonlin_type="tanh",
                                     nonlin_strength=4)
    X_t <- target_dat$X
    y_t <- as.vector(target_dat$Y)
    
    ## ----- Fit all methods -----
    ## our method
    # result <- transfer_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    # eta_hat    <- as.vector(result$eta)
    # beta_t_hat <- as.vector(result$beta_t)
    # beta_s_hat <- as.vector(result$beta_s)  # available if needed
    
    ## keep seed
    
    result = with_preserve_seed(
      transfer_estimation(X_s, y_s, X_t, y_t, true_beta_s = beta_s_true,
                          trim_form = trim_form)
    )
    eta_hat    <- as.vector(result$eta)
    beta_t_hat <- as.vector(result$beta_t)
    beta_s_hat <- as.vector(result$beta_s)
    burn_rng_like_original_ours(X_s, y_s, X_t, y_t, trim_form = trim_form, rho = 0.5)
    
    
    ## baseline 1: (e.g., no transfer / target-only per your implementation)
    baseline1_result <- baseline1_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    eta_hat_baseline1    <- as.vector(baseline1_result$eta)
    beta_t_hat_baseline1 <- as.vector(baseline1_result$beta_t)
    
    ## baseline 2: transfer β (JRSSB 2021 approach)
    baseline2_result <- baseline2_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    eta_hat_baseline2    <- as.vector(baseline2_result$eta)
    beta_t_hat_baseline2 <- as.vector(baseline2_result$beta_t)
    
    ## baseline 3: transfer X-structure
    baseline3_result <- baseline3_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    eta_hat_baseline3    <- as.vector(baseline3_result$eta)
    beta_t_hat_baseline3 <- as.vector(baseline3_result$beta_t)
    
    ## baseline 4: transfer both X-structure and β
    baseline4_result <- baseline4_estimation(X_s, y_s, X_t, y_t, trim_form = trim_form)
    eta_hat_baseline4    <- as.vector(baseline4_result$eta)
    beta_t_hat_baseline4 <- as.vector(baseline4_result$beta_t)
    
    ## -----------------------------------------
    ## L1 = sum |.|, L2 = sum of squared errors
    ## -----------------------------------------
    
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
  
  ## =========================================
  ## Summarize mean, SD, 95% CI for this p_cur
  ## =========================================
  to_row <- function(method, metric_label, vec, repeat_times) {
    m  <- mean(vec)
    std  <- stats::sd(vec)
    se <- std / sqrt(repeat_times)
    data.frame(
      p      = p_cur,
      method = method,
      metric = metric_label,              # one of "eta_L1","eta_L2","beta_L1","beta_L2"
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
  stopifnot(all(c("p","method","metric","mean","sd","ci_lo","ci_hi") %in% names(df_results)))
  library(dplyr)
  library(tidyr)
  
  df_results |>
    dplyr::mutate(
      method = factor(method,
                      levels = c("our_method","baseline_1","baseline_2","baseline_3","baseline_4"),
                      labels = c("Our","B1: Target-only","B2: β-transfer",
                                 "B3: X-structure","B4: Both"))
    ) |>
    dplyr::group_by(p, metric, method) |>
    dplyr::summarise(
      mean = mean(mean),
      sd   = mean(sd),
      ci_lo = mean(ci_lo),
      ci_hi = mean(ci_hi),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols = c(p, metric),
      names_from = method,
      values_from = mean,
      names_glue = "{method}_mean"
    ) |>
    dplyr::arrange(metric, p)
}

#df_summary1 <- compare_methods_summary(df_results)
#View(df_summary1)

save(df_results, file = "varying_p_single_source_nonlinear_with_sd_ini2_tan.Rdata")

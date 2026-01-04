folder_path = '/home/yuanyuanbenben/project_transfer_learning/transfer_learning_confounder/realdata2'
setwd(folder_path)

# load functions and data``
source(paste(folder_path,"/realdata_estimation_functions_fixed_p.R",sep=""))
source(paste(folder_path,"/realdata_estimation_functions.R",sep=""))
load("NLSM_data_preprocessed_fixed_p.RData")

set.seed(123456)

# our method 
result = transfer_estimation_fixed_p(source_X, source_Y, target_X, target_Y, 7,4,method='approx')
# result = transfer_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')

# baselines
# result_baseline1 <- baseline1_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')
# result_baseline2 <- baseline2_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')
# result_baseline3 <- baseline3_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')
# result_baseline4 <- baseline4_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')
result_baseline1 <- baseline1_estimation_fixed_p(source_X, source_Y, target_X, target_Y, 7,7,4,method='approx')
result_baseline2 <- baseline2_estimation_fixed_p(source_X, source_Y, target_X, target_Y, 7,4,method='approx')

# deconfounded target
P_target_confounder <- target_confounder %*% solve(t(target_confounder) %*% target_confounder, t(target_confounder))
deconfounded_target_X <- target_X - P_target_confounder %*% target_X
deconfounded_target_Y <- target_Y - P_target_confounder %*% target_Y

# deconfounded source
P_source_confounder <- source_confounder %*% solve(t(source_confounder) %*% source_confounder, t(source_confounder))
deconfounded_source_X <- source_X - P_source_confounder %*% source_X
deconfounded_source_Y <- source_Y - P_source_confounder %*% source_Y

# deconfounded baseline
result_baseline_deconfounded <- baseline2_estimation(deconfounded_source_X, deconfounded_source_Y, deconfounded_target_X, deconfounded_target_Y, trim_form = 'no_trim')


beta_t_l2_error <- (mean((result$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
beta_t_l2_error_baseline1 <- (mean((result_baseline1$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
beta_t_l2_error_baseline2 <- (mean((result_baseline2$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
# beta_t_l2_error_baseline3 <- (mean((result_baseline3$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
# beta_t_l2_error_baseline4 <- (mean((result_baseline4$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5

eta_l2_error <- (mean((result$eta - result_baseline_deconfounded$eta)**2))**0.5
eta_l2_error_baseline1 <- (mean((result_baseline1$eta - result_baseline_deconfounded$eta)**2))**0.5
eta_l2_error_baseline2 <- (mean((result_baseline2$eta - result_baseline_deconfounded$eta)**2))**0.5
# eta_l2_error_baseline3 <- (mean((result_baseline3$eta - result_baseline_deconfounded$eta)**2))**0.5
# eta_l2_error_baseline4 <- (mean((result_baseline4$eta - result_baseline_deconfounded$eta)**2))**0.5


eta_cate_error <- abs(sum((result$eta  - result_baseline_deconfounded$eta) *(means_x/std_x)))
eta_cate_error_baseline1 <- abs(sum((result_baseline1$eta  - result_baseline_deconfounded$eta) *(means_x/std_x)))
eta_cate_error_baseline2 <- abs(sum((result_baseline2$eta  - result_baseline_deconfounded$eta) *(means_x/std_x)))
# eta_cate_error_baseline3 <- abs(sum((result_baseline3$eta  - result_baseline_deconfounded$eta) *(means_x/std_x)))
# eta_cate_error_baseline4 <- abs(sum((result_baseline4$eta  - result_baseline_deconfounded$eta) *(means_x/std_x)))

error_collect <- c(beta_t_l2_error,beta_t_l2_error_baseline1,beta_t_l2_error_baseline2,
                   eta_l2_error,eta_l2_error_baseline1,eta_l2_error_baseline2,eta_cate_error,
                   eta_cate_error_baseline1,eta_cate_error_baseline2)

print(error_collect)

write.csv(error_collect,'results_fixed_p.csv')

print(abs(means_x%*%result$eta - means_x%*%result_baseline_deconfounded$eta))
print(abs(means_x%*%result_baseline1$eta - means_x%*%result_baseline_deconfounded$eta))
print(abs(means_x%*%result_baseline2$eta - means_x%*%result_baseline_deconfounded$eta))


folder_path = '/home/yuanyuanbenben/project_transfer_learning/transfer_learning_confounder/realdata1'
setwd(folder_path)

# load functions and data``
source(paste(folder_path,"/realdata_estimation_functions.R",sep=""))
load("gene_data_preprocessed.RData")

# select target tissues that have more than 50 samples
target_available = c()
for (target_index in 1:length(target_resp_files)) {
  target_response_file <- target_resp_files[target_index]
  target_Y_head <- try(
    read.csv(target_response_file),
    silent = TRUE
  )
  if (inherits(target_Y_head, "try-error")) {
    target_available = c(target_available,FALSE)
  } 
  else{
    target_Y <- target_Y_head[,2:dim(target_Y_head)[2]]
    if (length(target_Y) < 50){
      target_available = c(target_available,FALSE)
    }
    else{
      target_available = c(target_available,TRUE)
    }
  }
}

# selected target files 
target_cov_files <- target_cov_files[target_available]
target_resp_files <- target_resp_files[target_available]
target_confounder_files <- target_confounder_files[target_available]

# each time use a tissue as target 
target_num = length(target_cov_files)
source_num = length(source_cov_files)

if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
library(foreach)
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
library(doParallel)

cl.cores = detectCores(logical = F)
cl <- makeCluster(target_num + 2)
registerDoParallel(cl)

results_out <- foreach(target_index = 1:target_num,.verbose=TRUE,.combine = rbind,.packages=c("glmnet","lpSolve"))%dopar%{
  set.seed(123456)
  # target covariate 
  target_covariate_file <- target_cov_files[target_index]
  target_X_head <- read.csv(target_covariate_file)
  target_X <- target_X_head[,2:dim(target_X_head)[2]]
  
  # target response 
  target_response_file <- paste(sub("\\covariates.csv$", "", target_covariate_file),'response.csv',sep="")
  target_Y_head <- read.csv(target_response_file)
  target_Y <- target_Y_head[,2:dim(target_X_head)[2]]
  
  # target confounder
  target_confounder_file <- target_confounder_files[target_index]
  target_confounder_head <- read.csv(target_confounder_file)
  target_confounder <- target_confounder_head[,2:dim(target_confounder_head)[2]]
  
  
  # target normalization
  target_X = (target_X - mean_X) / (std_X)
  target_Y = target_Y - mean_Y
  target_confounder = (target_confounder - mean_confounder) / (std_confounder)
  target_confounder = t(target_confounder)
  target_X = t(target_X)
  target_Y = t(target_Y)
  
  # deconfounded target 
  P_target_confounder <- target_confounder %*% solve(t(target_confounder) %*% target_confounder, t(target_confounder))
  deconfounded_target_X <- target_X - P_target_confounder %*% target_X
  deconfounded_target_Y <- target_Y - P_target_confounder %*% target_Y
  
  # our method (no selection)
  set.seed(123456)
  result = transfer_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')
  
  # baselines
  result_baseline1 <- baseline1_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')
  result_baseline2 <- baseline2_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')
  result_baseline3 <- baseline3_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')
  result_baseline4 <- baseline4_estimation(source_X, source_Y, target_X, target_Y,trim_form = 'fixed')
  
  # deconfounded baseline
  result_baseline_deconfounded <- baseline2_estimation(deconfounded_source_X, deconfounded_source_Y, deconfounded_target_X, deconfounded_target_Y, trim_form = 'no_trim')
  
  # results 
  beta_t_l2_error <- (mean((result$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
  beta_t_l2_error_baseline1 <- (mean((result_baseline1$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
  beta_t_l2_error_baseline2 <- (mean((result_baseline2$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
  beta_t_l2_error_baseline3 <- (mean((result_baseline3$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
  beta_t_l2_error_baseline4 <- (mean((result_baseline4$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
  
  eta_l2_error <- (mean((result$eta - result_baseline_deconfounded$eta)**2))**0.5
  eta_l2_error_baseline1 <- (mean((result_baseline1$eta - result_baseline_deconfounded$eta)**2))**0.5
  eta_l2_error_baseline2 <- (mean((result_baseline2$eta - result_baseline_deconfounded$eta)**2))**0.5
  eta_l2_error_baseline3 <- (mean((result_baseline3$eta - result_baseline_deconfounded$eta)**2))**0.5
  eta_l2_error_baseline4 <- (mean((result_baseline4$eta - result_baseline_deconfounded$eta)**2))**0.5
  
  # our method (with selection)
  v_value_list = rep(0,source_num)
  for (source_index in 1:source_num){
    source_covariate_file <- source_cov_files[source_index]
    source_X_head <- read.csv(source_covariate_file)
    source_X <- source_X_head[,2:dim(source_X_head)[2]]
    
    source_response_file <- paste(sub("\\covariates.csv$", "", source_covariate_file),'response.csv',sep="")
    source_Y_head <- read.csv(source_response_file)
    source_Y <- source_Y_head[,2:dim(source_X_head)[2]]
    
    source_X <- t(source_X)
    source_Y <- t(source_Y)
    set.seed(123456)
    result = transfer_estimation_selection(source_X,source_Y,target_X,target_Y, trim_form = 'fixed')$v_value
    v_value_list[source_index] = result
  }
  v_min = min(v_value_list)
  rho_index = c(1,1.01,1.02,1.05,1.1,1.2,1.5,2,3,5)
  beta_t_loss_l2 = rep(0,length(rho_index))
  eta_loss_l2 = rep(0,length(rho_index))
  for (rho_i in 1:length(rho_index)){
    rho = rho_index[rho_i]
    v_index = v_value_list <= rho * v_min
    # selected source
    init_index = 0
    for (source_index in (1:length(source_cov_files))[v_index]){
      if (init_index == 0){
        source_covariate_file <- source_cov_files[source_index]
        source_X_head <- read.csv(source_covariate_file)
        source_X <- source_X_head[,2:dim(source_X_head)[2]]
        
        source_response_file <- paste(sub("\\covariates.csv$", "", source_covariate_file),'response.csv',sep="")
        source_Y_head <- read.csv(source_response_file)
        source_Y <- source_Y_head[,2:dim(source_X_head)[2]]
        init_index = 1
      }
      else{
        source_covariate_file <- source_cov_files[source_index]
        source_X_head <- read.csv(source_covariate_file)
        source_X <- cbind(source_X,source_X_head[,2:dim(source_X_head)[2]])
        
        source_response_file <- paste(sub("\\covariates.csv$", "", source_covariate_file),'response.csv',sep="")
        source_Y_head <- read.csv(source_response_file)
        source_Y <- cbind(source_Y,source_Y_head[,2:dim(source_X_head)[2]])
      }
    }
    source_X <- t(source_X)
    source_Y <- t(source_Y)
    # our method 
    set.seed(123456)
    result = transfer_estimation(source_X,source_Y,target_X,target_Y,trim_form = 'fixed')
    beta_t_loss_l2[rho_i] <- (mean((result$beta_t - result_baseline_deconfounded$beta_t)**2))**0.5
    eta_loss_l2[rho_i] <- (mean((result$eta - result_baseline_deconfounded$eta)**2))**0.5
  }
  
  error_collect <- c(beta_t_l2_error,beta_t_l2_error_baseline1,beta_t_l2_error_baseline2,beta_t_l2_error_baseline3,beta_t_l2_error_baseline4,beta_t_loss_l2,
                     eta_l2_error,eta_l2_error_baseline1,eta_l2_error_baseline2,eta_l2_error_baseline3,eta_l2_error_baseline4,eta_loss_l2)
  
  write.csv(error_collect,target_confounder_file |>
              sub(
                pattern = "data_preprocessed/", 
                replacement = "output/", 
                fixed = TRUE
              ) |>
              sub(
                pattern = ".v10_confounder.csv", 
                replacement = ".csv", 
                fixed = TRUE
              ))
  error_collect
}

stopCluster(cl)

output_files <- list.files(
  path = paste(folder_path,'/output',sep=""),
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

error_list <-rep(0,30)
for (i in 1:length(output_files)) {
  target_response_file <- target_resp_files[i]
  target_Y_head <- read.csv(target_response_file)
  target_Y <- target_Y_head[,2:dim(target_Y_head)[2]]
  error_list <- error_list + read.csv(output_files[i])$x
  print(target_response_file)
  print(length(target_Y))
  print(read.csv(output_files[i])$x)
}
print(error_list/length(output_files))
 
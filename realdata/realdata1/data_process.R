# using preprocessed data
rm(list = ls())
gc()

folder_path = '/home/yuanyuanbenben/project_transfer_learning/transfer_learning_confounder/realdata1'
setwd(folder_path)

# list all data files
covariates_csv_files <- list.files(
  path = paste(folder_path,'/data_preprocessed',sep=""),       
  pattern = "^gene.*covariates\\.csv$",         
  full.names = TRUE,           
  recursive = FALSE,           
  ignore.case = TRUE          
)

response_csv_files <- list.files(
  path = paste(folder_path,'/data_preprocessed',sep=""),       
  pattern = "^gene.*response\\.csv$",         
  full.names = TRUE,           
  recursive = FALSE,           
  ignore.case = TRUE          
)

confounder_csv_files <- list.files(
  path = paste(folder_path,'/data_preprocessed',sep=""),       
  pattern = ".*confounder\\.csv$",         
  full.names = TRUE,           
  recursive = FALSE,           
  ignore.case = TRUE          
)

# split source and target dataset
source_confounder_files <- confounder_csv_files[grepl("Brain", confounder_csv_files, ignore.case = TRUE)]
target_confounder_files <- confounder_csv_files[!grepl("Brain", confounder_csv_files, ignore.case = TRUE)]

source_cov_files <- covariates_csv_files[grepl("brain", covariates_csv_files, ignore.case = TRUE)]
source_resp_files <- response_csv_files[!grepl("brain", response_csv_files, ignore.case = TRUE)]
target_cov_files <- covariates_csv_files[!grepl("brain", covariates_csv_files, ignore.case = TRUE)]
target_resp_files <- response_csv_files[!grepl("brain", response_csv_files, ignore.case = TRUE)]


# data normalization
init_index = 0
for (source_index in 1:length(source_cov_files)){
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

mean_X <- rowMeans(source_X)
mean_Y <- rowMeans(source_Y)
std_X <- (diag(var(t(source_X))))**0.5
source_Y = source_Y - mean_Y
source_X = (source_X - mean_X) / (std_X)
source_X = t(source_X)
source_Y = t(source_Y)

init_index = 0
for (source_index in 1:length(source_confounder_files)){
  if (init_index == 0){
    source_confounder_file <- source_confounder_files[source_index]
    source_confounder_head <- read.csv(source_confounder_file)
    source_confounder <- source_confounder_head[,2:dim(source_confounder_head)[2]]
    init_index = 1
  }
  else{
    source_confounder_file <- source_confounder_files[source_index]
    source_confounder_head <- read.csv(source_confounder_file)
    source_confounder <- cbind(source_confounder,source_confounder_head[,2:dim(source_confounder_head)[2]])
  }
}

mean_confounder <- rowMeans(source_confounder)
std_confounder <- (diag(var(t(source_confounder))))**0.5
source_confounder = (source_confounder - mean_confounder) / (std_confounder)
source_confounder = t(source_confounder)


P_source_confounder <- source_confounder %*% solve(t(source_confounder) %*% source_confounder, t(source_confounder))
deconfounded_source_X <- source_X - P_source_confounder %*% source_X
deconfounded_source_Y <- source_Y - P_source_confounder %*% source_Y
# save as RData
save.image(file = "gene_data_preprocessed.RData")

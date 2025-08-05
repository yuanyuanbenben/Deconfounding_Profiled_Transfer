# preprocessing raw data

# set path
folder_path = '/home/yuanyuanbenben/project_transfer_learning/transfer_learning_confounder/realdata1'
setwd(folder_path)

# read raw gct files
read_gct <- function(file_path) {
  
  # title
  con <- file(file_path, "r")
  version <- readLines(con, n = 1)
  dim_line <- readLines(con, n = 1)
  close(con)
  
  dims <- scan(text = dim_line, quiet = TRUE)
  num_rows <- dims[1]
  num_cols <- dims[2]
  
  # raw data
  data <- read.delim(file_path, skip = 2, header = TRUE, check.names = FALSE)
  
  stopifnot(nrow(data) == num_rows)
  stopifnot(ncol(data) == num_cols + 2) 
  
  # expression matrix
  expr_matrix <- as.matrix(data[, -c(1, 2)])
  rownames(expr_matrix) <- data[, 1]
  
  list(
    version = version,
    dimensions = c(num_rows, num_cols),
    annotations = data[, 1:2],
    matrix = expr_matrix
  )
}

# matching the gene ID in CNS genes, enriched genes and the raw data
CNS_gene <- fromJSON("dataset/MODULE_137.v2024.1.Hs.json")$MODULE_137$geneSymbols

# enriched genes
enriched_gene <- fromJSON("dataset/enriched_module_137_ensembl_ids.json")
error_index <- is.na(enriched_gene$error)
enriched_gene_ensembl_ids <- enriched_gene$ensembl_id[error_index]
enriched_gene_entrez_ids <- enriched_gene$entrez_id[error_index]
# enriched_gene <- data.frame('ensembl_id' = enriched_gene_ensembl_ids, 'entrez_id' = enriched_gene_entrez_ids)

# gene covariates
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
library(jsonlite)

gct_result <- read_gct("dataset/gene_tpm_2022-06-06_v10_brain_amygdala.gct")
expression_data <- gct_result$matrix

# find genes ID in expression_data
expression_data_index <- c()
for (i in 1:length(CNS_gene)){
  if (CNS_gene[i] %in% gct_result$annotations$Description){
    expression_data_index <- c(expression_data_index, which(gct_result$annotations$Description == CNS_gene[i]))
  }
}
ensembl_clean <- sub("\\.\\d+$", "", gct_result$annotations$Name)
for (i in 1:length((enriched_gene_ensembl_ids))){
  if (enriched_gene_ensembl_ids[i] %in% ensembl_clean){
    expression_data_index <- c(expression_data_index, which(ensembl_clean == enriched_gene_ensembl_ids[i]))
  }
}

expression_data_index_unique <- unique(expression_data_index)
# write.csv(expression_data_index_unique,"dataset/selected_CNS_related_gene_index.csv")

# select the target JAM2 gene ID
JAM2_gene_index <- which(gct_result$annotations$Description == 'JAM2')
expression_data_index_unique <- expression_data_index_unique[expression_data_index_unique != JAM2_gene_index]

# select gene and construct expression data and response data for different tissues
gct_files <- list.files(
  path = paste(folder_path,'/dataset',sep=""),       
  pattern = "\\.gct$",         
  full.names = TRUE,           
  recursive = FALSE,           
  ignore.case = TRUE          
)

if (length(gct_files) == 0) {
  stop("no '.gct' files founded, please check the folder path!")
}

for (file in gct_files) {
  cat("waiting for", basename(file), "\n")
  gct_result <- read_gct(file)
  expression_data <- gct_result$matrix
  expression_data_selected <- expression_data[expression_data_index_unique,]
  # delete all zero expression genes
  expression_data_selected <- expression_data_selected[-c(202,592),]
  write.csv(expression_data_selected,paste(sub("\\.gct$", "", file),'covariates.csv',sep="_"))
  expression_data_y_selected <- matrix(expression_data[JAM2_gene_index,],1)
  write.csv(expression_data_y_selected,paste(sub("\\.gct$", "", file),'response.csv',sep="_"))
}

# list all covariates, responses and confounders files
covariates_csv_files <- list.files(
  path = paste(folder_path,'/dataset',sep=""),       
  pattern = "^gene.*covariates\\.csv$",         
  full.names = TRUE,           
  recursive = FALSE,           
  ignore.case = TRUE          
)

response_csv_files <- list.files(
  path = paste(folder_path,'/dataset',sep=""),       
  pattern = "^gene.*response\\.csv$",         
  full.names = TRUE,           
  recursive = FALSE,           
  ignore.case = TRUE          
)

confounder_txt_files <- list.files(
  path = paste(folder_path,'/dataset/GTEx_Analysis_v10_eQTL_covariates',sep=""),       
  pattern = "\\.txt$",         
  full.names = TRUE,           
  recursive = FALSE,           
  ignore.case = TRUE          
)

# delete tissues have no confounder data
covariates_csv_files <- covariates_csv_files[-c(24,25,31,35)]
response_csv_files <- response_csv_files[-c(24,25,31,35)]

# delete samples mismatched in covariates and confounders in each tissue
for (index in 1:length(confounder_txt_files)) {
  cat("waiting for", basename(covariates_csv_files[index]), "\n")
  
  # confounder data
  # using the first 5 principle components and first 15 inferred covariates
  confounder_file <- confounder_txt_files[index]
  confounder_head <- read.table(confounder_file)
  confounder <- confounder_head[c(2:21),2:dim(confounder_head)[2]]
  
  # covariate data
  covariate_file <- covariates_csv_files[index]
  X_head <- read.csv(covariate_file)
  X <- X_head[,2:dim(X_head)[2]]
  
  # response data
  response_file <- paste(sub("\\covariates.csv$", "", covariate_file),'response.csv',sep="")
  Y_head <- read.csv(response_file)
  Y <- Y_head[,2:dim(X_head)[2]]
  
  # ID_1: confounder data ID
  ID_index_1 <- sub("^GTEX-", "", confounder_head[1,2:dim(confounder_head)[2]])
  
  # ID_2: covariate data ID
  ID_index_2 <-  sub("^GTEX\\.([^.]+)\\..*", "\\1", names(X))
  
  # matching the ID
  ID_intersect <- intersect(ID_index_1, ID_index_2)
  confounder_index <- match(ID_intersect, ID_index_1)
  covariate_index <- match(ID_intersect, ID_index_2) 
  
  # select matched ID
  X <- X[,covariate_index]
  Y <- Y[,covariate_index]
  confounder <- confounder[,confounder_index]
  
  # save preprocessed data to the data_preprocessed folder
  write.csv(X,sub("dataset/", "data_preprocessed/", covariate_file, fixed = TRUE))
  write.csv(Y,sub("dataset/", "data_preprocessed/", response_file, fixed = TRUE))
  write.csv(confounder,confounder_file |>
              sub(
                pattern = "/dataset/GTEx_Analysis_v10_eQTL_covariates/", 
                replacement = "/data_preprocessed/", 
                fixed = TRUE
              ) |>
              sub(
                pattern = ".covariates.txt", 
                replacement = "_confounder.csv", 
                fixed = TRUE
              ))
}

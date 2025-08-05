# Codes for "Deconfounding via Profiled Transfer Learning"


## Overview

- ### Directory ***simulation***  
  Contains codes and outputs related to the simulations in *Section 5: Simulation Study*. 
  - ***single_source***: Simulation under the single-source scenario.  
  - ***multi_source***: Simulation under the multi-source scenario. 
    - ***output*** folder stores simulation results.  
    - ***varying_p.R*** and ***varying_phi.R***: Main R scripts for simulations under different $p$ and $\phi$ in the linear scenario.  
    - ***varying_p_nonlinear.R*** and ***varying_phi_nonlinear.R***: Main R scripts for simulations under different $p$ and $\phi$ in the nonlinear scenario.  
    - ***output_xxx.R***: Functions for saving results to the output folder.  
    - ***estimation_function.R*** and ***estimation_function_nonlinear.R***: Estimation methods proposed in our paper, including four baselines.  
    - ***utils.R*** and ***utils_nonlinear.R***: Data generation functions.  
    - ***varying_xxx.Rdata***: Simulation results.  
    - ***multi_source_selection.R*** (in ***multi_source***): Implements the source selection procedure.

- ### Directory ***realdata***  
  Contains datasets, codes, and outputs related to the two real-data experiments described in *Section 6: Real Data Examples*.  
  - **Subdirectory: realdata1 (Section 6.2: Gene Expression Data Analysis)**  
    - ***dataset***: Gene IDs for module 137 and enriched module 137.  
    - ***data_preprocessed***: Preprocessed dataset.  
    - ***output***: Results for different tissues.  
    - ***data_preprocess.R***: Functions for preprocessing (matching gene IDs, identifying covariates, responses, and confounders).  
    - ***data_process.R***: Functions for defining sources and targets with normalization.  
    - ***data_test_selection.R***: Main R script for running our method and baselines.  
    - ***data_estimation_function.R***: Estimation methods proposed in our paper, including four baselines.  

  - **Subdirectory: realdata2 (Section 6.1: Treatment Effect of Education on Earnings)**  
    - ***dataset***: Raw dataset.  
    - ***data_preprocess.R***: Functions for preprocessing (splitting control and treatment groups, imputing missing data, normalization).  
    - ***data_test.R***: Main R script for running our method and baselines.  
    - ***realdata_estimation_functions.R***: Estimation methods proposed in our paper, including four baselines.  
    - ***results.csv***: Output file.  

## Workflows

- ### Simulation: Performance comparison between our method and baselines
  1. Update the working directory in the `setwd(...)` function of the R scripts.  
  2. Run `varying_p.R` and `varying_p_nonlinear.R` (***single_source*** and ***multi_source***) for simulations under different $p$.  
  3. Run `output_p_linear.R` and `output_p_nonlinear.R` to save results in the ***output*** folder.  
  4. Run `varying_phi.R` and `varying_phi_nonlinear.R` for simulations under different $\phi$.  
  5. Run `output_phi_linear.R` and `output_phi_nonlinear.R` to save results.  

- ### Simulation: Source selection procedure  
  1. Update the working directory in `multi_source_selection.R` (***multi_source*** folder).  
  2. Run `multi_source_selection.R`.  

- ### Real-data Experiment 1: Treatment effect of education on earnings  
  1. Navigate to the ***realdata2*** folder.  
  2. Raw data is available [here](https://www.nlsinfo.org/investigator/pages/search?s=NLSM) or in the ***dataset*** folder.  
  3. Run `data_preprocess.R` for preprocessing.  
  4. Run `data_test.R` for applying our method and baselines.  

- ### Real-data Experiment 2: Gene expression data analysis  
  1. Navigate to the ***realdata1*** folder.  
  2. Download the raw data from [here](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression) or use the data in the ***dataset*** folder.  
  3. Run `data_preprocess.R` for preprocessing. Alternatively, use preprocessed data in ***data_preprocessed***.  
  4. Run `data_process.R` to prepare source and target data.  
  5. Run `data_test_selection.R` for results across different tissues.  

# Codes for "Deconfounding via Profiled Transfer Learning"


## Overview

- ### Directory ***simulation***  
  Contains codes and outputs related to the simulations in *Section 5: Simulation Study*. 
  - ***single_source***: Simulation under the single-source scenario.  
  - ***multi_source***: Simulation under the multi-source scenario. 
    - ***output***: Folder stores simulation results.  
    - ***varying_p.R*** and ***varying_phi.R***: Main R scripts for simulations under different $p$ and $\phi$ in the linear scenario.  
    - ***varying_p_nonlinear.R*** and ***varying_phi_nonlinear.R***: Main R scripts for simulations under different $p$ and $\phi$ in the nonlinear scenario.  
    - ***output_xxx.R***: Functions for saving results to the output folder.  
    - ***estimation_function.R*** and ***estimation_function_nonlinear.R***: Estimation methods proposed in our paper, including four baselines.  
    - ***utils.R*** and ***utils_nonlinear.R***: Data generation functions.  
    - ***varying_xxx.Rdata***: Simulation results.  
    - ***multi_source_selection.R*** (in ***multi_source***): Implements the source selection procedure.
  - ***additional_settings***: Simulation under more settings in supplementary materials. 
    - ***fixed_p***: Folder about the simulation in fixed-p scenario.
    - ***nonlinear_setting_3***: Folder about the third kind of nonlinear confounder structure in simulation
    - ***sensitivity_analysis***: Folder about the sensitivity analysis of the initial estimator $\hat{\beta}_{\text{init}}$ used in the model shift estimation. 

- ### Directory ***realdata***  
  Contains datasets, codes, and outputs related to the two real-data experiments described in *Section 6: Real Data Examples*.  
  - **Subdirectory: realdata1 (Section 6.1: Gene Expression Data Analysis)**  
    - ***dataset***: Gene IDs for module 137 and enriched module 137.  
    - ***data_preprocessed***: Preprocessed dataset.  
    - ***output***: Results for different tissues.  
    - ***data_preprocess.R***: Functions for preprocessing (matching gene IDs, identifying covariates, responses, and confounders).  
    - ***data_process.R***: Functions for defining sources and targets with normalization.  
    - ***data_test_selection.R***: Main R script for running our method and baselines.  
    - ***data_estimation_function.R***: Estimation methods proposed in our paper, including four baselines.  

  - **Subdirectory: realdata2 (Section 6.2, Section S.4.3: Treatment Effect of Education on Earnings)**  
    - ***dataset***: Raw dataset.  
    - ***data_preprocess.R***, ***data_preprocess_fixed_p.R***: Functions for preprocessing (splitting control and treatment groups, imputing missing data, normalization).  
    - ***data_test.R***, ***data_test_fixed_p.R***: Main R scripts for running our method and baselines in high-diemensional and fixed-p settings.  
    - ***realdata_estimation_functions.R***, ***realdata_estimation_functions_fixed_p.R***: Estimation methods proposed in our paper, including four baselines in the high-dimensional setting and two baselines in the fixed-p setting.  
    - ***results.csv***, ***results_fixed_p.csv***: Output files.  

## Workflows

- ### Simulation: Performance comparison between our method and baselines
  1. Update the working directory in the `setwd(...)` function of the R scripts.  
  2. Run `varying_p.R` and `varying_p_nonlinear.R` (in ***single_source*** and ***multi_source***) for simulations under different $p$.  
  3. Run `output_p_linear.R` and `output_p_nonlinear.R` to save results in the ***output*** folder.  
  4. Run `varying_phi.R` and `varying_phi_nonlinear.R` for simulations under different $\phi$.  
  5. Run `output_phi_linear.R` and `output_phi_nonlinear.R` to save results.  
  6. For the fixed-p simulation and more confounder structures, turn to the directory `~/simulation/fixed_p` and `~/simulation/nonlinear_setting_3` and the others are similar to steps 2-5. 

- ### Simulation: Source selection procedure  
  1. Update the working directory in `multi_source_selection.R` (***multi_source*** folder).  
  2. Run `multi_source_selection.R`.  

- ### Simulation: Sensitivity analysis
  1. Turn to the path `~/simulation/additional_settings/sensitivity_analysis/...`, the `single_source` and `multi_source` path are for single source and multi source respectively. 
  2. Run `varying_p.R` and `varying_p_nonlinear.R` for simulations under different $p$.  
  3. Run `varying_phi.R` and `varying_phi_nonlinear.R` for simulations under different $\phi$. 

- ### Real-data Experiment 1: Gene expression data analysis  
  1. Navigate to the ***realdata1*** folder.  
  2. Download the raw data from [here](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression) or use the data in the ***dataset*** folder.  
  3. Run `data_preprocess.R` for preprocessing. Alternatively, use preprocessed data in ***data_preprocessed***.  
  4. Run `data_process.R` to prepare source and target data.  
  5. Run `data_test_selection.R` for results across different tissues.  

- ### Real-data Experiment 2: Treatment effect of education on earnings  
  1. Navigate to the ***realdata2*** folder.  
  2. Raw data is available [here](https://www.nlsinfo.org/investigator/pages/search?s=NLSM) or in the ***dataset*** folder.  
  3. Run `data_preprocess.R`, `data_preprocess_fixed_p.R` for preprocessing.  
  4. Run `data_test.R`, `data_test_fixed_p.R` for applying our method and baselines.  

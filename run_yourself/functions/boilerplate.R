
# Libraries ----------------------------------------
# remotes::install_github("burchill/zplyr")
# remotes::install_github("burchill/cs")

library(tidyverse)
library(future)
library(furrr)
library(cs)
library(beepr)
library(broom.mixed) # need to load this for weird S3 reasons


# Constants ------------------------------------

data_path     = paste0(main_path, "og_data/")
models_path   = paste0(main_path, "models/")
model_files_path = models_path
bb_files_path = paste0(main_path, "bb_files/")
mm_files_path = paste0(main_path, "mm_files/")
complete_path = paste0(main_path, "collated_files/")
code_path     = paste0(main_path, "functions/")
script_path   = paste0(main_path, "runs_and_scripts/")


nsc_file_path = data_path
nsc_df_file = "nsc_sim_df.RDS"
nsc_noexcl_df_file = "nsc_noexcl_fillers.RDS"

setal_file_path = data_path
setal_df_file = "setal_sim_df.RDS"
setal_noexcl_df_file = "setal_noexcl_fillers.RDS"

fetal_file_path = data_path
fetal_df_file = "fetal_sim_df.RDS"
fetal_noexcl_df_file = "fetal_noexcl_fillers.RDS"


# Loading code -------------------------------------------

source(paste0(code_path, "utils.R"))
source(paste0(code_path, "amlap_constants.R"))
source(paste0(code_path, "amlap_modeling_data_functions.R"))
source(paste0(code_path, "modeling_data_functions.R"))
source(paste0(code_path, "amlap_bb_making_functions.R"))

done <- cs::done



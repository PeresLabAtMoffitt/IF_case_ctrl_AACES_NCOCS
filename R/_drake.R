# For help with drake
# https://ropenscilabs.github.io/drake-manual/projects.html#safer-interactivity or https://docs.ropensci.org/drake
# and
# https://ropensci.github.io/drake/reference/r_make.html or https://docs.ropensci.org/drake

# Load your packages from packages.R
source("R/packages_AACES1_NCOCS.R")
# Load the code as function that drake need to run
source("R/function_AACES1_NCOCS.R") 
# Load the plan that drake has to execute
source("R/plan_AACES1_NCOCS.R")    

# End drake
config <- drake_config(plan, parallelism = "future", jobs = 4, verbose = 1)
if (!interactive()) config

# make(plan)
loadd(#clinical_data, ancestry_data, tx_data,
      ROI_tumor ,ROI_stroma ,ROI_total,
      TMA_tumor ,TMA_stroma, TMA_total,
      TMAcases_remove, ROIcases_remove,
      common_ROITMA_IDs, cases_match)

# Cleaning
rm(fct_name_repair, var_names, # roir_import,
   data_import,ancestry_import, tx_import, 
   roit_import, rois_import, roi_import, 
   tmat_import, tmas_import, tma17_import, tma18_import,
   tma2t_import, tma2s_import, common_ROITMA_import, # tmar_import, 
   case_remove_import, binding, match_cases_import,
   plan, config)


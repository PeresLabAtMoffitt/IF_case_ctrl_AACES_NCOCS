# Import library

library(tidyverse); library(janitor)
library(haven)


############################################################################## I ### Load new ROIs data----
path <- fs::path("","Volumes","Peres_Research")
fct_name_repair <- function(colnms) {
  tolower(
    gsub("\\-", "minus", 
         (gsub("\\+", "plus", colnms))
         ))
}
#-----------------------------------------------------------------------------------------------------------------
aaces_clinical <- 
  read_csv(
    paste0(path,
           "/K99_R00/Image analysis data/AACES MCC18207 mIF Data/aaces_all.csv"))
ncocs_clinical <- 
  read_csv(
    paste0(path,
           "/K99_R00/Image analysis data/AACES MCC18207 mIF Data/ncocs_all.csv"))
path2 <- fs::path("", "Volumes", "Peres_Research", "AACES2", "Analgesic medications and survival")
survival_time_Nov2023 <-
  read_sas(paste0(path2, 
                  "/data/raw data/analgesics_with_phase2.sas7bdat")) %>% 
  select(suid, days_int_to_event, vital_status_fin) %>% 
  rename(vital_status_nov2023 = vital_status_fin) %>% 
  mutate(suid = as.character(suid))
#-----------------------------------------------------------------------------------------------------------------
ROI_tumor_2022jan <-
  readxl::read_xlsx(
    paste0(path,
           "/K99_R00/Image analysis data/AACES MCC18207 mIF Data/L.Peres_P1_AACES MCC18207_ROI_Results-2022-1.xlsx"
  ), sheet = "Tumor (PCK+)",
  .name_repair = fct_name_repair) %>% 
  janitor::clean_names(replace=janitor:::mu_to_u)
ROI_stroma_2022jan <-
  readxl::read_xlsx(
    paste0(path,
           "/K99_R00/Image analysis data/AACES MCC18207 mIF Data/L.Peres_P1_AACES MCC18207_ROI_Results-2022-1.xlsx"
  ), sheet = "Stroma (PCK-)",
  .name_repair = fct_name_repair) %>% 
  janitor::clean_names(replace=janitor:::mu_to_u)
ROI_total_2022jan <-
  readxl::read_xlsx(
    paste0(path,
           "/K99_R00/Image analysis data/AACES MCC18207 mIF Data/L.Peres_P1_AACES MCC18207_ROI_Results-2022-1.xlsx"
    ), sheet = "AACES MCC18207 ROI Counts",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names(replace=janitor:::mu_to_u) %>% 
  rename_at(vars(starts_with("fox") | 
                   starts_with("cd") |
                   starts_with("percent_fox") | 
                   starts_with("percent_cd") |
                   starts_with("area")
  ), ~ paste0("total_", .)) 


# End loading

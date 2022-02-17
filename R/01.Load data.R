# Import library

library(tidyverse); library(janitor)


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
           "/K99_R00/Image analysis data//AACES MCC18207 mIF Data/aaces_all.csv"))
ncocs_clinical <- 
  read_csv(
    paste0(path,
           "/K99_R00/Image analysis data//AACES MCC18207 mIF Data/ncocs_all.csv"))
#-----------------------------------------------------------------------------------------------------------------
ROI_tumor <-
  readxl::read_xlsx(
    paste0(path,
           "/K99_R00/Image analysis data//AACES MCC18207 mIF Data/L.Peres_P1_AACES MCC18207_ROI_Results-2022-1.xlsx"
  ), sheet = "Tumor (PCK+)",
  .name_repair = fct_name_repair) %>% 
  janitor::clean_names()
ROI_stroma <-
  readxl::read_xlsx(
    paste0(path,
           "/K99_R00/Image analysis data//AACES MCC18207 mIF Data/L.Peres_P1_AACES MCC18207_ROI_Results-2022-1.xlsx"
  ), sheet = "Stroma (PCK-)",
  .name_repair = fct_name_repair) %>% 
  janitor::clean_names()
ROI_total <-
  readxl::read_xlsx(
    paste0(path,
           "/K99_R00/Image analysis data//AACES MCC18207 mIF Data/L.Peres_P1_AACES MCC18207_ROI_Results-2022-1.xlsx"
    ), sheet = "AACES MCC18207 ROI Counts",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names() %>% 
  `colnames<-`(c(paste0("total_", colnames(.))))


# End loading

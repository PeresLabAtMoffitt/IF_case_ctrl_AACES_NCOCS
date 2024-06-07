# Import library

library(tidyverse); library(janitor)


############################################################################## I ### Load new ROIs data----
path <- fs::path("","Volumes","Peres_Research", "K99_R00",
                 "Image analysis data", "Panel3_MOTIF_Data_2024")
fct_name_repair <- function(colnms) {
  tolower(
    gsub("\\-", "minus", 
         (gsub("\\+", "plus", colnms))
    ))
}
#-----------------------------------------------------------------------------------------------------------------
TMA2017_total_2024june <-
  readxl::read_xlsx(
    paste0(path,
           "/TMA_Data/AACES 2017A (June 2024)/Peres_P3_MOTIF_AACES 2017A TMA_June2024_Results.xlsx"
    ), sheet = "Peres_P3_MOTIF_AACES 2017A TMA_",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names()
TMA2017_tumor_2024june <-
  readxl::read_xlsx(
    paste0(path,
           "/TMA_Data/AACES 2017A (June 2024)/Peres_P3_MOTIF_AACES 2017A TMA_June2024_Results.xlsx"
    ), sheet = "Tumor",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names()
TMA2017_stroma_2024june <-
  readxl::read_xlsx(
    paste0(path,
           "/TMA_Data/AACES 2017A (June 2024)/Peres_P3_MOTIF_AACES 2017A TMA_June2024_Results.xlsx"
    ), sheet = "Stroma",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names()

TMA2018_total_2024june <-
  readxl::read_xlsx(
    paste0(path,
           "/TMA_Data/AACES 2018A (June 2024)/Peres_P3_MOTIF_AACES 2018A TMA_June2024_Results.xlsx"
    ), sheet = "Peres_P3_MOTIF_AACES 2018A TMA_",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names()
TMA2018_tumor_2024june <-
  readxl::read_xlsx(
    paste0(path,
           "/TMA_Data/AACES 2018A (June 2024)/Peres_P3_MOTIF_AACES 2018A TMA_June2024_Results.xlsx"
    ), sheet = "Tumor",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names()
TMA2018_stroma_2024june <-
  readxl::read_xlsx(
    paste0(path,
           "/TMA_Data/AACES 2018A (June 2024)/Peres_P3_MOTIF_AACES 2018A TMA_June2024_Results.xlsx"
    ), sheet = "Stroma",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names()

TMAaaces2_total_2024june <-
  readxl::read_xlsx(
    paste0(path,
           "/TMA_Data/TMA 1 (June 2024)/Peres_P3_MOTIF_TMA1_June2024_Results.xlsx"
    ), sheet = "Total_Summary_Results",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names()
TMAaaces2_tumor_2024june <-
  readxl::read_xlsx(
    paste0(path,
           "/TMA_Data/TMA 1 (June 2024)/Peres_P3_MOTIF_TMA1_June2024_Results.xlsx"
    ), sheet = "Tumor",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names()
TMAaaces2_stroma_2024june <-
  readxl::read_xlsx(
    paste0(path,
           "/TMA_Data/TMA 1 (June 2024)/Peres_P3_MOTIF_TMA1_June2024_Results.xlsx"
    ), sheet = "Stroma",
    .name_repair = fct_name_repair) %>% 
  janitor::clean_names()
#-----------------------------------------------------------------------------------------------------------------
path1 <- fs::path("","Volumes","Peres_Research", "K99_R00",
                 "Image analysis data", "TMA manifests")
TMA2017_manifest <-
  readxl::read_xlsx(
    paste0(path1,
           "/TMA manifest exhaustion panel/TMAmanifest_2017.xlsx"
    ), sheet = "suid_core_match") %>% 
  janitor::clean_names()
TMA2018_manifest <-
  readxl::read_xlsx(
    paste0(path1,
           "/TMA manifest exhaustion panel/TMAmanifest_2018.xlsx"
    ), sheet = "suid_core_match") %>% 
  janitor::clean_names()
#-----------------------------------------------------------------------------------------------------------------
# Use list.files as data was saved over multiple files
file_list <- list.files(
  path =
    "/Volumes/Peres_Research/K99_R00/Image analysis data/Panel3_MOTIF_Data_2024/WTS_Data",
  pattern = "*.xlsx",
  recursive=FALSE,
  full.names = TRUE)

WTS_total_2024june <- do.call("rbind",lapply(Sys.glob(file_list), readxl::read_xlsx,
                                          sheet = 1,
                              .name_repair = fct_name_repair)) %>% 
  janitor::clean_names()
WTS_tumor_2024june <- do.call("rbind",lapply(Sys.glob(file_list), readxl::read_xlsx,
                                             sheet = "Tumor",
                              .name_repair = fct_name_repair)) %>% 
  janitor::clean_names()
WTS_stroma_2024june <- do.call("rbind",lapply(Sys.glob(file_list), readxl::read_xlsx,
                                             sheet = "Stroma",
                               .name_repair = fct_name_repair)) %>% 
  janitor::clean_names()

# End data load

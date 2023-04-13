
# I  ### Load data

path <- fs::path("","Volumes","Peres_Research")

fct_name_repair <- function(colnms) {
  tolower(gsub("[ ():]", "_", colnms))
}
data_import <- function(data_path){
  read_csv(paste0(data_path,
              "/K99_R00/Image analysis data/AACES and NCOCS data/aaces_ncocs_08142020.csv"))
}
ancestry_import <- function(data_path){
  read_csv(paste0(data_path,
                  "/K99_R00/Image analysis data/AACES and NCOCS data/abund_ancestry_brca.csv"))
}
tx_import <- function(data_path){
  readxl::read_xlsx(paste0(data_path,
                  "/Christelle Colin-Leitzinger/IF_AACES_NCOCS/data/AACES_tx_12082020.xlsx"))
}
# tx_import2 <- function(data_path){
#   read_csv(paste0(data_path,
#                   "/K99_R00/Image analysis data/AACES and NCOCS data/ncocstreat.csv"))
# }
#-----------------------------------------------------------------------------------------------------------------
roit_import <- function(path){
  readxl::read_xlsx(paste0(path,
              "/K99_R00/Image analysis data/Immune marker count data/Peres P1 ROI Set 1 Index 2 July 2021.xlsx"
  ), sheet = "Tumor", .name_repair = fct_name_repair)
}
rois_import <- function(path){
  readxl::read_xlsx(paste0(path,
              "/K99_R00/Image analysis data/Immune marker count data/Peres P1 ROI Set 1 Index 2 July 2021.xlsx"
  ), sheet = "Stroma", .name_repair = fct_name_repair)
}
roi_import <- function(path){
  readxl::read_xlsx(paste0(path,
                           "/K99_R00/Image analysis data/Immune marker count data/Peres P1 ROI Set 1 Index 2 July 2021.xlsx"
  ), sheet = "Peres P1 ROI Set 1 Index 2 July", .name_repair = fct_name_repair)
}
# roir_import <- function(path){
#   readxl::read_xlsx(paste0(path,
#               "/K99_R00/Image analysis data/Immune marker count data/Peres P1 ROI Analysis with location_updated 1-24-2020.xlsx"
#   ), sheet = "Analysis Info", skip = 24) %>% 
#     select(c(1, 5)) %>% 
#     `colnames<-`(c("image_tag", "REMOVED")) %>% 
#     drop_na("image_tag")
# }
#-----------------------------------------------------------------------------------------------------------------
tmat_import <- function(path){
  readxl::read_xlsx(paste0(path,
                "/K99_R00/Image analysis data/Immune marker count data/Peres P1 AACES 2017 TMA Summary.xlsx"
  ), sheet = "Tumor", .name_repair = fct_name_repair)
}
tmas_import <- function(path){
  readxl::read_xlsx(paste0(path,
                "/K99_R00/Image analysis data/Immune marker count data/Peres P1 AACES 2017 TMA Summary.xlsx"
  ), sheet = "Stroma", .name_repair = fct_name_repair)
}
tma17_import <- function(path){
  readxl::read_xlsx(paste0(path,
                "/K99_R00/Image analysis data/Immune marker count data/Peres P1 AACES 2017 TMA Summary.xlsx"
  ), sheet = "Peres P1 AACES 2017 TMA", .name_repair = fct_name_repair)
}
tma2t_import <- function(path){
  readxl::read_xlsx(paste0(path,
                "/K99_R00/Image analysis data/Immune marker count data/Peres P1 AACES 2018 TMA Summary.xlsx"
  ), sheet = "Tumor", .name_repair = fct_name_repair)
}
tma2s_import <- function(path){
  readxl::read_xlsx(paste0(path,
                "/K99_R00/Image analysis data/Immune marker count data/Peres P1 AACES 2018 TMA Summary.xlsx"
  ), sheet = "Stroma", .name_repair = fct_name_repair)
}
tma18_import <- function(path){
  readxl::read_xlsx(paste0(path,
                "/K99_R00/Image analysis data/Immune marker count data/Peres P1 AACES 2018 TMA Summary.xlsx"
  ), sheet = "Peres P1 AACES 2018 TMA", .name_repair = fct_name_repair)
}
#-----------------------------------------------------------------------------------------------------------------
case_remove_import <- function(path){
  read_csv(paste0(path,
                "/Christelle Colin-Leitzinger/IF_AACES_NCOCS/data/Subject_IDs to remove from TMA.csv"))
}
caseROI_remove_import <- function(path){
  read_csv(paste0(path,
                  "/Christelle Colin-Leitzinger/IF_AACES_NCOCS/data/Subject_IDs to remove from ROI.csv"))
}
imageROI_remove_import <- function(path){
  read_csv(paste0(path,
                  "/Christelle Colin-Leitzinger/IF_AACES_NCOCS/data/Remove images that were not stained properly in batch1.csv"))
}
#-----------------------------------------------------------------------------------------------------------------
common_ROITMA_import <- function(path){
  read_csv(paste0(path,
                  "/Christelle Colin-Leitzinger/IF_AACES_NCOCS/data/Subject_IDs common TMA ROI.csv"))
}
#-----------------------------------------------------------------------------------------------------------------
match_cases_import <- function(path){
  readxl::read_xlsx(paste0(path,
                           "/K99_R00/Image analysis data/AACES and NCOCS data/Case matches_12312019.xlsx"))
}
#-----------------------------------------------------------------------------------------------------------------

# Start Data Cleaning
binding <- function(data1, data2){
  bind_rows(data1, data2, .id = "TMA")
}

# Update variable names
# var_names <- function(data){
#   rename_all(data %>% gsub("__", "_", .) %>% gsub("%", "percent", .))
# }
var_names <- function(data) {
  names <- names(data) %<>%
    tolower() %>%
    str_replace_all("__", "_") %>%
    gsub("%", "percent", .) %>%
    gsub("²", "2", .) %>%
    gsub("\\+", "plus", .)
  data %>% `colnames<-`(names)
}





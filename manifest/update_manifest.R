# Import library

library(tidyverse); library(janitor)


############################################################################## I ### Load data----
path <- fs::path("","Volumes","Peres_Research", "AACES2", "Immune panel manifest")
path1 <- fs::path("","Volumes","Peres_Research")

# Manifest
manifest <- 
  readxl::read_xlsx(
    paste0(path,"/AACES_IF panel manifest.xlsx"), na = "NA")

# markers data
markers <- readRDS("/Users/colinccm/Documents/GitHub/Peres/IF_case_ctrl_AACES_NCOCS/markers_AACES_NCOCS_batch1_2_07132023.rds")

# Emory manifest
manifest_emory <- 
  readxl::read_xlsx(
    paste0(path,"/Panel 2_Emory to Moffitt Sept 2023.xlsx")) %>% 
  janitor::clean_names()

# Dukes manifest
manifest_duke <- 
  readxl::read_xlsx(
    paste0(path,"/Panel 2_Batch 1_Slides from Duke.xlsx")) %>% 
  janitor::clean_names()

tma_ids <- 
  readxl::read_xlsx(
    paste0(path1,"/K99_R00/Image analysis data/TMA manifests/IDs on TMAs.xlsx")) %>% 
  janitor::clean_names()


ROI_failed_2022jan <-
  readxl::read_xlsx(
    paste0(path1,
           "/K99_R00/Image analysis data/AACES MCC18207 mIF Data/L.Peres_P1_AACES MCC18207_ROI_Results-2022-1.xlsx"
    ), sheet = "Notes", skip = 38)

ROI_failed_2021july <- 
  readxl::read_xlsx(paste0(path1,
                           "/K99_R00/Image analysis data/Immune marker count data/Peres P1 ROI Set 1 Index 2 July 2021.xlsx"
  ), sheet = "Info", skip = 53)

roi_failed_excluded_patients <- 
  read_csv(paste0(path1,
                  "/Christelle Colin-Leitzinger/IF_AACES_NCOCS/data/Subject_IDs to remove from ROI.csv"
  ))
roi_failed_not_properly_stained <- 
  read_csv(paste0(path1,
                  "/Christelle Colin-Leitzinger/IF_AACES_NCOCS/data/Remove images that were not stained properly in batch1.csv"
  ))

tma_failed_2017 <- 
  readxl::read_xlsx(paste0(path1,
                           "/K99_R00/Image analysis data/Immune marker count data/Peres P1 AACES 2017 TMA Summary.xlsx"
  ), sheet = "Analysis Information")

tma_failed_2018 <- 
  readxl::read_xlsx(paste0(path1,
                           "/K99_R00/Image analysis data/Immune marker count data/Peres P1 AACES 2018 TMA Summary.xlsx"
  ), sheet = "Analysis Information", skip = 18)

tma_failed_Lauren_list <- 
  read_csv(paste0(path1,
                  "/Christelle Colin-Leitzinger/IF_AACES_NCOCS/data/Subject_IDs to remove from TMA.csv"
  ))



############################################################################## II ### Clean data----
# Duke
manifest_duke <- manifest_duke %>% 
  # Clean and add var
  mutate(suid = as.character(suid)) %>% 
  select(suid) %>%
  mutate(site = "Duke") %>% 
  mutate(panel2_slide_type = "Whole section")

# Emory
# Fix Emory TMA ids
tma_ids <- tma_ids %>% 
  mutate(suid = as.character(suid)) %>% 
  # pivot longer for patient who has 2 TMA for both year
  separate_wider_delim(cols = tma, delim = " and ",
                       names = c("year1", "year2"), 
                       too_few = "align_start", too_many = "error", 
                       cols_remove = TRUE) %>% 
  pivot_longer(cols = -suid, names_to = NULL, 
               values_drop_na = TRUE) %>% 
  mutate(tma_year = value) %>% 
  mutate(panel2_slide_type = "TMA")

manifest_emory <- manifest_emory %>% 
  mutate(suid = as.character(suid)) %>% 
  # Add the TMA ida
  full_join(., tma_ids,
            by = c("suid" = "value")) %>% 
  mutate(suid = coalesce(suid.y, suid), .keep = "unused") %>% 
  # Clean and add var
  select(suid, phase : panel2_slide_type) %>% 
  mutate(panel2_slide_type = case_when(
    panel2_slide_type == "TMA"           ~ "TMA",
    suid == "TMA"                        ~ "TMA",
    is.na(panel2_slide_type)             ~ "Whole section"
  )) %>% 
  mutate(phase = paste("AACES", phase)) %>% 
  mutate(site = "Emory")

non_mcc <- bind_rows(manifest_emory,
                     manifest_duke)

markers <- markers %>% 
  mutate(data_version = str_replace(data_version, "_v", "_batch")) %>% 
  # include AACES only - ids have 6 digits
  mutate(str_id = str_count(suid), .after = suid) %>% 
  filter(str_id == 6) %>% 
  # create var of interest
  mutate(has_intratumoral = case_when(
    annotation == "Intratumoral"       ~ "Yes"
  )) %>% 
  mutate(has_peripheral = case_when(
    annotation == "Peripheral"         ~ "Yes"
  )) %>% 
  mutate(has_ROI = case_when(
    slide_type == "ROI"                ~ "Yes"
  )) %>% 
  mutate(has_TMA = case_when(
    slide_type == "TMA"                ~ "Yes"
  )) %>% 
  select(suid, data_version, has_intratumoral, has_peripheral,
         has_ROI, has_TMA) %>% 
  group_by(suid, data_version) %>% 
  fill(has_intratumoral, has_peripheral, has_ROI, has_TMA, 
       .direction = "updown") %>% 
  distinct() %>% 
  ungroup() %>% 
  arrange(data_version, suid)
############################################################################## III ### Merge data----
# 2 patient are "Complete - Consented but no/insufficient tissue" + 
# "Yes" in panel 1 status but we have them in AACES 1 so considered as completed

manifest_ <- manifest %>% 
  mutate(suid = as.character(suid)) %>% 
  full_join(non_mcc, ., by = "suid") %>% 
    mutate(site = case_when(
      !is.na(site)                                ~ site,
      is.na(site)                                 ~ "Moffitt"
    ))

manifest_1 <- manifest_ %>%
  full_join(., markers, by = "suid") %>% 
  select(suid, data_version, has_ROI, has_TMA, everything()) %>% 
  # group_by(suid) %>%
  # mutate(imaged_sequence = row_number(suid)) %>% 
  # ungroup() %>% 
  mutate(phase = case_when(
    !is.na(phase)                              ~ phase,
    is.na(phase)                               ~ "AACES 1"
  ))

manifest_2 <- manifest_1 %>% 
  mutate(panel_1_data_from_whole_sections = case_when(
    (data_version == "AACES_batch1_NCOCS" | 
      data_version == "AACES_batch2") &
      has_ROI == "Yes"                         ~ "yes",
    (`AACES tissue status` == "Consented but never received tissue or problem with tissue" |
      `AACES tissue status` == "Did not consent for tissue procurement") &
      is.na(data_version)                      ~ "no"
  ), .after = `Panel 1 Data from whole sections`) %>% 
  # mutate(panel_2_data_from_whole_sections = case_when(
  #   (data_version == "AACES_batch1_NCOCS" | 
  #      data_version == "AACES_batch2") &
  #     panel == "panel 2" &
  #     has_ROI == "Yes"                         ~ "yes"
  # ), .after = `Panel 2 Data from whole sections`) %>% 
  
  mutate(panel_1_data_from_tma = case_when(
    (data_version == "AACES_batch1_NCOCS" | 
      data_version == "AACES_batch2") &
      has_TMA == "Yes"                         ~ "yes",
    is.na(data_version) &
      is.na(has_TMA) &
      !is.na(`Panel 1 TMA received at MCC`)    ~ "no/failed",
    (`AACES tissue status` == "Consented but never received tissue or problem with tissue" |
       `AACES tissue status` == "Did not consent for tissue procurement") &
      is.na(data_version)                      ~ "no"
  ), .after = `Panel 1 Data from TMA`) %>% 
  # mutate(panel_2_data_from_tma = case_when(
  #   (data_version == "AACES_batch1_NCOCS" | 
  #      data_version == "AACES_batch2") &
  #     panel == "panel 2" &
  #     has_TMA == "Yes"                         ~ "yes"
  # ), .after = `Panel 2 Data from TMA`) %>% 
  
  
  mutate(panel_1_status = case_when(
    `Panel 1 status` == "Did not consent"       ~ "Did not consent",
      data_version == "AACES_batch1_NCOCS" | 
      data_version == "AACES_batch2"            ~ "Completed",
      (`Panel 1 Whole sections received at MCC` == "yes" |
      `Panel 1 TMA received at MCC` == "yes") &
      is.na(data_version)                       ~ "Section received but no data",
      is.na(data_version) &
      `Panel 1 status` == "Failed"              ~ "Failed",
      is.na(data_version) &
      `Panel 1 status` == "Complete"            ~ "Weird, need to go back"
  ), .after = `Panel 1 status`) #%>% 
  
  # mutate(panel_2_status = case_when(
  #   panel == "panel 2" &
  #     data_version == "AACES_batch1_NCOCS" | 
  #     data_version == "AACES_batch2"            ~ "Completed",
  #   panel == "panel 2" &
  #     (`Panel 2 Whole sections received at MCC` == "yes" |
  #        `Panel 2 TMA received at MCC` == "yes") &
  #     is.na(data_version)                       ~ "Section received but no data",
  #   panel == "panel 2" &
  #     is.na(data_version) &
  #     `Panel 1 status` == "Failed"              ~ "Failed",
  #   panel == "panel 2" &
  #     is.na(data_version) &
  #     `Panel 1 status` == "Complete"            ~ "Weird, need to go back"
  #   
  # ), .after = `Panel 2 status`)
  

# panel 1 or 2??? 150133


############################################################################## IV ### Add reason for failing ----
# ROI_failed_2022jan
ROI_failed_2022jan_removed <- ROI_failed_2022jan %>% 
  slice(1 : 43) %>% purrr::keep(~!all(is.na(.))) %>% 
  mutate(`Removed cores:` = str_replace(`Removed ROI`, "16-", "16"),
         suid = str_match(`Removed cores:`,
                          "(L.Peres_P1_OV|L.Peres_P1_)([:digit:]*)")[,3],
         .before = 1) %>% 
  mutate(cause_ROI_core_remove = 
           "Removed at first analysis by Core") %>% 
  select(suid, cause_ROI_core_remove) %>% 
  group_by(suid) %>% 
  mutate(number_of_roi_core_removed = n()) %>% 
  ungroup() %>% 
  distinct()

ROI_failed_2022jan_cc <- ROI_failed_2022jan %>% 
  slice(48 : 53) %>% 
  select(`Removed ROI`) %>% 
  mutate(`Removed cores:` = str_replace(`Removed ROI`, "16-", "16"),
                  suid = str_match(`Removed cores:`,
                                   "(L.Peres_P1_OV|L.Peres_P1_)([:digit:]*)")[,3],
                  .before = 1) %>% 
  mutate(cause_ROI_core_remove = 
           'Removed by Christelle after Jonathan said it is "out of focus" and he forgot to remove it (email 02.17.2022)') %>% 
  select(suid, cause_ROI_core_remove) %>% 
  filter(!is.na(suid)) %>% 
  group_by(suid) %>% 
  mutate(number_of_roi_core_removed = n()) %>% 
  ungroup() %>% 
  distinct()


# ROI_failed_2021july
ROI_failed_2021july_PKCremoved <- ROI_failed_2021july %>% 
  slice(1 : 13) %>% 
  mutate(`Removed cores:` = 
           str_replace(`Discarded 1/24/2020 - Poor PCK Staining = no tumor cells counting`, 
                       "16-", "16"),
         suid = str_match(`Removed cores:`,
                          "(Peres_P1_AACES |Peres_P1_)([:digit:]*)")[,3],
         .before = 1) %>% 
  mutate(cause_ROI_core_remove = 
           'Discarded 1/24/2020 - Poor PCK Staining = no tumor cells counting') %>% 
  select(suid, cause_ROI_core_remove) %>% 
  group_by(suid) %>% 
  mutate(number_of_roi_core_removed = n()) %>% 
  ungroup() %>% 
  distinct()

ROI_failed_2021july_removed <- ROI_failed_2021july %>% 
  slice(15 : 69) %>% 
  mutate(`Removed cores:` = 
           str_replace(`Discarded 1/24/2020 - Poor PCK Staining = no tumor cells counting`, 
                       "16-", "16"),
         suid = str_match(`Removed cores:`,
                          "(Peres_P1_AACES |Peres_P1_|Peres_P3_)([:digit:]*)")[,3],
         .before = 1) %>% 
  mutate(cause_ROI_core_remove = 
           'ROIs removed during initial analysis in Dec/Jan') %>% 
  select(suid, cause_ROI_core_remove) %>% 
  filter(!is.na(suid)) %>% 
  group_by(suid) %>% 
  mutate(number_of_roi_core_removed = n()) %>% 
  ungroup() %>% 
  distinct()

roi_failed_excluded_patients <- roi_failed_excluded_patients %>% 
  mutate(number_of_roi_core_removed = "Patient excluded") %>% 
  mutate(number_of_roi_core_removed = "All") %>% 
  rename(suid = Subject_IDs) %>% 
  mutate(suid = as.character(suid))

roi_failed_not_properly_stained <- roi_failed_not_properly_stained %>% 
  mutate(`Removed cores:` = 
           str_replace(remove_images_that_were_not_stained_properly, 
                       "16-", "16"),
         suid = str_match(`Removed cores:`,
                          "(Peres_P1_AACES |Peres_P1_|Peres_P3_)([:digit:]*)")[,3],
         .before = 1) %>% 
  mutate(cause_ROI_core_remove = paste("Image stil in data but need to be removed for each analysis,",
                                       "not properly stained,",
                                       "Jonathan forgot to remove them during the 2021 halo run")) %>% 
  select(suid, cause_ROI_core_remove) %>% 
  group_by(suid) %>% 
  mutate(number_of_roi_core_removed = n()) %>% 
  ungroup() %>% 
  distinct()
  
# tma_failed_2017 <- 
#   readxl::read_xlsx(paste0(path1,
#                            "/K99_R00/Image analysis data/Immune marker count data/Peres P1 AACES 2017 TMA Summary.xlsx"
#   ), sheet = "Analysis Information", skip = 18)
# 
# tma_failed_2017 <- tma_failed_2017 %>% 
#   slice(1:54) %>% 
#   mutate(`Removed cores:` = str_replace(`Removed cores:`, "16-", "16"),
#          suid = str_match(`Removed cores:`,
#                           "(L.Peres_P1_OV|L.Peres_P1_)([:digit:]*)")[,3],
#          .before = 1)

tma_failed_2018 <- tma_failed_2018 %>% 
  rename(suid = ...5) %>% 
  mutate(cause_tma2018_core_remove = "Removed at first analysis by Core") %>% 
  select(suid, cause_tma2018_core_remove) %>% 
  group_by(suid) %>% 
  mutate(number_of_tma_core_removed = n()) %>% 
  ungroup() %>% 
  distinct()

tma_failed_Lauren_list <- tma_failed_Lauren_list %>% 
  mutate(cause_tma2018_core_remove = "Patient removed initially from Lauren's list")%>% 
  mutate(number_of_tma_core_removed = "All") %>% 
  rename(suid = Subject_IDs) %>% 
  mutate(suid = as.character(suid))


cause_core_remove <- bind_rows(ROI_failed_2021july_removed, ROI_failed_2021july_PKCremoved,
                               ROI_failed_2022jan_removed, ROI_failed_2022jan_cc,
                               roi_failed_not_properly_stained,
                               tma_failed_2018) %>% 
  mutate(number_of_tma_core_removed = as.character(number_of_tma_core_removed)) %>% 
  mutate(number_of_roi_core_removed = as.character(number_of_roi_core_removed)) %>% 
  bind_rows(., tma_failed_Lauren_list, roi_failed_excluded_patients) %>% 
  group_by(suid) %>% 
  summarise_at(vars(cause_ROI_core_remove, number_of_roi_core_removed,
                    cause_tma2018_core_remove, number_of_tma_core_removed), 
               str_c, collapse = " + ")


############################################################################## IV ### Merge manifest with reason for failing ----
manifest_final <- manifest_2 %>% 
  left_join(., cause_core_remove, 
            by = "suid") %>% 
  select(suid, site, phase, data_version, 
         has_ROI_data_panel1 = has_ROI, has_TMA_data_panel1 = has_TMA, tma_year,
         "AACES tissue status",
         "Panel 1 status" : panel_1_data_from_tma,
         panel2_slide_type,
         "Panel 2 status" : "Panel 2 Data from TMA",
         
         cause_ROI_core_remove, number_of_roi_core_removed,
         has_intratumoral, has_peripheral,
         cause_tma2018_core_remove, number_of_tma_core_removed,
         everything()
         ) %>% 
  mutate(panel_1_status = case_when(
    str_detect(panel_1_status, "Weird") &
      `Panel 1 status` == "Complete" &
      str_detect(cause_tma2018_core_remove,
                 "Patient remove")            ~ "Completed but exclude patient after the fact",
    is.na(panel_1_status) &
      str_detect(`Panel 1 status`, 
                 "Consented but")             ~ `Panel 1 status`,
    TRUE                                      ~ panel_1_status,
  )) %>% 
  mutate(flag = case_when(
    `AACES tissue status` == "Consented but never received tissue or problem with tissue" &
      panel_1_status == "Completed"           ~ "Why data/completed if never received tissues?",
    panel_1_status == "Failed" &
    panel_1_data_from_tma == "no/failed"      ~ "Reason TMA failed?",
    panel_1_status == "Section received but no data" &
      `Panel 1 status` == "Failed"            ~ "Reason ROI failed?",
    # `Panel 1 Whole sections received at MCC` == "yes" &
    #   !is.na(has_ROI_data)                     ~ "Why ROI failed?",
    has_TMA_data_panel1 == "Yes" &
    is.na(`Panel 1 Data from TMA`) &
      is.na(`Panel 1 TMA received at MCC`) &
      panel_1_status == "Completed"           ~ "Why TMA not received?"
  ), .after = panel_1_status) %>% 
  mutate(`Panel 2 Whole sections received at MCC` = case_when(
    panel2_slide_type == "Whole section"      ~ "Yes"
  )) %>% 
  mutate(`Panel 2 TMA received at MCC` = case_when(
    panel2_slide_type == "TMA"                ~ "Yes"
  )) %>% 
  select(-`Panel 1 status`)


write_csv(manifest_final, "manifest/AACES_IF panel manifest updated 02_15_2024.csv")


# End 
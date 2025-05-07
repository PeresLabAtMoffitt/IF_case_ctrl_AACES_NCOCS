# Import library

library(tidyverse); library(janitor)


############################################################################## I ### Load data----
path <- fs::path("","Volumes","Peres_Research", "AACES2", "Immune panel manifest")
path1 <- fs::path("","Volumes","Peres_Research", "AACES", "Panel 2 Immune Markers/Marker Data")

# Manifest
manifest_2024 <- 
  read_csv(paste0(path,
           "/AACES_IF panel manifest updated 02_15_2024.csv"))

core_inventory <- 
  readxl::read_xlsx(paste0(path,
                           "/Peres Inventory with Final Status Comments_JVN.xlsx")) %>% 
  clean_names()

# Maps
wts_map <- 
  read_csv(paste0(path1,
                  "/AACES_Panel2_WTS_ID_Map_23Oct2024.csv"))
tma_map <- 
  read_csv(paste0(path1,
                  "/AACES_Panel2_TMA_ID_Map_23Oct2024.csv"))

# WTS marker
markers_wts <- 
  read_csv(paste0(path1,
                  "/AACES_Panel2_WTS_immune_markers_23Oct2024.csv")) %>% 
  mutate(suid = as.character(suid))
# TMA marker
markers_tma <- 
  read_csv(paste0(path1,
                  "/AACES_Panel2_TMA_core_immune_markers_23Oct2024.csv")) %>% 
  mutate(suid = as.character(suid))

############################################################################## II ### Clean data----
core_inventory1 <- core_inventory %>% 
  filter(slide_id == "L.Peres" | slide_id == "PERES") %>% 
  mutate(site_panel2 = case_when(
    slide_id == "L.Peres"                ~ "Duke",
    slide_id == "PERES"                  ~ "Emory"
  )) %>% 
  mutate(suid = as.character(suid)) %>% 
  mutate(suid_ = str_replace(slide_id_2, "16-", "16"),
         suid2 = str_remove(suid_, "AACES |ACCESS |OV"),
         suid_ = str_match(suid2, "(SP.*?OR |)([:digit:]*)")[, 3],
         suid2 = str_remove(suid_, "^2017$|^2018$"),
         suid_ = na_if(suid2, "")) %>% 
  mutate(suid_ = coalesce(suid, suid_)) %>% 
  mutate(tma_year = case_when(
    is.na(suid_) &
      slide_id_2 == "TMA 1"        ~ 2024,
    is.na(suid_)                    ~ as.numeric(str_extract(slide_id_2, "[:digit:]+"))
  )) %>% 
  mutate(`AACES tissue status panel2` = "Consented") %>% # Fix here "Consented and received tissue"-------------
  mutate("Panel 2 Whole sections received at MCC pan2" = case_when(
    type == "WTS"                  ~ "yes"
  )) %>% 
  mutate("Panel 2 TMA received at MCC pan2" = case_when(
    type == "TMA"                  ~ "yes"
  )) %>% 
  select(suid_, tma_year, site_panel2, `AACES tissue status panel2`, type) %>% 
  distinct()

markers_wts1 <- markers_wts %>% 
  arrange(suid) %>% 
  mutate(removed_wts_reason = case_when(
    low_tumor == 1              ~ 1
  )) %>% 
  filter(remove_dup == 0) %>% 
  mutate(has_WTS_data_panel2 = case_when(
    is.na(removed_wts_reason)              ~ "Yes"
  )) %>% 
  select(suid, has_WTS_data_panel2, 
         `number of WTS removed for low tumor content panel2` = removed_wts_reason, 
         wts_panel2_batch_set = set)
  
markers_tma1 <- markers_tma %>% 
  arrange(suid) %>% 
  mutate(removed_tma_reason = case_when(
    low_tumor == 1              ~ "number of TMA removed for low tumor content panel2",
    TRUE                        ~ "number of good core panel 2"
  )) %>% 
  filter(remove_dup == 0) %>% 
  select(suid, removed_tma_reason, tma_year = tma) %>% 
  group_by(suid, tma_year, removed_tma_reason) %>% 
  mutate(number_adequate_core = n()) %>% 
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(id_cols = c(suid, tma_year), 
              names_from = removed_tma_reason, 
              values_from = number_adequate_core) %>% 
  mutate(has_TMA_data_panel2 = case_when(
    !is.na(`number of good core panel 2`)              ~ "Yes"
  )) %>% 
  mutate(tma_year = as.numeric(str_remove(tma_year, "TMA ")))
  
core_inventory1_ <- core_inventory1 %>% 
  full_join(., markers_tma1 %>% 
              select(suid, tma_year), by = "tma_year") %>% 
  mutate(suid = coalesce(suid_, suid)) %>% 
  select(-suid_)

manifest_2025 <- manifest_2024 %>% 
  mutate(tma_year = case_when(
    !is.na(tma_year)                ~ tma_year,
    suid == "TMA"                   ~ 2024
  )) %>% 
  full_join(markers_tma1 %>% 
              select(new_suid = suid, tma_year) %>% 
              filter(tma_year == 2024),
            ., 
            by = "tma_year") %>% 
  mutate(suid = case_when(
    suid == "TMA"                   ~ new_suid,
    !is.na(suid)                    ~ suid
  )) %>% 
  select(-new_suid) %>% 
  full_join(., core_inventory1_, by = c("suid", "tma_year")) %>%  # need tms_year as some were run in 2017 AND 2028
  full_join(., markers_wts1, by = "suid") %>% 
  full_join(., markers_tma1, by = c("suid", "tma_year")) %>% 
  mutate(`Panel 2 status` = case_when(
    has_WTS_data_panel2 == "Yes"          ~ "Completed",
    has_TMA_data_panel2 == "Yes"          ~ "Completed",
    `number of WTS removed for low tumor content panel2` == 1
                                          ~ "Received but WTS removed",
    is.na(`number of good core panel 2`) &
      !is.na(`number of TMA removed for low tumor content panel2`)
                                          ~ "Received but all TMA cores removed",
    cause_tma2018_core_remove ==
      "Patient removed initially from Lauren's list"
                                          ~ "Patient removed initially from Lauren's list (TMA 2018)",
    is.na(has_WTS_data_panel2) &
    is.na(has_TMA_data_panel2)            ~ "No core received for panel 2"
  )) %>% 
  mutate("Panel 2 Data from whole sections"= case_when(
    `Panel 2 Whole sections received at MCC` == "Yes" &
      has_WTS_data_panel2 == "Yes"        ~ "Yes",
    `Panel 2 Whole sections received at MCC` == "Yes" &
      is.na(has_WTS_data_panel2)          ~ "No"
  )) %>% 
  mutate("Panel 2 Data from TMA"= case_when(
    `Panel 2 TMA received at MCC` == "Yes" &
      has_TMA_data_panel2 == "Yes"        ~ "Yes",
    `Panel 2 TMA received at MCC` == "Yes" &
      is.na(has_TMA_data_panel2)          ~ "No"
  )) %>% 
  select(-has_WTS_data_panel2, -has_TMA_data_panel2, -type) %>% 
  select(suid : has_TMA_data_panel1, tma_year,
         "AACES tissue status" : "panel_1_data_from_tma",
         "cause_ROI_core_remove" : "Panel 1 Block ID",
         "AACES tissue status panel2",
         panel2_slide_type : `Panel 2 Data from TMA`,
         "wts_panel2_batch_set",
         "number of WTS removed for low tumor content panel2",
         "number of good core panel 2",
         "number of TMA removed for low tumor content panel2") %>% 
  rename(data_version_panel1 = data_version,
         "AACES tissue status panel 1" = "AACES tissue status")


write_rds(manifest_2025, 
          paste0("manifest/AACES_IF panel manifest updated_", today(), ".rds"))
write_csv(manifest_2025, 
          paste0("manifest/AACES_IF panel manifest updated_", today(), ".csv"))

write_rds(manifest_2025, 
          paste0(path, "/AACES_IF panel manifest updated_", today(), ".rds"))

write_csv(manifest_2025, 
          paste0(path, "/AACES_IF panel manifest updated_", today(), ".csv"))








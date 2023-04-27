library(tidyverse)

clinical_markers_AACES_NCOCS_batch1_2_04132023 <- 
  readRDS("~/Documents/GitHub/Peres/IF_case_ctrl_AACES_NCOCS/clinical_markers_AACES_NCOCS_batch1_2_04132023.rds")

a <- clinical_markers_AACES_NCOCS_batch1_2_04132023 %>% 
  select(suid, data_version, slide_type, site) %>%
  filter(site == "AAS") %>% 
  distinct()

b <- clinical_markers_AACES_NCOCS_batch1_2_04132023 %>% 
  select(suid, data_version, slide_type, site) %>%
  filter(site == "AAS") %>% 
  arrange(slide_type) %>% 
  distinct(suid, .keep_all = TRUE)

table(a$data_version)
table(a$slide_type)

table(b$data_version)
table(b$slide_type)

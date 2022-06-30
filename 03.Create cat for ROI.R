# Run only for doing the preliminary analysis on ROIs only
library(tidyverse)

complete_AACES_NCOCS_global <- 
  read_rds(paste0(here::here(), 
                  "/allmarkers_AACES_NCOCS_global.rds"))

ROI_global <- complete_AACES_NCOCS_global %>% filter(slide_type == "ROI")


############################################################################## VII ### Create cat with immune cells----
# Create average by patients by annotation
markers_full <- ROI_global %>% 
  select(-image_tag) %>% 
  group_by(suid, version, annotation) %>% 
  mutate(across(where(is.numeric), .fns = ~ mean(.))) %>% 
  ungroup() %>% 
  distinct()

# Create wide format to separate each markers
markers_ROI <-
  markers_full %>% 
  mutate_at(("annotation"), ~ case_when(
    annotation == "Peripheral"            ~ "p",
    annotation == "Intratumoral"          ~ "i",
    annotation == "Stromal"               ~ "s"
  )) %>% 
  pivot_wider(names_from = annotation,
              values_from = -c(suid, annotation, version),
              names_glue = "{.value}_{annotation}"
  )

# Create cat
markers_ROI <- markers_ROI %>% 
  mutate(cd3_tumor_i = case_when(
    tumor_percent_cd3_i <= 1      ~ "low",
    tumor_percent_cd3_i > 1       ~ "high"
  ), cd3_tumor_i = factor(cd3_tumor_i, levels = c("low","high"))) %>%
  mutate(cd3plus_cd8plus_tumor_i = case_when(
    tumor_percent_cd3plus_cd8plus_i <= 1      ~ "low",
    tumor_percent_cd3plus_cd8plus_i > 1       ~ "high"
  ), cd3plus_cd8plus_tumor_i = factor(cd3plus_cd8plus_tumor_i, levels = c("low","high"))) %>%
  mutate(cd3plus_foxp3plus_tumor_i = case_when(
    tumor_percent_cd3plus_foxp3plus_i <= 1      ~ "low",
    tumor_percent_cd3plus_foxp3plus_i > 1       ~ "high"
  ), cd3plus_foxp3plus_tumor_i = factor(cd3plus_foxp3plus_tumor_i, levels = c("low","high"))) %>%
  mutate(cd11b_tumor_i = case_when(
    tumor_percent_cd11b_i <= 1      ~ "low",
    tumor_percent_cd11b_i > 1       ~ "high"
  ), cd11b_tumor_i = factor(cd11b_tumor_i, levels = c("low","high"))) %>%
  mutate(cd11bplus_cd15plus_tumor_i = case_when(
    tumor_percent_cd11bplus_cd15plus_i == 0      ~ "Absence",
    tumor_percent_cd11bplus_cd15plus_i > 0       ~ "Presence"
  ), cd11bplus_cd15plus_tumor_i = factor(cd11bplus_cd15plus_tumor_i, levels = c("Absence","Presence"))) %>% 
  mutate(cd3_stroma_i = case_when(
    stroma_percent_cd3_i <= 1      ~ "low",
    stroma_percent_cd3_i > 1       ~ "high"
  ), cd3_stroma_i = factor(cd3_stroma_i, levels = c("low","high"))) %>%
  mutate(cd3plus_cd8plus_stroma_i = case_when(
    stroma_percent_cd3plus_cd8plus_i <= 1      ~ "low",
    stroma_percent_cd3plus_cd8plus_i > 1       ~ "high"
  ), cd3plus_cd8plus_stroma_i = factor(cd3plus_cd8plus_stroma_i, levels = c("low","high"))) %>%
  mutate(cd3plus_foxp3plus_stroma_i = case_when(
    stroma_percent_cd3plus_foxp3plus_i <= 1      ~ "low",
    stroma_percent_cd3plus_foxp3plus_i > 1       ~ "high"
  ), cd3plus_foxp3plus_stroma_i = factor(cd3plus_foxp3plus_stroma_i, levels = c("low","high"))) %>%
  mutate(cd11b_stroma_i = case_when(
    stroma_percent_cd11b_i <= 1      ~ "low",
    stroma_percent_cd11b_i > 1       ~ "high"
  ), cd11b_stroma_i = factor(cd11b_stroma_i, levels = c("low","high"))) %>%
  mutate(cd11bplus_cd15plus_stroma_i = case_when(
    stroma_percent_cd11bplus_cd15plus_i == 0      ~ "Absence",
    stroma_percent_cd11bplus_cd15plus_i > 0       ~ "Presence"
  ), cd11bplus_cd15plus_stroma_i = factor(cd11bplus_cd15plus_stroma_i, levels = c("Absence","Presence"))) %>% 
  mutate(cd3_total_i = case_when(
    total_percent_cd3_i <= 1      ~ "low",
    total_percent_cd3_i > 1       ~ "high"
  ), cd3_total_i = factor(cd3_total_i, levels = c("low","high"))) %>%
  mutate(cd3plus_cd8plus_total_i = case_when(
    total_percent_cd3plus_cd8plus_i <= 1      ~ "low",
    total_percent_cd3plus_cd8plus_i > 1       ~ "high"
  ), cd3plus_cd8plus_total_i = factor(cd3plus_cd8plus_total_i, levels = c("low","high"))) %>%
  mutate(cd3plus_foxp3plus_total_i = case_when(
    total_percent_cd3plus_foxp3plus_i <= 1      ~ "low",
    total_percent_cd3plus_foxp3plus_i > 1       ~ "high"
  ), cd3plus_foxp3plus_total_i = factor(cd3plus_foxp3plus_total_i, levels = c("low","high"))) %>%
  mutate(cd11b_total_i = case_when(
    total_percent_cd11b_i <= 1      ~ "low",
    total_percent_cd11b_i > 1       ~ "high"
  ), cd11b_total_i = factor(cd11b_total_i, levels = c("low","high"))) %>%
  mutate(cd11bplus_cd15plus_total_i = case_when(
    total_percent_cd11bplus_cd15plus_i == 0      ~ "Absence",
    total_percent_cd11bplus_cd15plus_i > 0       ~ "Presence"
  ), cd11bplus_cd15plus_total_i = factor(cd11bplus_cd15plus_total_i, levels = c("Absence","Presence"))) %>% 
  
  mutate(cd3_tumor_p = case_when(
    tumor_percent_cd3_p <= 1      ~ "low",
    tumor_percent_cd3_p > 1       ~ "high"
  ), cd3_tumor_p = factor(cd3_tumor_p, levels = c("low","high"))) %>%
  mutate(cd3plus_cd8plus_tumor_p = case_when(
    tumor_percent_cd3plus_cd8plus_p <= 1      ~ "low",
    tumor_percent_cd3plus_cd8plus_p > 1       ~ "high"
  ), cd3plus_cd8plus_tumor_p = factor(cd3plus_cd8plus_tumor_p, levels = c("low","high"))) %>%
  mutate(cd3plus_foxp3plus_tumor_p = case_when(
    tumor_percent_cd3plus_foxp3plus_p <= 1      ~ "low",
    tumor_percent_cd3plus_foxp3plus_p > 1       ~ "high"
  ), cd3plus_foxp3plus_tumor_p = factor(cd3plus_foxp3plus_tumor_p, levels = c("low","high"))) %>%
  mutate(cd11b_tumor_p = case_when(
    tumor_percent_cd11b_p <= 1      ~ "low",
    tumor_percent_cd11b_p > 1       ~ "high"
  ), cd11b_tumor_p = factor(cd11b_tumor_p, levels = c("low","high"))) %>%
  mutate(cd11bplus_cd15plus_tumor_p = case_when(
    tumor_percent_cd11bplus_cd15plus_p == 0      ~ "Absence",
    tumor_percent_cd11bplus_cd15plus_p > 0       ~ "Presence"
  ), cd11bplus_cd15plus_tumor_p = factor(cd11bplus_cd15plus_tumor_p, levels = c("Absence","Presence"))) %>% 
  mutate(cd3_stroma_p = case_when(
    stroma_percent_cd3_p <= 1      ~ "low",
    stroma_percent_cd3_p > 1       ~ "high"
  ), cd3_stroma_p = factor(cd3_stroma_p, levels = c("low","high"))) %>%
  mutate(cd3plus_cd8plus_stroma_p = case_when(
    stroma_percent_cd3plus_cd8plus_p <= 1      ~ "low",
    stroma_percent_cd3plus_cd8plus_p > 1       ~ "high"
  ), cd3plus_cd8plus_stroma_p = factor(cd3plus_cd8plus_stroma_p, levels = c("low","high"))) %>%
  mutate(cd3plus_foxp3plus_stroma_p = case_when(
    stroma_percent_cd3plus_foxp3plus_p <= 1      ~ "low",
    stroma_percent_cd3plus_foxp3plus_p > 1       ~ "high"
  ), cd3plus_foxp3plus_stroma_p = factor(cd3plus_foxp3plus_stroma_p, levels = c("low","high"))) %>%
  mutate(cd11b_stroma_p = case_when(
    stroma_percent_cd11b_p <= 1      ~ "low",
    stroma_percent_cd11b_p > 1       ~ "high"
  ), cd11b_stroma_p = factor(cd11b_stroma_p, levels = c("low","high"))) %>%
  mutate(cd11bplus_cd15plus_stroma_p = case_when(
    stroma_percent_cd11bplus_cd15plus_p == 0      ~ "Absence",
    stroma_percent_cd11bplus_cd15plus_p > 0       ~ "Presence"
  ), cd11bplus_cd15plus_stroma_p = factor(cd11bplus_cd15plus_stroma_p, levels = c("Absence","Presence"))) %>% 
  mutate(cd3_total_p = case_when(
    total_percent_cd3_p <= 1      ~ "low",
    total_percent_cd3_p > 1       ~ "high"
  ), cd3_total_p = factor(cd3_total_p, levels = c("low","high"))) %>%
  mutate(cd3plus_cd8plus_total_p = case_when(
    total_percent_cd3plus_cd8plus_p <= 1      ~ "low",
    total_percent_cd3plus_cd8plus_p > 1       ~ "high"
  ), cd3plus_cd8plus_total_p = factor(cd3plus_cd8plus_total_p, levels = c("low","high"))) %>%
  mutate(cd3plus_foxp3plus_total_p = case_when(
    total_percent_cd3plus_foxp3plus_p <= 1      ~ "low",
    total_percent_cd3plus_foxp3plus_p > 1       ~ "high"
  ), cd3plus_foxp3plus_total_p = factor(cd3plus_foxp3plus_total_p, levels = c("low","high"))) %>%
  mutate(cd11b_total_p = case_when(
    total_percent_cd11b_p <= 1      ~ "low",
    total_percent_cd11b_p > 1       ~ "high"
  ), cd11b_total_p = factor(cd11b_total_p, levels = c("low","high"))) %>%
  mutate(cd11bplus_cd15plus_total_p = case_when(
    total_percent_cd11bplus_cd15plus_p == 0      ~ "Absence",
    total_percent_cd11bplus_cd15plus_p > 0       ~ "Presence"
  ), cd11bplus_cd15plus_total_p = factor(cd11bplus_cd15plus_total_p, levels = c("Absence","Presence")))


######################################################################################## V ### Create immunoscore----
markers_ROI <- markers_ROI %>% 
  mutate(percentile_score_CD3_i = ntile(total_percent_cd3_i, 100) ) %>% 
  mutate(percentile_score_CD3_p = ntile(total_percent_cd3_p, 100) ) %>% 
  mutate(percentile_score_CD8_i = ntile(total_percent_cd8_i, 100) ) %>% 
  mutate(percentile_score_CD8_p = ntile(total_percent_cd8_p, 100) ) %>%
  mutate(percentile_score_mean = 
           rowMeans(select(.,`percentile_score_CD3_i`:`percentile_score_CD8_p`), 
                    na.rm = TRUE)
  ) %>%
  mutate(immunoscore_patients = case_when(
    percentile_score_mean <= 10        ~ 0,
    percentile_score_mean <= 25        ~ 1,
    percentile_score_mean <= 70        ~ 2,
    percentile_score_mean <= 95        ~ 3,
    percentile_score_mean > 95         ~ 4 
  )) %>% 
  mutate(immunoscore_patients = factor(immunoscore_patients)) %>% 
  mutate(immunoscore_2018lancet_patients = case_when(
    percentile_score_mean <= 25        ~ "Low",
    percentile_score_mean <= 70        ~ "Intermediate",
    percentile_score_mean > 70         ~ "High"
  )) %>% 
  mutate(immunoscore_2018lancet_patients = 
           factor(immunoscore_2018lancet_patients, 
                  levels = c("Low", "Intermediate", "High"))) 


######################################################################################## VI ### Join data----
markers_ROI <- right_join(case_ctrl_data, markers_ROI, 
                          by = "suid")
ROI_global <- right_join(case_ctrl_data, ROI_global,
                         by = "suid")

saveRDS(ROI_global, file = "ROI_global.rds")
saveRDS(markers_ROI, file = "markers_ROI.rds")

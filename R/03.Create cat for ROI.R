# Run only for doing the preliminary analysis on ROIs only
library(tidyverse)


############################################################################## I ### Load data----
marker_AACES_NCOCS <- 
  read_rds(paste0(here::here(), 
                  "/markers_AACES_NCOCS_batch1_2_07132023.rds"))
case_ctrl_data <- readRDS(paste0(here::here(), "/case_ctrl_data_07132023.rds"))


############################################################################## II ### Summarize----
# Create average by patients by annotation
markers_summ <- marker_AACES_NCOCS %>% 
  select(-image_tag) %>% 
  group_by(suid, data_version, annotation) %>% 
  mutate(across(where(is.numeric), .fns = ~ mean(.))) %>% 
  ungroup() %>% 
  distinct()

write_rds(markers_summ, "summarized_markers_clinical_AACES_NCOCS_batch1_2_07132023.rds")


ROI_summ <- markers_summ %>% filter(slide_type == "ROI")
ROI_global <- marker_AACES_NCOCS %>% filter(slide_type == "ROI")


# Create wide format to separate each markers
markers_ROI <-
  ROI_summ %>% 
  mutate_at(("annotation"), ~ case_when(
    annotation == "Peripheral"            ~ "p",
    annotation == "Intratumoral"          ~ "i",
    annotation == "Stromal"               ~ "s"
  )) %>% 
  pivot_wider(names_from = annotation,
              values_from = -c(suid, annotation, data_version, slide_type),
              names_glue = "{.value}_{annotation}"
  )

# ROI_global1 <- ROI_global %>% 
#   # select(-image_tag) %>% 
#   mutate_at(("annotation"), ~ case_when(
#     annotation == "Peripheral"            ~ "p",
#     annotation == "Intratumoral"          ~ "i",
#     annotation == "Stromal"               ~ "s"
#   )) %>% 
#   pivot_wider(names_from = annotation,
#               values_from = -c(suid, annotation, data_version, slide_type),
#               names_glue = "{.value}_{annotation}"
#   )


############################################################################## III ### Create cat with immune cells----
markers_ROI <- markers_ROI %>% 
    mutate(cd3_tumor_cat.i = case_when(
      tumor_percent_cd3_i <= 1      ~ "low",
      tumor_percent_cd3_i > 1       ~ "high"
    ), cd3_tumor_cat.i = factor(cd3_tumor_cat.i, levels = c("low","high"))) %>%
    mutate(cd3plus_cd8plus_tumor_cat.i = case_when(
      tumor_percent_cd3plus_cd8plus_i <= 1      ~ "low",
      tumor_percent_cd3plus_cd8plus_i > 1       ~ "high"
    ), cd3plus_cd8plus_tumor_cat.i = factor(cd3plus_cd8plus_tumor_cat.i, levels = c("low","high"))) %>%
    mutate(cd3plus_foxp3plus_tumor_cat.i = case_when(
      tumor_percent_cd3plus_foxp3plus_i <= 1      ~ "low",
      tumor_percent_cd3plus_foxp3plus_i > 1       ~ "high"
    ), cd3plus_foxp3plus_tumor_cat.i = factor(cd3plus_foxp3plus_tumor_cat.i, levels = c("low","high"))) %>%
    mutate(cd11b_tumor_cat.i = case_when(
      tumor_percent_cd11b_i == 0      ~ "Absence",
      tumor_percent_cd11b_i > 0       ~ "Presence"
    ), cd11b_tumor_cat.i = factor(cd11b_tumor_cat.i, levels = c("Absence","Presence"))) %>%
    mutate(cd11bplus_cd15plus_tumor_cat.i = case_when(
      tumor_percent_cd11bplus_cd15plus_i == 0      ~ "Absence",
      tumor_percent_cd11bplus_cd15plus_i > 0       ~ "Presence"
    ), cd11bplus_cd15plus_tumor_cat.i = factor(cd11bplus_cd15plus_tumor_cat.i, levels = c("Absence","Presence"))) %>% 
    mutate(cd3_stroma_cat.i = case_when(
      stroma_percent_cd3_i <= 1      ~ "low",
      stroma_percent_cd3_i > 1       ~ "high"
    ), cd3_stroma_cat.i = factor(cd3_stroma_cat.i, levels = c("low","high"))) %>%
    mutate(cd3plus_cd8plus_stroma_cat.i = case_when(
      stroma_percent_cd3plus_cd8plus_i <= 1      ~ "low",
      stroma_percent_cd3plus_cd8plus_i > 1       ~ "high"
    ), cd3plus_cd8plus_stroma_cat.i = factor(cd3plus_cd8plus_stroma_cat.i, levels = c("low","high"))) %>%
    mutate(cd3plus_foxp3plus_stroma_cat.i = case_when(
      stroma_percent_cd3plus_foxp3plus_i <= 1      ~ "low",
      stroma_percent_cd3plus_foxp3plus_i > 1       ~ "high"
    ), cd3plus_foxp3plus_stroma_cat.i = factor(cd3plus_foxp3plus_stroma_cat.i, levels = c("low","high"))) %>%
    mutate(cd11b_stroma_cat.i = case_when(
      stroma_percent_cd11b_i == 0      ~ "Absence",
      stroma_percent_cd11b_i > 0       ~ "Presence"
    ), cd11b_stroma_cat.i = factor(cd11b_stroma_cat.i, levels = c("Absence","Presence"))) %>%
    mutate(cd11bplus_cd15plus_stroma_cat.i = case_when(
      stroma_percent_cd11bplus_cd15plus_i == 0      ~ "Absence",
      stroma_percent_cd11bplus_cd15plus_i > 0       ~ "Presence"
    ), cd11bplus_cd15plus_stroma_cat.i = factor(cd11bplus_cd15plus_stroma_cat.i, levels = c("Absence","Presence"))) %>% 
    mutate(cd3_total_cat.i = case_when(
      total_percent_cd3_i <= 1      ~ "low",
      total_percent_cd3_i > 1       ~ "high"
    ), cd3_total_cat.i = factor(cd3_total_cat.i, levels = c("low","high"))) %>%
    mutate(cd3plus_cd8plus_total_cat.i = case_when(
      total_percent_cd3plus_cd8plus_i <= 1      ~ "low",
      total_percent_cd3plus_cd8plus_i > 1       ~ "high"
    ), cd3plus_cd8plus_total_cat.i = factor(cd3plus_cd8plus_total_cat.i, levels = c("low","high"))) %>%
    mutate(cd3plus_foxp3plus_total_cat.i = case_when(
      total_percent_cd3plus_foxp3plus_i <= 1      ~ "low",
      total_percent_cd3plus_foxp3plus_i > 1       ~ "high"
    ), cd3plus_foxp3plus_total_cat.i = factor(cd3plus_foxp3plus_total_cat.i, levels = c("low","high"))) %>%
    mutate(cd11b_total_cat.i = case_when(
      total_percent_cd11b_i == 0      ~ "Absence",
      total_percent_cd11b_i > 0       ~ "Presence"
    ), cd11b_total_cat.i = factor(cd11b_total_cat.i, levels = c("Absence","Presence"))) %>%
    mutate(cd11bplus_cd15plus_total_cat.i = case_when(
      total_percent_cd11bplus_cd15plus_i == 0      ~ "Absence",
      total_percent_cd11bplus_cd15plus_i > 0       ~ "Presence"
    ), cd11bplus_cd15plus_total_cat.i = factor(cd11bplus_cd15plus_total_cat.i, levels = c("Absence","Presence"))) %>% 
    
    mutate(cd3_tumor_cat.p = case_when(
      tumor_percent_cd3_p <= 1      ~ "low",
      tumor_percent_cd3_p > 1       ~ "high"
    ), cd3_tumor_cat.p = factor(cd3_tumor_cat.p, levels = c("low","high"))) %>%
    mutate(cd3plus_cd8plus_tumor_cat.p = case_when(
      tumor_percent_cd3plus_cd8plus_p <= 1      ~ "low",
      tumor_percent_cd3plus_cd8plus_p > 1       ~ "high"
    ), cd3plus_cd8plus_tumor_cat.p = factor(cd3plus_cd8plus_tumor_cat.p, levels = c("low","high"))) %>%
    mutate(cd3plus_foxp3plus_tumor_cat.p = case_when(
      tumor_percent_cd3plus_foxp3plus_p <= 1      ~ "low",
      tumor_percent_cd3plus_foxp3plus_p > 1       ~ "high"
    ), cd3plus_foxp3plus_tumor_cat.p = factor(cd3plus_foxp3plus_tumor_cat.p, levels = c("low","high"))) %>%
    mutate(cd11b_tumor_cat.p = case_when(
      tumor_percent_cd11b_p == 0      ~ "Absence",
      tumor_percent_cd11b_p > 0       ~ "Presence"
    ), cd11b_tumor_cat.p = factor(cd11b_tumor_cat.p, levels = c("Absence","Presence"))) %>%
    mutate(cd11bplus_cd15plus_tumor_cat.p = case_when(
      tumor_percent_cd11bplus_cd15plus_p == 0      ~ "Absence",
      tumor_percent_cd11bplus_cd15plus_p > 0       ~ "Presence"
    ), cd11bplus_cd15plus_tumor_cat.p = factor(cd11bplus_cd15plus_tumor_cat.p, levels = c("Absence","Presence"))) %>% 
    mutate(cd3_stroma_cat.p = case_when(
      stroma_percent_cd3_p <= 1      ~ "low",
      stroma_percent_cd3_p > 1       ~ "high"
    ), cd3_stroma_cat.p = factor(cd3_stroma_cat.p, levels = c("low","high"))) %>%
    mutate(cd3plus_cd8plus_stroma_cat.p = case_when(
      stroma_percent_cd3plus_cd8plus_p <= 1      ~ "low",
      stroma_percent_cd3plus_cd8plus_p > 1       ~ "high"
    ), cd3plus_cd8plus_stroma_cat.p = factor(cd3plus_cd8plus_stroma_cat.p, levels = c("low","high"))) %>%
    mutate(cd3plus_foxp3plus_stroma_cat.p = case_when(
      stroma_percent_cd3plus_foxp3plus_p <= 1      ~ "low",
      stroma_percent_cd3plus_foxp3plus_p > 1       ~ "high"
    ), cd3plus_foxp3plus_stroma_cat.p = factor(cd3plus_foxp3plus_stroma_cat.p, levels = c("low","high"))) %>%
    mutate(cd11b_stroma_cat.p = case_when(
      stroma_percent_cd11b_p == 0      ~ "Absence",
      stroma_percent_cd11b_p > 0       ~ "Presence"
    ), cd11b_stroma_cat.p = factor(cd11b_stroma_cat.p, levels = c("Absence","Presence"))) %>%
    mutate(cd11bplus_cd15plus_stroma_cat.p = case_when(
      stroma_percent_cd11bplus_cd15plus_p == 0      ~ "Absence",
      stroma_percent_cd11bplus_cd15plus_p > 0       ~ "Presence"
    ), cd11bplus_cd15plus_stroma_cat.p = factor(cd11bplus_cd15plus_stroma_cat.p, levels = c("Absence","Presence"))) %>% 
    mutate(cd3_total_cat.p = case_when(
      total_percent_cd3_p <= 1      ~ "low",
      total_percent_cd3_p > 1       ~ "high"
    ), cd3_total_cat.p = factor(cd3_total_cat.p, levels = c("low","high"))) %>%
    mutate(cd3plus_cd8plus_total_cat.p = case_when(
      total_percent_cd3plus_cd8plus_p <= 1      ~ "low",
      total_percent_cd3plus_cd8plus_p > 1       ~ "high"
    ), cd3plus_cd8plus_total_cat.p = factor(cd3plus_cd8plus_total_cat.p, levels = c("low","high"))) %>%
    mutate(cd3plus_foxp3plus_total_cat.p = case_when(
      total_percent_cd3plus_foxp3plus_p <= 1      ~ "low",
      total_percent_cd3plus_foxp3plus_p > 1       ~ "high"
    ), cd3plus_foxp3plus_total_cat.p = factor(cd3plus_foxp3plus_total_cat.p, levels = c("low","high"))) %>%
    mutate(cd11b_total_cat.p = case_when(
      total_percent_cd11b_p == 0      ~ "Absence",
      total_percent_cd11b_p  > 0       ~ "Presence"
    ), cd11b_total_cat.p = factor(cd11b_total_cat.p, levels = c("Absence","Presence"))) %>%
    mutate(cd11bplus_cd15plus_total_cat.p = case_when(
      total_percent_cd11bplus_cd15plus_p == 0      ~ "Absence",
      total_percent_cd11bplus_cd15plus_p > 0       ~ "Presence"
    ), cd11bplus_cd15plus_total_cat.p = factor(cd11bplus_cd15plus_total_cat.p, levels = c("Absence","Presence")))

ROI_global <- ROI_global %>% 
  mutate(cd3_tumor_cat = case_when(
    tumor_percent_cd3 <= 1      ~ "low",
    tumor_percent_cd3 > 1       ~ "high"
  ), cd3_tumor_cat = factor(cd3_tumor_cat, levels = c("low","high"))) %>%
  mutate(cd3plus_cd8plus_tumor_cat = case_when(
    tumor_percent_cd3plus_cd8plus <= 1      ~ "low",
    tumor_percent_cd3plus_cd8plus > 1       ~ "high"
  ), cd3plus_cd8plus_tumor_cat = factor(cd3plus_cd8plus_tumor_cat, levels = c("low","high"))) %>%
  mutate(cd3plus_foxp3plus_tumor_cat = case_when(
    tumor_percent_cd3plus_foxp3plus <= 1      ~ "low",
    tumor_percent_cd3plus_foxp3plus > 1       ~ "high"
  ), cd3plus_foxp3plus_tumor_cat = factor(cd3plus_foxp3plus_tumor_cat, levels = c("low","high"))) %>%
  mutate(cd11b_tumor_cat = case_when(
    tumor_percent_cd11b == 0      ~ "Absence",
    tumor_percent_cd11b > 0       ~ "Presence"
  ), cd11b_tumor_cat = factor(cd11b_tumor_cat, levels = c("Absence","Presence"))) %>%
  mutate(cd11bplus_cd15plus_tumor_cat = case_when(
    tumor_percent_cd11bplus_cd15plus == 0      ~ "Absence",
    tumor_percent_cd11bplus_cd15plus > 0       ~ "Presence"
  ), cd11bplus_cd15plus_tumor_cat = factor(cd11bplus_cd15plus_tumor_cat, levels = c("Absence","Presence"))) %>% 
  mutate(cd3_stroma_cat = case_when(
    stroma_percent_cd3 <= 1      ~ "low",
    stroma_percent_cd3 > 1       ~ "high"
  ), cd3_stroma_cat = factor(cd3_stroma_cat, levels = c("low","high"))) %>%
  mutate(cd3plus_cd8plus_stroma_cat = case_when(
    stroma_percent_cd3plus_cd8plus <= 1      ~ "low",
    stroma_percent_cd3plus_cd8plus > 1       ~ "high"
  ), cd3plus_cd8plus_stroma_cat = factor(cd3plus_cd8plus_stroma_cat, levels = c("low","high"))) %>%
  mutate(cd3plus_foxp3plus_stroma_cat = case_when(
    stroma_percent_cd3plus_foxp3plus <= 1      ~ "low",
    stroma_percent_cd3plus_foxp3plus > 1       ~ "high"
  ), cd3plus_foxp3plus_stroma_cat = factor(cd3plus_foxp3plus_stroma_cat, levels = c("low","high"))) %>%
  mutate(cd11b_stroma_cat = case_when(
    stroma_percent_cd11b == 0      ~ "Absence",
    stroma_percent_cd11b > 0       ~ "Presence"
  ), cd11b_stroma_cat = factor(cd11b_stroma_cat, levels = c("Absence","Presence"))) %>%
  mutate(cd11bplus_cd15plus_stroma_cat = case_when(
    stroma_percent_cd11bplus_cd15plus == 0      ~ "Absence",
    stroma_percent_cd11bplus_cd15plus > 0       ~ "Presence"
  ), cd11bplus_cd15plus_stroma_cat = factor(cd11bplus_cd15plus_stroma_cat, levels = c("Absence","Presence"))) %>% 
  mutate(cd3_total_cat = case_when(
    total_percent_cd3 <= 1      ~ "low",
    total_percent_cd3 > 1       ~ "high"
  ), cd3_total_cat = factor(cd3_total_cat, levels = c("low","high"))) %>%
  mutate(cd3plus_cd8plus_total_cat = case_when(
    total_percent_cd3plus_cd8plus <= 1      ~ "low",
    total_percent_cd3plus_cd8plus > 1       ~ "high"
  ), cd3plus_cd8plus_total_cat = factor(cd3plus_cd8plus_total_cat, levels = c("low","high"))) %>%
  mutate(cd3plus_foxp3plus_total_cat = case_when(
    total_percent_cd3plus_foxp3plus <= 1      ~ "low",
    total_percent_cd3plus_foxp3plus > 1       ~ "high"
  ), cd3plus_foxp3plus_total_cat = factor(cd3plus_foxp3plus_total_cat, levels = c("low","high"))) %>%
  mutate(cd11b_total_cat = case_when(
    total_percent_cd11b == 0      ~ "Absence",
    total_percent_cd11b > 0       ~ "Presence"
  ), cd11b_total_cat = factor(cd11b_total_cat, levels = c("Absence","Presence"))) %>%
  mutate(cd11bplus_cd15plus_total_cat = case_when(
    total_percent_cd11bplus_cd15plus == 0      ~ "Absence",
    total_percent_cd11bplus_cd15plus > 0       ~ "Presence"
  ), cd11bplus_cd15plus_total_cat = factor(cd11bplus_cd15plus_total_cat, levels = c("Absence","Presence")))



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
saveRDS(markers_ROI, file = "cat_summarized_markers_clinical_ROI_07132023.rds")
ROI_global <- right_join(case_ctrl_data, ROI_global,
                         by = "suid")
saveRDS(ROI_global, file = "cat_long_format_markers_clinical_ROI_07132023.rds")

# END create cat for ROI

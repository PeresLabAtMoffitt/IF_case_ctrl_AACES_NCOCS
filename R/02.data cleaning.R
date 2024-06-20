# Import library

library(tidyverse)

############################################################################## I III Load data----
# Run _drake for AACES1 batch1 and NCOCS data
# Run "012.Load data AACES1batch2.R" for AACES1 batch2 data
# Run "013.Load data exhaustion marker_June2024.R" for AACES2 batch1 data

############################################################################## II ### Cleaning clinical data----
cases_match <- cases_match %>% mutate(suid = as.character(suid))

case_ctrl_data <- bind_rows(aaces_clinical, ncocs_clinical) %>% 
  # remove variable not for this data
  select(suid, everything(), -casematch) %>%
  mutate(suid = factor(suid)) %>%
  # Add matched paired id AA-White
  full_join(cases_match, .,
            by= "suid") %>% 
  mutate(casecon = case_when(
    casecon == 1                                       ~ "Case",
    casecon == 2                                       ~ "Control"
  )) %>% 
  # Rename older vital status but keep in data
  rename(vitalstatus_2022 = vitalstatus,
         timelastfu_2022 = timelastfu) %>% 
  # Recode older vital status for NCOCS cases
  mutate(vitalstatus_2022 = na_if(vitalstatus_2022, 88)) %>% 
  mutate(vitalstatus_2022 = case_when(
    vitalstatus_2022 == 1         ~ 0,
    vitalstatus_2022 == 2         ~ 1
  )) %>% 
  # Calculate follow up time as time from interview for NCOCS cases
  mutate(timelastfu_2022 = na_if(timelastfu_2022, 8888)) %>% 
  mutate(os_time_days_NCOCS = timelastfu_2022 - timeint) %>%
  # Create survival var
  left_join(., survival_time_Nov2023, by = "suid") %>% 
  mutate(os_event = coalesce(vital_status_nov2023, vitalstatus_2022)) %>%
  mutate(vitalstatus = case_when(
    os_event == 0                                      ~ "Alive",
    os_event == 1                                      ~ "Deceased",
    TRUE                                               ~ NA_character_
  )) %>% 
  mutate(os_time_days = coalesce(days_int_to_event, os_time_days_NCOCS)) %>% 
  mutate(os_time_months = round(os_time_days/30.417, digit=1)) %>% 
  # Recode other variables
  mutate(cancersite = case_when(
    cancersite == 1                                    ~ "Ovarian",
    cancersite == 2                                    ~ "Tubal",
    cancersite == 3                                    ~ "Peritoneal",
    cancersite == 4                                    ~ "Ovarian or tubal, can't distinguish",
    cancersite == 5                                    ~ "Ovarian, tubal or peritoneal, can't distinguish",
    TRUE                                               ~ NA_character_
  )) %>% 
  mutate_at(c("morphology", "hysteryear", "oophoryear", "tubeligyear",
              "anyfhdur", "eonlydur", "epdur"), 
            ~ case_when(
              . %in% c("8888","9998", "9999")          ~ NA_real_,
              TRUE                                     ~ as.numeric(.)
            )) %>% 
  mutate_at(c("refage", "pregmos", "agefbirth", "agelbirth", 
              "height", "wt_recent", "wt_YA", "wtgain", "BMI_recent", "BMI_YA",
              "hysterage", "oophorage", "tubeligage", "ocdur", "ocstart", "ocstop",
              "breastfedtimes", "breastfeddur", "menarch_age", "menopause_age", "talcgenfreq", 
              "talcgendur", "talcnongenfreq", "talcnongendur", "talcagegen", "talcagenongen",
              "endomet_age", "fibroids_age", "pid_age", "pcos_age", "ovcyst_age", "anyfhstart",
              "anyfhstop", "eonlystart", "eonlystop", "epstart", "epstop", "smokstart", "smokstop",
              "cigday", "smokyrs", "packyrs", "diabage", "hrtdisage", "hbpage", "hcholage", "osteoage",
              "thyrdage", "prvcanage", "prbreastage", "prcolage", "prcervage", "prlungage", 
              "prmelage", "prutage"), 
            ~ case_when(
              . %in% c("888","998", "999")             ~ NA_real_,
              TRUE                                     ~ as.numeric(.)
            )) %>% 
  mutate_at(c("menopause_age"), 
            ~ case_when(
              . %in% c("777")                          ~ NA_real_,
              TRUE                                     ~ as.numeric(.)
            )) %>% 
  mutate_at(c("pregnum", "fullpregnum", "numsis", "numbro", "numsibs", "strenpa"), 
            ~ case_when(
              . %in% c("88","98", "99")                ~ NA_real_,
              TRUE                                     ~ as.numeric(.)
            )) %>% 
  mutate(histology = case_when(
    histology == 1                                     ~ "Serous",
    histology == 2                                     ~ "Endometrioid",
    histology == 3                                     ~ "Clear cell",
    histology == 4                                     ~ "Mucinous",
    histology == 5                                     ~ "Carcinosarcoma",
    histology == 6                                     ~ "Carcinoma, NOS",
    histology == 7                                     ~ "Other specified epithelial ovarian cancer \n(e.g. Malignant Brenner, mixed)",
    histology == 8                                     ~ "Epithelial, NOS",
    histology == 9                                     ~ "Synchronous",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(behavior = case_when(
    behavior == 1                                      ~ "Borderline",
    behavior == 2                                      ~ "Invasive",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(stage_cat = case_when(
    stage == 1 |
    stage == 2                                         ~ "Early",
    stage == 3                                         ~ "Late",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(stage = case_when(
    stage == 1                                         ~ "Localized",
    stage == 2                                         ~ "Regional",
    stage == 3                                         ~ "Distant",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(grade = case_when(
    grade == 1                                         ~ "well differentiated",
    grade == 2                                         ~ "moderately differentiated",
    grade == 3                                         ~ "poorly differentiated",
    grade == 4                                         ~ "undifferentiated",
    TRUE                                               ~ NA_character_
  )) %>% 
  mutate(histotype1 = case_when(
    histotype == 1                                     ~ "high-grade serous",
    histotype == 2                                     ~ "low-grade serous",
    histotype == 5 |
      histotype == 10                                  ~ "mucinous",
    histotype == 3                                     ~ "endometrioid",
    histotype == 4                                     ~ "clear cell",
    histotype %in% (6:13)                              ~ "other epithelial", # Will not take the 10
    TRUE                                               ~ NA_character_
  )) %>% 
  mutate(histotype2 = case_when(
    histotype == 1                                     ~ "high-grade serous",
    histotype %in% (2:13)                              ~ "non-high-grade serous",
    TRUE                                               ~ NA_character_
  )) %>%
  mutate(histotype_cat = case_when(
    histotype == 1 |
      histotype == 6                                   ~ "Type II",
    histotype == 2 |
    histotype == 5 |
      histotype == 10 |
    histotype == 3 |
    histotype == 4                                     ~ "Type I",
    histotype %in% (7:13)                              ~ "other epithelial", # Will not take the 10
    TRUE                                               ~ NA_character_
  )) %>% 
  mutate(histotype = case_when(
    histotype == 1                                     ~ "high-grade serous",
    histotype == 2                                     ~ "low-grade serous",
    histotype == 3                                     ~ "endometrioid",
    histotype == 4                                     ~ "clear cell",
    histotype == 5                                     ~ "mucinous",
    histotype == 6                                     ~ "carcinosarcoma",
    histotype == 7                                     ~ "other epithelial ovarian cancer \n(e.g. Malignant Brenner, mixed, carcinoma, NOS)",
    histotype == 9                                     ~ "serous borderline",
    histotype == 10                                    ~ "mucinous borderline",
    histotype == 11                                    ~ "other epithelial borderline",
    histotype == 13                                    ~ "synchronous",
    TRUE                                               ~ NA_character_
  )) %>% 
  mutate(race = case_when(
    race == 1                                          ~ "White",
    race == 2                                          ~ "Black",
    race == 3                                          ~ "Biracial",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(hispanic = case_when(
    hispanic == 1                                      ~ "Hispanic",
    hispanic == 2                                      ~ "Non-hispanic",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(birthplace = case_when(
    birthplace == 1                                    ~ "born in United States",
    birthplace == 2                                    ~ "born outside of United States",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(education = case_when(
    education == 1                                      ~ "high school graduate/GED or less",
    education == 2                                      ~ "some college",
    education == 3                                      ~ "college graduate",
    education == 4                                      ~ "graduate/professional school",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(married = case_when(
    married == 1                                        ~ "single/never married",
    married == 2                                        ~ "married/living as married",
    married == 3                                        ~ "divorced/separated",
    married == 4                                        ~ "widowed"
  )) %>% 
  mutate_at(c("pregever", "nulliparous", "oophor", 
              "oophor1yr", "ocever", "breastfedever", 
              "talcever", "talcgen", "talcnongen", "talcpartner", "talc_occ", "brcancer",
              "brcancermom", "brcancersis", "brcancergrandma", "brcanceraunt", "brcancerdau",
              "brcancerdad", "brcancerbro", "brcancerson", "brcancerchildren", "ovcancermom",
              "ovcancersis", "ovcancergrandma", "ovcanceraunt", "ovcancerdaughter", "famhxbr",
              "famhxov", "pcos", "ovcyst", "eonlyever",
              "epever", "eonlyever", "hrtdis", "hbp", 
              "hchol", "osteo", "prvcan", "prbreast",  "prcol", "prcerv", "prlung", "prmel", 
              "prut", "infert", "trypreg1yr", "mdvisit", "fertmed", "infertsurgivf",
              "neoadj_treat", "adj_treat", "dblkstat_treat"), 
            ~ case_when(
              . == 1                                              ~ "yes",
              . == 2                                              ~ "no",
              TRUE                                                ~ NA_character_
            )) %>%
  mutate_at(c("diab", "aspirin", "NSAID", "aceta",
              "hyster", "hyster1yr", "hyster2yr",
              "tubelig", "tubelig1yr", "anyfhever", 
              "endomet", "fibroids", "pid"),
            ~ case_when(
              . == 0                                              ~ "no",
              . == 1                                              ~ "yes",
              TRUE                                                ~ NA_character_
            )) %>% 
  mutate_at(c("diab", "aspirin", "NSAID", "aceta",
              "hyster", "hyster1yr", "hyster2yr",
              "tubelig", "tubelig1yr", "anyfhever", 
              "endomet", "fibroids", "pid"), 
            ~ factor(., levels = c("no", "yes"))) %>% 

  mutate(aspirin_use = factor(aspirin, levels = c("no", "yes"))) %>% 
  mutate(hysterreason = case_when(
    hysterreason == 1                                   ~ "Ovarian/FT/peritoneal cancer diagnosis",
    hysterreason == 2                                   ~ "Any reason not due to ovarian/FT/peritoneal cancer diagnosis",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate_at(c("hystertype", "menopause"), ~ case_when(
    . == 1                                              ~ "Premenopausal",
    . == 2                                              ~ "Postmenopausal",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(oophortype = case_when(
    oophortype == 1                                     ~ "Unilateral",
    oophortype == 2                                     ~ "Bilateral",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(ocindication = case_when(
    ocever == 1                                         ~ "birth control",
    ocever == 2                                         ~ "regulate menstrual periods",
    ocever == 3                                         ~ "acne or skin problems",
    ocever == 4                                         ~ "endometriosis",
    ocever == 5                                         ~ "premenstrual syndrome (PMS)",
    ocever == 6                                         ~ "menopausal symptoms",
    ocever == 7                                         ~ "vaginal bleeding",
    ocever == 8                                         ~ "other",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(menarch_agecat = case_when(
    menarch_agecat == 1                                 ~ "<11 years",
    menarch_agecat == 2                                 ~ "11-12 years",
    menarch_agecat == 3                                 ~ "13-14 years",
    menarch_agecat == 4                                 ~ "15-16 years",
    menarch_agecat == 5                                 ~ "17 or older",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(regperiod_agecat = case_when(
    regperiod_agecat == 1                               ~ "<11 years",
    regperiod_agecat == 2                               ~ "11-12 years",
    regperiod_agecat == 3                               ~ "13-14 years",
    regperiod_agecat == 4                               ~ "15-16 years",
    regperiod_agecat == 5                               ~ "17 or older",
    regperiod_agecat == 6                               ~ "never became regular",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(periodstopreason = case_when(
    periodstopreason == 1                               ~ "stopped naturally",
    periodstopreason == 2                               ~ "stopped due to surgery",
    periodstopreason == 3                               ~ "stopped due to radiation or chemotherapy",
    periodstopreason == 4                               ~ "stopped due to hormonal medications (OC or menopausal hormone therapy)",
    periodstopreason == 5                               ~ "stopped due to other reason",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate_at(c("brcancermom_age", "brcancersis_age", "brcancergrandma_age", 
              "brcancerdau_age", "ovcancermom_age", "ovcancersis_age", 
              "ovcancergrandma_age", "ovcancerdaughter_age"), ~ case_when(
                . == 1                                              ~ "<50 years",
                . == 2                                              ~ "50+ years",
                TRUE                                                ~ NA_character_
              )) %>% 
  mutate(smokever = case_when(
    smokever == 1                                       ~ "ever",
    smokever == 0                                       ~ "never",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(smokcurrent = case_when(
    smokcurrent == 1                                    ~ "current",
    smokcurrent == 2                                    ~ "never",
    smokcurrent == 3                                    ~ "former",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(smoker = factor(smokcurrent, levels = c("never", "former", "current"))) %>% 
  mutate(thyrd = case_when(
    thyrd == 1                                          ~ "yes unspecified",
    thyrd == 2                                          ~ "no",
    thyrd == 3                                          ~ "yes hyperthyroid/overreactive",
    thyrd == 4                                          ~ "yes hypothyroid/underreactive",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(mdvisitrsn = case_when(
    mdvisitrsn == 1                                     ~ "Female problem",
    mdvisitrsn == 2                                     ~ "Partner problem",
    mdvisitrsn == 3                                     ~ "Both female and partner problem",
    mdvisitrsn == 4                                     ~ "No problem was found",
    TRUE                                                ~ NA_character_
  )) %>% 
  mutate(diab = factor(diab, levels = c("no", "yes")))

case_ctrl_data <- case_ctrl_data %>% 
  mutate(refage_cat = case_when(
    refage < 50                      ~ "<50",
    refage >= 50 &
      refage < 60                    ~ "50-59",
    refage >= 60 &
      refage < 70                    ~ "60-69",
    refage >= 70 &
      refage < 80                    ~ "70-79",
    refage >= 80                     ~ "≥80"
  )) %>% 
  mutate(BMI_classification = case_when(
    BMI_recent < 18.5	~ "underweight",
    BMI_recent >=18.5 & BMI_recent <25 ~ "normal",
    BMI_recent >=25.0 & BMI_recent <30 ~ "overweight",
    BMI_recent >=30.0 & BMI_recent <35 ~ "obesity I",
    BMI_recent >=35.0 & BMI_recent <40 ~ "obesity II",
    BMI_recent >= 40.0 ~ "obesity III"
  )) %>% 
  mutate(BMI_recent_grp = case_when(
    BMI_recent < 25                  ~ "<25",
    BMI_recent >= 25 &
      BMI_recent < 30                ~ "25-29",
    BMI_recent >= 30 &
      BMI_recent < 35                ~ "30-35",
    BMI_recent >= 35                 ~ "≥35"
  )) %>% 
  mutate(BMI_YA_grp = case_when(
    BMI_YA < 20                      ~ "<20",
    BMI_YA >= 20 &
      BMI_YA < 25                    ~ "20-24",
    BMI_YA >= 25                     ~ "≥25"
  ))%>% 
  mutate(smokever = factor(smokever, levels = c("never", "ever"))) %>% 
  mutate(smokcurrent = factor(smokcurrent, levels = c("never", "former", "current"))) %>% 
  mutate(packyrs_cat = case_when(
    packyrs <= 15                        ~ "≤ 15",
    packyrs > 15                         ~ "> 15",
    smokcurrent == "never"               ~ "none"
  )) %>%
  mutate(packyrs_cat = factor(packyrs_cat, levels = c("none", "≤ 15", "> 15"))) %>% 
  mutate(education = case_when(
    education == "high school graduate/GED or less"|
      education == "some college"                     ~ "didn't graduated",
    education == "college graduate"|
      education == "graduate/professional"            ~ "graduated"
  )) %>% 
  mutate(education = factor(education, levels = c("graduated", "didn't graduated"))) %>% 
  mutate(anyfhdur_cat = case_when(
    anyfhdur <= 100                      ~ "≤ 100",
    anyfhdur > 100                       ~ "> 100",
    anyfhever == "no"                    ~ "none"
  )) %>%
  mutate(anyfhdur_cat = factor(anyfhdur_cat, levels = c("none", "≤ 100", "> 100"))) %>% 
  mutate(BMI_recent_grp = case_when(
    BMI_recent <25                       ~ "<25 kg/m2",
    BMI_recent >=25 & BMI_recent <= 30   ~ "25-30 kg/m2",
    BMI_recent > 30.0                    ~ ">30 kg/m2"
  )) %>% 
  mutate(BMI_YA_grp = case_when(
    BMI_YA <25                           ~ "<25 kg/m2",
    BMI_YA >=25 & BMI_YA <= 30           ~ "25-30 kg/m2",
    BMI_YA > 30.0                        ~ ">30 kg/m2"
  )) %>% 
  mutate(stage = case_when(
    stage == "Localized"                 ~ "Regional",
    TRUE                                 ~ stage
  )) %>% 
  mutate(married = case_when(
    married == "married/living as married"      ~ "married",
    TRUE                                        ~ "others"
  ))

save_path <- fs::path("","Volumes","Peres_Research", "AACES",
                      "AACES_NCOCS_immune_data", "processed_data")
write_rds(case_ctrl_data, "case_ctrl_data_06142024.rds")
write_rds(case_ctrl_data, paste0(save_path, "/case_ctrl_data_06142024.rds"))
write_csv(case_ctrl_data, paste0(save_path, "/case_ctrl_data_06142024.csv"))

survival_time_Nov2023 <- case_ctrl_data %>% 
  select(suid, timeint, os_event, os_time_days, os_time_months)

write_csv(survival_time_Nov2023, paste0("survival_time_Nov2023_", Sys.Date(),".csv"))
write_csv(survival_time_Nov2023, paste0(save_path, "/survival_time_Nov2023_", Sys.Date(),".csv"))

rm(aaces_clinical, ncocs_clinical, cases_match, 
   survival_time_Nov2023)


############################################################################## III ### Join New ROIs data----
ROI_global_2022jan <- 
  full_join(ROI_tumor_2022jan, ROI_stroma_2022jan,
            by = "image_tag") %>% 
  full_join(., ROI_total_2022jan,
            by = "image_tag") %>% 
  mutate(image_tag = str_replace(image_tag, "16-", "16"),
         suid = str_match(image_tag,
                          "(L.Peres_P1_OV|L.Peres_P1_)([:digit:]*)")[,3],
         .before = 1) %>%
  select(image_tag, suid, roi_region = annotation, everything()) %>%
  arrange(suid) %>% 
  # Remove patients who were excluded from the study probably after path review
  filter(!str_detect(suid, paste0(ROIcases_remove$Subject_IDs, collapse = "|"))) %>% 
  mutate_at(("roi_region"), ~ case_when(
    roi_region == "P"                     ~ "Peripheral",
    roi_region == "I"                     ~ "Intratumoral",
    roi_region == "S"                     ~ "Stromal")
  ) %>%
  mutate(aaces_phase = "AACES_1", .after = image_tag) %>%
  mutate(panel = "1 (Abundance)", .after = aaces_phase) %>% 
  mutate(panel1_batch = "roi 2022", .after = panel) %>% 
  mutate(slide_type = "ROI", .after = panel1_batch) %>% 
  mutate(data_version = "AACES_1_v2", .after = slide_type) %>% 
  # rename(tumor_area_analyzed_μm2 = tumor_area_analyzed_mm2,
  #        stroma_area_analyzed_μm2 = stroma_area_analyzed_mm2,
  #        total_area_analyzed_μm2 = total_area_analyzed_mm2,
  #        tumor_area_analyzed_mm2 = tumor_area_analyzed_mm2_2,
  #        stroma_area_analyzed_mm2 = stroma_area_analyzed_mm2_2,
  #        total_area_analyzed_mm2 = total_area_analyzed_mm2_2
  #        )

rm(ROI_stroma_2022jan, ROI_total_2022jan, ROI_tumor_2022jan)



############################################################################## IV ### Join R00 TMA / ROIs data----
# uid <- 
TMA_tumor <- TMA_tumor %>% 
  # Plus remove TMA with no IDs = controls images
  filter(!is.na(suid))
TMA_stroma <- TMA_stroma %>% 
  filter(!is.na(suid))
TMA_total <- TMA_total %>% 
  filter(!is.na(suid)) %>% 
  # Add "total" to each marker variable
  rename_at(vars(starts_with("fox") | 
                   starts_with("cd") |
                   starts_with("percent_fox") | 
                   starts_with("percent_cd") |
                   starts_with("area")
  ), ~ paste0("total_", .)) 

TMA_global <- 
  full_join(TMA_tumor, TMA_stroma %>% select(-suid),
            by = "image_tag") %>% 
  full_join(., TMA_total %>% select(-suid),
            by = "image_tag") %>% 
  mutate(suid = as.character(suid)) %>% 
  # Remove the TMA IDs of excluded patient from the study
  filter(!str_detect(
    paste(unique(TMAcases_remove$Subject_IDs), collapse = "|"), 
    suid)) %>% 
  mutate(aaces_phase = "AACES_1", .after = image_tag) %>%
  mutate(panel = "1 (Abundance)", .after = aaces_phase) %>% 
  mutate(panel1_batch = "tma 2017-18", .after = panel) %>% 
  mutate(slide_type = "TMA", .after = panel1_batch) %>% 
  mutate(data_version = "AACES_1_NCOCS_v1", .after = slide_type)
  
  rename_at(vars(ends_with("_")
  ), ~ str_replace(., "_$", "")) %>% 
  select(image_tag, suid, everything(),
         -c(starts_with("tma"), starts_with("other_tissue"), 
            starts_with("block_id")))

rm(TMA_tumor, TMA_total, TMA_stroma)


# 3.2.Create suid for ROIs----
ROI_tumor$suid <- str_match(ROI_tumor$image_tag, 
                            "(Peres_P1_AACES.|Peres_P1_AACEES.|Peres_P1_OV|Peres_P1.|)([:digit:]*)")[,3]
ROI_stroma$suid <- str_match(ROI_stroma$image_tag, 
                             "(Peres_P1_AACES.|Peres_P1_AACEES.|Peres_P1_OV|Peres_P1.|)([:digit:]*)")[,3]
ROI_total$suid <- str_match(ROI_total$image_tag, 
                            "(Peres_P1_AACES.|Peres_P1_AACEES.|Peres_P1_OV|Peres_P1.|)([:digit:]*)")[,3]

# 3.3.Merging stroma and tumor for TMAs and ROIs
ROI_global_2021 <- 
  full_join(ROI_tumor, ROI_stroma %>% select(-c("intratumoral_i_vs_peripheral_p_", "suid")),
            by = "image_tag") %>% 
  full_join(., ROI_total %>% select(-c("intratumoral_i_vs_peripheral_p_", "suid")),
            by = "image_tag") %>% 
  # Remove Ctrl image
  filter(!str_detect(image_tag, "Ctrl")) %>%
  # Remove images that were not stained properly
  filter(!image_tag %in% c(ROIimage_remove$remove_images_that_were_not_stained_properly)) %>% 
  # homogenize / AACES2
  rename(roi_region = intratumoral_i_vs_peripheral_p_) %>% 
  mutate(roi_region = case_when(
    roi_region == "p" ~ "Peripheral",
    roi_region == "i" ~ "Intratumoral")
  ) %>% 
  mutate(suid = as.character(suid)) %>% 
  
  mutate(aaces_phase = "AACES_1", .after = image_tag) %>%
  mutate(panel = "1 (Abundance)", .after = aaces_phase) %>% 
  mutate(panel1_batch = "roi 2021", .after = panel) %>% 
  mutate(slide_type = "ROI", .after = panel1_batch) %>% 
  mutate(data_version = "AACES_1_NCOCS_v1", .after = slide_type)
  select(image_tag, suid, everything()) %>% 
    
  rename_at(vars(ends_with("_")
  ), ~ str_replace(., "_$", "")) %>% 
    # Add "total" to each marker variable
    rename_at(vars(starts_with("fox") | 
                   starts_with("cd") |
                   starts_with("percent_fox") | 
                   starts_with("percent_cd") |
                   starts_with("area")
                 ), ~ paste0("total_", .)) 

rm(ROI_total, ROI_tumor, ROI_stroma, ROIimage_remove)


############################################################################## VI ### Clean TMA/WTS exhaustion data 2024----
# TMA----
TMA_manifest <- bind_rows(TMA2017_manifest %>% 
                            mutate(panel2_batch = "tma 2017"),
                          TMA2018_manifest %>% 
                            mutate(panel2_batch = "tma 2018"),
                          TMA2024_manifest %>% 
                            mutate(panel2_batch = "tma 2024"))

TMA_exhaustion_aaces1 <- bind_rows(
  full_join(TMA2017_tumor_2024june, TMA2017_stroma_2024june,
            by = "image_tag") %>% 
  full_join(., TMA2017_total_2024june,
            by = "image_tag") %>% 
    mutate(panel2_batch = "tma 2017", .after = 1),
  full_join(TMA2018_tumor_2024june, TMA2018_stroma_2024june,
            by = "image_tag") %>% 
    full_join(., TMA2018_total_2024june,
              by = "image_tag") %>% 
    mutate(panel2_batch = "tma 2018", .after = 1)
  ) %>% 
  mutate(aaces_phase = "AACES_1", .after = 1) %>% 
  # Add "total" to each marker variable
  rename_at(vars(starts_with("tim") | 
                   starts_with("pd") |
                   starts_with("cd") |
                   starts_with("percent_tim") | 
                   starts_with("percent_pd") | 
                   starts_with("percent_cd") |
                   starts_with("area")
  ), ~ paste0("total_", .)) 

TMA_exhaustion_aaces2 <- 
  # Join
  full_join(TMAaaces2_tumor_2024june, TMAaaces2_stroma_2024june,
            by = "image_tag") %>% 
    full_join(., TMAaaces2_total_2024june,
              by = "image_tag") %>% 
  # Add "total" to each marker variable
  rename_at(vars(starts_with("tim") | 
                   starts_with("pd") |
                   starts_with("cd") |
                   starts_with("percent_tim") | 
                   starts_with("percent_pd") | 
                   starts_with("percent_cd") |
                   starts_with("area")
  ), ~ paste0("total_", .)) 
  
  
TMA_exhaustion2024 <- bind_rows(TMA_exhaustion_aaces1, 
                                TMA_exhaustion_aaces2) %>% 
  mutate(aaces_phase = "AACES_2", .after = 1) %>% 
  mutate(panel = "2 (Exhaustion)", .after = 3) %>% 
  mutate(panel2_batch = "tma 2024", .after = panel) %>% 
  mutate(slide_type = "TMA", .after = panel2_batch) %>% 
  mutate(data_version = "AACES_2_v1", .after = slide_type)

TMA_exhaustion2024 <- TMA_exhaustion2024 %>% 
  # Extract suid
  mutate(core = str_match(image_tag, 
                          ".*(\\[.*\\]).tif")[,2], .after = image_tag) %>% 
  full_join(., TMA_manifest, by = c("core", "panel2_batch")) %>% 
  select(suid, image_tag, core, everything()) %>% 
  mutate(suid = as.character(suid)) %>% 
  # Remove ctrl (no suid) and TMA that were remove from core (no image_tag)
  filter(!is.na(image_tag) & !is.na(suid)) %>% 
  # Remove the TMA IDs of excluded patient from the study
  filter(!str_detect(
    paste(unique(TMAcases_remove$Subject_IDs), collapse = "|"),
    suid))

# Note :
# 2017 core 4A was removed as is broken pieces of core 
# 2017 cores row 1 and 28 are controls
# 2017 - 4 additional cores [5,G], [5,J], [3,O], [11,Q] were removed compared to the manifest
# because of no data
# 2018 cores row 1 and 20 are controls
# 2024 cores row 1 and (9?) are controls
write_rds(TMA_exhaustion2024, "raw_data_all_TMAs_exhaustion_panel_2024.rds")

rm(TMA2017_tumor_2024june, TMA2017_stroma_2024june, 
   TMA2017_total_2024june, TMA2018_tumor_2024june, 
   TMA2018_stroma_2024june, TMA2018_total_2024june,
   TMAaaces2_tumor_2024june, TMAaaces2_stroma_2024june,
   TMAaaces2_total_2024june,
   TMA_exhaustion_aaces1, TMA_exhaustion_aaces2,
   TMA_manifest, TMA2017_manifest, TMA2018_manifest,
   TMAcases_remove)


# WTS----
# We have 6 duplicate scanned slides between set but with different counts
# Looks like the slides were scanned twice
# Using set 1 for now
# Peres_P3_MOTIF_110393
# Peres_P3_MOTIF_110388
# Peres_P3_MOTIF_190047
# Peres_P3_MOTIF_110393
# Peres_P3_MOTIF_110388
# Peres_P3_MOTIF_190047

WTS_exhaustion_aaces2 <- 
  # Join
  full_join(WTS_tumor_2024june, WTS_stroma_2024june,
            by = c("image_tag", "set")) %>% 
    full_join(., WTS_total_2024june,
              by = c("image_tag", "set")) %>% 
  # Add "total" to each marker variable
  rename_at(vars(starts_with("tim") | 
                   starts_with("pd") |
                   starts_with("cd") |
                   starts_with("percent_tim") | 
                   starts_with("percent_pd") | 
                   starts_with("percent_cd") |
                   starts_with("area")
  ), ~ paste0("total_", .)) %>% 
  mutate(aaces_phase = "AACES_2", .after = image_tag) %>% 
  mutate(panel = "2 (Exhaustion)", .after = aaces_phase) %>% 
  mutate(panel2_batch = "wts 2024", .after = panel) %>% 
  mutate(slide_type = "WTS", .after = panel2_batch) %>% 
  mutate(data_version = "AACES_2_v1", .after = slide_type)

WTS_exhaustion_aaces2 <- WTS_exhaustion_aaces2 %>% 
  # extract suid
  mutate(suid = str_replace(image_tag, "16-", "16"),
         suid = str_match(suid, 
                          "(Peres_P3_MOTIF_AACES |Peres_P3_MOTIF_ACCESS |Peres_P3_MOTIF_OV|Peres_P3_MOTIF_)(\\d+).*.tif")[,3],
         .after = image_tag) %>% 
  select(suid, image_tag, everything()) %>% 
  # Remove patients who were excluded from the study probably after path review
  filter(!str_detect(
    suid, paste0(ROIcases_remove$Subject_IDs, collapse = "|")))


############################################################################## VII ### Bind ROIs data----
markers_AACES_NCOCS <- bind_rows(ROI_global_2021, TMA_global, ROI_global_2022jan) %>% 
  `colnames<-`(str_remove_all(colnames(.), "_positive_cells|_cells|_opal_..._positive_cells")) %>% 
  select(data_version, image_tag, suid, roi_region, slide_type, everything()) %>% 
  rename(total_total = total)

saveRDS(markers_AACES_NCOCS, file = "markers_AACES_NCOCS_batch1_2_07132023.rds")
write_csv(markers_AACES_NCOCS, "markers_AACES_NCOCS_batch1_2_07132023.csv")

complete_AACES_NCOCS_global <- left_join(markers_AACES_NCOCS,
                                    case_ctrl_data,
                                    by= "suid") %>% 
  select(image_tag, suid, roi_region, slide_type, site, data_version, everything())

complete_AACES_NCOCS_global %>% 
  distinct(suid, .keep_all = TRUE) %>% 
  select(refage, race, histotype_cat, stage_cat, vitalstatus, timeint_fu, 
         site) %>% 
  tbl_summary(by = site)

complete_AACES_NCOCS_global %>% 
  distinct(suid, slide_type, .keep_all = TRUE) %>% 
  select(refage, race, histotype_cat, stage_cat, vitalstatus, timeint_fu, 
         site, slide_type) %>% 
  tbl_strata(
    strata = site,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(by = slide_type, missing = "no") %>%
      add_n(),
    .header = "**{strata}**, N = {n}"
  )

complete_AACES_NCOCS_global %>% 
  select(refage,
         site, slide_type) %>% 
  tbl_strata(
    strata = site,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(by = slide_type, missing = "no") %>%
      add_n(),
    .header = "**{strata}**, N = {n}"
  )


saveRDS(complete_AACES_NCOCS_global, file = "clinical_markers_AACES_NCOCS_batch1_2_07132023.rds")
write_csv(complete_AACES_NCOCS_global, "clinical_markers_AACES_NCOCS_batch1_2_07132023.csv")

# End cleaning
# YOU ARE DONE CREATING THE MAIN DATA USED FOR ALL AACES BASED PROJECTS!


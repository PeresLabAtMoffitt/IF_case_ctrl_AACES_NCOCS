# Import library

library(tidyverse)


############################################################################## II ### Cleaning clinical data----
case_ctrl_data <- bind_rows(aaces_clinical, ncocs_clinical) %>% 
  # remove variable not for thidata
  select(suid, everything(), -casematch) %>% 
  mutate(suid = factor(suid)) %>% 
  mutate(casecon = case_when(
    casecon == 1                                       ~ "Case",
    casecon == 2                                       ~ "Control"
  )) %>% 
  mutate(vitalstatus = case_when(
    vitalstatus == 1                                   ~ "Alive",
    vitalstatus == 2                                   ~ "Deceased",
    TRUE                                               ~ NA_character_
  )) %>% 
  mutate(os_event = case_when(
    vitalstatus == "Alive"                             ~ 0,
    vitalstatus == "Deceased"                          ~ 1,
    TRUE                                               ~ NA_real_
  )) %>% 
  mutate(cancersite = case_when(
    cancersite == 1                                    ~ "Ovarian",
    cancersite == 2                                    ~ "Tubal",
    cancersite == 3                                    ~ "Peritoneal",
    cancersite == 4                                    ~ "Ovarian or tubal, can't distinguish",
    cancersite == 5                                    ~ "Ovarian, tubal or peritoneal, can't distinguish",
    TRUE                                               ~ NA_character_
  )) %>% 
  mutate_at(c("timelastfu", "morphology", "hysteryear", "oophoryear", "tubeligyear",
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
  mutate(BMI_classification = case_when(
    BMI_recent < 18.5	~ "underweight",
    BMI_recent >=18.5 & BMI_recent <25 ~ "normal",
    BMI_recent >=25.0 & BMI_recent <30 ~ "overweight",
    BMI_recent >=30.0 & BMI_recent <35 ~ "obesity I",
    BMI_recent >=35.0 & BMI_recent <40 ~ "obesity II",
    BMI_recent >= 40.0 ~ "obesity III"
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
  mutate_at(c("pregever", "nulliparous", "hyster", "hyster1yr", "hyster2yr", "oophor", 
              "oophor1yr", "tubelig", "tubelig1yr", "ocever", "breastfedever", 
              "talcever", "talcgen", "talcnongen", "talcpartner", "talc_occ", "brcancer",
              "brcancermom", "brcancersis", "brcancergrandma", "brcanceraunt", "brcancerdau",
              "brcancerdad", "brcancerbro", "brcancerson", "brcancerchildren", "ovcancermom",
              "ovcancersis", "ovcancergrandma", "ovcanceraunt", "ovcancerdaughter", "famhxbr",
              "famhxov", "endomet", "fibroids", "pid", "pcos", "ovcyst", "anyfhever", "eonlyever",
              "epever", "eonlyever", "aspirin", "NSAID", "aceta", "diab", "hrtdis", "hbp", 
              "hchol", "osteo", "prvcan", "prbreast",  "prcol", "prcerv", "prlung", "prmel", 
              "prut", "infert", "trypreg1yr", "mdvisit", "fertmed", "infertsurgivf"), 
            ~ case_when(
              . == 1                                              ~ "yes",
              . == 2                                              ~ "no",
              TRUE                                                ~ NA_character_
            )) %>%
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
    smokever == 2                                       ~ "never",
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
  mutate(diab = factor(diab, levels = c("no", "yes"))) %>% 
  # Calculate follow up time as time from interview to follow up
  mutate(os_time = timelastfu - timeint) %>% 
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
  ))

saveRDS(case_ctrl_data, file = "case_ctrl_data.rds")


############################################################################## III ### Join ROIs data----
ROI_global <- 
  full_join(ROI_tumor, ROI_stroma,
            by = "image_tag") %>% 
  full_join(., ROI_total,
            by = c("image_tag" = "total_image_tag")) %>% 
  mutate(image_tag = str_replace(image_tag, "16-", "16"),
         suid = str_match(image_tag,
                          "(L.Peres_P1_OV|L.Peres_P1_)([:digit:]*)")[,3],
         .before = 1) %>%
  `colnames<-`(str_remove(colnames(.), "positive_")) %>% 
  select(image_tag, suid, annotation = total_annotation,
         total_cells = total_total_cells, everything()) %>%
  arrange(suid) %>% 
  # mutate(suid = as.character(suid)) %>% 
  # Calculate percent of tumor cell
  mutate(percent_tumor = round((tumor_total_cells / total_cells)*100, 2) 
  ) %>% 
  # Calculate percent of stromal cell
  mutate(percent_stroma = round((stroma_total_cells / total_cells)*100, 2) 
  ) %>% 
  # Calculate percent of stromal cell
  mutate(percent_total = round((total_cells / total_cells)*100, 2) 
  ) %>% 
  mutate_at(("annotation"), ~ case_when(
    annotation == "P"                     ~ "Peripheral",
    annotation == "I"                     ~ "Intratumoral",
    annotation == "S"                     ~ "Stromal")
  ) %>% 
  mutate(slide_type = "ROI")


######################################################################################## IV ### Create tertile cat with immune cells----
# Create average by patients by annotation
markers_full <- ROI_global %>% 
  select(-image_tag) %>% 
  group_by(suid, annotation) %>% 
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
              values_from = -c(suid, annotation),
              names_glue = "{.value}.{annotation}_"
  )

for (i in colnames(markers_ROI %>% 
                   select(c(contains("percent_cd3"), 
                            starts_with("percent_cd8"), 
                            contains("percent_cd11b"))))) {
  
  cell_type <-
    str_match(i, "percent_(.*?)(_opal|_cells)")[, 2]
  compartment <-
    str_match(i, "^(.*?)_")[, 2]
  region <-
    str_match(i, "_cells.(.*?)$")[, 2]

  markers_ROI[[paste0(cell_type, "_", compartment, "_", region, "grp")]] <-  case_when(
    i <= 1      ~ "low",
    i > 1       ~ "high"
  )
}
colnames(markers_ROI)

markers_ROI <- markers_ROI %>% 
  mutate(across(ends_with("_grp"), ~ factor(.)))


######################################################################################## V ### Create immunoscore and cluster----

######################################################################################## VI ### Join data----
markers <- right_join(case_ctrl_data, markers_ROI, ############# WRONG suid for NCOCS
                     by = "suid")
ROI_global <- full_join(case_ctrl_data, ROI_global,
                     by = "suid")


saveRDS(ROI_global, file = "ROI_global.rds") ############# WRONG suid for NCOCS
saveRDS(markers, file = "markers_ROI.rds")


# End cleaning

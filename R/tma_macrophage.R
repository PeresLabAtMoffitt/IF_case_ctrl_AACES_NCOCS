# Import library

library(tidyverse); library(janitor)


############################################################################## I ### Load data----
path <- fs::path("","Volumes","Peres_Research", "K99_R00", "Image analysis data")

# Total
tma_2017_total <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/Tworoger_P4_AACES_2017A.xlsx"),
    sheet = "P4_AACES 2017A") %>% 
  janitor::clean_names()
tma_2018_total <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/Tworoger_P4_AACES_2018A.xlsx"),
    sheet = "P4_AACES 2018A") %>% 
  janitor::clean_names()
macrophage_total <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/AACES_intensity_simple_macrophages_only.xlsx"),
    sheet = "Overall") %>% 
  janitor::clean_names()
# Tumor
tma_2017_tumor <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/Tworoger_P4_AACES_2017A.xlsx"),
    sheet = "Tumor (PCK+)") %>% 
  janitor::clean_names()
tma_2018_tumor <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/Tworoger_P4_AACES_2018A.xlsx"),
    sheet = "Tumor (PCK+)") %>% 
  janitor::clean_names()
macrophage_tumor <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/AACES_intensity_simple_macrophages_only.xlsx"),
    sheet = "Tumor") %>% 
  janitor::clean_names()
# Stroma
tma_2017_stroma <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/Tworoger_P4_AACES_2017A.xlsx"),
    sheet = "Stroma (PCK-)") %>% 
  janitor::clean_names()
tma_2018_stroma <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/Tworoger_P4_AACES_2018A.xlsx"),
    sheet = "Stroma (PCK-)") %>% 
  janitor::clean_names()
macrophage_stroma <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/AACES_intensity_simple_macrophages_only.xlsx"),
    sheet = "Stroma") %>% 
  janitor::clean_names()

tma_map2017 <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/TAM TMA Map.xlsx"),
    sheet = "TMA 2017")
tma_map2018 <- 
  readxl::read_xlsx(
    paste0(path,
           "/TAM panel/TAM TMA Map.xlsx"),
    sheet = "TMA 2018")

############################################################################## I ### Clean data----
# map data
tma_map <- bind_rows(tma_map2017,
                     tma_map2018, .id = "version") %>% 
  mutate(control = case_when(
    !is.na(control)             ~ paste0("ctrl=", control)
  )) %>% 
  unite(suid, suid:control, na.rm = TRUE) %>% 
  mutate(version = case_when(
    version == 1                ~ "2017A_TMA",
    version == 2                ~ "2018A_TMA"
  )) %>% 
  unite(core, core:version, na.rm = TRUE)

rm(tma_map2017, tma_map2018)


# TMA data
total_tma <- bind_rows(tma_2017_total, tma_2018_total)
tumor_tma <- bind_rows(tma_2017_tumor, tma_2018_tumor)
stroma_tma <- bind_rows(tma_2017_stroma, tma_2018_stroma)
rm(tma_2017_total, tma_2018_total,
   tma_2017_tumor, tma_2018_tumor,
   tma_2017_stroma, tma_2018_stroma)

# Merge data
total_data <- 
  inner_join(macrophage_total,
             total_tma,
             by = "image_tag") %>% 
  mutate(core = str_match(image_tag, "\\[(.*)\\]")[,2]) %>% 
  unite(core, c(core,tma), remove = FALSE) %>% 
  full_join(tma_map,
            .,
            by = "core") %>% 
  filter(!str_detect(suid, "ctrl")) %>% 
  drop_na(image_location) %>% 
  rename(area_analyzed_µm2 = area_analyzed_mm2,
         area_analyzed_mm2 = area_analyzed_mm2_2)

tumor_data <- 
  inner_join(macrophage_tumor %>% 
               select(image_location, image_tag, tma,
                      m1_simple_cells, m2_simple_cells, m0_simple_cells,
                      m1_simple_cells_percent, m2_simple_cells_percent, 
                      m0_simple_cells_percent, m1_70p, m2_30p,
                      neutral_3070p, m1_60p, m2_40p, neutral_4060p),
             tumor_tma,
             by = "image_tag") %>% 
  mutate(core = str_match(image_tag, "\\[(.*)\\]")[,2]) %>% 
  unite(core, c(core,tma), remove = FALSE) %>% 
  full_join(tma_map,
            .,
            by = "core") %>% 
  filter(!str_detect(suid, "ctrl")) %>% 
  drop_na(image_location) %>% 
  rename(tumor_area_analyzed_µm2 = tumor_area_analyzed_mm2,
         tumor_area_analyzed_mm2 = tumor_area_analyzed_mm2_2)

stroma_data <- 
  inner_join(macrophage_stroma %>% 
               select(image_location, image_tag, tma,
                      m1_simple_cells, m2_simple_cells, m0_simple_cells,
                      m1_simple_cells_percent, m2_simple_cells_percent, 
                      m0_simple_cells_percent, m1_70p, m2_30p,
                      neutral_3070p, m1_60p, m2_40p, neutral_4060p),
             stroma_tma,
             by = "image_tag") %>% 
  mutate(core = str_match(image_tag, "\\[(.*)\\]")[,2]) %>% 
  unite(core, c(core,tma), remove = FALSE) %>% 
  full_join(tma_map,
            .,
            by = "core") %>% 
  filter(!str_detect(suid, "ctrl")) %>% 
  drop_na(image_location) %>% 
  rename(stroma_area_analyzed_µm2 = stroma_area_analyzed_mm2,
         stroma_area_analyzed_mm2 = stroma_area_analyzed_mm2_2)

write_csv(total_data, "Mapped_total_macrophage_raw_data_07192023.csv")
write_csv(tumor_data, "Mapped_tumor_macrophage_raw_data_07192023.csv")
write_csv(stroma_data, "Mapped_stroma_macrophage_raw_data_07192023.csv")


# Summarize
data_summ <- function(data){
  data <- data %>% 
    select(-c(core, image_tag, image_location)) %>% 
    group_by(suid) %>% 
    mutate(across(where(is.numeric), .fns = ~ mean(., na.rm = TRUE))) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric), .fns = ~ na_if(., NaN))) %>% 
    arrange(suid, tma) %>% 
    distinct(suid, .keep_all = TRUE)
}

total_data <- data_summ(total_data) %>% 
  mutate(compartment = "Total") %>% 
  select(suid, compartment, everything(), -tma, tma)
tumor_data <- data_summ(tumor_data) %>% 
  mutate(compartment = "Tumor") %>% 
  select(-tma)
stroma_data <- data_summ(stroma_data) %>% 
  mutate(compartment = "Stroma") %>% 
  select(-tma)


# Bind / Join data
long_data <- bind_rows(total_data,
                       tumor_data %>% 
                         `colnames<-`(str_remove(colnames(.), "tumor_")),
                       stroma_data %>% 
                         `colnames<-`(str_remove(colnames(.), "stroma_"))) %>% 
  arrange(suid) %>% 
  mutate(across(where(is.numeric), .fns = ~ round(., 3))) %>% 
  
  group_by(suid) %>% 
  fill(tma, .direction = "updown")
  
write_csv(long_data, "Clean_mapped_macrophage_data_long_format_07192023.csv")


wide_data <- 
  full_join(total_data %>% 
              select(suid, tma, everything(), -compartment) %>% 
              rename_at(vars(-"suid", -"tma",
              ), ~ paste0("total_", .)),
            tumor_data %>% 
              select(-compartment) %>% 
              rename_at(vars(starts_with("m") |
                                starts_with("neutral") | 
                                starts_with("area")
              ), ~ paste0("tumor_", .)),
            by = "suid") %>% 
  full_join(.,
            stroma_data %>% 
              select(-compartment) %>% 
              rename_at(vars(starts_with("m") |
                               starts_with("neutral") | 
                               starts_with("area")
              ), ~ paste0("stroma_", .)),
            by = "suid") %>% 
  # `colnames<-`(str_replace(colnames(.), "area_analyzed_mm2", "area_analyzed_µm2")) %>% 
  # `colnames<-`(str_replace(colnames(.), "area_analyzed_µm2_2", "area_analyzed_mm2")) %>% 
  # drop_na(total_m1_simple_cells) %>% 
  arrange(suid) %>% 
  mutate(across(where(is.numeric), .fns = ~ round(., 3))) 

write_csv(wide_data, "Clean_mapped_macrophage_data_wide_format_07192023.csv")









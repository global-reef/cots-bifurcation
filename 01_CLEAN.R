### 01_CLEAN.R
### Clean CoTS, CPCe substrate, and SST data for bifurcation simulation

######### 1. Paths #######

data_raw_dir <- file.path(getwd(), "data", "raw")
data_processed_dir <- file.path(getwd(), "data", "processed")

cots_path <- file.path(data_raw_dir, paste0(analysis_date, "_COTS.csv"))
cpce_path <- file.path(data_raw_dir, paste0(analysis_date, "_cpce.csv"))
sst_path1  <- file.path(data_raw_dir, "copernicus_sst.nc") # has 2020 to 2024
sst_path2  <- file.path(data_raw_dir, "copernicus_thetao.nc") # has 06.01.2022 to 06.24.2022


######### 2. Load Raw Data #######

cots_raw <- read_csv(cots_path, show_col_types = FALSE)
cpce_raw <- read_csv(cpce_path, show_col_types = FALSE)

str(cots_raw)
str(cpce_raw)

sst_nc <- terra::rast(sst_path1)
tao_nc <- terra::rast(sst_path2)


sst_nc
tao_nc

######### 3. Standardise Site Names #######

site_code_lookup <- function(site) {
  case_when(
    str_detect(site, regex("^aow\\s*leuk|^ao\\s*leuk", ignore_case = TRUE)) ~ "AL",
    str_detect(site, regex("^shark\\s*island", ignore_case = TRUE)) ~ "SI",
    str_detect(site, regex("^tanote", ignore_case = TRUE)) ~ "TB",
    str_detect(site, regex("^red\\s*rock", ignore_case = TRUE)) ~ "RR",
    str_detect(site, regex("^green\\s*rock", ignore_case = TRUE)) ~ "GR",
    str_detect(site, regex("^twins", ignore_case = TRUE)) ~ "TW",
    TRUE ~ NA_character_
  )
}

site_name_lookup <- function(site_code) {
  case_when(
    site_code == "AL" ~ "Aow Leuk",
    site_code == "SI" ~ "Shark Island",
    site_code == "TB" ~ "Tanote Bay",
    site_code == "GR" ~ "Green Rock",
    site_code == "RR" ~ "Red Rock",
    site_code == "TW" ~ "Twins",
    TRUE ~ NA_character_
  )
}

site_type_lookup <- function(site_code) {
  case_when(
    site_code %in% c("AL", "SI", "TB") ~ "Control",
    site_code %in% c("GR", "RR", "TW") ~ "High Density",
    TRUE ~ NA_character_
  )
}
# check 
cots_raw %>%
  distinct(Site) %>%
  mutate(site_code = site_code_lookup(Site))

cpce_raw %>%
  distinct(Site) %>%
  mutate(site_code = site_code_lookup(Site))
### 04. CLEAN COTS #######
#### 4A. Clean CoTS Individual Data #######

cots_indiv <- cots_raw %>%
  clean_names() %>%
  select(-starts_with("x")) %>%
  rename(
    date_raw = date_mm_dd_yyyy
  ) %>%
  filter(!is.na(site)) %>%
  mutate(
    date = lubridate::mdy(date_raw),
    
    site_code = site_code_lookup(site),
    site = site_name_lookup(site_code),
    site_type = site_type_lookup(site_code),
    
    survey_id = str_squish(survey_id),
    active = str_to_upper(str_squish(active)),
    substrate = str_to_upper(str_squish(substrate)),
    growth_form = str_to_upper(str_squish(growth_form)),
    removal = str_to_upper(str_squish(removal)),
    
    size_cm = parse_number(size_cm),
    depth_m = parse_number(depth_m),
    
    site_code = factor(site_code, levels = c("AL", "SI", "TB", "GR", "RR", "TW")),
    site = factor(site, levels = c("Aow Leuk", "Shark Island", "Tanote Bay", "Green Rock", "Red Rock", "Twins")),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    survey_id = factor(survey_id),
    active = factor(active, levels = c("N", "Y")),
    substrate = factor(substrate),
    growth_form = factor(growth_form),
    removal = factor(removal)
  )
# check 
str(cots_indiv)

cots_indiv %>%
  count(site, site_code, site_type)

cots_indiv %>%
  filter(is.na(site_code)) %>%
  distinct(site)
#### 4B. Clean CoTS Survey Data #######

cots_survey <- cots_indiv %>%
  group_by(site, site_code, site_type, survey_id, date) %>%
  summarise(
    cots_count = n(),
    mean_size_cm = mean(size_cm, na.rm = TRUE),
    mean_depth_m = mean(depth_m, na.rm = TRUE),
    avg_survey_depth_m = first(avg_depth_m),
    duration_min = first(duration_min),
    vis_m = first(vis_m),
    .groups = "drop"
  ) %>%
  mutate(
    area_ha = 0.2,
    cots_ha = cots_count / area_ha
  )

str(cots_survey)
summary(cots_survey$cots_ha)

# Check for duplicated specimen numbers within survey
cots_indiv %>%
  count(site_code, survey_id, date, specimen) %>%
  filter(n > 1, !is.na(specimen))

### 05. CLEAN SUBSTRATE ####
#### 5A. Clean CPCe Point Data #######

valid_substrate <- c("HC", "TUR", "MAC", "OB", "AB", "AN", "UKN")
valid_gf <- c("MA", "SMA", "CAE", "CB", "TB", "ARB", "DI", "SOL", "FOL", "LM", "ENC")

substrate <- cpce_raw %>%
  clean_names() %>%
  mutate(
    across(c(site, transect, substrate, gf, scar), ~ str_to_upper(str_squish(.x))),
    
    site_code = site_code_lookup(site),
    site = site_name_lookup(site_code),
    site_type = site_type_lookup(site_code),
    
    quadrat = as.integer(quadrat),
    cpce_point = as.integer(cpce_point),
    
    substrate = case_when(
      substrate %in% valid_substrate ~ substrate,
      substrate %in% c("UK", "??", "@@@@", "0", "NUMBER NOT IN IMAGE", "UKN/OB") ~ "UKN",
      substrate %in% c("TUIR", "TURB", "TURF", "TURTUR", "ATUR", "UR") ~ "TUR",
      substrate %in% c("A B", "ABA", "RAB", "SAB", "HAB") ~ "AB",
      substrate %in% c("OB?", "OB (SP)", "O", "SP", "SC") ~ "OB",
      substrate %in% c("HHC", "HCC", "HCH", "HCMA", "HCENC", "HCTUR") ~ "HC",
      TRUE ~ NA_character_
    ),
    
    gf = case_when(
      gf %in% valid_gf ~ gf,
      gf %in% c("MAS", "M", "MAMA") ~ "MA",
      gf %in% c("SUB", "SM", "MS") ~ "SMA",
      gf == "DIG" ~ "DI",
      gf == "LAM" ~ "LM",
      gf == "TAB" ~ "TB",
      gf %in% c("FOLF", "FO") ~ "FOL",
      gf %in% c("SM/ENC", "ENC/SM") ~ "ENC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(site_code), !is.na(substrate)) %>%
  select(site_code, site, site_type, transect, quadrat, cpce_point, substrate, gf, scar) %>%
  mutate(
    site_code = factor(site_code, levels = c("AL", "SI", "TB", "GR", "RR", "TW")),
    site = factor(site, levels = c("Aow Leuk", "Shark Island", "Tanote Bay", "Green Rock", "Red Rock", "Twins")),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    transect = factor(transect),
    substrate = factor(substrate, levels = valid_substrate),
    gf = factor(gf, levels = valid_gf),
    scar = factor(scar)
  )

str(substrate)

#### 5B. Substrate Transect Cover #######

substrate_transect <- substrate %>%
  count(site_code, site, site_type, transect, substrate, name = "n_points") %>%
  group_by(site_code, site, site_type, transect) %>%
  mutate(
    total_points = sum(n_points),
    prop = n_points / total_points,
    pct = prop * 100
  ) %>%
  ungroup() %>%
  mutate(
    site_code = factor(site_code, levels = c("AL", "SI", "TB", "GR", "RR", "TW")),
    site = factor(site, levels = c("Aow Leuk", "Shark Island", "Tanote Bay", "Green Rock", "Red Rock", "Twins")),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    transect = factor(transect),
    substrate = factor(substrate, levels = valid_substrate)
  )

str(substrate_transect)

#### 5C. Substrate Wide #######

substrate_wide <- substrate_transect %>%
  select(site_code, site, site_type, transect, substrate, pct) %>%
  pivot_wider(
    names_from = substrate,
    values_from = pct,
    values_fill = 0
  ) %>%
  mutate(
    site_code = factor(site_code, levels = c("AL", "SI", "TB", "GR", "RR", "TW")),
    site = factor(site, levels = c("Aow Leuk", "Shark Island", "Tanote Bay", "Green Rock", "Red Rock", "Twins")),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    transect = factor(transect)
  )

str(substrate_wide)
# checks 
substrate %>%
  count(substrate, sort = TRUE)

substrate_transect %>%
  filter(substrate == "HC") %>%
  arrange(site_code, transect)

substrate_wide %>%
  summarise(across(c(HC, TUR, MAC, OB, AB, UKN), ~ mean(.x, na.rm = TRUE)))


### 06 TEMPERATURE DATA ####
### Extract and combine near-surface ocean temperature from:
### 1) MULTIOBS `to` file, 2020-2024
### 2) ANALYSISFORECAST `thetao` file, 2022-2026
###
### Unified model variable = temp_c

#### 6A. Define Koh Tao extraction point ####

kny_lat <- 10 + 7/60 + 2.42/3600
kny_lon <- 99 + 48/60 + 51.84/3600

kny_point <- terra::vect(
  data.frame(lon = kny_lon, lat = kny_lat),
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

# https://doi.org/10.48670/moi-00052
#### 6B. Prepare MULTIOBS temperature data, 2020-2024 ####
#### variable: to
#### source file: copernicus_sst.nc

sst_depths <- terra::depth(sst_nc)
surface_depth <- sst_depths[which.min(abs(sst_depths))]

surface_depth

sst_surface <- sst_nc[[sst_depths == surface_depth]]

sst_dates <- as.Date(terra::time(sst_surface))

range(sst_dates)
terra::nlyr(sst_surface)
terra::depth(sst_surface)[1:10]

extract_multobs_temp <- function(temp_date) {
  
  if (is.na(temp_date)) return(NA_real_)
  
  if (temp_date < min(sst_dates) | temp_date > max(sst_dates)) {
    return(NA_real_)
  }
  
  nearest_layer <- which.min(abs(sst_dates - temp_date))
  
  temp_value <- terra::extract(
    sst_surface[[nearest_layer]],
    kny_point
  )
  
  as.numeric(temp_value[1, 2])
}

multobs_daily <- tibble(
  date = sst_dates,
  temp_multobs_c = map_dbl(sst_dates, extract_multobs_temp)
)

summary(multobs_daily$temp_multobs_c)



#### 6C. Prepare thetao temperature data, 2022-2026 ####
#### variable: thetao
#### source file: copernicus_thetao.nc - https://doi.org/10.48670/moi-00016
#### note: this file is already a 1-cell Koh Tao extraction

thetao_dates <- as.Date(terra::time(tao_nc))

range(thetao_dates)
terra::depth(tao_nc)

thetao_daily <- tibble(
  date = thetao_dates,
  temp_thetao_c = as.numeric(terra::values(tao_nc))
)

summary(thetao_daily$temp_thetao_c)


#### 6D. Combine temperature sources ####
#### Use MULTIOBS through 2024-12-31
#### Use thetao from 2025-01-01 onward

temp_daily <- full_join(multobs_daily, thetao_daily, by = "date") %>%
  arrange(date) %>%
  mutate(
    temp_c = case_when(
      date <= as.Date("2024-12-31") ~ temp_multobs_c,
      date >= as.Date("2025-01-01") ~ temp_thetao_c,
      TRUE ~ NA_real_
    ),
    temp_source = case_when(
      date <= as.Date("2024-12-31") & !is.na(temp_multobs_c) ~ "MULTIOBS_to_surface",
      date >= as.Date("2025-01-01") & !is.na(temp_thetao_c) ~ "ANALYSISFORECAST_thetao_0.49m",
      TRUE ~ NA_character_
    )
  )

str(temp_daily)
summary(temp_daily$temp_c)

temp_daily %>%
  count(temp_source)


#### 6E. Join temperature to CoTS survey table ####

cots_survey <- cots_survey %>%
  left_join(
    temp_daily %>%
      select(date, temp_c, temp_source),
    by = "date"
  )

str(cots_survey)

cots_survey %>%
  summarise(
    n_surveys = n(),
    n_missing_temp = sum(is.na(temp_c)),
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    min_temp_c = min(temp_c, na.rm = TRUE),
    max_temp_c = max(temp_c, na.rm = TRUE)
  )

cots_survey %>%
  count(temp_source)


#### 6F. Check surveys with missing temperature ####

cots_survey %>%
  filter(is.na(temp_c)) %>%
  select(site, site_code, survey_id, date, cots_count, cots_ha) %>%
  arrange(date)


#### 6G. Save temperature table ####

write_csv(
  temp_daily,
  file.path(data_processed_dir, paste0(analysis_date, "_temp_daily.csv"))
)

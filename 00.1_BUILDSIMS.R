# ============================================================
# CoTS bifurcation simulation model
# Plug-and-play draft script
# Global Reef | Koh Tao
# ============================================================

# Purpose:
# Build and test the full modelling workflow using real CoTS data
# and dummy substrate data. When real substrate data is available,
# replace the dummy substrate input with the real file path and keep
# the rest of the workflow unchanged.

# ============================================================
# 0. Setup
# ============================================================

library(tidyverse)
library(lubridate)
library(terra)

analysis_date <-"2026.06.04" # update each run 



# ---- Folder paths ----

data_raw_dir       <- "data/raw"
data_processed_dir <- "data/processed"
outputs_dir        <- "outputs"
figures_dir        <- file.path(outputs_dir, "figures")
tables_dir         <- file.path(outputs_dir, "tables")

dir.create(data_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Input paths ----

cots_path <- file.path(data_raw_dir, paste0(analysis_date, "_COTS.csv"))
cpce_path <- file.path(data_raw_dir, paste0(analysis_date, "_cpce.csv"))
# This is where the real substrate file will eventually go:
# real_substrate_path <- file.path(data_raw_dir, "mario_substrate_raw.csv")

use_dummy_substrate <- TRUE # change this when we have real substrate data 

# ---- Survey area assumption ----
# Current density logic:
# cots_density_ha = cots / 0.2
# This assumes each survey covers 0.2 ha.

survey_area_ha <- 0.2


# ============================================================
# 1. Helper functions
# ============================================================

clean_names_simple <- function(df) {
  names(df) <- names(df) |>
    str_trim() |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("_$", "")
  
  df
}

standardise_site_names <- function(site) {
  site |> str_trim()
}


cots_growth <- function(N, r_T, State_1, K_t, K_c) {
  r_T * (N - State_1) * (1 - N / K_c) * (N / K_t - 1)
}


simulate_cots <- function(N0, r_T, State_1, K_t, K_c,
                          time_steps = 100,
                          dt = 1) {
  
  N <- numeric(time_steps)
  N[1] <- N0
  
  for (t in 2:time_steps) {
    
    growth <- cots_growth(
      N = N[t - 1],
      r_T = r_T,
      State_1 = State_1,
      K_t = K_t,
      K_c = K_c
    )
    
    N[t] <- N[t - 1] + growth * dt
    
    # Prevent impossible values
    N[t] <- max(N[t], 0)
    
    # Prevent numerical blow-up during early testing
    N[t] <- min(N[t], K_c * 1.5)
  }
  
  tibble(
    time = 1:time_steps,
    cots_density_ha = N
  )
}


calc_r_T <- function(SST, r_max = 0.005, T_opt = 28.7, thermal_width = 3) {
  r_max * exp(-((SST - T_opt)^2) / (2 * thermal_width^2))
}

calc_r_T_ashna <- function(SST) {
  
  case_when(
    is.na(SST) ~ NA_real_,
    
    # negligible growth below 26°C
    SST < 26 ~ 0,
    
    # rising from 26 to 27°C
    SST >= 26 & SST < 27 ~ scales::rescale(SST, to = c(0, 0.1), from = c(26, 27)),
    
    # rising from 27 to 28°C
    SST >= 27 & SST < 28 ~ scales::rescale(SST, to = c(0.1, 0.25), from = c(27, 28)),
    
    # optimum plateau from 28 to 30°C
    SST >= 28 & SST <= 30 ~ 0.25,
    
    # decline from 30 to 32°C
    SST > 30 & SST <= 32 ~ scales::rescale(SST, to = c(0.25, 0.05), from = c(30, 32)),
    
    # negligible/low growth above 32°C
    SST > 32 ~ 0.05
  )
}

classify_threshold_position <- function(N0, K_t, buffer = 0.10) {
  case_when(
    N0 < K_t * (1 - buffer) ~ "below_threshold",
    N0 >= K_t * (1 - buffer) & N0 <= K_t * (1 + buffer) ~ "near_threshold",
    N0 > K_t * (1 + buffer) ~ "above_threshold",
    TRUE ~ NA_character_
  )
}


classify_sim_outcome <- function(final_density, State_1, K_t, K_c) {
  case_when(
    final_density <= K_t ~ "returns_to_low_density",
    final_density > K_t & final_density < K_c ~ "intermediate",
    final_density >= K_c ~ "outbreak_domain",
    TRUE ~ NA_character_
  )
}


# ============================================================
# 2. Create dummy substrate data
# ============================================================

# This dummy dataset mimics the expected raw substrate format:
# Site | Location | Date | Quadrat | Point | Substrate | Genus | MORPH
#
# Dummy hard coral cover is based on previous mean estimates:
# Twins Wall = 31.92%
# Green Rock Wall = 39.52%
# Red Rock Wall = 48.40%
#
# Later, Mario's real substrate data should be cleaned into this same logic.

set.seed(42)

dummy_coral_targets <- tibble(
  Site = c("Twins Wall", "Green Rock Wall", "Red Rock Wall"),
  target_hc_prop = c(0.3192, 0.3952, 0.4840)
)

dummy_substrate <- expand_grid(
  Site = dummy_coral_targets$Site,
  Location = c("Wall"),
  Date = as.Date(c("2024-05-01", "2024-11-01")),
  Quadrat = 1:10,
  Point = 1:20
) %>%
  left_join(dummy_coral_targets, by = "Site") %>%
  rowwise() %>%
  mutate(
    non_hc_remaining = 1 - target_hc_prop,
    
    Substrate = sample(
      c("HC", "RK", "RU", "SP", "SC", "UNKN"),
      size = 1,
      prob = c(
        target_hc_prop,
        non_hc_remaining * 0.50, # RK
        non_hc_remaining * 0.20, # RU
        non_hc_remaining * 0.12, # SP
        non_hc_remaining * 0.10, # SC
        non_hc_remaining * 0.08  # UNKN
      )
    ),
    
    Genus = if_else(
      Substrate == "HC",
      sample(c("Acropora", "Porites", "Pavona", "Favites", NA_character_), 1),
      NA_character_
    ),
    
    MORPH = if_else(
      Substrate == "HC",
      sample(c("ENC", "SUBM", "BR", "MASS", NA_character_), 1),
      NA_character_
    )
  ) %>%
  ungroup() %>%
  select(-target_hc_prop, -non_hc_remaining)

write_csv(
  dummy_substrate,
  file.path(data_raw_dir, "dummy_substrate.csv")
)


# ============================================================
# 3. Read and clean CoTS data
# ============================================================

# Expected raw CoTS format:
# Site | Survey_ID | Date_mm.dd.yyyy | Avg_depth_m | Time | Duration_min |
# Vis_m | Researcher | Specimen | Active | Substrate | Growth Form |
# Size_cm | Depth_m | Site_Type

 
cots_raw <- read_csv(cots_path, show_col_types = FALSE) %>%
  clean_names_simple() %>%
  select(!starts_with("x")) %>%      # removes x17/x18/x19 if clean_names renamed them that way
  select(!starts_with("_")) %>%
  select(!starts_with("..."))        # removes ...17/...18/...19 if they remain

target_sites <- c(
  "Green Rock Wall",
  "Red Rock Wall",
  "Twins Wall"
)

cots_clean <- cots_raw %>%
  mutate(
    site = standardise_site_names(site),
    survey_id = as.character(survey_id),
    date = mdy(date_mm_dd_yyyy),
    specimen = readr::parse_number(as.character(specimen)),
    avg_depth_m = readr::parse_number(as.character(avg_depth_m)),
    duration_min = readr::parse_number(as.character(duration_min)),
    vis_m = readr::parse_number(as.character(vis_m)),
    size_cm = readr::parse_number(as.character(size_cm)),
    depth_m = readr::parse_number(as.character(depth_m))
  ) %>%
  filter(site %in% target_sites) %>%
  select(!starts_with("_"))


# CoTS raw data are one row per individual.
# Summarise to survey level.

cots_summary <- cots_clean %>%
  group_by(site, survey_id, date, site_type) %>%
  summarise(
    cots = n(),
    mean_size_cm = mean(size_cm, na.rm = TRUE),
    mean_depth_m = mean(depth_m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    survey_area_ha = survey_area_ha,
    cots_density_ha = cots / survey_area_ha
  )

write_csv(
  cots_clean,
  file.path(data_processed_dir, paste0(analysis_date, "_cots_clean.csv"))
)




# extracting SST from copernicus marine service information ### 
sst_nc <- rast("data/raw/copernicus_sst.nc")
## reference: DOI (product): https://doi.org/10.48670/moi-00052

# References:
#   
#   Guinehut S., A.-L. Dhomps, G. Larnicol and P.-Y. Le Traon, 2012: High resolution 3D temperature and salinity fields derived from in situ and satellite observations. Ocean Sci., 8(5):845–857.
# Mulet, S., M.-H. Rio, A. Mignot, S. Guinehut and R. Morrow, 2012: A new estimate of the global 3D geostrophic ocean circulation based on satellite data and in-situ measurements. Deep Sea Research Part II : Topical Studies in Oceanography, 77–80(0):70–81.
#  We recommend this convention to cite a specific product:
# "Product Title. E.U. Copernicus Marine Service Information (CMEMS). Marine Data Store (MDS). DOI: 10.48670/moi-xxxxx (Accessed on DD MMM YYYY)"
sst_nc

kny_lat <- 10 + 7/60 + 2.42/3600
kny_lon <- 99 + 48/60 + 51.84/3600

kny_point <- terra::vect(
  data.frame(lon = kny_lon, lat = kny_lat),
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)
sst_dates <- as.Date(terra::time(sst_nc))
sst_dates[1:10]

terra::depth(sst_nc)
sst_surface <- sst_nc[[terra::depth(sst_nc) == 0]]
sst_dates <- as.Date(terra::time(sst_surface))

extract_sst_by_date <- function(survey_date) {
  
  nearest_layer <- which.min(abs(sst_dates - survey_date))
  
  sst_value <- terra::extract(
    sst_surface[[nearest_layer]],
    kny_point
  )
  
  as.numeric(sst_value[1, 2])
}
# add sst to cots 
cots_summary <- cots_summary %>%
  mutate(
    SST = purrr::map_dbl(date, extract_sst_by_date)
  )
write_csv(
  cots_summary,
  file.path(data_processed_dir, paste0(analysis_date, "_cots_summary.csv"))
)

# ============================================================
# 4. Read and clean substrate data
# ============================================================

if (use_dummy_substrate) {
  
  substrate_raw <- read_csv(
    file.path(data_raw_dir, "dummy_substrate.csv"),
    show_col_types = FALSE
  )
  
} else {
  
  # Uncomment this once real substrate data are available:
  # substrate_raw <- read_csv(real_substrate_path, show_col_types = FALSE)
  
}

substrate_clean <- substrate_raw %>%
  clean_names_simple() %>%
  mutate(
    site = standardise_site_names(site),
    date = as.Date(date),
    quadrat = as.character(quadrat),
    point = as.integer(point),
    substrate = str_to_upper(substrate),
    is_live_hard_coral = substrate == "HC",
    quadrat_id = paste(site, location, date, quadrat, sep = "_")
  )

# Summarise live coral cover at site-date level.
# This is the object the model needs.

live_coral_summary <- substrate_clean %>%
  group_by(site, date) %>%
  summarise(
    n_points = n(),
    n_live_hard_coral = sum(is_live_hard_coral, na.rm = TRUE),
    live_coral_prop = n_live_hard_coral / n_points,
    live_coral_percent = live_coral_prop * 100,
    n_quadrats = n_distinct(quadrat_id),
    .groups = "drop"
  )

write_csv(
  substrate_clean,
  file.path(data_processed_dir, paste0(analysis_date, "_substrate_clean.csv"))
)

write_csv(
  live_coral_summary,
  file.path(data_processed_dir, paste0(analysis_date, "_live_coral_summary.csv"))
)


# ============================================================
# 5. Create observed site-level CoTS starting conditions
# ============================================================

# For the first version, use site-level maximum and mean observed density.
# Max density is useful as a stress-test starting condition.
# Mean density is useful as a typical observed condition.
site_cots_inputs <- cots_summary %>%
  group_by(site) %>%
  summarise(
    median_density = median(cots_density_ha, na.rm = TRUE),
    mean_density   = mean(cots_density_ha, na.rm = TRUE),
    q75_density    = quantile(cots_density_ha, 0.75, na.rm = TRUE),
    q90_density    = quantile(cots_density_ha, 0.90, na.rm = TRUE),
    max_density    = max(cots_density_ha, na.rm = TRUE),
    n_surveys      = n(),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(
      median_density,
      mean_density,
      q75_density,
      q90_density,
      max_density
    ),
    names_to = "N0_type",
    values_to = "N0"
  ) %>%
  mutate(
    N0_type = recode(
      N0_type,
      median_density = "median_observed_density",
      mean_density   = "mean_observed_density",
      q75_density    = "q75_observed_density",
      q90_density    = "q90_observed_density",
      max_density    = "max_observed_density"
    )
  )
site_n0_sst <- site_cots_inputs %>%
  left_join(
    cots_summary %>%
      select(site, date, survey_id, cots_density_ha, SST),
    by = "site"
  ) %>%
  group_by(site, N0_type, N0) %>%
  slice_min(
    order_by = abs(cots_density_ha - N0),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    site,
    N0_type,
    N0,
    representative_survey_id = survey_id,
    representative_date = date,
    SST,
    r_T = calc_r_T_ashna(SST)
  )

site_cots_inputs <- site_n0_sst

# ============================================================
# 6. Build parameter grid
# ============================================================

# ---- Endemic density scenarios ----
# State_1 should be around but below 10 CoTS/ha.

state1_scenarios <- tibble(
  state1_scenario = c("low_endemic", "base_endemic", "high_endemic"),
  State_1 = c(3, 8, 10)
)

# ---- Carrying capacity scenarios ----
# These are placeholder values for now.
# They should eventually be replaced with literature-cited scenario values.

kc_scenarios <- tibble(
  kc_scenario = c(
    "low_coral_low_K",
    "moderate_coral_moderate_K",
    "high_coral_high_K",
    "disturbed_reef_reduced_K",
    "recovery_reef_increasing_K"
  ),
  K_c = c(25, 50, 100, 20, 75)
)


# ---- Use latest live coral value per site for first plug-and-play version ----
# Later this can be changed to site-date matching or nearest-date matching.

live_coral_site <- live_coral_summary %>%
  group_by(site) %>%
  arrange(date) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  select(site, substrate_date = date, live_coral_prop, live_coral_percent)

# ---- Full parameter grid ----
parameter_grid <- site_cots_inputs %>%
  left_join(live_coral_site, by = "site") %>%
  crossing(state1_scenarios, kc_scenarios) %>%
  mutate(
    K_t = K_c * live_coral_prop,
    
    valid_state_order = State_1 < K_t & K_t < K_c,
    
    observed_threshold_position = case_when(
      N0 < State_1 ~ "below_endemic_state",
      N0 >= State_1 & N0 < K_t ~ "allee_basin",
      N0 >= K_t & N0 < K_c ~ "above_tipping_threshold",
      N0 >= K_c ~ "outbreak_overshoot",
      TRUE ~ NA_character_
    ),
    
    sim_id = row_number()
  )
# Keep invalid scenarios in the table for transparency,
# but only valid scenarios will be simulated.

write_csv(
  parameter_grid,
  file.path(tables_dir, paste0(analysis_date, "_parameter_grid.csv"))
)

parameter_grid_valid <- parameter_grid %>%
  filter(valid_state_order)


# ============================================================
# 7. Run simulations
# ============================================================

sim_results <- parameter_grid_valid %>%
  select(
    sim_id,
    site,
    N0_type,
    N0,
    representative_survey_id,
    representative_date,
    SST,
    r_T,
    substrate_date,
    live_coral_prop,
    live_coral_percent,
    state1_scenario,
    State_1,
    kc_scenario,
    K_c,
    K_t,
    observed_threshold_position
  ) %>%
  pmap_dfr(function(sim_id, site, N0_type, N0,
                    representative_survey_id, representative_date,
                    SST, r_T,
                    substrate_date,
                    live_coral_prop, live_coral_percent,
                    state1_scenario, State_1,
                    kc_scenario, K_c, K_t,
                    observed_threshold_position) {
    
    simulate_cots(
      N0 = N0,
      r_T = r_T,
      State_1 = State_1,
      K_t = K_t,
      K_c = K_c,
      time_steps = 100,
      dt = 1
    ) %>%
      mutate(
        sim_id = sim_id,
        site = site,
        N0_type = N0_type,
        N0 = N0,
        representative_survey_id = representative_survey_id,
        representative_date = representative_date,
        SST = SST,
        r_T = r_T,
        substrate_date = substrate_date,
        live_coral_prop = live_coral_prop,
        live_coral_percent = live_coral_percent,
        state1_scenario = state1_scenario,
        State_1 = State_1,
        kc_scenario = kc_scenario,
        K_c = K_c,
        K_t = K_t,
        observed_threshold_position = observed_threshold_position
      )
  })

write_csv(
  sim_results,
  file.path(tables_dir, paste0(analysis_date, "_simulation_results.csv"))
)


# ============================================================
# 8. Summarise simulation outcomes
# ============================================================

sim_summary <- sim_results %>%
  group_by(
    sim_id,
    site,
    N0_type,
    N0,
    live_coral_prop,
    live_coral_percent,
    state1_scenario,
    State_1,
    kc_scenario,
    K_c,
    K_t,
    SST,
    r_T,
    observed_threshold_position
  ) %>%
  summarise(
    final_density = last(cots_density_ha),
    max_density = max(cots_density_ha),
    outcome = classify_sim_outcome(final_density, State_1, K_t, K_c),
    .groups = "drop"
  )

write_csv(
  sim_summary,
  file.path(tables_dir, paste0(analysis_date, "_simulation_summary.csv"))
)

# Outcome probabilities by site

outcome_by_site <- sim_summary %>%
  count(site, N0_type, outcome) %>%
  group_by(site, N0_type) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()

write_csv(
  outcome_by_site,
  file.path(tables_dir, paste0(analysis_date, "_outcome_by_site.csv"))
)


# ============================================================
# 9. Basic plots
# ============================================================

# ---- Plot 1: Observed density vs derived threshold ----

p_threshold <- parameter_grid_valid %>%
  ggplot(aes(x = K_t, y = N0)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(aes(shape = N0_type), alpha = 0.7) +
  facet_wrap(~ site, scales = "free") +
  labs(
    x = "Derived tipping threshold, K_t (CoTS/ha)",
    y = "Observed starting density, N0 (CoTS/ha)",
    shape = "Observed density type",
    title = "Observed CoTS density relative to derived tipping thresholds"
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, paste0(analysis_date, "_observed_density_vs_threshold.png")),
  p_threshold,
  width = 10,
  height = 7,
  dpi = 300
)
p_threshold


# ---- Plot 2: Example simulation trajectories ----

p_trajectories <- sim_results %>%
  filter(N0_type == "max_observed_density") %>%
  ggplot(aes(x = time, y = cots_density_ha, group = sim_id)) +
  geom_line(alpha = 0.15) +
  facet_wrap(~ site, scales = "free_y") +
  labs(
    x = "Time step",
    y = expression("CoTS density"~ha^{-1}),
    title = "Simulated CoTS trajectories across parameter scenarios"
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, paste0(analysis_date, "_simulation_trajectories.png")),
  p_trajectories,
  width = 10,
  height = 7,
  dpi = 300
)

p_trajectories
# ---- Plot 3: Outcome proportions by site ----

p_outcomes <- outcome_by_site %>%
  ggplot(aes(x = site, y = prop, fill = outcome)) +
  geom_col(position = "stack") +
  facet_wrap(~ N0_type) +
  coord_flip() +
  labs(
    x = "Site",
    y = "Proportion of simulations",
    fill = "Outcome",
    title = "Simulation outcomes by site"
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, paste0(analysis_date, "_outcome_proportions_by_site.png")),
  p_outcomes,
  width = 10,
  height = 7,
  dpi = 300
)

p_outcomes

# check parameter grid # 
parameter_grid %>%
  count(site, observed_threshold_position)
# ============================================================
# 10. Quick console summaries
# ============================================================

cat("\nCoTS bifurcation simulation complete.\n")
cat("Number of CoTS surveys:", nrow(cots_summary), "\n")
cat("Number of substrate site-date summaries:", nrow(live_coral_summary), "\n")
cat("Number of parameter combinations:", nrow(parameter_grid), "\n")
cat("Number of valid parameter combinations:", nrow(parameter_grid_valid), "\n")
cat("Number of simulations:", n_distinct(sim_results$sim_id), "\n")

cat("\nOutputs saved to:\n")
cat("- Processed data:", data_processed_dir, "\n")
cat("- Tables:", tables_dir, "\n")
cat("- Figures:", figures_dir, "\n")


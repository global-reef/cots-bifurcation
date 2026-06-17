### 03_SIMULATION.R
### Monte Carlo bifurcation simulation for Koh Tao CoTS dynamics

### Goal:
### Use observed CoTS density, observed coral cover, and observed temperature
### to simulate a range of possible CoTS population trajectories.
###
### This is a Monte Carlo simulation, not full MCMC.
### We are sampling plausible ecological parameter values and propagating
### uncertainty through the bifurcation model.


######### 1. Prepare coral cover input data #######

### We keep coral cover at the transect level.
### This preserves within-site variation instead of reducing each site to one mean value.

coral_cover_transect <- substrate_wide %>%
  select(site_code, site, site_type, transect, HC) %>%
  mutate(
    HC_prop = HC / 100
  ) %>%
  filter(
    !is.na(HC_prop),
    HC_prop > 0
  )

### Check spread of coral cover values
summary(coral_cover_transect$HC_prop)

coral_cover_transect %>%
  group_by(site_code, site, site_type) %>%
  summarise(
    n_transects = n(),
    mean_HC = mean(HC, na.rm = TRUE),
    min_HC = min(HC, na.rm = TRUE),
    max_HC = max(HC, na.rm = TRUE),
    .groups = "drop"
  )


######### 2. Prepare CoTS starting density input data #######

### We use survey-level CoTS density as possible starting densities.
### This preserves survey-to-survey variation within each site.

cots_starting_density <- cots_survey %>%
  filter(
    !is.na(cots_ha),
    cots_ha >= 0,
    !is.na(site_code)
  ) %>%
  select(site_code, site, site_type, survey_id, date, cots_ha)

### Check spread of observed starting densities
summary(cots_starting_density$cots_ha)

cots_starting_density %>%
  group_by(site_code, site, site_type) %>%
  summarise(
    n_surveys = n(),
    mean_cots_ha = mean(cots_ha, na.rm = TRUE),
    median_cots_ha = median(cots_ha, na.rm = TRUE),
    max_cots_ha = max(cots_ha, na.rm = TRUE),
    .groups = "drop"
  )


######### 3. Prepare temperature forcing #######

### Use the cleaned daily temperature time series from 01_CLEAN.R.
### The model will run one time step per day.

temp_forcing <- temp_daily %>%
  filter(!is.na(temp_c)) %>%
  arrange(date) %>%
  select(date, temp_c, temp_source)

summary(temp_forcing$temp_c)

temp_forcing %>%
  count(temp_source)


######### 4. Define temperature response function #######

### r_temp() converts daily temperature into a growth / recruitment multiplier.
###
### Current version follows the broad Ashna / Kamya-style logic:
### - below 26 C = negligible growth
### - 28 to 30 C = highest growth
### - above 30 C = reduced growth
###
### Later, this can be replaced with a smoother continuous function.

r_temp <- function(temp_c) {
  case_when(
    temp_c < 26 ~ 0,
    temp_c >= 26 & temp_c < 27 ~ 0.05,
    temp_c >= 27 & temp_c < 28 ~ 0.10,
    temp_c >= 28 & temp_c <= 30 ~ 0.25,
    temp_c > 30 & temp_c <= 31 ~ 0.15,
    temp_c > 31 ~ 0.05,
    TRUE ~ NA_real_
  )
}

### Check temperature response values
temp_forcing <- temp_forcing %>%
  mutate(r_t = r_temp(temp_c))

summary(temp_forcing$r_t)


######### 5. Define bifurcation change function #######

### This function calculates the daily change in CoTS density.
###
### N      = current CoTS density, individuals ha-1
### r      = temperature-dependent growth/recruitment value
### state1 = lower stable endemic state
### K_T    = unstable tipping threshold
### K_C    = upper carrying capacity / outbreak state
### scale  = tuning parameter to keep daily changes biologically reasonable
###
### The scale parameter is important because the cubic equation can otherwise
### produce changes that are too large for daily time steps.

dN_dt <- function(N, r, state1, K_T, K_C, scale = 0.01) {
  scale * r * (N - state1) * (1 - N / K_C) * (N / K_T - 1)
}


######### 6. Define single simulation function #######

### This function simulates one CoTS trajectory through time.
###
### It starts at N0, then updates density each day using:
###
### N[t + 1] = N[t] + dN_dt
###
### Density is bounded between 0 and K_C.
### This prevents negative CoTS density and prevents unlimited growth above carrying capacity.

simulate_cots <- function(N0, temp_series, state1, K_T, K_C, scale = 0.01) {
  
  N <- numeric(length(temp_series))
  N[1] <- N0
  
  for (i in 2:length(temp_series)) {
    
    r_i <- r_temp(temp_series[i - 1])
    
    change_i <- dN_dt(
      N = N[i - 1],
      r = r_i,
      state1 = state1,
      K_T = K_T,
      K_C = K_C,
      scale = scale
    )
    
    N[i] <- N[i - 1] + change_i
    
    ### Keep density within biologically possible bounds
    N[i] <- max(0, min(N[i], K_C))
  }
  
  N
}


######### 7. Sample one parameter set #######

### This function creates one random simulation setup for one site.
###
### For each simulation, it samples:
### - one observed starting CoTS density from that site
### - one observed coral-cover transect from that site
### - one plausible State 1 value
### - one plausible K_T value
###
### Then it calculates:
###
### K_C = K_T / HC_prop

sample_parameter_set <- function(site_i) {
  
  N0_i <- cots_starting_density %>%
    filter(site_code == site_i) %>%
    slice_sample(n = 1) %>%
    pull(cots_ha)
  
  HC_i <- coral_cover_transect %>%
    filter(site_code == site_i) %>%
    slice_sample(n = 1) %>%
    pull(HC_prop)
  
  ### State 1 is sampled around 8 CoTS ha-1,
  ### but constrained between 3 and 10.
  state1_i <- rnorm(1, mean = 8, sd = 2)
  state1_i <- pmin(pmax(state1_i, 3), 10)
  
  ### K_T is sampled around 15 CoTS ha-1,
  ### but constrained between 10 and 25.
  K_T_i <- rnorm(1, mean = 15, sd = 3)
  K_T_i <- pmin(pmax(K_T_i, 10), 25)
  
  ### Carrying capacity is inferred from tipping threshold and coral cover.
  K_C_i <- K_T_i / HC_i
  
  tibble(
    site_code = site_i,
    N0 = N0_i,
    HC_prop = HC_i,
    state1 = state1_i,
    K_T = K_T_i,
    K_C = K_C_i
  )
}


######### 8. Generate Monte Carlo parameter draws #######

### Set seed so the simulation can be exactly repeated.
set.seed(42)

### Number of simulations.
### Start small while testing. Increase later.
n_sims <- 1000

### Sites to simulate
sites_to_simulate <- levels(cots_starting_density$site_code)
sites_to_simulate <- sites_to_simulate[!is.na(sites_to_simulate)]

### Generate parameter sets.
### Each row is one simulation.
params <- map_dfr(1:n_sims, function(i) {
  
  site_i <- sample(sites_to_simulate, size = 1)
  
  sample_parameter_set(site_i)
  
}) %>%
  filter(
    K_T > state1,
    K_C > K_T
  ) %>%
  mutate(
    sim_id = row_number()
  )

### Check sampled parameter distributions
summary(params)

params %>%
  group_by(site_code) %>%
  summarise(
    n_sims = n(),
    mean_N0 = mean(N0, na.rm = TRUE),
    mean_HC = mean(HC_prop, na.rm = TRUE),
    mean_state1 = mean(state1, na.rm = TRUE),
    mean_K_T = mean(K_T, na.rm = TRUE),
    mean_K_C = mean(K_C, na.rm = TRUE),
    .groups = "drop"
  )


######### 9. Run simulations #######

### This runs the bifurcation model for every parameter set.
###
### Output is a long table with:
### - sim_id
### - site_code
### - date
### - simulated CoTS density

sim_results <- params %>%
  split(.$sim_id) %>%
  map_dfr(function(p) {
    
    N_sim <- simulate_cots(
      N0 = p$N0,
      temp_series = temp_forcing$temp_c,
      state1 = p$state1,
      K_T = p$K_T,
      K_C = p$K_C,
      scale = 0.01
    )
    
    tibble(
      sim_id = p$sim_id,
      site_code = p$site_code,
      date = temp_forcing$date,
      N = N_sim,
      N0 = p$N0,
      HC_prop = p$HC_prop,
      state1 = p$state1,
      K_T = p$K_T,
      K_C = p$K_C
    )
  })

str(sim_results)
summary(sim_results$N)
# inspect 
simulation_summary

# check if trajectories are exploding or collapsing 
sim_results %>%
  group_by(sim_id) %>%
  summarise(
    min_N = min(N, na.rm = TRUE),
    max_N = max(N, na.rm = TRUE),
    final_N = last(N),
    .groups = "drop"
  ) %>%
  summary()





######### 10. Summarise final simulation states #######

### Take the final simulated density from each simulation.
### Then classify whether each simulation ends near endemic or outbreak state.

sim_final <- sim_results %>%
  group_by(sim_id, site_code) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  mutate(
    final_ratio_KC = N / K_C,
    final_ratio_KT = N / K_T,
    ended_above_KT = N >= K_T,
    ended_near_KC = final_ratio_KC >= 0.95,
    ended_near_state1 = N <= state1,
    outcome = case_when(
      ended_near_KC ~ "Carrying capacity",
      ended_above_KT ~ "Above threshold",
      ended_near_state1 ~ "Endemic",
      TRUE ~ "Intermediate"
    )
  )
sim_final %>%
  count(site_code, outcome) %>%
  group_by(site_code) %>%
  mutate(prop = n / sum(n))

simulation_summary <- sim_final %>%
  group_by(site_code) %>%
  summarise(
    n_sims = n(),
    median_final_N = median(N, na.rm = TRUE),
    lower_final_N = quantile(N, 0.025, na.rm = TRUE),
    upper_final_N = quantile(N, 0.975, na.rm = TRUE),
    p_above_KT = mean(ended_above_KT, na.rm = TRUE),
    p_near_state1 = mean(ended_near_state1, na.rm = TRUE),
    median_K_C = median(K_C, na.rm = TRUE),
    .groups = "drop"
  )

simulation_summary


######### 11. Save simulation outputs #######

write_csv(
  params,
  file.path(data_processed_dir, paste0(analysis_date, "_simulation_params.csv"))
)

write_csv(
  sim_results,
  file.path(data_processed_dir, paste0(analysis_date, "_simulation_results.csv"))
)

write_csv(
  sim_final,
  file.path(data_processed_dir, paste0(analysis_date, "_simulation_final.csv"))
)

write_csv(
  simulation_summary,
  file.path(output_dir, "tables", paste0(analysis_date, "_simulation_summary.csv"))
)
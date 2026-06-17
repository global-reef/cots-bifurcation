### 03_SIM.R
### Monte Carlo bifurcation simulation for Koh Tao CoTS dynamics
### Assumes 00_SETUP.R, 01_CLEAN.R, and 02_EXPLORE.R have already been run


######### 1. Simulation Settings #######

set.seed(42)

n_sims <- 10000 # for testing 100 

# lower endemic equilibrium range
KE_min <- 3
KE_max <- 10

# tipping threshold range
KT_min <- 15
KT_max <- 45

# daily scale range, sampled log-uniformly
scale_min <- 0.001
scale_max <- 0.02

# KC must be at least this much above observed N0
KC_buffer <- 1.10


######### 2. Prepare Coral Cover Input Data #######

coral_cover_transect <- substrate_wide %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order)
  ) %>%
  select(site_code, site, site_type, transect, HC) %>%
  mutate(HC_prop = HC / 100) %>%
  filter(
    !is.na(HC_prop),
    HC_prop > 0
  )

coral_cover_summary <- coral_cover_transect %>%
  group_by(site_code, site, site_type) %>%
  summarise(
    n_transects = n(),
    mean_HC = mean(HC, na.rm = TRUE),
    min_HC = min(HC, na.rm = TRUE),
    max_HC = max(HC, na.rm = TRUE),
    .groups = "drop"
  )

save_table(coral_cover_summary, "S1_coral_cover_input_summary.csv")


######### 3. Prepare CoTS Starting Density Input Data #######

cots_starting_density <- cots_survey %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order)
  ) %>%
  filter(
    !is.na(cots_ha),
    cots_ha >= 0,
    !is.na(site_code)
  ) %>%
  select(site_code, site, site_type, survey_id, date, cots_ha)

cots_starting_summary <- cots_starting_density %>%
  group_by(site_code, site, site_type) %>%
  summarise(
    n_surveys = n(),
    mean_cots_ha = mean(cots_ha, na.rm = TRUE),
    median_cots_ha = median(cots_ha, na.rm = TRUE),
    max_cots_ha = max(cots_ha, na.rm = TRUE),
    .groups = "drop"
  )

save_table(cots_starting_summary, "S2_cots_starting_density_summary.csv")


######### 4. Prepare Temperature Forcing #######

r_temp <- function(temp_c) {
  case_when(
    is.na(temp_c) ~ NA_real_,
    temp_c < 26 ~ 0,
    temp_c >= 26 & temp_c < 27 ~ 0.05,
    temp_c >= 27 & temp_c < 28 ~ 0.10,
    temp_c >= 28 & temp_c <= 30 ~ 0.25,
    temp_c > 30 & temp_c <= 31 ~ 0.15,
    temp_c > 31 ~ 0.05,
    TRUE ~ NA_real_
  )
}

temp_forcing <- temp_daily %>%
  filter(!is.na(temp_c)) %>%
  arrange(date) %>%
  select(date, temp_c, temp_source) %>%
  mutate(r_t = r_temp(temp_c))

temp_forcing_summary <- temp_forcing %>%
  summarise(
    n_days = n(),
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    mean_temp_c = mean(temp_c, na.rm = TRUE),
    min_temp_c = min(temp_c, na.rm = TRUE),
    max_temp_c = max(temp_c, na.rm = TRUE),
    mean_r_t = mean(r_t, na.rm = TRUE)
  )

save_table(temp_forcing_summary, "S3_temperature_forcing_summary.csv")


######### 5. Define Bifurcation Functions #######

dN_dt <- function(N, r, KE, KT, KC, scale = 0.01) {
  scale * r * (N - KE) * (1 - N / KC) * (N / KT - 1)
}

simulate_cots <- function(N0, temp_forcing, KE, KT, KC, scale = 0.01) {
  
  N <- numeric(nrow(temp_forcing))
  N[1] <- N0
  
  for (i in 2:nrow(temp_forcing)) {
    
    change_i <- dN_dt(
      N = N[i - 1],
      r = temp_forcing$r_t[i - 1],
      KE = KE,
      KT = KT,
      KC = KC,
      scale = scale
    )
    
    N[i] <- N[i - 1] + change_i
    
    # Keep density within bounded ecological range
    N[i] <- max(0, min(N[i], KC))
  }
  
  N
}


######### 6. Sample Monte Carlo Parameter Sets #######

sample_parameter_set <- function(site_i) {
  
  N0_i <- cots_starting_density %>%
    filter(site_code == site_i) %>%
    slice_sample(n = 1) %>%
    pull(cots_ha)
  
  HC_i <- coral_cover_transect %>%
    filter(site_code == site_i) %>%
    slice_sample(n = 1) %>%
    pull(HC_prop)
  
  # KE = lower endemic equilibrium
  # sampled from a bounded plausible range, not a fake normal mean
  KE_i <- runif(1, min = KE_min, max = KE_max)
  
  # KT = tipping threshold
  # sampled from a bounded plausible range
  KT_i <- runif(1, min = KT_min, max = KT_max)
  
  # scale sampled log-uniformly so small and large rates are both explored
  scale_i <- 10^runif(
    1,
    min = log10(scale_min),
    max = log10(scale_max)
  )
  
  # raw implied upper equilibrium from coral cover relationship
  KC_raw_i <- KT_i / HC_i
  
  # constrain KC so it cannot be lower than observed starting density
  KC_i <- max(KC_raw_i, N0_i * KC_buffer)
  
  tibble(
    site_code = site_i,
    N0 = N0_i,
    HC_prop = HC_i,
    KE = KE_i,
    KT = KT_i,
    KC_raw = KC_raw_i,
    KC = KC_i,
    scale = scale_i,
    KC_constrained = KC_i > KC_raw_i
  )
}
n_sims_per_site <- 2000

sites_to_simulate <- cots_starting_density %>%
  distinct(site_code) %>%
  filter(!is.na(site_code)) %>%
  pull(site_code) %>%
  as.character()

params <- map_dfr(sites_to_simulate, function(site_i) {
  
  map_dfr(seq_len(n_sims_per_site), function(i) {
    sample_parameter_set(site_i)
  })
  
}) %>%
  filter(
    KT > KE,
    KC > KT,
    KC >= N0
  ) %>%
  mutate(
    sim_id = row_number(),
    site_code = factor(site_code, levels = site_code_order),
    basin_start = case_when(
      N0 < KT ~ "Endemic",
      N0 == KT ~ "Threshold",
      N0 > KT ~ "Outbreak basin",
      TRUE ~ NA_character_
    ),
    basin_start = factor(
      basin_start,
      levels = c("Endemic", "Threshold", "Outbreak basin")
    )
  )

param_summary <- params %>%
  group_by(site_code) %>%
  summarise(
    n_sims = n(),
    mean_N0 = mean(N0, na.rm = TRUE),
    max_N0 = max(N0, na.rm = TRUE),
    mean_HC_prop = mean(HC_prop, na.rm = TRUE),
    mean_KE = mean(KE, na.rm = TRUE),
    min_KT = min(KT, na.rm = TRUE),
    median_KT = median(KT, na.rm = TRUE),
    max_KT = max(KT, na.rm = TRUE),
    median_KC_raw = median(KC_raw, na.rm = TRUE),
    median_KC = median(KC, na.rm = TRUE),
    prop_KC_constrained = mean(KC_constrained, na.rm = TRUE),
    median_scale = median(scale, na.rm = TRUE),
    .groups = "drop"
  )

save_table(params, "S4_simulation_params.csv")
save_table(param_summary, "S5_simulation_param_summary.csv")


######### 7. Plot S1: Implied KC by Hard Coral Cover #######

KC_plot_data <- crossing(
  KT = c(10, 15, 20, 30, 40, 45),
  HC_prop = seq(0.2, 0.85, by = 0.01)
) %>%
  mutate(KC = KT / HC_prop)

KC_observed_data <- crossing(
  KT = c(10, 15, 20, 30, 40, 45),
  coral_cover_transect
) %>%
  mutate(KC = KT / HC_prop)

p_KC <- ggplot() +
  geom_line(
    data = KC_plot_data,
    aes(x = HC_prop * 100, y = KC, colour = factor(KT)),
    linewidth = 0.8
  ) +
  geom_point(
    data = KC_observed_data,
    aes(x = HC_prop * 100, y = KC, colour = factor(KT)),
    size = 2,
    alpha = 0.45
  ) +
  labs(
    x = "Hard coral cover (%)",
    y = expression("Implied K"[C]*" (CoTS ha"^-1*")"),
    colour = expression(K[T])
  ) +
  theme_clean

save_plot(p_KC, "S1_implied_KC_by_hard_coral.png", width = 7, height = 5)


######### 8. Plot S2: Starting Basin Classification by Site #######

threshold_site_summary <- params %>%
  count(site_code, basin_start) %>%
  group_by(site_code) %>%
  mutate(prop_sims = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order),
    basin_start = factor(
      basin_start,
      levels = c("Endemic", "Threshold", "Outbreak basin")
    )
  )

save_table(threshold_site_summary, "S6_starting_basin_classification_summary.csv")

p_threshold <- threshold_site_summary %>%
  ggplot(aes(x = site_code, y = prop_sims, fill = basin_start)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = basin_cols, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Site",
    y = "Simulation proportion",
    fill = NULL
  ) +
  theme_clean

save_plot(p_threshold, "S2_starting_basin_classification_by_site.png", width = 7, height = 5)


######### 9. Run Monte Carlo Simulations #######

sim_results <- params %>%
  split(.$sim_id) %>%
  map_dfr(function(p) {
    
    N_sim <- simulate_cots(
      N0 = p$N0,
      temp_forcing = temp_forcing,
      KE = p$KE,
      KT = p$KT,
      KC = p$KC,
      scale = p$scale
    )
    
    tibble(
      sim_id = p$sim_id,
      site_code = p$site_code,
      date = temp_forcing$date,
      temp_c = temp_forcing$temp_c,
      r_t = temp_forcing$r_t,
      N = N_sim,
      N0 = p$N0,
      HC_prop = p$HC_prop,
      KE = p$KE,
      KT = p$KT,
      KC_raw = p$KC_raw,
      KC = p$KC,
      scale = p$scale,
      KC_constrained = p$KC_constrained,
      basin_start = p$basin_start
    )
  }) %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order)
  )

save_table(sim_results, "S7_simulation_results.csv")


######### 10. Summarise Final Simulation States #######
sim_final <- sim_results %>%
  group_by(sim_id, site_code) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  mutate(
    final_ratio_KC = N / KC,
    final_ratio_KT = N / KT,
    ended_above_KT = N >= KT,
    ended_near_KC = final_ratio_KC >= 0.95,
    ended_near_KE = abs(N - KE) <= 1,
    outcome = case_when(
      ended_near_KC ~ "Upper outbreak equilibrium",
      ended_above_KT ~ "Outbreak basin",
      ended_near_KE ~ "Endemic equilibrium",
      N < KT ~ "Endemic basin",
      TRUE ~ NA_character_
    ),
    outcome = factor(
      outcome,
      levels = c(
        "Endemic equilibrium",
        "Endemic basin",
        "Outbreak basin",
        "Upper outbreak equilibrium"
      )
    )
  )

simulation_summary <- sim_final %>%
  group_by(site_code) %>%
  summarise(
    n_sims = n(),
    median_final_N = median(N, na.rm = TRUE),
    lower_final_N = quantile(N, 0.025, na.rm = TRUE),
    upper_final_N = quantile(N, 0.975, na.rm = TRUE),
    p_above_KT = mean(ended_above_KT, na.rm = TRUE),
    p_near_KE = mean(ended_near_KE, na.rm = TRUE),
    p_near_KC = mean(ended_near_KC, na.rm = TRUE),
    median_KC = median(KC, na.rm = TRUE),
    .groups = "drop"
  )

outcome_summary <- sim_final %>%
  count(site_code, outcome) %>%
  group_by(site_code) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

trajectory_check <- sim_results %>%
  group_by(sim_id, site_code) %>%
  summarise(
    min_N = min(N, na.rm = TRUE),
    max_N = max(N, na.rm = TRUE),
    final_N = last(N),
    .groups = "drop"
  )

save_table(sim_final, "S8_simulation_final_states.csv")
save_table(simulation_summary, "S9_simulation_summary.csv")
save_table(outcome_summary, "S10_simulation_outcome_summary.csv")
save_table(trajectory_check, "S11_trajectory_check.csv")


######### 11. Plot S3: Bifurcation Curve Diagnostics #######

bifurcation_scenarios <- crossing(
  KE = c(8, 10),
  KT = c(15, 30),
  HC_prop = c(0.35, 0.50, 0.70),
  scale = 0.01
) %>%
  filter(KT > KE) %>%
  mutate(
    KC = KT / HC_prop,
    scenario = paste0("KE = ", KE, ", KT = ", KT),
    HC_label = paste0("HC = ", round(HC_prop * 100), "%")
  )

bifurcation_curve_data <- bifurcation_scenarios %>%
  pmap_dfr(function(KE, KT, HC_prop, scale, KC, scenario, HC_label) {
    
    tibble(
      KE = KE,
      KT = KT,
      HC_prop = HC_prop,
      scale = scale,
      KC = KC,
      scenario = scenario,
      HC_label = HC_label,
      N = seq(0, KC * 1.2, length.out = 500)
    ) %>%
      mutate(
        dN = dN_dt(
          N = N,
          r = 0.25,
          KE = KE,
          KT = KT,
          KC = KC,
          scale = scale
        )
      )
  })

bifurcation_equilibria <- bifurcation_scenarios %>%
  select(scenario, HC_label, KE, KT, KC) %>%
  pivot_longer(
    cols = c(KE, KT, KC),
    names_to = "equilibrium",
    values_to = "N_eq"
  )

save_table(bifurcation_curve_data, "S12_bifurcation_curve_data.csv")

p_bifurcation <- bifurcation_curve_data %>%
  ggplot(aes(x = N, y = dN, colour = HC_label)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_vline(
    data = bifurcation_equilibria,
    aes(xintercept = N_eq, linetype = equilibrium),
    colour = "grey40",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ scenario, scales = "free_x") +
  labs(
    x = expression("CoTS density, N (individuals ha"^-1*")"),
    y = expression("dN / dt"),
    colour = NULL,
    linetype = "Equilibrium"
  ) +
  theme_clean

save_plot(p_bifurcation, "S3_bifurcation_curve_diagnostics.png", width = 9, height = 6)


######### 12. Plot S4: Monte Carlo Individual Trajectories #######


trajectory_ids_sample <- sim_results %>%
  mutate(site_code = factor(site_code, levels = site_code_order)) %>%
  distinct(site_code, sim_id, N0, KT, KC) %>%
  group_by(site_code) %>%
  slice_sample(prop = 1) %>%
  slice_head(n = 500) %>%
  ungroup() %>%
  mutate(
    starting_basin = case_when(
      N0 < KT ~ "Started below KT",
      N0 >= KT ~ "Started above KT",
      TRUE ~ NA_character_
    ),
    starting_basin = factor(
      starting_basin,
      levels = c("Started below KT", "Started above KT")
    ),
    starts_near_KC = N0 / KC >= 0.95,
    starts_above_KC = N0 > KC
  )

trajectory_plot_sample <- sim_results %>%
  mutate(site_code = factor(site_code, levels = site_code_order)) %>%
  semi_join(
    trajectory_ids_sample %>% select(site_code, sim_id),
    by = c("site_code", "sim_id")
  ) %>%
  left_join(
    trajectory_ids_sample %>%
      select(site_code, sim_id, starting_basin, starts_near_KC, starts_above_KC),
    by = c("site_code", "sim_id")
  ) %>%
  group_by(site_code, sim_id) %>%
  mutate(
    sim_year = as.numeric(date - min(date, na.rm = TRUE)) / 365
  ) %>%
  ungroup()

trajectory_refs <- sim_results %>%
  mutate(site_code = factor(site_code, levels = site_code_order)) %>%
  distinct(sim_id, site_code, KT) %>%
  group_by(site_code) %>%
  summarise(
    median_KT = median(KT, na.rm = TRUE),
    .groups = "drop"
  )

trajectory_labels <- trajectory_refs %>%
  mutate(
    label = paste0("median KT = ", round(median_KT, 1)),
    sim_year = max(trajectory_plot_sample$sim_year, na.rm = TRUE)
  )

trajectory_basin_cols <- c(
  "Started below KT" = "#95B971",
  "Started above KT" = "#FF9683"
)

save_table(trajectory_plot_sample, "S13_sampled_individual_trajectories.csv")
save_table(trajectory_refs, "S14_trajectory_reference_values.csv")

p_trajectory_summary <- ggplot(
  trajectory_plot_sample,
  aes(
    x = sim_year,
    y = N,
    group = sim_id,
    colour = starting_basin
  )
) +
  geom_line(
    alpha = 0.10,
    linewidth = 0.35,
    key_glyph = draw_key_rect
  ) +
  geom_hline(
    data = trajectory_refs,
    aes(yintercept = median_KT),
    inherit.aes = FALSE,
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_label(
    data = trajectory_labels,
    aes(
      x = sim_year,
      y = median_KT,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = 2.8,
    label.size = 0,
    fill = "white",
    alpha = 0.85
  ) +
  facet_wrap(~ site_code, scales = "free_y") +
  scale_colour_manual(
    values = trajectory_basin_cols,
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(
        alpha = 1,
        linewidth = 6
      )
    )
  ) +
  labs(
    x = "Simulation time (years)",
    y = expression("Simulated CoTS density (individuals ha"^-1*")"),
    colour = NULL
  ) +
  theme_clean +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.45, "cm")
  )

save_plot(
  p_trajectory_summary,
  "S4_monte_carlo_individual_trajectories.png",
  width = 9,
  height = 6
)

p_trajectory_summary

######### 13. Plot S5: Final Simulation Outcomes by Site #######

p_outcome <- outcome_summary %>%
  ggplot(aes(x = site_code, y = prop, fill = outcome)) +
  geom_col(
    position = "stack",
    colour = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = outcome_cols, drop = FALSE) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Site",
    y = "Simulation proportion",
    fill = "Final outcome"
  ) +
  theme_clean

p_outcome

save_plot(
  p_outcome,
  "S5_final_simulation_outcomes_by_site.png",
  width = 7,
  height = 5
)

######### 14. Print Main Checks #######

param_summary
threshold_site_summary
simulation_summary
outcome_summary
trajectory_check %>% summary()
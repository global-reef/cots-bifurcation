# ============================================================
# Figure 1: Saddle-node bifurcation curves by site
# ============================================================

# Choose representative parameters
plot_kc_value <- 50
plot_state1 <- 8
plot_sst <- 29

plot_r_T <- calc_r_T(plot_sst)

# Site-level inputs
site_bifurcation_inputs <- live_coral_site %>%
  filter(site %in% target_sites) %>%
  mutate(
    State_1 = plot_state1,
    K_c = plot_kc_value,
    K_t = K_c * live_coral_prop,
    r_T = plot_r_T
  )

# Build density sequence
bifurcation_curves <- site_bifurcation_inputs %>%
  group_by(site) %>%
  summarise(
    State_1 = first(State_1),
    K_t = first(K_t),
    K_c = first(K_c),
    r_T = first(r_T),
    .groups = "drop"
  ) %>%
  crossing(
    N = seq(0, 90, by = 0.1)
  ) %>%
  mutate(
    dN_dt = cots_growth(
      N = N,
      r_T = r_T,
      State_1 = State_1,
      K_t = K_t,
      K_c = K_c
    )
  )

# Equilibrium points
equilibrium_points <- site_bifurcation_inputs %>%
  select(site, State_1, K_t, K_c) %>%
  pivot_longer(
    cols = c(State_1, K_t, K_c),
    names_to = "equilibrium",
    values_to = "N"
  ) %>%
  mutate(
    dN_dt = 0,
    stability = case_when(
      equilibrium %in% c("State_1", "K_c") ~ "stable",
      equilibrium == "K_t" ~ "unstable"
    ),
    equilibrium_label = recode(
      equilibrium,
      State_1 = "State I\nendemic",
      K_t = "State II\ntipping",
      K_c = "State III\noutbreak"
    )
  )

p_saddle_node <- ggplot(bifurcation_curves, aes(x = N, y = dN_dt)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(
    data = equilibrium_points,
    aes(x = N, y = dN_dt, shape = stability),
    size = 3,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = equilibrium_points,
    aes(x = N, y = dN_dt, label = equilibrium_label),
    vjust = -0.8,
    size = 3,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ site, scales = "free_y") +
  labs(
    x = expression("CoTS density, N (ha"^{-1}*")"),
    y = expression("Population growth rate, dN/dt"),
    shape = "Equilibrium stability",
    title = "Saddle-node state-transition model for CoTS outbreak dynamics",
    subtitle = paste0(
      "Curves shown for Kc = ", plot_kc_value,
      " CoTS/ha and SST = ", plot_sst, "°C"
    )
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, paste0(analysis_date, "_fig1_saddle_node_curves.png")),
  p_saddle_node,
  width = 10,
  height = 6,
  dpi = 300
)
p_saddle_node

# ============================================================
# Figure 2: Coral cover shifts the tipping threshold
# ============================================================

coral_threshold_curves <- expand_grid(
  live_coral_prop = seq(0.05, 0.80, by = 0.01),
  K_c = c(25, 50, 75, 100)
) %>%
  mutate(
    K_t = K_c * live_coral_prop
  )

observed_coral_points <- live_coral_site %>%
  filter(site %in% target_sites) %>%
  crossing(K_c = c(25, 50, 75, 100)) %>%
  mutate(
    K_t = K_c * live_coral_prop
  )

p_coral_threshold <- ggplot(coral_threshold_curves,
                            aes(x = live_coral_prop * 100, y = K_t,
                                group = K_c)) +
  geom_line(aes(linetype = factor(K_c)), linewidth = 1) +
  geom_point(
    data = observed_coral_points,
    aes(x = live_coral_prop * 100, y = K_t),
    size = 2.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ K_c) +
  labs(
    x = "Live hard coral cover (%)",
    y = expression("Derived tipping threshold,"~K[t]~"(CoTS"~ha^{-1}*")"),
    title = "Live coral cover shifts the derived CoTS tipping threshold",
    subtitle = expression(K[t] == K[c] %.% "live coral proportion")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(
  file.path(figures_dir, paste0(analysis_date, "_fig2_coral_cover_tipping_threshold.png")),
  p_coral_threshold,
  width = 9,
  height = 6,
  dpi = 300
)
p_coral_threshold

# ============================================================
# Figure 3: SST modifies growth rate around the bifurcation
# ============================================================

# Choose one site or make facets by site
sst_bifurcation_curves <- site_bifurcation_inputs %>%
  select(site, State_1, K_t, K_c) %>%
  crossing(
    SST = c(27, 29, 31, 32),
    N = seq(0, 90, by = 0.1)
  ) %>%
  mutate(
    r_T = calc_r_T(SST),
    dN_dt = cots_growth(
      N = N,
      r_T = r_T,
      State_1 = State_1,
      K_t = K_t,
      K_c = K_c
    )
  )

p_sst_bifurcation <- ggplot(
  sst_bifurcation_curves,
  aes(x = N, y = dN_dt, group = factor(SST), colour = SST)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  facet_wrap(~ site, scales = "free_y") +
  scale_colour_gradient(
    low = "gray",
    high = "red",
    name = "SST (°C)"
  ) +
  labs(
    x = expression("CoTS density, N (ha"^{-1}*")"),
    y = expression("Population growth rate, dN/dt"),
    title = "SST modifies the rate of CoTS outbreak dynamics",
    subtitle = "In this formulation, SST changes the magnitude of growth, not the equilibrium densities"
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, paste0(analysis_date, "_fig3_sst_bifurcation_curves.png")),
  p_sst_bifurcation,
  width = 10,
  height = 6,
  dpi = 300
)
p_sst_bifurcation

# ============================================================
# Bifurcation diagram: coral cover vs CoTS density
# ============================================================

plot_kc_value <- 50
plot_state1 <- 8

bifurcation_diagram <- tibble(
  live_coral_prop = seq(0.05, 0.80, by = 0.005)
) %>%
  mutate(
    live_coral_percent = live_coral_prop * 100,
    State_1 = plot_state1,
    K_t = plot_kc_value * live_coral_prop,
    K_c = plot_kc_value
  ) %>%
  pivot_longer(
    cols = c(State_1, K_t, K_c),
    names_to = "branch",
    values_to = "cots_density_ha"
  ) %>%
  mutate(
    stability = case_when(
      branch %in% c("State_1", "K_c") ~ "stable",
      branch == "K_t" ~ "unstable"
    ),
    branch_label = recode(
      branch,
      State_1 = "State I: endemic",
      K_t = "State II: tipping threshold",
      K_c = "State III: outbreak capacity"
    )
  )

observed_site_points <- site_cots_inputs %>%
  left_join(live_coral_site, by = "site") %>%
  mutate(
    live_coral_percent = live_coral_prop * 100
  )

p_bifurcation_coral <- ggplot(
  bifurcation_diagram,
  aes(
    x = live_coral_percent,
    y = cots_density_ha,
    linetype = stability
  )
) +
  geom_line(aes(group = branch_label), linewidth = 1) +
  geom_point(
    data = observed_site_points,
    aes(
      x = live_coral_percent,
      y = N0,
      shape = N0_type
    ),
    inherit.aes = FALSE,
    size = 2.5,
    alpha = 0.75
  ) +
  labs(
    x = "Live hard coral cover (%)",
    y = expression("CoTS density"~ha^{-1}),
    linetype = "Equilibrium stability",
    shape = "Observed density",
    title = "Bifurcation diagram for CoTS outbreak dynamics",
    subtitle = paste0(
      "Shown for Kc = ", plot_kc_value,
      " CoTS/ha and State I = ", plot_state1,
      " CoTS/ha"
    )
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, paste0(analysis_date, "_bifurcation_diagram_coral_cover.png")),
  p_bifurcation_coral,
  width = 9,
  height = 6,
  dpi = 300
)
p_bifurcation_coral


# ============================================================
# Ashna-style simulation envelope plot
# x = SST, y = CoTS density
# ============================================================

# SST axis for plotting
sst_axis <- tibble(
  SST = seq(26, 33, by = 0.05)
)

# Scenario grid
# These values should eventually be literature-justified.
threshold_scenarios <- crossing(
  site = target_sites,
  K_c = c(25, 50, 75, 100),
  State_1 = c(3, 8, 10),
  sst_sensitivity = c(0.08, 0.12, 0.16, 0.22, 0.28)
) %>%
  left_join(
    live_coral_site %>%
      select(site, live_coral_prop, live_coral_percent),
    by = "site"
  ) %>%
  mutate(
    K_t_base = K_c * live_coral_prop
  ) %>%
  filter(
    State_1 < K_t_base,
    K_t_base < K_c
  )

# Generate one SST-dependent tipping curve per scenario
threshold_curves_sim <- threshold_scenarios %>%
  crossing(sst_axis) %>%
  group_by(site, K_c, State_1, sst_sensitivity, live_coral_prop) %>%
  mutate(
    K_t_SST = K_t_base * exp(-sst_sensitivity * (SST - min(SST))),
    K_t_SST = pmax(K_t_SST, State_1)
  ) %>%
  ungroup()

# Summarise simulation envelope by site and SST
threshold_envelope <- threshold_curves_sim %>%
  group_by(site, SST) %>%
  summarise(
    kt_median = median(K_t_SST, na.rm = TRUE),
    kt_q05 = quantile(K_t_SST, 0.05, na.rm = TRUE),
    kt_q25 = quantile(K_t_SST, 0.25, na.rm = TRUE),
    kt_q75 = quantile(K_t_SST, 0.75, na.rm = TRUE),
    kt_q95 = quantile(K_t_SST, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# Observed density summaries
observed_density_summaries <- site_cots_inputs %>%
  filter(
    N0_type %in% c(
      "median_observed_density",
      "mean_observed_density",
      "q75_observed_density",
      "q90_observed_density",
      "max_observed_density"
    )
  ) %>%
  mutate(
    N0_label = recode(
      N0_type,
      median_observed_density = "Median",
      mean_observed_density = "Mean",
      q75_observed_density = "75th percentile",
      q90_observed_density = "90th percentile",
      max_observed_density = "Maximum"
    ),
    N0_label = factor(
      N0_label,
      levels = c("Median", "Mean", "75th percentile", "90th percentile", "Maximum")
    )
  )

# Optional: observed density band from median to q90
observed_density_band <- site_cots_inputs %>%
  filter(N0_type %in% c("median_observed_density", "q90_observed_density")) %>%
  select(site, N0_type, N0) %>%
  pivot_wider(names_from = N0_type, values_from = N0) %>%
  rename(
    observed_median = median_observed_density,
    observed_q90 = q90_observed_density
  )

# Plot
p_threshold_envelope <- ggplot() +
  
  # Broad simulation uncertainty envelope
  geom_ribbon(
    data = threshold_envelope,
    aes(x = SST, ymin = kt_q05, ymax = kt_q95),
    alpha = 0.12
  ) +
  
  # Central 50% simulation envelope
  geom_ribbon(
    data = threshold_envelope,
    aes(x = SST, ymin = kt_q25, ymax = kt_q75),
    alpha = 0.25
  ) +
  
  # Median simulated tipping threshold
  geom_line(
    data = threshold_envelope,
    aes(x = SST, y = kt_median),
    linewidth = 1.2
  ) +
  
  # Observed median to q90 density band
  geom_rect(
    data = observed_density_band,
    aes(
      xmin = 26,
      xmax = 33,
      ymin = observed_median,
      ymax = observed_q90
    ),
    inherit.aes = FALSE,
    alpha = 0.08
  ) +
  
  # Observed density summaries
  geom_point(
    data = observed_density_summaries,
    aes(x = 33.15, y = N0, shape = N0_label),
    size = 2.8,
    alpha = 0.85
  ) +
  
  facet_wrap(~ site) +
  
  coord_cartesian(
    xlim = c(26, 33.4),
    ylim = c(0, max(observed_density_summaries$N0, threshold_envelope$kt_q95, na.rm = TRUE) * 1.1),
    clip = "off"
  ) +
  
  labs(
    x = "Sea surface temperature (°C)",
    y = expression("CoTS density"~ha^{-1}),
    shape = "Observed density",
    title = "SST-dependent tipping thresholds across simulated parameter scenarios",
    subtitle = "Line = median simulated tipping threshold; ribbons = 50% and 90% simulation envelopes"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.margin = margin(10, 35, 10, 10)
  )

ggsave(
  file.path(figures_dir, paste0(analysis_date, "_sst_threshold_simulation_envelope.png")),
  p_threshold_envelope,
  width = 10,
  height = 6,
  dpi = 300
)
p_threshold_envelope


# ============================================================
# Plot: observed SST vs CoTS density
# ============================================================

p_sst_density <- cots_summary %>%
  filter(site %in% target_sites) %>%
  ggplot(aes(x = SST, y = cots_density_ha)) +
  geom_point(aes(shape = site), size = 3, alpha = 0.75) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  facet_wrap(~ site) +
  labs(
    x = "Sea surface temperature (°C)",
    y = expression("CoTS density"~ha^{-1}),
    shape = "Site",
    title = "Observed CoTS density across sea surface temperature",
    subtitle = "SST values matched to survey dates from Copernicus NetCDF data"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(figures_dir, paste0(analysis_date, "_observed_sst_vs_cots_density.png")),
  p_sst_density,
  width = 10,
  height = 6,
  dpi = 300
)
p_sst_density

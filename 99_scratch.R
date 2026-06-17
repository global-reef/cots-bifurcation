######### 12. Plot S4: Monte Carlo Individual Trajectories #######

### Plot sampled individual trajectories.
### Colour indicates whether the simulation started below or above its own K_T.

set.seed(42)

trajectory_ids_sample <- sim_results %>%
  mutate(site_code = factor(site_code, levels = site_code_order)) %>%
  distinct(site_code, sim_id, N0, K_T, K_C) %>%
  group_by(site_code) %>%
  slice_sample(prop = 1) %>%
  slice_head(n = 175) %>%
  ungroup() %>%
  mutate(
    starting_basin = case_when(
      N0 < K_T ~ "Started below K[T]",
      N0 >= K_T ~ "Started above K[T]",
      TRUE ~ NA_character_
    ),
    starting_basin = factor(
      starting_basin,
      levels = c("Started below K[T]", "Started above K[T]")
    ),
    starts_near_KC = N0 / K_C >= 0.95,
    starts_above_KC = N0 > K_C
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
  distinct(sim_id, site_code, K_T) %>%
  group_by(site_code) %>%
  summarise(
    median_KT = median(K_T, na.rm = TRUE),
    .groups = "drop"
  )

trajectory_labels <- trajectory_refs %>%
  mutate(
    label = paste0("median K[T] = ", round(median_KT, 1)),
    sim_year = max(trajectory_plot_sample$sim_year, na.rm = TRUE)
  )

trajectory_basin_cols <- c(
  "Started below K[T]" = "#95B971",
  "Started above K[T]" = "#FF9683"
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
    linewidth = 0.35
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
  scale_colour_manual(values = trajectory_basin_cols, drop = FALSE) +
  labs(
    x = "Simulation time (years)",
    y = expression("Simulated CoTS density (individuals ha"^-1*")"),
    colour = NULL
  ) +
  theme_clean

p_trajectory_summary

save_plot(
  p_trajectory_summary,
  "S4_monte_carlo_individual_trajectories.png",
  width = 9,
  height = 6
)
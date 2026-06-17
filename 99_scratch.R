######### 12. Plot S4: Monte Carlo Individual Trajectories #######

### Plot sampled individual trajectories only.
### This shows the spread of possible pathways rather than a site-level average.

trajectory_ids_sample <- sim_results %>%
  mutate(site_code = factor(site_code, levels = site_code_order)) %>%
  distinct(site_code, sim_id) %>%
  group_by(site_code) %>%
  slice_sample(prop = 1) %>%
  slice_head(n = 75) %>%
  ungroup()

trajectory_plot_sample <- sim_results %>%
  mutate(site_code = factor(site_code, levels = site_code_order)) %>%
  semi_join(
    trajectory_ids_sample,
    by = c("site_code", "sim_id")
  )

trajectory_refs <- sim_results %>%
  mutate(site_code = factor(site_code, levels = site_code_order)) %>%
  distinct(sim_id, site_code, state1, K_T, K_C) %>%
  group_by(site_code) %>%
  summarise(
    median_KE = median(state1, na.rm = TRUE),
    median_KT = median(K_T, na.rm = TRUE),
    median_KC = median(K_C, na.rm = TRUE),
    .groups = "drop"
  )

trajectory_ref_lines <- trajectory_refs %>%
  pivot_longer(
    cols = c(median_KE, median_KT, median_KC),
    names_to = "ref_type",
    values_to = "ref_value"
  ) %>%
  mutate(
    ref_type = recode(
      ref_type,
      "median_KE" = "K[E]",
      "median_KT" = "K[T]",
      "median_KC" = "K[C]"
    )
  )

trajectory_ref_labels <- trajectory_ref_lines %>%
  group_by(site_code) %>%
  mutate(
    label_x = max(trajectory_plot_sample$N[trajectory_plot_sample$site_code == first(site_code)], na.rm = TRUE) * 0.98,
    label_y = ref_value,
    label = paste0(
      ref_type,
      " = ",
      round(ref_value, 1)
    )
  ) %>%
  ungroup()

save_table(trajectory_plot_sample, "S13_sampled_individual_trajectories.csv")
save_table(trajectory_refs, "S14_trajectory_reference_values.csv")

p_trajectory_summary <- ggplot() +
  geom_line(
    data = trajectory_plot_sample,
    aes(
      x = date,
      y = N,
      group = sim_id
    ),
    alpha = 0.18,
    linewidth = 0.35,
    colour = "#2C7FB8"
  ) +
  geom_hline(
    data = trajectory_ref_lines,
    aes(yintercept = ref_value),
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_label(
    data = trajectory_ref_labels,
    aes(
      x = max(trajectory_plot_sample$date, na.rm = TRUE),
      y = label_y,
      label = label
    ),
    hjust = 1,
    size = 2.8,
    label.size = 0,
    fill = "white",
    alpha = 0.85
  ) +
  facet_wrap(~ site_code, scales = "free_y") +
  labs(
    x = "Date",
    y = expression("Simulated CoTS density (individuals ha"^-1*")")
  ) +
  theme_clean

p_trajectory_summary

save_plot(
  p_trajectory_summary,
  "S4_monte_carlo_individual_trajectories.png",
  width = 9,
  height = 6
)
#### 99_scratch.R ####

# this area is for trying out stuff that may or may not make it into the final analysis. 
sim_results %>%
  mutate(at_KC = abs(N - K_C) < 0.001) %>%
  group_by(sim_id) %>%
  summarise(
    final_N = last(N),
    K_C = first(K_C),
    max_N = max(N, na.rm = TRUE),
    ever_hit_KC = any(at_KC),
    final_at_KC = abs(final_N - K_C) < 0.001,
    .groups = "drop"
  ) %>%
  summarise(
    n_sims = n(),
    p_ever_hit_KC = mean(ever_hit_KC),
    p_final_at_KC = mean(final_at_KC),
    median_final_N = median(final_N),
    median_K_C = median(K_C)
  )
sim_final %>%
  mutate(
    HC_pct = HC_prop * 100
  ) %>%
  arrange(desc(N)) %>%
  select(sim_id, site_code, N, N0, HC_pct, state1, K_T, K_C) %>%
  head(20)



ggplot(
  sim_results %>% filter(sim_id %in% sample(unique(sim_results$sim_id), 100)),
  aes(x = date, y = N, group = sim_id)
) +
  geom_line(alpha = 0.2) +
  facet_wrap(~ site_code) +
  theme_clean +
  labs(
    x = "Date",
    y = "Simulated CoTS density (ha-1)"
  )

sim_ribbon <- sim_results %>%
  group_by(site_code, date) %>%
  summarise(
    median_N = median(N, na.rm = TRUE),
    lower_N = quantile(N, 0.025, na.rm = TRUE),
    upper_N = quantile(N, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(sim_ribbon, aes(x = date, y = median_N)) +
  geom_ribbon(aes(ymin = lower_N, ymax = upper_N), alpha = 0.2) +
  geom_line() +
  facet_wrap(~ site_code) +
  theme_clean +
  labs(
    x = "Date",
    y = "Simulated CoTS density (ha-1)"
  )




# --------------------
  cots_survey %>%
  group_by(site_code) %>%
  summarise(
    n_surveys = n(),
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    mean_cots = mean(cots_ha, na.rm = TRUE),
    sd_cots = sd(cots_ha, na.rm = TRUE),
    min_cots = min(cots_ha, na.rm = TRUE),
    max_cots = max(cots_ha, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(cots_survey, aes(x = date, y = cots_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ site_code) +
  theme_clean

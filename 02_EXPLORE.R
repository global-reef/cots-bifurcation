### 02_EXPLORE.R
### Explore cleaned CoTS, substrate, and SST inputs for bifurcation simulation
### Assumes 00_SETUP.R and 01_CLEAN.R have already been run


######### 1. Check Required Objects #######

required_objects <- c(
  "cots_indiv",
  "cots_survey",
  "substrate",
  "substrate_transect",
  "substrate_wide",
  "temp_daily",
  "theme_clean",
  "save_plot",
  "save_table",
  "site_code_order",
  "site_order",
  "site_type_order",
  "substrate_order",
  "key_substrates",
  "site_type_cols",
  "threshold_cols",
  "density_class_cols"
)

missing_objects <- required_objects[!required_objects %in% ls()]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", "),
    "\nRun 00_SETUP.R and 01_CLEAN.R first."
  )
}


######### 2. Enforce Plot Orders #######

cots_indiv <- cots_indiv %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order),
    site = factor(site, levels = site_order),
    site_type = factor(site_type, levels = site_type_order),
    substrate = factor(as.character(substrate), levels = substrate_order)
  )

cots_survey <- cots_survey %>%
  mutate(
    date = as.Date(date),
    site_code = factor(site_code, levels = site_code_order),
    site = factor(site, levels = site_order),
    site_type = factor(site_type, levels = site_type_order)
  )

substrate <- substrate %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order),
    site = factor(site, levels = site_order),
    site_type = factor(site_type, levels = site_type_order),
    substrate = factor(as.character(substrate), levels = substrate_order)
  )

substrate_transect <- substrate_transect %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order),
    site = factor(site, levels = site_order),
    site_type = factor(site_type, levels = site_type_order),
    substrate = factor(as.character(substrate), levels = substrate_order)
  )

substrate_wide <- substrate_wide %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order),
    site = factor(site, levels = site_order),
    site_type = factor(site_type, levels = site_type_order)
  )

temp_daily <- temp_daily %>%
  mutate(date = as.Date(date))


######### 3. CoTS Survey Summaries #######

cots_site_summary <- cots_survey %>%
  group_by(site_code, site, site_type) %>%
  summarise(
    n_surveys = n(),
    total_cots = sum(cots_count, na.rm = TRUE),
    mean_cots_ha = mean(cots_ha, na.rm = TRUE),
    sd_cots_ha = sd(cots_ha, na.rm = TRUE),
    median_cots_ha = median(cots_ha, na.rm = TRUE),
    min_cots_ha = min(cots_ha, na.rm = TRUE),
    max_cots_ha = max(cots_ha, na.rm = TRUE),
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    .groups = "drop"
  )

save_table(cots_site_summary, "E1_cots_site_summary.csv")

cots_threshold_summary <- cots_survey %>%
  mutate(
    density_class = case_when(
      cots_ha < 10 ~ "Below 10",
      cots_ha >= 10 & cots_ha < 15 ~ "10 to <15",
      cots_ha >= 15 & cots_ha < 30 ~ "15 to <30",
      cots_ha >= 30 & cots_ha < 40 ~ "30 to <40",
      cots_ha >= 40 ~ "40+",
      TRUE ~ NA_character_
    ),
    density_class = factor(
      density_class,
      levels = c("Below 10", "10 to <15", "15 to <30", "30 to <40", "40+")
    )
  ) %>%
  count(site_code, site, site_type, density_class) %>%
  group_by(site_code, site, site_type) %>%
  mutate(prop_surveys = n / sum(n)) %>%
  ungroup()

save_table(cots_threshold_summary, "E2_cots_threshold_summary.csv")


######### 4. Substrate Summaries #######

hc_transect <- substrate_wide %>%
  mutate(HC_prop = HC / 100) %>%
  select(site_code, site, site_type, transect, HC, HC_prop)

save_table(hc_transect, "E3_hard_coral_transect_with_prop.csv")

substrate_site_summary <- substrate_wide %>%
  group_by(site_code, site, site_type) %>%
  summarise(
    n_transects = n(),
    mean_HC = mean(HC, na.rm = TRUE),
    sd_HC = sd(HC, na.rm = TRUE),
    mean_TUR = mean(TUR, na.rm = TRUE),
    sd_TUR = sd(TUR, na.rm = TRUE),
    mean_MAC = mean(MAC, na.rm = TRUE),
    sd_MAC = sd(MAC, na.rm = TRUE),
    mean_AB = mean(AB, na.rm = TRUE),
    sd_AB = sd(AB, na.rm = TRUE),
    mean_OB = mean(OB, na.rm = TRUE),
    sd_OB = sd(OB, na.rm = TRUE),
    .groups = "drop"
  )

save_table(substrate_site_summary, "E4_substrate_site_summary.csv")


######### 5. Temperature Growth Function #######

r_temp <- function(temp_c) {
  case_when(
    is.na(temp_c) ~ NA_real_,
    temp_c < 26 ~ 0,
    temp_c >= 26 & temp_c < 27 ~ 0.05,
    temp_c >= 27 & temp_c < 28 ~ 0.10,
    temp_c >= 28 & temp_c <= 30 ~ 0.25,
    temp_c > 30 & temp_c <= 31 ~ 0.10,
    temp_c > 31 ~ 0,
    TRUE ~ NA_real_
  )
}

temp_daily <- temp_daily %>%
  mutate(r_t = r_temp(temp_c))

temp_summary <- temp_daily %>%
  summarise(
    n_days = n(),
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    min_temp_c = min(temp_c, na.rm = TRUE),
    mean_temp_c = mean(temp_c, na.rm = TRUE),
    max_temp_c = max(temp_c, na.rm = TRUE),
    mean_r_t = mean(r_t, na.rm = TRUE),
    prop_peak_growth = mean(r_t == 0.25, na.rm = TRUE),
    prop_zero_growth = mean(r_t == 0, na.rm = TRUE)
  )

save_table(temp_summary, "E5_temperature_summary.csv")

temp_by_source_summary <- temp_daily %>%
  group_by(temp_source) %>%
  summarise(
    n_days = n(),
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    mean_temp_c = mean(temp_c, na.rm = TRUE),
    sd_temp_c = sd(temp_c, na.rm = TRUE),
    min_temp_c = min(temp_c, na.rm = TRUE),
    max_temp_c = max(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

save_table(temp_by_source_summary, "E6_temperature_by_source_summary.csv")


######### 6. Plot E1: CoTS Density by Site #######

p_cots_site <- ggplot(cots_survey, aes(x = site_code, y = cots_ha, fill = site_type)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(aes(colour = site_type), width = 0.15, height = 0, size = 2, alpha = 0.7) +
  geom_hline(
    aes(yintercept = 15, colour = "15 CoTS ha-1 \n Green Fins (2021)"),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = site_type_cols, drop = FALSE) +
  scale_colour_manual(
    values = c(site_type_cols, threshold_cols),
    drop = FALSE
  ) +
  labs(
    x = "Site",
    y = expression("CoTS density (individuals ha"^-1*")"),
    fill = NULL,
    colour = NULL
  ) +
  theme_clean

save_plot(p_cots_site, "E1_cots_density_by_site.png", width = 7, height = 5)


######### 7. Plot E2: Threshold Classes by Site #######

p_cots_threshold <- cots_threshold_summary %>%
  ggplot(aes(x = site_code, y = prop_surveys, fill = density_class)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = density_class_cols, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Site",
    y = "Survey proportion",
    fill = expression("CoTS density (ha"^-1*")")
  ) +
  theme_clean

save_plot(p_cots_threshold, "E2_cots_threshold_classes_by_site.png", width = 7, height = 5)


######### 8. Plot E3: Hard Coral Cover by Site #######

p_hc_site <- ggplot(substrate_wide, aes(x = site_code, y = HC, fill = site_type)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(aes(colour = site_type), width = 0.15, height = 0, size = 2, alpha = 0.7) +
  scale_fill_manual(values = site_type_cols, drop = FALSE) +
  scale_colour_manual(values = site_type_cols, drop = FALSE) +
  labs(
    x = "Site",
    y = "Hard coral cover (%)",
    fill = NULL,
    colour = NULL
  ) +
  theme_clean

save_plot(p_hc_site, "E3_hard_coral_cover_by_site.png", width = 7, height = 5)


######### 9. Plot E4: Focal Substrate Cover by Site #######

p_substrate <- substrate_transect %>%
  filter(substrate %in% key_substrates) %>%
  mutate(substrate = factor(as.character(substrate), levels = key_substrates)) %>%
  ggplot(aes(x = site_code, y = pct, fill = site_type)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(aes(colour = site_type), width = 0.15, height = 0, size = 1.8, alpha = 0.65) +
  facet_wrap(~ substrate, scales = "free_y") +
  scale_fill_manual(values = site_type_cols, drop = FALSE) +
  scale_colour_manual(values = site_type_cols, drop = FALSE) +
  labs(
    x = "Site",
    y = "Cover (%)",
    fill = NULL,
    colour = NULL
  ) +
  theme_clean

save_plot(p_substrate, "E4_focal_substrate_cover_by_site.png", width = 8, height = 5)


######### 10. Plot E5: Temperature Time Series #######

p_temp <- ggplot(temp_daily, aes(x = date, y = temp_c)) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = c(26, 28, 30), linetype = "dashed", linewidth = 0.4) +
  labs(
    x = "Date",
    y = "Temperature (°C)"
  ) +
  theme_clean

save_plot(p_temp, "E5_temperature_time_series.png", width = 8, height = 5)


######### 11. Save Exploration Objects for Simulation #######

save_table(temp_daily, "E7_temp_daily_with_rt.csv")


######### 12. Print Main Checks #######

cots_site_summary
cots_threshold_summary
substrate_site_summary
temp_summary
temp_by_source_summary
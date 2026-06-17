######### Conceptual Figure: Bifurcation Pathways #######

concept_cols <- c(
  "Endemic basin" = "#95B971",
  "Outbreak basin" = "#FF9683",
  "Threshold" = "#9E2F3D",
  "Reference" = "grey35"
)

KE_y <- 8
KT_y <- 15
KC_y <- 38

concept_levels <- tibble(
  state = c("KE", "KT", "KC"),
  y = c(KE_y, KT_y, KC_y),
  label = c(
    "endemic equilibrium (KE)",
    "unstable tipping threshold (KT)",
    "upper outbreak equilibrium (KC)"
  )
)

concept_curves <- bind_rows(
  
  tibble(
    pathway = "Endemic basin",
    trajectory = "below_1",
    x = seq(0, 6, length.out = 200),
    y = KE_y + (12 - KE_y) * exp(-0.75 * x)
  ),
  
  tibble(
    pathway = "Endemic basin",
    trajectory = "below_2",
    x = seq(0, 6, length.out = 200),
    y = KE_y - (KE_y - 5) * exp(-0.65 * x)
  ),
  
  tibble(
    pathway = "Endemic basin",
    trajectory = "below_3",
    x = seq(0, 6, length.out = 200),
    y = KE_y + (12 - KE_y) * exp(-0.70 * x) -
      1.3 * exp(-0.35 * x) * sin(1.5 * x)
  ),
  
  tibble(
    pathway = "Outbreak basin",
    trajectory = "above_1",
    x = seq(0, 6, length.out = 200),
    y = KC_y - (KC_y - 20) * exp(-0.75 * x)
  ),
  
  tibble(
    pathway = "Outbreak basin",
    trajectory = "above_2",
    x = seq(0, 6, length.out = 200),
    y = KC_y - (KC_y - 28) * exp(-0.55 * x)
  ),
  
  tibble(
    pathway = "Outbreak basin",
    trajectory = "above_3",
    x = seq(0, 6, length.out = 200),
    y = KC_y - (KC_y - 21) * exp(-0.65 * x) +
      3.0 * exp(-0.35 * x) * sin(1.4 * x)
  )
)

concept_points <- concept_curves %>%
  group_by(pathway, trajectory) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    label = case_when(
      pathway == "Endemic basin" ~ "starting density below KT",
      pathway == "Outbreak basin" ~ "starting density above KT"
    )
  )

concept_start_labels <- concept_points %>%
  group_by(pathway) %>%
  slice_head(n = 1) %>%
  ungroup()

concept_domain_labels <- tibble(
  label = c(
    "OUTBREAK BASIN\nN > KT\ntrajectories return toward KC",
    "ENDEMIC BASIN\nN < KT\ntrajectories return toward KE"
  ),
  x = c(3.4, 3.4),
  y = c(43.0, 4.0)
)

concept_curve_labels <- tibble(
  label = c(
    "density above KT trends\ntoward upper outbreak equilibrium",
    "density below KT trends\ntoward endemic equilibrium"
  ),
  x = c(3.1, 3.1),
  y = c(30.5, 11.3)
)

p_concept <- ggplot() +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = KT_y,
    ymax = Inf,
    fill = concept_cols[["Outbreak basin"]],
    alpha = 0.10
  ) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = KT_y,
    fill = concept_cols[["Endemic basin"]],
    alpha = 0.12
  ) +
  geom_hline(
    data = concept_levels,
    aes(yintercept = y),
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey40"
  ) +
  geom_text(
    data = concept_levels %>% filter(state != "KT"),
    aes(x = 6.2, y = y, label = label),
    hjust = 1,
    vjust = -0.45,
    size = 3.0,
    colour = "grey35"
  ) +
  geom_text(
    data = concept_levels %>% filter(state == "KT"),
    aes(x = 6.2, y = y, label = label),
    hjust = 1,
    vjust = -0.45,
    size = 3.1,
    fontface = "bold",
    colour = concept_cols[["Threshold"]]
  ) +
  geom_text(
    data = concept_domain_labels,
    aes(x = x, y = y, label = label),
    size = 3.5,
    fontface = "bold",
    lineheight = 0.95,
    colour = "grey20"
  ) +
  geom_path(
    data = concept_curves,
    aes(
      x = x,
      y = y,
      colour = pathway,
      group = trajectory
    ),
    linewidth = 1.05,
    alpha = 0.85,
    lineend = "round"
  ) +
  geom_point(
    data = concept_points,
    aes(x = x, y = y),
    size = 1,
    colour = concept_cols[["Reference"]]
  ) +
  geom_text(
    data = concept_start_labels,
    aes(x = x, y = y, label = label),
    nudge_x = 0.25,
    size = 3.0,
    hjust = 0
  ) +
  geom_text(
    data = concept_curve_labels,
    aes(x = x, y = y, label = label),
    size = 3.0,
    hjust = 0.5
  ) +
  scale_colour_manual(values = concept_cols[c("Endemic basin", "Outbreak basin")]) +
  scale_x_continuous(
    limits = c(-0.15, 6.6),
    breaks = c(0, 6),
    labels = c("start", "time")
  ) +
  scale_y_continuous(
    limits = c(0, 46),
    breaks = NULL
  ) +
  labs(
    x = "Simulation time",
    y = expression("CoTS density (individuals ha"^-1*")"),
    colour = NULL
  ) +
  theme_clean +
  theme(
    legend.position = "none",
    axis.line.x = element_line(colour = "grey35"),
    axis.line.y = element_line(colour = "grey35"),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

p_concept

save_plot(
  p_concept,
  "S6_conceptual_bifurcation_pathways.png",
  width = 7,
  height = 5
)

### CoTS Bifurcation Project
### Run all scripts to update analysis

######### 0. Set Analysis Date & Create Project Folders #######

analysis_date <- "2026.06.15"  # Update for each analysis run

output_dir <- file.path(getwd(), paste0("Analysis_", analysis_date))

project_dirs <- c(
  "R",
  "data",
  "data/raw",
  "data/processed",
  "docs",
  output_dir,
  file.path(output_dir, "figures"),
  file.path(output_dir, "tables"),
  file.path(output_dir, "model_objects")
)

for (dir in project_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}


######### 1. Load Packages #######
suppressPackageStartupMessages({
library(dplyr);
library(tidyr);
library(ggplot2);
library(readr);
library(stringr);
library(purrr);
library(tibble);
library(lubridate);
library(here);
library(terra);
library(janitor)
})

# call citations with:
citation(package="base")
######### 2. Custom Theme & Colour Palettes #######

theme_clean <- theme_minimal(base_family = "Arial") +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank()
  )

# CoTS / outbreak palette
cots_cols <- c(
  "Endemic" = "#95B971",
  "Elevated" = "#d8c39f",
  "Outbreak" = "#FF9683",
  "Severe outbreak" = "#9E2F3D"
)


# Site-type palette, adjust later if needed
site_cols <- c(
  "Control" = "#95B971",
  "Outbreak" = "#FF9683"
)
threshold_cols <- c(
  "15 CoTS ha-1" = "#B8464F"
)


######### 3. Helpers #######
save_plot <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(output_dir, "figures", filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

save_table <- function(x, filename) {
  write_csv(x, file.path(output_dir, "tables", filename))
}
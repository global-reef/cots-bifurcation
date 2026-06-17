### 00_SETUP.R ###

######### 0. Set Analysis Date & Create Project Folders #######

analysis_date <- "2026.06.17"  # Update for each analysis run

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

######### 2. Custom Theme & Colour Palettes & Forced Orders  #######

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
outcome_cols <- c(
  "Endemic equilibrium" = "#95B971",
  "Endemic basin" = "#d8c39f",
  "Outbreak basin" = "#FF9683",
  "Upper outbreak equilibrium" = "#9E2F3D"
)

basin_cols <- c( 
  "Endemic" = "#95B971",
  "Threshold" = "#d8c39f",
  "Outbreak basin" = "#9E2F3D"
  )
  

# Site-type palette, adjust later if needed
site_type_cols <- c(
  "Control" = "#95B971",
  "High Density" = "#FF9683"
)
threshold_cols <- c(
  "15 CoTS ha-1 \n Green Fins (2021)" = "red"
)

density_class_cols <- c(
  "Below 10" = "#95B971",
  "10 to <15" = "#d8c39f",
  "15 to <30" = "#FF9683",
  "30 to <40" = "#D95F5F",
  "40+" = "#9E2F3D"
)

# Site Orders 
site_code_order <- c("GR", "RR", "TW", "AL", "SI", "TB")

site_order <- c(
  "Green Rock",
  "Red Rock",
  "Twins",
  "Aow Leuk",
  "Shark Island",
  "Tanote Bay"
)

site_type_order <- c("High Density", "Control")

# Substrate orders 
substrate_order <- c("HC", "AB", "TUR", "MAC", "OB", "AN", "UKN")

key_substrates <- c("HC", "TUR", "MAC", "AB", "OB")


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
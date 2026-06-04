### CoTS Bifurcation Project
### Run all scripts to update analysis

######### 0. Set Analysis Date & Create Project Folders #######

analysis_date <- "2026.06.04"  # Update for each analysis run

output_dir <- file.path(getwd(), paste0("Analysis_", analysis_date))

project_dirs <- c(
  "R",
  "data",
  "data/raw",
  "data/processed",
  "docs",
  "outputs",
  "outputs/figures",
  "outputs/tables",
  "outputs/model_objects",
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

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(lubridate)
library(here)


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
  "Endemic" = "#66c2a4",
  "Elevated" = "#41b6c4",
  "Outbreak" = "#2c7fb8",
  "Severe outbreak" = "#253494"
)

# Site-type palette, adjust later if needed
site_cols <- c(
  "Control" = "#66BFA6",
  "Outbreak" = "#007A87"
)

'''
######### 3. Source Helper Functions #######

source(here("R", "functions_cleaning.R"))
source(here("R", "functions_bifurcation.R"))
source(here("R", "functions_simulation.R"))
source(here("R", "functions_plotting.R"))


######### 4. Run Analysis Scripts #######

# 01 clean
source(here("01_CLEAN.R"))

# 01.1 explore
source(here("01.1_EXPLORE.R"))

# 02 parameter setup
source(here("02_PARAMETERS.R"))

# 03 bifurcation model
source(here("03_BIFURCATION_MODEL.R"))

# 04 simulation model
source(here("04_SIMULATION_MODEL.R"))

# 05 sensitivity analysis
source(here("05_SENSITIVITY.R"))

# 06 final plots and tables
source(here("06_PLOTS.R"))'''
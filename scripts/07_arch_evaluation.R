# ============================================================
# 07_arch_evaluation.R
# Archaeobotanical evaluation of EcoCrop suitability
# ============================================================

# This script compares archaeobotanical crop occurrences from southern
# Africa with EcoCrop suitability values extracted from matching crop-
# and time-specific suitability rasters.
#
# Data access note:
# The raw archaeobotanical workbook is not included in this repository
# because the dataset has not yet been openly published/licensed.
# It can be made available upon reasonable request to the corresponding
# author of the paper. To reproduce this script, place the restricted
# workbook at:
#   data/restricted/archaeobotanical_raw_SA.xlsx

# ---------- packages ----------
req_pkgs <- c(
  "readxl", "writexl", "dplyr", "tidyr", "terra",
  "sf", "ggplot2", "patchwork", "grid"
)

missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Install required packages before running this script: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(terra)
library(sf)
library(ggplot2)
library(patchwork)
library(grid)

# ---------- paths ----------
# Run this script from the repository root.

project_dir <- getwd()

raw_arch_file <- file.path(
  project_dir,
  "data", "restricted", "archaeobotanical_raw_SA.xlsx"
)

region_dir <- file.path(project_dir, "data", "raw", "regions")
sa_shape_file <- file.path(region_dir, "southern_africa.gpkg")

if (!file.exists(sa_shape_file)) {
  sa_shape_file <- file.path(region_dir, "Southern Africa Merged.shp")
}

raster_dir <- file.path(project_dir, "data", "processed", "ecocrop_rasters")

processed_dir <- file.path(project_dir, "data", "processed", "archaeobotanical")
fig_dir <- file.path(project_dir, "outputs", "figures")
supp_fig_dir <- file.path(project_dir, "outputs", "figures", "supplement")
table_dir <- file.path(project_dir, "outputs", "tables")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_arch_file)) {
  stop(
    "Restricted archaeobotanical workbook not found.\n",
    "This raw dataset is not distributed with the repository.\n",
    "Place the file at data/restricted/archaeobotanical_raw_SA.xlsx ",
    "after obtaining it from the corresponding author."
  )
}

if (!file.exists(sa_shape_file)) {
  stop(
    "Southern Africa boundary file not found. Expected either:\n",
    file.path(region_dir, "southern_africa.gpkg"), "\n",
    "or\n",
    file.path(region_dir, "Southern Africa Merged.shp")
  )
}

# ---------- settings ----------
timesteps <- seq(0, 10000, 100)
set.seed(123)
n_background <- 300

crop_cols <- c(
  "Sorghum" = "#2166AC",
  "Pearl millet" = "#1B9E77",
  "Finger millet" = "#D95F02"
)

crop_labs <- c(
  sorghum = "Sorghum",
  pearl_millet = "Pearl millet",
  finger_millet = "Finger millet"
)

crop_codes <- c(
  sorghum = "SG",
  pearl_millet = "PM",
  finger_millet = "FM"
)

# ============================================================
# 1. Prepare archaeobotanical occurrence table
# ============================================================

raw_arch <- read_excel(raw_arch_file)

arch_df <- raw_arch %>%
  mutate(
    lat = as.numeric(gsub("−", "-", `x-coordinates`)),
    lon = as.numeric(`y-coordinates`),
    year_bp = suppressWarnings(as.numeric(`Date cal BP`)),
    sorghum = ifelse(
      !is.na(`Sorghum bicolor`) | !is.na(`Sorghum caffra`),
      1, NA
    ),
    finger_millet = ifelse(
      !is.na(`Eleusine corocana (finger millet)`),
      1, NA
    ),
    pearl_millet = ifelse(
      !is.na(`Cenchrus americanus (pearl millet)`) | !is.na(`Pennisetum glaucum`),
      1, NA
    )
  ) %>%
  select(
    `Site Name`, Country, lon, lat, year_bp,
    sorghum, finger_millet, pearl_millet
  ) %>%
  pivot_longer(
    cols = c(sorghum, finger_millet, pearl_millet),
    names_to = "crop",
    values_to = "presence"
  ) %>%
  filter(
    !is.na(presence),
    !is.na(year_bp),
    !is.na(lon),
    !is.na(lat)
  ) %>%
  rowwise() %>%
  mutate(model_timestep = timesteps[which.min(abs(timesteps - year_bp))]) %>%
  ungroup()

write_xlsx(
  arch_df,
  file.path(processed_dir, "arch_occurrences.xlsx")
)

# ============================================================
# 2. Extract EcoCrop suitability at archaeobotanical sites
# ============================================================

make_raster_file <- function(crop, model_timestep) {
  crop_code <- crop_codes[[crop]]
  raster_index <- model_timestep / 100

  file.path(
    raster_dir,
    crop_code,
    paste0(crop_code, "_suit_t", raster_index, ".tif")
  )
}

add_raster_links <- function(df) {
  df %>%
    mutate(
      record_id = paste0("rec_", seq_len(n())),
      site_id = as.integer(factor(`Site Name`)),
      point_type = "archaeo",
      crop_code = unname(crop_codes[crop]),
      raster_index = model_timestep / 100,
      raster_file = mapply(make_raster_file, crop, model_timestep)
    )
}

extract_suitability <- function(df) {

  out <- df
  out$suitability <- NA_real_

  for (i in seq_len(nrow(out))) {

    if (!file.exists(out$raster_file[i])) {
      warning("Missing raster: ", out$raster_file[i])
      next
    }

    r <- rast(out$raster_file[i])

    p <- vect(
      data.frame(x = out$lon[i], y = out$lat[i]),
      geom = c("x", "y"),
      crs = "EPSG:4326"
    )

    out$suitability[i] <- extract(r, p)[1, 2]
  }

  out
}

add_suitability_classes <- function(df) {
  df %>%
    mutate(
      suit_class = case_when(
        suitability < 0.3 ~ "low",
        suitability <= 0.6 ~ "moderate",
        suitability > 0.6 ~ "high",
        TRUE ~ NA_character_
      ),
      is_low = ifelse(!is.na(suitability) & suitability < 0.3, 1, 0),
      is_moderate = ifelse(!is.na(suitability) & suitability >= 0.3 & suitability <= 0.6, 1, 0),
      is_high = ifelse(!is.na(suitability) & suitability > 0.6, 1, 0),
      dist_from_optimum = ifelse(!is.na(suitability), 1 - suitability, NA_real_)
    )
}

arch_modelled <- arch_df %>%
  add_raster_links() %>%
  extract_suitability() %>%
  add_suitability_classes()

write_xlsx(
  arch_modelled,
  file.path(processed_dir, "arch_modelled.xlsx")
)

# ============================================================
# 3. Create crop- and time-matched background sample
# ============================================================

sa_shape <- st_read(sa_shape_file, quiet = TRUE)

if (is.na(st_crs(sa_shape))) {
  stop("Southern Africa boundary has no CRS.")
}

sa_shape <- st_transform(sa_shape, 4326)

bg_points <- st_sample(sa_shape, size = n_background)

bg_df <- st_as_sf(bg_points) %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  mutate(
    crop = sample(arch_modelled$crop, n(), replace = TRUE),
    model_timestep = sample(arch_modelled$model_timestep, n(), replace = TRUE),
    point_type = "background",
    `Site Name` = NA_character_,
    Country = NA_character_,
    year_bp = model_timestep,
    presence = NA_real_,
    record_id = paste0("bg_", seq_len(n())),
    site_id = NA_integer_,
    crop_code = unname(crop_codes[crop]),
    raster_index = model_timestep / 100,
    raster_file = mapply(make_raster_file, crop, model_timestep)
  ) %>%
  extract_suitability() %>%
  add_suitability_classes()

bg_df <- bg_df %>%
  select(names(arch_modelled))

master_df <- bind_rows(arch_modelled, bg_df)

write_xlsx(
  master_df,
  file.path(processed_dir, "arch_ecocrop_master.xlsx")
)

# ============================================================
# 4. Summary statistics and effect-size table
# ============================================================

plot_df <- master_df %>%
  mutate(
    crop = recode(crop, !!!crop_labs),
    crop = factor(crop, levels = names(crop_cols)),
    point_type = recode(
      point_type,
      archaeo = "Archaeobotanical sites",
      background = "Potential niche"
    ),
    point_type = factor(
      point_type,
      levels = c("Archaeobotanical sites", "Potential niche")
    )
  )

tab_main <- plot_df %>%
  group_by(point_type, crop) %>%
  summarise(
    n = n(),
    mean_suitability = mean(suitability, na.rm = TRUE),
    median_suitability = median(suitability, na.rm = TRUE),
    min_suitability = min(suitability, na.rm = TRUE),
    max_suitability = max(suitability, na.rm = TRUE),
    low_n = sum(suit_class == "low", na.rm = TRUE),
    moderate_n = sum(suit_class == "moderate", na.rm = TRUE),
    high_n = sum(suit_class == "high", na.rm = TRUE),
    low_pct = 100 * mean(suit_class == "low", na.rm = TRUE),
    moderate_pct = 100 * mean(suit_class == "moderate", na.rm = TRUE),
    high_pct = 100 * mean(suit_class == "high", na.rm = TRUE),
    .groups = "drop"
  )

boot_diff <- function(dat, nboot = 2000) {

  arch <- dat$suitability[dat$point_type == "Archaeobotanical sites"]
  pot <- dat$suitability[dat$point_type == "Potential niche"]

  sims <- replicate(
    nboot,
    mean(sample(arch, replace = TRUE), na.rm = TRUE) -
      mean(sample(pot, replace = TRUE), na.rm = TRUE)
  )

  c(
    archaeo_mean = mean(arch, na.rm = TRUE),
    potential_mean = mean(pot, na.rm = TRUE),
    diff = mean(arch, na.rm = TRUE) - mean(pot, na.rm = TRUE),
    ci_low = unname(quantile(sims, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(sims, 0.975, na.rm = TRUE)),
    p_value = wilcox.test(suitability ~ point_type, data = dat)$p.value
  )
}

set.seed(123)

tab_effect <- plot_df %>%
  group_by(crop) %>%
  group_modify(~ {
    out <- boot_diff(.x)

    data.frame(
      archaeo_mean = out["archaeo_mean"],
      potential_mean = out["potential_mean"],
      mean_difference = out["diff"],
      ci_low = out["ci_low"],
      ci_high = out["ci_high"],
      p_value = out["p_value"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_value_3dp = sprintf("%.3f", p_value),
    archaeo_mean = round(archaeo_mean, 3),
    potential_mean = round(potential_mean, 3),
    mean_difference = round(mean_difference, 3),
    ci_low = round(ci_low, 3),
    ci_high = round(ci_high, 3)
  )

tab_low_records <- plot_df %>%
  filter(point_type == "Archaeobotanical sites", suit_class == "low") %>%
  select(
    `Site Name`, Country, crop, year_bp,
    model_timestep, suitability, raster_file
  )

write_xlsx(
  list(
    main_summary = tab_main,
    effect_sizes = tab_effect,
    low_suitability_records = tab_low_records
  ),
  file.path(table_dir, "archaeobotanical_evaluation.xlsx")
)

write.csv(
  tab_main,
  file.path(table_dir, "arch_summary.csv"),
  row.names = FALSE
)

write.csv(
  tab_effect,
  file.path(table_dir, "arch_effects.csv"),
  row.names = FALSE
)

# ============================================================
# 5. Figures
# ============================================================

map_df <- master_df %>%
  filter(point_type == "archaeo") %>%
  arrange(desc(suitability)) %>%
  mutate(
    crop = factor(
      recode(crop, !!!crop_labs),
      levels = names(crop_cols)
    ),
    lon_j = jitter(lon, amount = 0.35),
    lat_j = jitter(lat, amount = 0.35)
  )

map_sf <- st_as_sf(map_df, coords = c("lon_j", "lat_j"), crs = 4326)

map_labs <- c(
  "Sorghum" = paste0(
    "Sorghum (n = ",
    sum(map_df$crop == "Sorghum"),
    ")"
  ),
  "Pearl millet" = paste0(
    "Pearl millet (n = ",
    sum(map_df$crop == "Pearl millet"),
    ")"
  ),
  "Finger millet" = paste0(
    "Finger millet (n = ",
    sum(map_df$crop == "Finger millet"),
    ")"
  )
)

p_map <- ggplot() +
  geom_sf(
    data = sa_shape,
    fill = "#F7F4EA",
    color = "#B8B8B8",
    linewidth = 0.3
  ) +
  geom_sf(
    data = map_sf,
    aes(color = crop, size = suitability),
    alpha = 0.9
  ) +
  scale_color_manual(values = crop_cols, labels = map_labs) +
  scale_size_continuous(
    name = "Modelled suitability",
    range = c(2.5, 8.5)
  ) +
  coord_sf(xlim = c(12, 41), ylim = c(-36, -8), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", color = "Crop") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

fig_df <- master_df %>%
  mutate(
    crop = recode(crop, !!!crop_labs),
    crop = factor(crop, levels = names(crop_cols)),
    site_type = ifelse(
      point_type == "archaeo",
      "Archaeobotanical sites",
      "Regional climate availability"
    ),
    x_pos = ifelse(site_type == "Archaeobotanical sites", 1, 2)
  )

p_df <- fig_df %>%
  group_by(crop) %>%
  summarise(
    p = wilcox.test(suitability ~ point_type)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_lab = paste0("p = ", formatC(p, format = "f", digits = 3)))

p_compare <- ggplot(fig_df, aes(x = x_pos, y = suitability)) +
  geom_violin(
    data = subset(fig_df, point_type == "background"),
    fill = "#F7F4EA",
    color = NA,
    alpha = 1,
    width = 0.7
  ) +
  geom_violin(
    data = subset(fig_df, point_type == "archaeo"),
    aes(fill = crop),
    color = NA,
    alpha = 0.50,
    width = 0.7
  ) +
  geom_jitter(
    data = subset(fig_df, point_type == "archaeo"),
    aes(color = crop),
    width = 0.06,
    size = 2.6,
    alpha = 0.95
  ) +
  stat_summary(
    data = subset(fig_df, point_type == "archaeo"),
    aes(color = crop),
    fun = median,
    geom = "crossbar",
    width = 0.24,
    fatten = 0,
    linewidth = 1.3
  ) +
  stat_summary(
    data = subset(fig_df, point_type == "background"),
    fun = median,
    geom = "crossbar",
    width = 0.24,
    fatten = 0,
    linewidth = 1.15,
    color = "#5A5A5A"
  ) +
  geom_text(
    data = p_df,
    aes(x = 1.5, y = 1.02, label = p_lab),
    inherit.aes = FALSE,
    size = 3.6,
    color = "black"
  ) +
  facet_wrap(~ crop, nrow = 1) +
  scale_fill_manual(values = crop_cols) +
  scale_color_manual(values = crop_cols) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Realized niche", "Potential niche")
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(x = NULL, y = "Modelled suitability") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.25),
    strip.text = element_text(face = "bold", color = "black", size = 12),
    strip.background = element_blank(),
    panel.spacing.x = grid::unit(1.2, "lines"),
    panel.border = element_rect(color = "#D9D9D9", fill = NA, linewidth = 0.5),
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

fig_arch <- p_map / p_compare +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(fig_dir, "Fig_archaeobotanical_evaluation.png"),
  fig_arch,
  width = 11,
  height = 11,
  dpi = 500,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "Fig_archaeobotanical_evaluation.pdf"),
  fig_arch,
  width = 11,
  height = 11,
  bg = "white"
)

p_ecdf <- ggplot(plot_df, aes(x = suitability, color = point_type)) +
  stat_ecdf(linewidth = 1) +
  facet_wrap(~ crop, nrow = 1) +
  scale_color_manual(values = c(
    "Archaeobotanical sites" = "#2166AC",
    "Potential niche" = "#9E9E9E"
  )) +
  labs(
    x = "Modelled suitability",
    y = "Cumulative proportion",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.25),
    strip.text = element_text(face = "bold", color = "black"),
    legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(supp_fig_dir, "SM_arch_ECDF.png"),
  p_ecdf,
  width = 8.5,
  height = 3.6,
  dpi = 500,
  bg = "white"
)

effect_plot_df <- tab_effect %>%
  mutate(p_lab = paste0("p = ", p_value_3dp))

p_effect <- ggplot(
  effect_plot_df,
  aes(x = crop, y = mean_difference, color = crop)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#7F7F7F",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.08,
    linewidth = 0.8
  ) +
  geom_point(size = 4) +
  geom_text(
    aes(y = ci_high + 0.04, label = p_lab),
    color = "black",
    size = 3.6
  ) +
  scale_color_manual(values = crop_cols) +
  labs(
    x = NULL,
    y = "Difference in mean suitability\n(archaeobotanical - potential niche)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.25),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(supp_fig_dir, "SM_arch_effect_size.png"),
  p_effect,
  width = 6.5,
  height = 4.5,
  dpi = 500,
  bg = "white"
)

# ============================================================
# 08_opportunity_space_figures.R
# Climate opportunity-space figures from EcoCrop suitability rasters
# ============================================================

# This script creates the main Late Holocene opportunity-space figure,
# supplementary Early/Mid-Holocene opportunity-space maps, a quantitative
# opportunity-space summary figure, and optional regional trajectory figure.

# ---------- packages ----------
req_pkgs <- c(
  "terra", "sf", "dplyr", "tidyr", "ggplot2", "patchwork",
  "tibble", "forcats", "scales", "readxl", "cowplot"
)

missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(forcats)
library(scales)
library(readxl)
library(cowplot)

# ---------- paths ----------
# Run this script from the repository root.

project_dir <- getwd()

raster_dir <- file.path(project_dir, "data", "processed", "ecocrop_rasters")
region_file <- file.path(project_dir, "data", "raw", "regions", "africa_UN_regions.shp")

suitability_table <- file.path(
  project_dir,
  "data", "processed", "ecocrop_suitability.xlsx"
)

fig_dir <- file.path(project_dir, "outputs", "figures")
supp_fig_dir <- file.path(project_dir, "outputs", "figures", "supplement")
table_dir <- file.path(project_dir, "outputs", "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(region_file)) {
  stop("Africa region shapefile not found: ", region_file)
}

if (!dir.exists(raster_dir)) {
  stop("EcoCrop raster folder not found: ", raster_dir)
}

# ---------- settings ----------
crop_codes <- c("SG", "PM", "FM")

crop_cols <- c(
  SG = "#2166AC",
  PM = "#1B9E77",
  FM = "#D95F02"
)

crop_cols_named <- c(
  "Sorghum" = "#2166AC",
  "Pearl millet" = "#1B9E77",
  "Finger millet" = "#D95F02"
)

breadth_cols <- c(
  `0` = "#f6e4e1",
  `1` = "#d6604d",
  `2` = "#f4a582",
  `3` = "#fee08b"
)

agg_fact <- 3
min_floor <- 0.30
plausible_thr <- 0.50

phase_order <- c(
  "Early Holocene (11.7–8.2 ka BP)",
  "Mid-Holocene (8.2–4.2 ka BP)",
  "Late Holocene (4.2–0 ka BP)"
)

region_order <- c("Northern", "Western", "Central", "Eastern", "Southern")

# ============================================================
# 1. Read regions and available timesteps
# ============================================================

regions_sf <- st_read(region_file, quiet = TRUE)

if ("region5" %in% names(regions_sf)) {
  regions_sf <- regions_sf |>
    mutate(region = recode(
      region5,
      "North" = "Northern",
      "West" = "Western",
      "East" = "Eastern",
      "Central" = "Central",
      "South" = "Southern",
      .default = as.character(region5)
    ))
} else if (!"region" %in% names(regions_sf)) {
  names(regions_sf)[1] <- "region"
}

africa_outline <- st_as_sf(st_union(regions_sf))

get_ts <- function(cr) {
  fs <- list.files(
    file.path(raster_dir, cr),
    pattern = paste0("^", cr, "_suit_t[0-9]+\\.tif$"),
    full.names = FALSE
  )

  as.integer(
    sub(paste0("^", cr, "_suit_t([0-9]+)\\.tif$"), "\\1", fs)
  )
}

ts <- Reduce(intersect, lapply(crop_codes, get_ts))
ts <- sort(ts, decreasing = TRUE)

if (length(ts) == 0) {
  stop("No common EcoCrop suitability timesteps found for SG, PM, and FM.")
}

ts_df <- tibble(
  timestep = ts,
  kaBP = ts / 10
)

# ============================================================
# 2. Helper functions
# ============================================================

fmt_lon <- function(x) {
  ifelse(
    x < 0, paste0(abs(x), "°W"),
    ifelse(x > 0, paste0(x, "°E"), "0°")
  )
}

fmt_lat <- function(y) {
  ifelse(
    y < 0, paste0(abs(y), "°S"),
    ifelse(y > 0, paste0(y, "°N"), "0°")
  )
}

load_phase_mean_stack <- function(ka_min, ka_max, agg_fact = 1) {

  ts_phase <- ts_df |>
    filter(kaBP >= ka_min, kaBP <= ka_max) |>
    arrange(desc(kaBP))

  if (nrow(ts_phase) == 0) stop("No timesteps found for phase.")

  out <- vector("list", length(crop_codes))
  names(out) <- crop_codes

  for (cr in crop_codes) {

    fs <- file.path(
      raster_dir,
      cr,
      paste0(cr, "_suit_t", ts_phase$timestep, ".tif")
    )

    fs <- fs[file.exists(fs)]

    if (length(fs) == 0) {
      stop("No raster files found for crop ", cr, " in selected phase.")
    }

    rs <- lapply(fs, function(f) {
      r <- rast(f)
      if (agg_fact > 1) {
        r <- aggregate(r, fact = agg_fact, fun = mean, na.rm = TRUE)
      }
      r
    })

    r_mean <- app(rast(rs), mean, na.rm = TRUE)
    names(r_mean) <- cr
    out[[cr]] <- r_mean
  }

  s <- c(out$SG, out$PM, out$FM)
  names(s) <- crop_codes
  s
}

compute_metric_rasters <- function(r_stack, min_floor = 0.30, plausible_thr = 0.50) {

  max_val <- app(r_stack, max)

  second_val <- app(r_stack, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sort(x, decreasing = TRUE)[2]
  })

  top_idx <- app(r_stack, function(x) {
    if (all(is.na(x))) return(NA_real_)
    which.max(x)
  })

  strength <- (max_val - second_val) / max_val

  n_options <- (r_stack[[1]] >= plausible_thr) +
    (r_stack[[2]] >= plausible_thr) +
    (r_stack[[3]] >= plausible_thr)

  top_idx <- ifel(max_val < min_floor, NA, top_idx)
  strength <- ifel(max_val < min_floor, NA, strength)
  n_options <- ifel(is.na(max_val), NA, n_options)
  n_options <- ifel(max_val < min_floor, 0, n_options)

  list(
    top_idx = top_idx,
    strength = strength,
    max_val = max_val,
    second_val = second_val,
    n_options = n_options
  )
}

make_phase_data <- function(label, ka_min, ka_max) {

  r_stack <- load_phase_mean_stack(ka_min, ka_max, agg_fact = agg_fact)

  regions_use <- regions_sf
  outline_use <- africa_outline

  if (!is.na(crs(r_stack)) && !is.na(st_crs(regions_sf)$wkt)) {
    if (crs(r_stack) != st_crs(regions_sf)$wkt) {
      regions_use <- st_transform(regions_sf, crs(r_stack))
      outline_use <- st_transform(africa_outline, crs(r_stack))
    }
  }

  met <- compute_metric_rasters(
    r_stack,
    min_floor = min_floor,
    plausible_thr = plausible_thr
  )

  lead_df <- as.data.frame(c(met$top_idx, met$strength), xy = TRUE, na.rm = TRUE)
  names(lead_df) <- c("x", "y", "top_idx", "strength")

  lead_df <- lead_df |>
    mutate(
      crop = crop_codes[top_idx],
      alpha_val = pmax(0.20, pmin(1, strength^0.85)),
      panel = "Leading climatic option",
      phase = label
    )

  breadth_df <- as.data.frame(met$n_options, xy = TRUE, na.rm = FALSE)
  names(breadth_df) <- c("x", "y", "n_options")

  breadth_df <- breadth_df |>
    filter(!is.na(n_options)) |>
    mutate(
      n_options = factor(n_options, levels = c(0, 1, 2, 3)),
      panel = "Climatic choice breadth",
      phase = label
    )

  region_r <- rasterize(vect(regions_use), r_stack[[1]], field = "region")
  area_r <- cellSize(r_stack[[1]], unit = "km")

  sum_df <- as.data.frame(
    c(region_r, area_r, met$top_idx, met$strength, met$n_options),
    xy = FALSE,
    na.rm = FALSE
  )

  names(sum_df) <- c("region", "cell_km2", "top_idx", "strength", "n_options")

  sum_df <- sum_df |>
    filter(!is.na(region), !is.na(cell_km2)) |>
    mutate(
      region = as.character(region),
      crop = case_when(
        top_idx == 1 ~ "Sorghum",
        top_idx == 2 ~ "Pearl millet",
        top_idx == 3 ~ "Finger millet",
        TRUE ~ NA_character_
      ),
      phase = label
    )

  crop_frac <- sum_df |>
    filter(!is.na(crop)) |>
    group_by(phase, region, crop) |>
    summarise(
      lead_area_km2 = sum(cell_km2, na.rm = TRUE),
      .groups = "drop"
    ) |>
    group_by(phase, region) |>
    mutate(frac_lead_area = lead_area_km2 / sum(lead_area_km2)) |>
    ungroup()

  region_stats <- sum_df |>
    group_by(phase, region) |>
    summarise(
      mean_strength = weighted.mean(strength, w = cell_km2, na.rm = TRUE),
      mean_n_options = weighted.mean(n_options, w = cell_km2, na.rm = TRUE),
      area_any_option_km2 = sum(cell_km2[n_options >= 1], na.rm = TRUE),
      area_two_plus_km2 = sum(cell_km2[n_options >= 2], na.rm = TRUE),
      area_three_km2 = sum(cell_km2[n_options >= 3], na.rm = TRUE),
      .groups = "drop"
    )

  list(
    r_stack = r_stack,
    regions = regions_use,
    outline = outline_use,
    lead_df = lead_df,
    breadth_df = breadth_df,
    crop_frac = crop_frac,
    region_stats = region_stats
  )
}

make_legend_bar <- function(fills, labels, title_text) {

  df <- tibble(
    x = seq_along(fills),
    fill = fills,
    lab = labels
  )

  ggplot(df, aes(x, 1, fill = fill)) +
    geom_tile(height = 0.38, width = 0.98, color = NA) +
    geom_text(aes(y = 0.60, label = lab), size = 3.0, fontface = "bold") +
    annotate(
      "text",
      x = mean(df$x),
      y = 1.40,
      label = title_text,
      fontface = "bold",
      size = 3.2
    ) +
    scale_fill_identity() +
    xlim(0.5, nrow(df) + 0.5) +
    ylim(0.45, 1.55) +
    theme_void() +
    theme(plot.margin = margin(0, 8, 0, 8))
}

theme_map_box <- function() {
  theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 8.5, face = "bold", colour = "black"),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5, colour = "black"),
      strip.text = element_text(size = 10, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.9),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.95),
      plot.margin = margin(4, 4, 4, 4)
    )
}

make_map_pair <- function(res, overall_title = NULL) {

  bb <- st_bbox(res$outline)
  xpad <- (bb$xmax - bb$xmin) * 0.03
  ypad <- (bb$ymax - bb$ymin) * 0.03

  xlim_use <- c(bb$xmin - xpad, bb$xmax + xpad)
  ylim_use <- c(bb$ymin - ypad, bb$ymax + ypad)

  p_lead <- ggplot() +
    geom_sf(data = res$outline, fill = "grey96", colour = NA) +
    geom_raster(
      data = res$lead_df,
      aes(x = x, y = y, fill = crop, alpha = alpha_val)
    ) +
    geom_sf(data = res$regions, fill = NA, colour = "grey35", linewidth = 0.82) +
    geom_sf(data = res$outline, fill = NA, colour = "black", linewidth = 1.05) +
    scale_fill_manual(values = crop_cols, guide = "none") +
    scale_alpha(range = c(0.20, 1), guide = "none") +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    facet_wrap(~panel) +
    scale_x_continuous(
      breaks = c(-10, 0, 10, 20, 30, 40, 50),
      labels = fmt_lon
    ) +
    scale_y_continuous(
      breaks = c(-30, -20, -10, 0, 10, 20, 30),
      labels = fmt_lat
    ) +
    theme_map_box()

  p_breadth <- ggplot() +
    geom_sf(data = res$outline, fill = "grey96", colour = NA) +
    geom_raster(
      data = res$breadth_df,
      aes(x = x, y = y, fill = n_options)
    ) +
    geom_sf(data = res$regions, fill = NA, colour = "grey35", linewidth = 0.82) +
    geom_sf(data = res$outline, fill = NA, colour = "black", linewidth = 1.05) +
    scale_fill_manual(values = breadth_cols, drop = FALSE, guide = "none") +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    facet_wrap(~panel) +
    scale_x_continuous(
      breaks = c(-10, 0, 10, 20, 30, 40, 50),
      labels = fmt_lon
    ) +
    scale_y_continuous(
      breaks = c(-30, -20, -10, 0, 10, 20, 30),
      labels = fmt_lat
    ) +
    theme_map_box()

  leg_lead <- make_legend_bar(
    fills = unname(crop_cols[c("SG", "PM", "FM")]),
    labels = c("Sorghum", "Pearl millet", "Finger millet"),
    title_text = "Leading climatic option"
  )

  leg_breadth <- make_legend_bar(
    fills = unname(breadth_cols[c("0", "1", "2", "3")]),
    labels = c("n = 0", "n = 1", "n = 2", "n = 3"),
    title_text = "Number of plausible crop options"
  )

  (p_lead | p_breadth) / (leg_lead | leg_breadth) +
    plot_layout(heights = c(1, 0.12), widths = c(1, 1)) +
    plot_annotation(
      title = overall_title,
      theme = theme(
        plot.title = element_text(
          size = 12,
          face = "bold",
          hjust = 0.5,
          colour = "black"
        )
      )
    )
}

save_both <- function(plot_obj, filename_base, width, height, folder) {
  ggsave(
    file.path(folder, paste0(filename_base, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(folder, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
}

# ============================================================
# 3. Opportunity-space maps and tables
# ============================================================

phases <- tribble(
  ~phase_key, ~phase_label, ~ka_min, ~ka_max,
  "late", "Late Holocene (4.2–0 ka BP)", 0.0, 4.2,
  "mid", "Mid-Holocene (8.2–4.2 ka BP)", 4.2, 8.2,
  "early", "Early Holocene (11.7–8.2 ka BP)", 8.2, 11.7
)

phase_results <- lapply(seq_len(nrow(phases)), function(i) {
  make_phase_data(
    label = phases$phase_label[i],
    ka_min = phases$ka_min[i],
    ka_max = phases$ka_max[i]
  )
})

names(phase_results) <- phases$phase_key

fig_late <- make_map_pair(
  phase_results$late,
  overall_title = "Late Holocene (4.2–0 ka BP)"
)

fig_mid <- make_map_pair(
  phase_results$mid,
  overall_title = "Mid-Holocene (8.2–4.2 ka BP)"
)

fig_early <- make_map_pair(
  phase_results$early,
  overall_title = "Early Holocene (11.7–8.2 ka BP)"
)

save_both(
  fig_late,
  filename_base = "Fig_opportunity_space",
  width = 12.0,
  height = 7.6,
  folder = fig_dir
)

save_both(
  fig_mid,
  filename_base = "SM_mid_holocene_opportunity_space",
  width = 12.0,
  height = 7.6,
  folder = supp_fig_dir
)

save_both(
  fig_early,
  filename_base = "SM_early_holocene_opportunity_space",
  width = 12.0,
  height = 7.6,
  folder = supp_fig_dir
)

regional_crop_fraction_tbl <- bind_rows(lapply(phase_results, `[[`, "crop_frac")) |>
  arrange(phase, region, crop)

regional_strength_tbl <- bind_rows(lapply(phase_results, `[[`, "region_stats")) |>
  arrange(phase, region)

opportunity_summary_tbl <- regional_crop_fraction_tbl |>
  left_join(regional_strength_tbl, by = c("phase", "region"))

write.csv(
  regional_crop_fraction_tbl,
  file.path(table_dir, "opportunity_crop_fraction.csv"),
  row.names = FALSE
)

write.csv(
  regional_strength_tbl,
  file.path(table_dir, "opportunity_strength_breadth.csv"),
  row.names = FALSE
)

write.csv(
  opportunity_summary_tbl,
  file.path(table_dir, "opportunity_space_summary_detailed.csv"),
  row.names = FALSE
)

# ============================================================
# 4. Quantitative supplementary summary figure/table
# ============================================================

summarise_phase_quant <- function(res, min_floor = 0.30, plausible_thr = 0.50) {

  met <- compute_metric_rasters(
    res$r_stack,
    min_floor = min_floor,
    plausible_thr = plausible_thr
  )

  region_r <- rasterize(vect(res$regions), res$r_stack[[1]], field = "region")
  area_r <- cellSize(res$r_stack[[1]], unit = "km")

  df <- as.data.frame(
    c(region_r, area_r, met$top_idx, met$strength, met$n_options),
    xy = FALSE,
    na.rm = FALSE
  )

  names(df) <- c("region", "cell_km2", "top_idx", "strength", "n_options")

  df <- df |>
    filter(!is.na(region), !is.na(cell_km2)) |>
    mutate(
      phase = unique(res$lead_df$phase),
      region = as.character(region),
      crop = case_when(
        top_idx == 1 ~ "Sorghum",
        top_idx == 2 ~ "Pearl millet",
        top_idx == 3 ~ "Finger millet",
        TRUE ~ NA_character_
      )
    )

  crop_frac <- df |>
    filter(!is.na(crop)) |>
    group_by(phase, region, crop) |>
    summarise(lead_area_km2 = sum(cell_km2, na.rm = TRUE), .groups = "drop_last") |>
    mutate(frac_lead_area = lead_area_km2 / sum(lead_area_km2, na.rm = TRUE)) |>
    ungroup()

  region_stats <- df |>
    group_by(phase, region) |>
    summarise(
      mean_strength_pct = 100 * weighted.mean(strength, w = cell_km2, na.rm = TRUE),
      mean_n_options = weighted.mean(n_options, w = cell_km2, na.rm = TRUE),
      .groups = "drop"
    )

  breadth_comp <- df |>
    mutate(n_options = factor(n_options, levels = c(0, 1, 2, 3))) |>
    group_by(phase, region, n_options) |>
    summarise(area_km2 = sum(cell_km2, na.rm = TRUE), .groups = "drop_last") |>
    mutate(area_share = area_km2 / sum(area_km2, na.rm = TRUE)) |>
    ungroup()

  list(
    crop_frac = crop_frac,
    region_stats = region_stats,
    breadth_comp = breadth_comp
  )
}

quant_list <- lapply(phase_results, summarise_phase_quant)

crop_frac_all <- bind_rows(lapply(quant_list, `[[`, "crop_frac")) |>
  mutate(
    phase = factor(phase, levels = phase_order),
    region = factor(region, levels = region_order),
    crop = factor(crop, levels = c("Sorghum", "Pearl millet", "Finger millet"))
  )

region_stats_all <- bind_rows(lapply(quant_list, `[[`, "region_stats")) |>
  mutate(
    phase = factor(phase, levels = phase_order),
    region = factor(region, levels = region_order)
  )

breadth_comp_all <- bind_rows(lapply(quant_list, `[[`, "breadth_comp")) |>
  mutate(
    phase = factor(phase, levels = phase_order),
    region = factor(region, levels = region_order),
    n_options = factor(n_options, levels = c(0, 1, 2, 3))
  )

crop_frac_wide <- crop_frac_all |>
  select(phase, region, crop, frac_lead_area) |>
  pivot_wider(
    names_from = crop,
    values_from = frac_lead_area,
    names_prefix = "frac_lead_"
  )

breadth_wide <- breadth_comp_all |>
  select(phase, region, n_options, area_share) |>
  pivot_wider(
    names_from = n_options,
    values_from = area_share,
    names_prefix = "share_n"
  )

supp_table <- region_stats_all |>
  left_join(crop_frac_wide, by = c("phase", "region")) |>
  left_join(breadth_wide, by = c("phase", "region")) |>
  arrange(phase, region)

write.csv(
  supp_table,
  file.path(table_dir, "opportunity_space_summary.csv"),
  row.names = FALSE
)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(
    supp_table,
    file.path(table_dir, "opportunity_space_summary.xlsx"),
    overwrite = TRUE
  )
}

theme_supp <- function() {
  theme_minimal(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = 9.5, colour = "black"),
      strip.background = element_rect(fill = "grey94", colour = "black", linewidth = 0.7),
      axis.title = element_text(face = "bold", size = 9.5, colour = "black"),
      axis.text = element_text(face = "bold", size = 8.5, colour = "black"),
      legend.title = element_text(face = "bold", size = 9, colour = "black"),
      legend.text = element_text(size = 8.2, colour = "black"),
      plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
      plot.margin = margin(4, 4, 4, 4)
    )
}

p_crop_fraction <- ggplot(
  crop_frac_all,
  aes(x = frac_lead_area, y = fct_rev(region), fill = crop)
) +
  geom_col(width = 0.72, colour = NA) +
  facet_wrap(~phase, ncol = 1) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = percent_format(accuracy = 1)
  ) +
  scale_fill_manual(values = crop_cols_named) +
  labs(
    title = "A. Regional fraction of leading-option area",
    x = "Share of leading-option area",
    y = NULL,
    fill = NULL
  ) +
  theme_supp() +
  theme(legend.position = "bottom")

p_breadth <- ggplot(
  breadth_comp_all,
  aes(x = area_share, y = fct_rev(region), fill = n_options)
) +
  geom_col(width = 0.72, colour = NA) +
  facet_wrap(~phase, ncol = 1) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = percent_format(accuracy = 1)
  ) +
  scale_fill_manual(values = breadth_cols) +
  labs(
    title = "B. Regional composition of climatic choice breadth",
    x = "Share of regional area",
    y = NULL,
    fill = "Number of plausible crop options"
  ) +
  theme_supp() +
  theme(legend.position = "bottom")

fig_summary <- p_crop_fraction / p_breadth +
  plot_layout(heights = c(1, 1))

save_both(
  fig_summary,
  filename_base = "SM_opportunity_space_summary",
  width = 8.3,
  height = 9.2,
  folder = supp_fig_dir
)

# ============================================================
# 5. Optional regional trajectory figure
# ============================================================
# This uses the suitability summary table from script 02. It is kept here
# because it was present in the loose figure script, but it is not required
# for the opportunity-space maps above.

if (file.exists(suitability_table)) {

  suit_tab <- read_excel(suitability_table) |>
    filter(region != "Africa") |>
    mutate(
      time_ka = time_bp / 1000,
      area_mkm2 = high_suit_area_km2 / 1e6,
      crop = recode(
        crop,
        "SG" = "Sorghum",
        "PM" = "Pearl millet",
        "FM" = "Finger millet"
      )
    )

  label_pts <- st_point_on_surface(regions_sf)

  p_map <- ggplot() +
    geom_sf(data = regions_sf, aes(fill = region), color = "grey35", linewidth = 0.55) +
    geom_sf(data = st_union(regions_sf), fill = NA, color = "black", linewidth = 0.9) +
    geom_sf_text(data = label_pts, aes(label = region), size = 3.2, fontface = "bold") +
    scale_fill_manual(values = c(
      "Northern" = "#f2f2f2",
      "Western" = "#e6e6e6",
      "Central" = "#d9d9d9",
      "Eastern" = "#cccccc",
      "Southern" = "#bdbdbd"
    )) +
    coord_sf(expand = FALSE) +
    theme_void() +
    theme(legend.position = "none")

  phase_cols <- c(
    "African Humid Period" = "#d9ead3",
    "Post-mid-Holocene aridification" = "#f4cccc"
  )

  make_reg_plot <- function(reg_name) {

    d <- suit_tab |> filter(region == reg_name)

    yr <- range(d$area_mkm2, na.rm = TRUE)
    pad <- diff(yr) * 0.18
    y0 <- yr[1] - pad
    y1 <- yr[1] - pad * 0.45

    phase_df <- tibble(
      xmin = c(10, 6),
      xmax = c(6, 0),
      ymin = y0,
      ymax = y1,
      phase = c(
        "African Humid Period",
        "Post-mid-Holocene aridification"
      )
    )

    ggplot(d, aes(time_ka, area_mkm2, color = crop)) +
      geom_rect(
        data = phase_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = phase),
        inherit.aes = FALSE,
        color = NA,
        alpha = 0.8
      ) +
      geom_line(alpha = 0.22, linewidth = 0.35) +
      geom_smooth(
        se = TRUE,
        method = "loess",
        span = 0.22,
        linewidth = 0.9,
        alpha = 0.12
      ) +
      scale_color_manual(values = crop_cols_named) +
      scale_fill_manual(values = phase_cols) +
      scale_x_reverse(limits = c(10, 0), breaks = c(10, 6, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0.14, 0.06))) +
      labs(
        title = paste0(reg_name, " Africa"),
        x = NULL,
        y = expression("Highly suitable area (10"^6*" km"^2*")")
      ) +
      theme_minimal(base_size = 8.5) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9.5),
        axis.title.y = element_text(face = "bold", size = 8),
        axis.text = element_text(face = "bold", size = 7.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.55),
        legend.position = "none",
        plot.margin = margin(3, 3, 3, 3)
      )
  }

  p_north <- make_reg_plot("Northern")
  p_west <- make_reg_plot("Western")
  p_south <- make_reg_plot("Southern")
  p_east <- make_reg_plot("Eastern")
  p_cent <- make_reg_plot("Central")

  leg_crop <- ggplot(suit_tab, aes(time_ka, area_mkm2, color = crop)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = crop_cols_named) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 9)
    ) +
    guides(color = guide_legend(nrow = 1))

  legend_crop <- cowplot::get_legend(leg_crop)

  fig_region_traj <-
    (
      (p_north / p_west / p_south) |
        p_map |
        (p_east / p_cent / plot_spacer())
    ) +
    plot_layout(widths = c(1.25, 1.05, 1.25)) +
    plot_annotation(
      title = "Regional trajectories of highly suitable crop area",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
      )
    )

  fig_region_traj <- fig_region_traj / wrap_elements(legend_crop) +
    plot_layout(heights = c(1, 0.06))

  save_both(
    fig_region_traj,
    filename_base = "SM_regional_suitability_trajectories",
    width = 11.5,
    height = 8.2,
    folder = supp_fig_dir
  )
}

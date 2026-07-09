# ============================================================
# 04_threshold_metrics.R
# Recalculate threshold-based suitable-area metrics
# ============================================================

# This script reads continuous EcoCrop suitability rasters and
# calculates area-weighted suitable-area metrics for African regions.
# The resulting table is the corrected source of suitable-area
# proportions used by the GAM and threshold-sensitivity analyses.

# ---------- packages ----------
req_pkgs <- c("terra", "dplyr", "purrr", "stringr", "openxlsx", "tibble")
missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install required packages before running this script: ", paste(missing_pkgs, collapse = ", "))
}

library(terra)
library(dplyr)
library(purrr)
library(stringr)
library(openxlsx)
library(tibble)

# ---------- paths ----------
# Run this script from the repository root.
project_dir <- getwd()

raster_dir <- file.path(project_dir, "data", "processed", "ecocrop_rasters")
region_dir <- file.path(project_dir, "data", "raw", "regions")
out_file <- file.path(project_dir, "data", "processed", "threshold_metrics.xlsx")

region_file <- file.path(region_dir, "africa_UN_regions.shp")
if (!file.exists(region_file)) {
  region_file <- file.path(region_dir, "africa_regions.gpkg")
}
if (!file.exists(region_file)) {
  stop("Africa region file not found in data/raw/regions/.")
}

# ---------- settings ----------
crops <- c("SG", "PM", "FM")
thresholds <- c(0.4, 0.6, 0.7, 0.8)
region_col <- "region5"

# ============================================================
# 1. Read regions and list suitability rasters
# ============================================================

regions <- vect(region_file)

if (!region_col %in% names(regions)) {
  if ("region" %in% names(regions)) {
    region_col <- "region"
  } else {
    stop("Region file must contain a region5 or region column.")
  }
}

files <- map_dfr(crops, function(cr) {
  crop_dir <- file.path(raster_dir, cr)
  if (!dir.exists(crop_dir)) stop("Raster folder not found: ", crop_dir)

  f <- list.files(
    crop_dir,
    pattern = paste0("^", cr, "_suit_t-?[0-9]+\\.tif$"),
    full.names = TRUE
  )

  tibble(
    crop = cr,
    file = f,
    filename = basename(f),
    time_id = str_extract(basename(f), "t-?[0-9]+"),
    time_num = as.numeric(str_remove(str_extract(basename(f), "t-?[0-9]+"), "t")),
    time_bp = time_num * 100
  )
}) %>%
  arrange(crop, time_num)

if (nrow(files) == 0) {
  stop("No continuous suitability rasters found in data/processed/ecocrop_rasters/.")
}

print(files %>% count(crop))

# ============================================================
# 2. Precompute regional area weights
# ============================================================

template <- rast(files$file[1])
regions_proj <- project(regions, crs(template))

cell_area <- cellSize(template, unit = "km", mask = TRUE)
names(cell_area) <- "cell_area_km2"

weights <- extract(cell_area, regions_proj, cells = TRUE, exact = TRUE)
weights$region <- as.data.frame(regions_proj)[[region_col]][weights$ID]

weights <- weights %>%
  filter(!is.na(region), !is.na(cell_area_km2), !is.na(fraction)) %>%
  transmute(
    cell,
    region,
    area_weighted_km2 = cell_area_km2 * fraction
  )

total_area <- weights %>%
  group_by(region) %>%
  summarise(
    total_area_km2 = sum(area_weighted_km2, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 3. Calculate metrics for one raster
# ============================================================

calc_one <- function(file, crop, time_id, time_bp, i, n_total) {

  pct <- round(100 * i / n_total, 1)
  cat("[", pct, "%] ", basename(file), "\n", sep = "")

  r <- rast(file)
  names(r) <- "suitability"

  vals <- as.data.frame(r, cells = TRUE, na.rm = TRUE)

  dat <- weights %>%
    inner_join(vals, by = "cell")

  out <- total_area

  suit_stats <- dat %>%
    group_by(region) %>%
    summarise(
      mean_suitability = weighted.mean(suitability, area_weighted_km2, na.rm = TRUE),
      median_suitability = median(suitability, na.rm = TRUE),
      min_suitability = min(suitability, na.rm = TRUE),
      max_suitability = max(suitability, na.rm = TRUE),
      sd_suitability = sd(suitability, na.rm = TRUE),
      .groups = "drop"
    )

  out <- out %>%
    left_join(suit_stats, by = "region")

  for (thr in thresholds) {
    area_thr <- dat %>%
      group_by(region) %>%
      summarise(
        area_km2 = sum(area_weighted_km2 * (suitability >= thr), na.rm = TRUE),
        .groups = "drop"
      )

    out <- out %>%
      left_join(area_thr, by = "region")

    names(out)[names(out) == "area_km2"] <- paste0("area_", thr, "_km2")
  }

  out %>%
    mutate(
      crop = crop,
      time_id = time_id,
      time_bp = time_bp,
      prop_0.4 = area_0.4_km2 / total_area_km2,
      prop_0.6 = area_0.6_km2 / total_area_km2,
      prop_0.7 = area_0.7_km2 / total_area_km2,
      prop_0.8 = area_0.8_km2 / total_area_km2
    ) %>%
    select(
      crop, region, time_id, time_bp,
      total_area_km2,
      mean_suitability, median_suitability,
      min_suitability, max_suitability, sd_suitability,
      area_0.4_km2, area_0.6_km2, area_0.7_km2, area_0.8_km2,
      prop_0.4, prop_0.6, prop_0.7, prop_0.8
    )
}

# ============================================================
# 4. Run all rasters and save output
# ============================================================

threshold_metrics <- pmap_dfr(
  files %>% mutate(i = row_number(), n_total = n()),
  function(crop, file, filename, time_id, time_num, time_bp, i, n_total) {
    calc_one(file, crop, time_id, time_bp, i, n_total)
  }
)

check_props <- threshold_metrics %>%
  mutate(
    diff_04 = prop_0.4 - area_0.4_km2 / total_area_km2,
    diff_06 = prop_0.6 - area_0.6_km2 / total_area_km2,
    diff_07 = prop_0.7 - area_0.7_km2 / total_area_km2,
    diff_08 = prop_0.8 - area_0.8_km2 / total_area_km2
  ) %>%
  summarise(
    max_diff_04 = max(abs(diff_04), na.rm = TRUE),
    max_diff_06 = max(abs(diff_06), na.rm = TRUE),
    max_diff_07 = max(abs(diff_07), na.rm = TRUE),
    max_diff_08 = max(abs(diff_08), na.rm = TRUE)
  )

print(threshold_metrics %>% count(crop, region))
print(check_props)

write.xlsx(threshold_metrics, out_file, overwrite = TRUE)

cat("\nSaved threshold metrics to:\n", out_file, "\n")

# ============================================================
# 01_climate_data_prep.R
# Prepare Africa climate rasters for EcoCrop modelling
# ============================================================

# ---------- packages ----------
req_pkgs <- c("terra", "sf", "dplyr", "rnaturalearth", "rnaturalearthdata")

to_install <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(terra)
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

# ---------- paths ----------
# Run this script from the repository root.
project_dir <- getwd()

# Raw CHELSA-TraCE21k source folder.
# This large external dataset is NOT uploaded to GitHub.
src_dir <- "path to the raw data folder/CHELSA-TraCE21k/Data"

data_raw_dir <- file.path(project_dir, "data", "raw")
clim_out_dir <- file.path(project_dir, "data", "processed", "climate")
tab_dir <- file.path(project_dir, "outputs", "tables")
log_dir <- file.path(project_dir, "outputs", "logs")

dir.create(data_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clim_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(clim_out_dir, "pr"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(clim_out_dir, "tasmin"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(clim_out_dir, "tasmax"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(clim_out_dir, "tmean"), recursive = TRUE, showWarnings = FALSE)

region_file <- file.path(data_raw_dir, "africa_regions.gpkg")
country_file <- file.path(data_raw_dir, "africa_countries.gpkg")

log_file <- file.path(log_dir, "climate_log.txt")
fail_file <- file.path(log_dir, "climate_failed.csv")

# ---------- settings ----------
time_ids <- -100:20
fact_agg <- 10

# time convention:
# label = 20 - time_id
# label 120 = 12 ka BP, label 0 = present

cat("CHELSA Africa climate processing\n", file = log_file)
cat(paste("Started:", Sys.time(), "\n\n"), file = log_file, append = TRUE)

# ============================================================
# 1. Create Africa regional polygons
# ============================================================

world <- ne_countries(scale = "medium", returnclass = "sf")

africa_countries <- world %>%
  filter(continent == "Africa") %>%
  mutate(region5 = case_when(
    subregion == "Northern Africa" ~ "North",
    subregion == "Western Africa" ~ "West",
    subregion == "Eastern Africa" ~ "East",
    subregion == "Middle Africa" ~ "Central",
    subregion == "Southern Africa" ~ "South",
    TRUE ~ NA_character_
  ))

bbox <- st_bbox(
  c(xmin = -20, xmax = 55, ymin = -35, ymax = 38),
  crs = st_crs(africa_countries)
)

africa_countries <- st_crop(africa_countries, bbox)

africa_regions <- africa_countries %>%
  filter(!is.na(region5)) %>%
  group_by(region5) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

if (file.exists(region_file)) file.remove(region_file)
if (file.exists(country_file)) file.remove(country_file)

st_write(africa_regions, region_file, quiet = TRUE)
st_write(africa_countries, country_file, quiet = TRUE)

# ============================================================
# 2. Helper functions
# ============================================================

load_stack <- function(var, tid, src_dir) {
  
  f <- list.files(
    file.path(src_dir, var),
    full.names = TRUE,
    pattern = paste0("^CHELSA_TraCE21k_", var, "_[0-9]+_", tid, "_V1\\.0\\.tif$"),
    ignore.case = TRUE
  )
  
  m <- as.integer(
    gsub(
      paste0("^CHELSA_TraCE21k_", var, "_([0-9]+)_.*$"),
      "\\1",
      basename(f)
    )
  )
  
  f <- f[order(m)]
  m <- m[order(m)]
  
  if (length(f) != 12) {
    stop(paste("Expected 12 files for", var, "time_id", tid, "but found", length(f)))
  }
  
  if (!all(m == 1:12)) {
    stop(paste("Month order problem for", var, "time_id", tid))
  }
  
  r <- rast(f)
  names(r) <- paste0("month_", 1:12)
  r
}

extract_region_stats <- function(r_pr, r_tm, regions_v, label_val) {
  
  regs <- unique(regions_v$region5)
  africa_mask <- aggregate(regions_v)
  
  pr_af <- mask(r_pr, africa_mask)
  tm_af <- mask(r_tm, africa_mask)
  
  pr_m_af <- as.numeric(global(pr_af, "mean", na.rm = TRUE)[1, ])
  tm_m_af <- as.numeric(global(tm_af, "mean", na.rm = TRUE)[1, ])
  
  annual_out <- list()
  monthly_out <- list()
  
  annual_out[[1]] <- data.frame(
    label = label_val,
    region = "Africa",
    annual_pr_mm = sum(pr_m_af, na.rm = TRUE),
    annual_tmean_c = mean(tm_m_af, na.rm = TRUE)
  )
  
  monthly_out[[1]] <- data.frame(
    label = label_val,
    region = "Africa",
    month = 1:12,
    pr_mm = pr_m_af,
    tmean_c = tm_m_af
  )
  
  k <- 2
  
  for (rg in regs) {
    
    poly <- regions_v[regions_v$region5 == rg, ]
    
    pr_rg <- mask(r_pr, poly)
    tm_rg <- mask(r_tm, poly)
    
    pr_m <- as.numeric(global(pr_rg, "mean", na.rm = TRUE)[1, ])
    tm_m <- as.numeric(global(tm_rg, "mean", na.rm = TRUE)[1, ])
    
    annual_out[[k]] <- data.frame(
      label = label_val,
      region = rg,
      annual_pr_mm = sum(pr_m, na.rm = TRUE),
      annual_tmean_c = mean(tm_m, na.rm = TRUE)
    )
    
    monthly_out[[k]] <- data.frame(
      label = label_val,
      region = rg,
      month = 1:12,
      pr_mm = pr_m,
      tmean_c = tm_m
    )
    
    k <- k + 1
  }
  
  list(
    annual = do.call(rbind, annual_out),
    monthly = do.call(rbind, monthly_out)
  )
}

# ============================================================
# 3. Process climate rasters
# ============================================================

africa_sf <- st_read(region_file, quiet = TRUE)

fail_log <- data.frame(
  time_id = integer(),
  label = integer(),
  error = character(),
  stringsAsFactors = FALSE
)

annual_all <- list()
monthly_all <- list()

for (i in seq_along(time_ids)) {
  
  tid <- time_ids[i]
  label_val <- 20 - tid
  
  pr_out <- file.path(clim_out_dir, "pr", paste0("CHELSA_pr_", label_val, ".tif"))
  tn_out <- file.path(clim_out_dir, "tasmin", paste0("CHELSA_tasmin_", label_val, ".tif"))
  tx_out <- file.path(clim_out_dir, "tasmax", paste0("CHELSA_tasmax_", label_val, ".tif"))
  tm_out <- file.path(clim_out_dir, "tmean", paste0("CHELSA_tmean_", label_val, ".tif"))
  
  cat(sprintf("\n[%d/%d] Processing time_id = %s | label = %s\n",
              i, length(time_ids), tid, label_val))
  
  cat(sprintf("[%d/%d] START time_id=%s label=%s time=%s\n",
              i, length(time_ids), tid, label_val, Sys.time()),
      file = log_file, append = TRUE)
  
  tryCatch({
    
    if (file.exists(pr_out) &&
        file.exists(tn_out) &&
        file.exists(tx_out) &&
        file.exists(tm_out)) {
      
      cat(sprintf("Skipping label %s; processed rasters already exist\n", label_val))
      
      pr <- rast(pr_out)
      tm <- rast(tm_out)
      
      africa_v <- vect(st_transform(africa_sf, crs(pr)))
      stats <- extract_region_stats(pr, tm, africa_v, label_val)
      
      annual_all[[length(annual_all) + 1]] <- stats$annual
      monthly_all[[length(monthly_all) + 1]] <- stats$monthly
      
      cat(sprintf("label=%s skipped; summaries refreshed at %s\n",
                  label_val, Sys.time()),
          file = log_file, append = TRUE)
      
    } else {
      
      pr <- load_stack("pr", tid, src_dir)
      tx <- load_stack("tasmax", tid, src_dir)
      tn <- load_stack("tasmin", tid, src_dir)
      
      africa_v <- vect(st_transform(africa_sf, crs(pr)))
      
      pr <- mask(crop(pr, africa_v), africa_v)
      tx <- mask(crop(tx, africa_v), africa_v)
      tn <- mask(crop(tn, africa_v), africa_v)
      
      pr <- aggregate(pr, fact = fact_agg, fun = mean, na.rm = TRUE)
      tx <- aggregate(tx, fact = fact_agg, fun = mean, na.rm = TRUE)
      tn <- aggregate(tn, fact = fact_agg, fun = mean, na.rm = TRUE)
      
      tx <- tx * 0.1 - 273.15
      tn <- tn * 0.1 - 273.15
      tm <- (tx + tn) / 2
      
      writeRaster(pr, pr_out, overwrite = TRUE)
      writeRaster(tn, tn_out, overwrite = TRUE)
      writeRaster(tx, tx_out, overwrite = TRUE)
      writeRaster(tm, tm_out, overwrite = TRUE)
      
      stats <- extract_region_stats(pr, tm, africa_v, label_val)
      
      annual_all[[length(annual_all) + 1]] <- stats$annual
      monthly_all[[length(monthly_all) + 1]] <- stats$monthly
      
      cat(sprintf("Finished label %s at %s\n", label_val, Sys.time()))
      cat(sprintf("label=%s finished at %s\n", label_val, Sys.time()),
          file = log_file, append = TRUE)
    }
    
  }, error = function(e) {
    
    fail_log <<- rbind(
      fail_log,
      data.frame(
        time_id = tid,
        label = label_val,
        error = e$message,
        stringsAsFactors = FALSE
      )
    )
    
    cat(sprintf("FAILED time_id = %s | label = %s\n", tid, label_val))
    cat(sprintf("FAILED time_id=%s label=%s error=%s time=%s\n",
                tid, label_val, e$message, Sys.time()),
        file = log_file, append = TRUE)
  })
}

# ============================================================
# 4. Save summary tables
# ============================================================

annual_df <- do.call(rbind, annual_all)
monthly_df <- do.call(rbind, monthly_all)

annual_df$ka_bp <- annual_df$label / 10
monthly_df$ka_bp <- monthly_df$label / 10

write.csv(
  annual_df,
  file.path(tab_dir, "climate_annual.csv"),
  row.names = FALSE
)

write.csv(
  monthly_df,
  file.path(tab_dir, "climate_monthly.csv"),
  row.names = FALSE
)

write.csv(
  fail_log,
  fail_file,
  row.names = FALSE
)

cat(sprintf("\nRun finished at %s\n", Sys.time()))
cat(sprintf("\nRun finished at %s\n", Sys.time()),
    file = log_file,
    append = TRUE)
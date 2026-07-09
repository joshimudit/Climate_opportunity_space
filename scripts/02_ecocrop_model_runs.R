# ============================================================
# 02_ecocrop_model_runs.R
# Run EcoCrop suitability models for African dryland cereals
# ============================================================

# This script uses monthly temperature and precipitation rasters
# prepared by 01_climate_data_prep.R and estimates EcoCrop
# suitability for sorghum, pearl millet, and finger millet.

# ---------- packages ----------
req_pkgs <- c("terra", "Recocrop", "openxlsx")

missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Install required packages before running this script: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(terra)
library(Recocrop)
library(openxlsx)

# ---------- paths ----------
# Run this script from the repository root.

project_dir <- getwd()

climate_dir <- file.path(project_dir, "data", "processed", "climate")
region_dir  <- file.path(project_dir, "data", "raw", "regions")

raster_dir <- file.path(project_dir, "data", "processed", "ecocrop_rasters")
table_dir  <- file.path(project_dir, "data", "processed")
log_dir    <- file.path(project_dir, "outputs", "logs")
qc_dir     <- file.path(project_dir, "outputs", "qc_maps")

dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

region_file <- file.path(region_dir, "africa_UN_regions.shp")

if (!file.exists(region_file)) {
  region_file <- file.path(region_dir, "africa_regions.gpkg")
}

if (!file.exists(region_file)) {
  stop("Africa region file not found in data/raw/regions/.")
}

log_file  <- file.path(log_dir, "ecocrop_log.txt")
fail_file <- file.path(log_dir, "ecocrop_failed.xlsx")

# ---------- settings ----------
time_ids <- 0:120
high_thr <- 0.8

save_binary_rasters <- FALSE
save_qc_maps <- FALSE

time_bp_fun <- function(x) x * 100

cat("EcoCrop suitability workflow\n", file = log_file)
cat(paste("Started:", Sys.time(), "\n\n"), file = log_file, append = TRUE)

# ============================================================
# 1. Crop parameter setup
# ============================================================

pm <- ecocropPars("Pearl millet")$parameters
fm <- ecocropPars("Finger millet")$parameters

sg1 <- ecocropPars("Sorghum (med. altitude)")$parameters
sg2 <- ecocropPars("Sorghum (low altitude)")$parameters
sg3 <- ecocropPars("Sorghum (high altitude)")$parameters

sg <- sg1
sg[] <- NA

sg[1, ] <- pmin(sg1[1, ], sg2[1, ], sg3[1, ], na.rm = TRUE)
sg[2, ] <- pmin(sg1[2, ], sg2[2, ], sg3[2, ], na.rm = TRUE)
sg[3, ] <- pmax(sg1[3, ], sg2[3, ], sg3[3, ], na.rm = TRUE)
sg[4, ] <- pmax(sg1[4, ], sg2[4, ], sg3[4, ], na.rm = TRUE)

crop_models <- list(
  SG = ecocrop(list(name = "SG", parameters = sg)),
  PM = ecocrop(list(name = "PM", parameters = pm)),
  FM = ecocrop(list(name = "FM", parameters = fm))
)

# ============================================================
# 2. Region file
# ============================================================

regions <- vect(region_file)

if ("region5" %in% names(regions)) {
  regions$region <- regions$region5
}

if (!"region" %in% names(regions)) {
  stop("Region file must contain a region or region5 column.")
}

# ============================================================
# 3. Helper functions
# ============================================================

log_msg <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  write(msg, file = log_file, append = TRUE)
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

africa_suitability_row <- function(suit, area_rast, crop_name, time_id, time_bp, high_thr) {
  
  sv <- values(suit, mat = FALSE)
  av <- values(area_rast, mat = FALSE)
  
  keep <- !is.na(sv)
  sv <- sv[keep]
  av <- av[keep]
  
  high <- sv >= high_thr
  
  data.frame(
    crop = crop_name,
    time_id = time_id,
    time_bp = time_bp,
    region = "Africa",
    mean_suitability = safe_mean(sv),
    median_suitability = safe_median(sv),
    min_suitability = safe_min(sv),
    max_suitability = safe_max(sv),
    high_suit_area_prop = safe_mean(high),
    high_suit_area_km2 = sum(av[high], na.rm = TRUE)
  )
}

region_suitability_rows <- function(suit, area_rast, regions, crop_name, time_id, time_bp, high_thr) {
  
  s_df <- extract(suit, regions)
  a_df <- extract(area_rast, regions)
  
  out <- lapply(seq_len(nrow(regions)), function(i) {
    
    sv <- s_df[s_df$ID == i, 2]
    av <- a_df[a_df$ID == i, 2]
    
    keep <- !is.na(sv)
    sv <- sv[keep]
    av <- av[keep]
    
    high <- sv >= high_thr
    
    data.frame(
      crop = crop_name,
      time_id = time_id,
      time_bp = time_bp,
      region = values(regions)$region[i],
      mean_suitability = safe_mean(sv),
      median_suitability = safe_median(sv),
      min_suitability = safe_min(sv),
      max_suitability = safe_max(sv),
      high_suit_area_prop = safe_mean(high),
      high_suit_area_km2 = sum(av[high], na.rm = TRUE)
    )
  })
  
  do.call(rbind, out)
}

region_mean_from_raster <- function(x, regions) {
  z <- extract(x, regions, fun = mean, na.rm = TRUE)
  data.frame(region = values(regions)$region[z$ID], value = z[[2]])
}

africa_mean_from_raster <- function(x) {
  global(x, "mean", na.rm = TRUE)[1, 1]
}

monthly_sd_raster <- function(x) {
  app(x, sd, na.rm = TRUE)
}

climate_summary_rows <- function(r_tavg, r_prec, regions, time_id, time_bp) {
  
  mat_r <- mean(r_tavg)
  map_r <- sum(r_prec)
  
  t_warmest_r <- app(r_tavg, max, na.rm = TRUE)
  t_coldest_r <- app(r_tavg, min, na.rm = TRUE)
  p_wettest_r <- app(r_prec, max, na.rm = TRUE)
  p_driest_r <- app(r_prec, min, na.rm = TRUE)
  
  temp_seasonality_r <- monthly_sd_raster(r_tavg)
  prec_seasonality_r <- monthly_sd_raster(r_prec)
  
  t_sorted <- sort(r_tavg, decreasing = TRUE)
  p_sorted <- sort(r_prec, decreasing = TRUE)
  
  gs_temp_r <- mean(t_sorted[[1:6]])
  gs_prec_r <- mean(p_sorted[[1:6]])
  
  reg_list <- list(
    mat = region_mean_from_raster(mat_r, regions),
    map = region_mean_from_raster(map_r, regions),
    t_warmest = region_mean_from_raster(t_warmest_r, regions),
    t_coldest = region_mean_from_raster(t_coldest_r, regions),
    p_wettest = region_mean_from_raster(p_wettest_r, regions),
    p_driest = region_mean_from_raster(p_driest_r, regions),
    temp_seasonality = region_mean_from_raster(temp_seasonality_r, regions),
    prec_seasonality = region_mean_from_raster(prec_seasonality_r, regions),
    gs_temp = region_mean_from_raster(gs_temp_r, regions),
    gs_prec = region_mean_from_raster(gs_prec_r, regions)
  )
  
  reg_tab <- data.frame(region = values(regions)$region)
  
  for (nm in names(reg_list)) {
    reg_tab[[nm]] <- reg_list[[nm]]$value[
      match(reg_tab$region, reg_list[[nm]]$region)
    ]
  }
  
  africa_row <- data.frame(
    region = "Africa",
    mat = africa_mean_from_raster(mat_r),
    map = africa_mean_from_raster(map_r),
    t_warmest = africa_mean_from_raster(t_warmest_r),
    t_coldest = africa_mean_from_raster(t_coldest_r),
    p_wettest = africa_mean_from_raster(p_wettest_r),
    p_driest = africa_mean_from_raster(p_driest_r),
    temp_seasonality = africa_mean_from_raster(temp_seasonality_r),
    prec_seasonality = africa_mean_from_raster(prec_seasonality_r),
    gs_temp = africa_mean_from_raster(gs_temp_r),
    gs_prec = africa_mean_from_raster(gs_prec_r)
  )
  
  out <- rbind(africa_row, reg_tab)
  out$time_id <- time_id
  out$time_bp <- time_bp
  
  out[, c(
    "time_id", "time_bp", "region",
    "mat", "map", "t_warmest", "t_coldest",
    "p_wettest", "p_driest",
    "temp_seasonality", "prec_seasonality",
    "gs_temp", "gs_prec"
  )]
}

save_png_map <- function(r, regions, file_out, main_txt) {
  png(filename = file_out, width = 1800, height = 1200, res = 200)
  plot(r, main = main_txt)
  plot(regions, add = TRUE, border = "black", lwd = 1)
  dev.off()
}

# ============================================================
# 4. Main EcoCrop loop
# ============================================================

table_A <- data.frame()
table_B <- data.frame()
failed_runs <- data.frame()

for (time_id in time_ids) {
  
  time_bp <- time_bp_fun(time_id)
  
  log_msg("Starting time step t", time_id, " (", time_bp, " BP)")
  
  f_tavg <- file.path(climate_dir, "tmean", paste0("CHELSA_tmean_", time_id, ".tif"))
  f_prec <- file.path(climate_dir, "pr", paste0("CHELSA_pr_", time_id, ".tif"))
  
  if (!file.exists(f_tavg) || !file.exists(f_prec)) {
    
    msg <- paste("Missing climate files for time_id", time_id)
    log_msg("FAILED: ", msg)
    
    failed_runs <- rbind(
      failed_runs,
      data.frame(
        crop = NA_character_,
        time_id = time_id,
        time_bp = time_bp,
        error = msg
      )
    )
    
    next
  }
  
  r_tavg <- rast(f_tavg)
  r_prec <- rast(f_prec)
  
  if (nlyr(r_tavg) != 12 || nlyr(r_prec) != 12) {
    
    msg <- paste("Climate stack does not have 12 layers for time_id", time_id)
    log_msg("FAILED: ", msg)
    
    failed_runs <- rbind(
      failed_runs,
      data.frame(
        crop = NA_character_,
        time_id = time_id,
        time_bp = time_bp,
        error = msg
      )
    )
    
    next
  }
  
  if (!compareGeom(r_tavg, r_prec, stopOnError = FALSE)) {
    
    msg <- paste("Temperature and precipitation geometry mismatch for time_id", time_id)
    log_msg("FAILED: ", msg)
    
    failed_runs <- rbind(
      failed_runs,
      data.frame(
        crop = NA_character_,
        time_id = time_id,
        time_bp = time_bp,
        error = msg
      )
    )
    
    next
  }
  
  if (!same.crs(regions, r_tavg)) {
    
    msg <- paste("Region polygons and raster CRS mismatch for time_id", time_id)
    log_msg("FAILED: ", msg)
    
    failed_runs <- rbind(
      failed_runs,
      data.frame(
        crop = NA_character_,
        time_id = time_id,
        time_bp = time_bp,
        error = msg
      )
    )
    
    next
  }
  
  area_rast <- cellSize(r_tavg[[1]], unit = "km")
  
  clim_rows <- climate_summary_rows(
    r_tavg = r_tavg,
    r_prec = r_prec,
    regions = regions,
    time_id = time_id,
    time_bp = time_bp
  )
  
  table_A <- rbind(table_A, clim_rows)
  
  for (crop_name in names(crop_models)) {
    
    log_msg("Running crop ", crop_name, " at t", time_id)
    
    tryCatch({
      
      m <- crop_models[[crop_name]]
      control(m, get_max = TRUE)
      
      suit <- predict(m, tavg = r_tavg, prec = r_prec)
      names(suit) <- "suitability"
      
      high_suit <- ifel(suit >= high_thr, 1, 0)
      names(high_suit) <- "high_suitability"
      
      crop_raster_dir <- file.path(raster_dir, crop_name)
      crop_qc_dir <- file.path(qc_dir, crop_name)
      
      dir.create(crop_raster_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(crop_qc_dir, recursive = TRUE, showWarnings = FALSE)
      
      suit_file <- file.path(
        crop_raster_dir,
        paste0(crop_name, "_suit_t", time_id, ".tif")
      )
      
      writeRaster(suit, suit_file, overwrite = TRUE)
      
      if (save_binary_rasters) {
        
        high_file <- file.path(
          crop_raster_dir,
          paste0(crop_name, "_highsuit_t", time_id, ".tif")
        )
        
        writeRaster(high_suit, high_file, overwrite = TRUE)
      }
      
      if (save_qc_maps) {
        
        map_file <- file.path(
          crop_qc_dir,
          paste0(crop_name, "_highsuit_t", time_id, ".png")
        )
        
        save_png_map(
          high_suit,
          regions,
          map_file,
          paste0(crop_name, " high suitability (>= ", high_thr, ") | t", time_id)
        )
      }
      
      africa_row <- africa_suitability_row(
        suit = suit,
        area_rast = area_rast,
        crop_name = crop_name,
        time_id = time_id,
        time_bp = time_bp,
        high_thr = high_thr
      )
      
      reg_rows <- region_suitability_rows(
        suit = suit,
        area_rast = area_rast,
        regions = regions,
        crop_name = crop_name,
        time_id = time_id,
        time_bp = time_bp,
        high_thr = high_thr
      )
      
      table_B <- rbind(table_B, africa_row, reg_rows)
      
      log_msg(
        "Finished crop ", crop_name, " at t", time_id,
        " | Africa mean suitability: ", round(africa_row$mean_suitability, 4),
        " | Africa high-suitability area km2: ", round(africa_row$high_suit_area_km2, 2)
      )
      
    }, error = function(e) {
      
      log_msg("FAILED crop ", crop_name, " at t", time_id, " :: ", e$message)
      
      failed_runs <<- rbind(
        failed_runs,
        data.frame(
          crop = crop_name,
          time_id = time_id,
          time_bp = time_bp,
          error = e$message
        )
      )
    })
  }
  
  log_msg("Completed all crops for t", time_id)
}

# ============================================================
# 5. Save summary tables
# ============================================================

write.xlsx(
  table_A,
  file.path(table_dir, "ecocrop_climate.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  table_B,
  file.path(table_dir, "ecocrop_suitability.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  failed_runs,
  fail_file,
  overwrite = TRUE
)

log_msg("EcoCrop workflow complete.")
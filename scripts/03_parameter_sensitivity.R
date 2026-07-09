# ============================================================
# 03_parameter_sensitivity.R
# EcoCrop parameter sensitivity analysis
# ============================================================

# This script evaluates whether Africa-level EcoCrop outputs are sensitive
# to uncertainty in crop climatic parameters. It uses two representative
# Holocene climate states and perturbs temperature, precipitation, and
# growing-season parameters using Monte Carlo sampling.

# ---------- packages ----------
required_pkgs <- c(
  "terra", "Recocrop", "dplyr", "tidyr", "ggplot2",
  "purrr", "writexl", "tibble"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Install required packages before running this script: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(terra)
library(Recocrop)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(writexl)
library(tibble)

# ---------- paths ----------
# Run this script from the repository root.

project_dir <- getwd()

climate_dir <- file.path(project_dir, "data", "processed", "climate")
pr_dir <- file.path(climate_dir, "pr")
tmean_dir <- file.path(climate_dir, "tmean")

sens_dir <- file.path(project_dir, "data", "processed", "parameter_sensitivity")
raw_dir <- file.path(sens_dir, "raw")
figdata_dir <- file.path(sens_dir, "figure_data")

table_dir <- file.path(project_dir, "outputs", "tables")
fig_dir <- file.path(project_dir, "outputs", "figures", "supplement")

invisible(lapply(
  c(sens_dir, raw_dir, figdata_dir, table_dir, fig_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# ---------- settings ----------
set.seed(123)

n_runs <- 500
time_ids <- c(60, 20)
high_thr <- 0.8
save_figures <- TRUE

# time_id 60 = 6 ka BP; time_id 20 = 2 ka BP

time_tbl <- tibble(
  time_id = time_ids,
  time_bp = time_id * 100,
  time_ka = time_bp / 1000,
  period_label = paste0(time_ka, " ka BP"),
  pr_file = file.path(pr_dir, paste0("CHELSA_pr_", time_id, ".tif")),
  tmean_file = file.path(tmean_dir, paste0("CHELSA_tmean_", time_id, ".tif"))
)

if (!all(file.exists(time_tbl$pr_file))) {
  stop(
    "Missing precipitation files:\n",
    paste(time_tbl$pr_file[!file.exists(time_tbl$pr_file)], collapse = "\n")
  )
}

if (!all(file.exists(time_tbl$tmean_file))) {
  stop(
    "Missing temperature files:\n",
    paste(time_tbl$tmean_file[!file.exists(time_tbl$tmean_file)], collapse = "\n")
  )
}

# ============================================================
# 1. EcoCrop parameters
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

crop_pars <- list(
  Sorghum = sg,
  Pearl_millet = pm,
  Finger_millet = fm
)

crop_labels <- c(
  Sorghum = "Sorghum",
  Pearl_millet = "Pearl millet",
  Finger_millet = "Finger millet"
)

param_audit <- purrr::imap_dfr(crop_pars, function(p, crop_name) {
  as.data.frame(p) |>
    tibble::rownames_to_column("limit_row") |>
    mutate(crop = crop_name, .before = 1)
})

# ============================================================
# 2. Monte Carlo parameter samples
# ============================================================

mc_params <- tidyr::expand_grid(
  crop = names(crop_pars),
  run = 1:n_runs
) |>
  group_by(crop) |>
  mutate(
    temp_shift = runif(n(), -2, 2),
    prec_mult = runif(n(), 0.8, 1.2),
    gs_mult = runif(n(), 0.85, 1.15)
  ) |>
  ungroup()

find_gs_cols <- function(pars) {
  grep(
    "grow|season|duration|cycle|length",
    colnames(pars),
    ignore.case = TRUE,
    value = TRUE
  )
}

perturb_pars <- function(pars, temp_shift, prec_mult, gs_mult) {
  p <- pars

  if ("tavg" %in% colnames(p)) {
    p[, "tavg"] <- p[, "tavg"] + temp_shift
  }

  if ("prec" %in% colnames(p)) {
    p[, "prec"] <- p[, "prec"] * prec_mult
  }

  gs_cols <- find_gs_cols(p)

  if (length(gs_cols) > 0) {
    p[, gs_cols] <- p[, gs_cols] * gs_mult
  }

  p
}

gs_cols_found <- lapply(crop_pars, find_gs_cols)

# ============================================================
# 3. Model functions
# ============================================================

load_climate <- function(time_id) {
  prec <- terra::rast(file.path(pr_dir, paste0("CHELSA_pr_", time_id, ".tif")))
  tavg <- terra::rast(file.path(tmean_dir, paste0("CHELSA_tmean_", time_id, ".tif")))

  names(prec) <- rep("prec", terra::nlyr(prec))
  names(tavg) <- rep("tavg", terra::nlyr(tavg))

  if (!terra::compareGeom(prec, tavg, stopOnError = FALSE)) {
    stop("Temperature and precipitation rasters do not match for time_id ", time_id)
  }

  list(prec = prec, tavg = tavg)
}

run_ecocrop <- function(pars, tavg, prec, crop_name) {
  m <- ecocrop(list(name = crop_name, parameters = pars))
  control(m, get_max = TRUE)
  predict(m, tavg = tavg, prec = prec)
}

summarise_africa <- function(r, area_vals, high_thr) {
  v <- terra::values(r, mat = FALSE)

  keep <- !is.na(v) & !is.na(area_vals)
  v <- v[keep]
  a <- area_vals[keep]

  total_area <- sum(a, na.rm = TRUE)
  high_area <- sum(a[v >= high_thr], na.rm = TRUE)

  tibble(
    mean_suitability = sum(v * a, na.rm = TRUE) / total_area,
    area_0.8_km2 = high_area,
    total_area_km2 = total_area,
    prop_0.8 = high_area / total_area
  )
}

# Area weights are derived from the climate raster grid used in the sensitivity test.
template <- terra::rast(time_tbl$pr_file[1])[[1]]
area_vals <- terra::values(terra::cellSize(template, unit = "km"), mat = FALSE)

# ============================================================
# 4. Baseline runs
# ============================================================

message("Running baseline EcoCrop models")

baseline_list <- list()

for (i in seq_len(nrow(time_tbl))) {
  tt <- time_tbl[i, ]
  clim <- load_climate(tt$time_id)

  for (crop_name in names(crop_pars)) {
    r <- run_ecocrop(crop_pars[[crop_name]], clim$tavg, clim$prec, crop_name)

    baseline_list[[length(baseline_list) + 1]] <-
      summarise_africa(r, area_vals, high_thr) |>
      mutate(
        model_type = "baseline",
        crop = crop_name,
        time_id = tt$time_id,
        time_bp = tt$time_bp,
        time_ka = tt$time_ka,
        period_label = tt$period_label,
        run = 0,
        temp_shift = 0,
        prec_mult = 1,
        gs_mult = 1
      )
  }

  rm(clim)
  gc()
}

baseline <- bind_rows(baseline_list)

# ============================================================
# 5. Monte Carlo runs
# ============================================================

message("Running Monte Carlo EcoCrop models")

mc_list <- list()
job <- 0
total_jobs <- nrow(time_tbl) * length(crop_pars) * n_runs
start_time <- Sys.time()

for (i in seq_len(nrow(time_tbl))) {
  tt <- time_tbl[i, ]
  clim <- load_climate(tt$time_id)

  for (crop_name in names(crop_pars)) {
    crop_mc <- mc_params |> filter(crop == crop_name)

    for (j in seq_len(nrow(crop_mc))) {
      job <- job + 1

      pars_j <- perturb_pars(
        crop_pars[[crop_name]],
        crop_mc$temp_shift[j],
        crop_mc$prec_mult[j],
        crop_mc$gs_mult[j]
      )

      r <- run_ecocrop(pars_j, clim$tavg, clim$prec, crop_name)

      mc_list[[length(mc_list) + 1]] <-
        summarise_africa(r, area_vals, high_thr) |>
        mutate(
          model_type = "perturbed",
          crop = crop_name,
          time_id = tt$time_id,
          time_bp = tt$time_bp,
          time_ka = tt$time_ka,
          period_label = tt$period_label,
          run = crop_mc$run[j],
          temp_shift = crop_mc$temp_shift[j],
          prec_mult = crop_mc$prec_mult[j],
          gs_mult = crop_mc$gs_mult[j]
        )

      if (j %% 25 == 0) {
        elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
        pct <- round(100 * job / total_jobs, 1)
        message("Completed ", pct, "% | elapsed ", elapsed, " min")
      }
    }
  }

  saveRDS(bind_rows(mc_list), file.path(raw_dir, "mc_partial.rds"))

  rm(clim)
  gc()
}

mc_runs <- bind_rows(mc_list)

# ============================================================
# 6. Save raw model outputs
# ============================================================

saveRDS(baseline, file.path(raw_dir, "baseline.rds"))
saveRDS(mc_runs, file.path(raw_dir, "mc_runs.rds"))
saveRDS(mc_params, file.path(raw_dir, "mc_params.rds"))

# ============================================================
# 7. Summary tables
# ============================================================

metric_labels <- c(
  mean_suitability = "Mean suitability",
  prop_0.8 = "Highly suitable area"
)

base_long <- baseline |>
  select(crop, time_id, time_bp, time_ka, period_label, mean_suitability, prop_0.8) |>
  pivot_longer(
    cols = c(mean_suitability, prop_0.8),
    names_to = "metric",
    values_to = "baseline"
  )

mc_long <- mc_runs |>
  select(
    crop, run, time_id, time_bp, time_ka, period_label,
    temp_shift, prec_mult, gs_mult, mean_suitability, prop_0.8
  ) |>
  pivot_longer(
    cols = c(mean_suitability, prop_0.8),
    names_to = "metric",
    values_to = "value"
  ) |>
  left_join(
    base_long,
    by = c("crop", "time_id", "time_bp", "time_ka", "period_label", "metric")
  ) |>
  mutate(
    delta = value - baseline,
    abs_delta = abs(delta),
    crop = factor(crop, levels = names(crop_labels), labels = crop_labels),
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels),
    period_label = factor(period_label, levels = c("6 ka BP", "2 ka BP"))
  )

summary_metric <- function(x) {
  tibble(
    mc_median = median(x, na.rm = TRUE),
    mc_mean = mean(x, na.rm = TRUE),
    mc_sd = sd(x, na.rm = TRUE),
    mc_cv = sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE),
    q05 = quantile(x, 0.05, na.rm = TRUE),
    q95 = quantile(x, 0.95, na.rm = TRUE)
  )
}

ST1_parameter_design <- tibble(
  parameter = c("Temperature limits", "Precipitation limits", "Growing-season length"),
  perturbation = c("Additive shift", "Multiplicative scaling", "Multiplicative scaling"),
  range = c("-2 to +2", "0.8 to 1.2", "0.85 to 1.15"),
  unit = c("°C", "proportion", "proportion"),
  interpretation = c(
    "Uncertainty in crop thermal tolerance",
    "Uncertainty in crop rainfall requirement",
    "Uncertainty in phenological duration, if present in the EcoCrop object"
  )
)

ST2_summary <- mc_long |>
  group_by(crop, period_label, metric) |>
  summarise(
    baseline = first(baseline),
    summary_metric(value),
    median_delta = median(delta, na.rm = TRUE),
    median_abs_delta = median(abs_delta, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    delta_percent_of_baseline = 100 * median_delta / baseline,
    abs_delta_percent_of_baseline = 100 * median_abs_delta / baseline
  )

ST3_effects <- mc_long |>
  group_by(crop, period_label, metric) |>
  summarise(
    cor_temperature = cor(temp_shift, value, use = "complete.obs"),
    cor_precipitation = cor(prec_mult, value, use = "complete.obs"),
    cor_growing_season = cor(gs_mult, value, use = "complete.obs"),
    std_slope_temperature = coef(lm(scale(value) ~ scale(temp_shift)))[2],
    std_slope_precipitation = coef(lm(scale(value) ~ scale(prec_mult)))[2],
    std_slope_growing_season = coef(lm(scale(value) ~ scale(gs_mult)))[2],
    .groups = "drop"
  )

ST4_claims <- ST2_summary |>
  select(crop, period_label, metric, mc_cv, median_abs_delta, abs_delta_percent_of_baseline) |>
  pivot_wider(
    names_from = metric,
    values_from = c(mc_cv, median_abs_delta, abs_delta_percent_of_baseline)
  ) |>
  mutate(
    cv_ratio_area_vs_mean = `mc_cv_Highly suitable area` / `mc_cv_Mean suitability`
  )

writexl::write_xlsx(
  list(
    ST1_parameter_design = ST1_parameter_design,
    ST2_sensitivity_summary = ST2_summary,
    ST3_parameter_effects = ST3_effects,
    ST4_key_claims = ST4_claims,
    original_EcoCrop_parameters = param_audit,
    MC_parameter_samples = mc_params
  ),
  file.path(table_dir, "parameter_sensitivity.xlsx")
)

saveRDS(mc_long, file.path(figdata_dir, "mc_long.rds"))
saveRDS(ST2_summary, file.path(figdata_dir, "summary.rds"))
saveRDS(ST3_effects, file.path(figdata_dir, "effects.rds"))
saveRDS(ST4_claims, file.path(figdata_dir, "claims.rds"))

# ============================================================
# 8. Supplementary figures
# ============================================================

save_plot <- function(p, name, width = 10, height = 7) {
  if (save_figures) {
    ggsave(
      file.path(fig_dir, paste0(name, ".png")),
      p,
      width = width,
      height = height,
      dpi = 400
    )
  }
}

SM1 <- mc_long |>
  filter(metric == "Mean suitability") |>
  ggplot(aes(temp_shift, prec_mult, colour = value)) +
  geom_point(alpha = 0.45, size = 1) +
  facet_grid(period_label ~ crop) +
  scale_colour_viridis_c(option = "C", name = "Mean\nsuitability") +
  theme_classic(base_size = 11) +
  labs(
    x = "Temperature shift (°C)",
    y = "Precipitation multiplier",
    title = "SM1. Monte Carlo parameter space coloured by mean suitability"
  )

save_plot(SM1, "SM1_parameter_space_mean_suitability", 11, 6)

SM2 <- mc_long |>
  filter(metric == "Highly suitable area") |>
  ggplot(aes(temp_shift, prec_mult, colour = value)) +
  geom_point(alpha = 0.45, size = 1) +
  facet_grid(period_label ~ crop) +
  scale_colour_viridis_c(option = "C", name = "Prop.\narea ≥0.8") +
  theme_classic(base_size = 11) +
  labs(
    x = "Temperature shift (°C)",
    y = "Precipitation multiplier",
    title = "SM2. Monte Carlo parameter space coloured by highly suitable area"
  )

save_plot(SM2, "SM2_parameter_space_high_suitability", 11, 6)

SM3 <- ST2_summary |>
  ggplot(aes(period_label, mc_median)) +
  geom_linerange(aes(ymin = q05, ymax = q95), linewidth = 0.8) +
  geom_point(size = 2) +
  geom_point(aes(y = baseline), shape = 1, stroke = 1, size = 2.8) +
  facet_grid(metric ~ crop, scales = "free_y") +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = NULL,
    title = "SM3. Baseline EcoCrop output against Monte Carlo uncertainty envelope",
    subtitle = "Open circles show baseline; filled circles show Monte Carlo median; vertical lines show 5–95% range"
  )

save_plot(SM3, "SM3_baseline_vs_MC_envelope", 11, 6.5)

SM4 <- ST2_summary |>
  ggplot(aes(metric, mc_cv)) +
  geom_col(width = 0.65) +
  facet_grid(period_label ~ crop, scales = "free_y") +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "Coefficient of variation",
    title = "SM4. Relative parameter sensitivity of continuous and threshold-derived outputs"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_plot(SM4, "SM4_relative_sensitivity_CV", 11, 6)

SM5 <- mc_long |>
  ggplot(aes(delta)) +
  geom_histogram(bins = 45, colour = "black", linewidth = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_grid(metric + period_label ~ crop, scales = "free") +
  theme_classic(base_size = 10) +
  labs(
    x = "Monte Carlo output − baseline output",
    y = "Frequency",
    title = "SM5. Deviations from baseline EcoCrop outputs"
  )

save_plot(SM5, "SM5_delta_from_baseline_distributions", 11, 8)

param_effect_long <- ST3_effects |>
  select(crop, period_label, metric, cor_temperature, cor_precipitation, cor_growing_season) |>
  pivot_longer(
    cols = starts_with("cor_"),
    names_to = "parameter",
    values_to = "correlation"
  ) |>
  mutate(
    parameter = recode(
      parameter,
      cor_temperature = "Temperature",
      cor_precipitation = "Precipitation",
      cor_growing_season = "Growing season"
    ),
    abs_correlation = abs(correlation)
  )

SM6 <- param_effect_long |>
  ggplot(aes(parameter, abs_correlation)) +
  geom_col(width = 0.65) +
  facet_grid(metric + period_label ~ crop) +
  theme_classic(base_size = 10) +
  labs(
    x = NULL,
    y = "|Correlation| with output",
    title = "SM6. Relative influence of perturbed EcoCrop parameter groups"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_plot(SM6, "SM6_parameter_importance", 11, 8)

response_df <- mc_long |>
  select(crop, period_label, metric, value, temp_shift, prec_mult, gs_mult) |>
  pivot_longer(
    cols = c(temp_shift, prec_mult, gs_mult),
    names_to = "parameter",
    values_to = "parameter_value"
  ) |>
  mutate(
    parameter = recode(
      parameter,
      temp_shift = "Temperature shift",
      prec_mult = "Precipitation multiplier",
      gs_mult = "Growing-season multiplier"
    )
  )

SM7 <- response_df |>
  ggplot(aes(parameter_value, value)) +
  geom_point(alpha = 0.12, size = 0.45) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.7) +
  facet_grid(metric + period_label ~ crop + parameter, scales = "free_x") +
  theme_classic(base_size = 9) +
  labs(
    x = NULL,
    y = "Output value",
    title = "SM7. Response curves of EcoCrop outputs to parameter perturbations"
  )

save_plot(SM7, "SM7_parameter_response_curves", 15, 8)

SM8 <- ST2_summary |>
  filter(metric == "Highly suitable area") |>
  ggplot(aes(period_label, mc_cv, group = crop)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~crop, nrow = 1, scales = "free_y") +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "CV of highly suitable area",
    title = "SM8. Parameter sensitivity of highly suitable area across climate states"
  )

save_plot(SM8, "SM8_climate_state_sensitivity_contrast", 9, 4)

message("Parameter sensitivity workflow complete.")

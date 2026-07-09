# ============================================================
# 05_gam_analysis.R
# Corrected GAM analysis of EcoCrop climate responses
# ============================================================

# This script estimates climate-response curves for EcoCrop
# suitable-area proportions using the corrected area-weighted
# threshold metrics. Result 3 uses the suitability >= 0.8 response.

# ---------- packages ----------
req_pkgs <- c(
  "readxl", "dplyr", "tidyr", "purrr", "ggplot2",
  "mgcv", "tibble", "openxlsx", "patchwork"
)
missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install required packages before running this script: ", paste(missing_pkgs, collapse = ", "))
}

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(mgcv)
library(tibble)
library(openxlsx)
library(patchwork)

# ---------- paths ----------
# Run this script from the repository root.
project_dir <- getwd()

threshold_file <- file.path(project_dir, "data", "processed", "threshold_metrics.xlsx")
climate_file <- file.path(project_dir, "data", "processed", "climate_for_gam.xlsx")

gam_dir <- file.path(project_dir, "data", "processed", "gam")
fig_dir <- file.path(project_dir, "outputs", "figures")
supp_fig_dir <- file.path(project_dir, "outputs", "figures", "supplement")
tab_dir <- file.path(project_dir, "outputs", "tables")

dir.create(gam_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(threshold_file)) stop("Missing input: ", threshold_file)
if (!file.exists(climate_file)) stop("Missing input: ", climate_file)

# ---------- settings ----------
region_order <- c("North", "West", "East", "Central", "South")

crop_labs <- c(
  SG = "Sorghum",
  PM = "Pearl millet",
  FM = "Finger millet"
)

crop_cols <- c(
  "Sorghum" = "#2166AC",
  "Pearl millet" = "#1B9E77",
  "Finger millet" = "#D95F02"
)

region_full_labs <- c(
  North = "Northern Africa",
  West = "Western Africa",
  East = "Eastern Africa",
  Central = "Central Africa",
  South = "Southern Africa"
)

# ============================================================
# 1. Read and prepare corrected GAM input
# ============================================================

threshold_raw <- read_excel(threshold_file)
climate_raw <- read_excel(climate_file)

threshold_dat <- threshold_raw %>%
  filter(region != "Africa") %>%
  mutate(time_num = as.numeric(gsub("t", "", as.character(time_id))))

climate_dat <- climate_raw %>%
  filter(region != "Africa") %>%
  mutate(time_num = as.numeric(gsub("t", "", as.character(time_id)))) %>%
  select(crop, region, time_num, mat, map) %>%
  distinct()

dat_r3 <- threshold_dat %>%
  left_join(climate_dat, by = c("crop", "region", "time_num")) %>%
  mutate(
    prop = prop_0.8,
    prop_clamped = pmin(pmax(prop, 0.001), 0.999),
    logit_prop = qlogis(prop_clamped),
    crop = factor(crop, levels = c("SG", "PM", "FM")),
    region = factor(region, levels = region_order)
  ) %>%
  filter(
    !is.na(mat),
    !is.na(map),
    !is.na(prop),
    !is.na(logit_prop)
  )

write.xlsx(dat_r3, file.path(gam_dir, "gam_input.xlsx"), overwrite = TRUE)

# ============================================================
# 2. Fit one GAM per region
# ============================================================

fit_region_gam <- function(d) {
  gam(
    logit_prop ~ crop +
      s(mat, by = crop, k = 5) +
      s(map, by = crop, k = 5),
    data = d,
    method = "REML"
  )
}

region_data <- split(dat_r3, dat_r3$region, drop = TRUE)
region_data <- region_data[region_order]
region_data <- region_data[!vapply(region_data, is.null, logical(1))]

r3_models <- map(region_data, fit_region_gam)
saveRDS(r3_models, file.path(gam_dir, "gam_models.rds"))

# ============================================================
# 3. Effect-size table
# ============================================================

extract_effects <- function(region_name, model, data) {

  d <- data %>% filter(region == region_name)

  base_dat <- d %>%
    group_by(crop) %>%
    summarise(
      mat = median(mat, na.rm = TRUE),
      map = median(map, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(crop = factor(crop, levels = levels(data$crop)))

  temp_dat <- base_dat %>% mutate(mat = mat + 1)
  prec_dat <- base_dat %>% mutate(map = map + 100)

  p_base <- plogis(predict(model, newdata = base_dat))
  p_temp <- plogis(predict(model, newdata = temp_dat))
  p_prec <- plogis(predict(model, newdata = prec_dat))

  tibble(
    region = region_name,
    crop = as.character(base_dat$crop),
    baseline_prop = as.numeric(p_base),
    delta_temp_pp = as.numeric(100 * (p_temp - p_base)),
    delta_prec_pp = as.numeric(100 * (p_prec - p_base))
  )
}

gam_effects <- map2_dfr(
  names(r3_models),
  r3_models,
  ~ extract_effects(.x, .y, dat_r3)
) %>%
  mutate(
    crop = recode(crop, !!!crop_labs),
    region = recode(region, !!!region_full_labs),
    dominant_driver = case_when(
      baseline_prop < 0.05 & abs(delta_temp_pp) < 1 & abs(delta_prec_pp) < 1 ~ "Weak",
      abs(delta_temp_pp) > abs(delta_prec_pp) ~ "Temperature",
      abs(delta_prec_pp) > abs(delta_temp_pp) ~ "Precipitation",
      TRUE ~ "Mixed"
    ),
    response_type = case_when(
      dominant_driver == "Weak" ~ "Weak",
      abs(delta_temp_pp) >= 1 & abs(delta_prec_pp) >= 1 &
        sign(delta_temp_pp) != sign(delta_prec_pp) ~ "Mixed / trade-off",
      dominant_driver == "Temperature" ~ "Temperature-dominated",
      dominant_driver == "Precipitation" ~ "Precipitation-dominated",
      TRUE ~ "Mixed"
    )
  ) %>%
  arrange(region, crop)

gam_table <- gam_effects %>%
  mutate(
    baseline_prop = round(baseline_prop, 3),
    delta_temp_pp = round(delta_temp_pp, 1),
    delta_prec_pp = round(delta_prec_pp, 1)
  )

write.xlsx(gam_table, file.path(tab_dir, "gam_effects.xlsx"), overwrite = TRUE)

# ============================================================
# 4. Smooth summaries, k-checks, and k sensitivity
# ============================================================

gam_smooths <- map_dfr(names(r3_models), function(r) {
  summary(r3_models[[r]])$s.table %>%
    as.data.frame() %>%
    rownames_to_column("smooth") %>%
    as_tibble() %>%
    mutate(region = r, .before = 1)
})

write.xlsx(gam_smooths, file.path(tab_dir, "gam_smooths.xlsx"), overwrite = TRUE)

gam_kcheck <- map_dfr(names(r3_models), function(r) {
  k.check(r3_models[[r]]) %>%
    as.data.frame() %>%
    rownames_to_column("smooth") %>%
    as_tibble() %>%
    mutate(region = r, .before = 1)
})

write.xlsx(gam_kcheck, file.path(tab_dir, "gam_kcheck.xlsx"), overwrite = TRUE)

fit_models_k <- function(k_value) {
  map(region_data, ~ gam(
    logit_prop ~ crop +
      s(mat, by = crop, k = k_value) +
      s(map, by = crop, k = k_value),
    data = .x,
    method = "REML"
  ))
}

extract_all_effects_k <- function(k_value) {
  models_k <- fit_models_k(k_value)
  map2_dfr(names(models_k), models_k, ~ extract_effects(.x, .y, dat_r3)) %>%
    mutate(k = k_value, .before = 1)
}

gam_k_sensitivity <- bind_rows(
  extract_all_effects_k(5),
  extract_all_effects_k(6),
  extract_all_effects_k(7)
) %>%
  mutate(
    crop = recode(crop, !!!crop_labs),
    region = recode(region, !!!region_full_labs)
  )

gam_k_compare <- gam_k_sensitivity %>%
  select(k, region, crop, baseline_prop, delta_temp_pp, delta_prec_pp) %>%
  pivot_wider(
    names_from = k,
    values_from = c(baseline_prop, delta_temp_pp, delta_prec_pp),
    names_prefix = "k"
  ) %>%
  mutate(
    temp_diff_k6 = delta_temp_pp_k6 - delta_temp_pp_k5,
    temp_diff_k7 = delta_temp_pp_k7 - delta_temp_pp_k5,
    prec_diff_k6 = delta_prec_pp_k6 - delta_prec_pp_k5,
    prec_diff_k7 = delta_prec_pp_k7 - delta_prec_pp_k5
  )

write.xlsx(gam_k_sensitivity, file.path(tab_dir, "gam_k_sensitivity.xlsx"), overwrite = TRUE)
write.xlsx(gam_k_compare, file.path(tab_dir, "gam_k_compare.xlsx"), overwrite = TRUE)

# ============================================================
# 5. Main GAM response figure
# ============================================================

dat_plot <- dat_r3 %>%
  mutate(
    crop_name = recode(as.character(crop), !!!crop_labs),
    crop_name = factor(crop_name, levels = names(crop_cols)),
    region = factor(as.character(region), levels = region_order)
  )

make_pred_grid <- function(region_name, xvar = "mat", n = 150) {

  d <- dat_plot %>% filter(region == region_name)
  model <- r3_models[[as.character(region_name)]]

  bind_rows(lapply(levels(d$crop_name), function(cr) {

    dc <- d %>% filter(crop_name == cr)

    if (xvar == "mat") {
      nd <- data.frame(
        mat = seq(min(dc$mat), max(dc$mat), length.out = n),
        map = median(dc$map, na.rm = TRUE),
        crop = factor(names(crop_labs)[crop_labs == cr], levels = levels(dat_r3$crop))
      )
    } else {
      nd <- data.frame(
        mat = median(dc$mat, na.rm = TRUE),
        map = seq(min(dc$map), max(dc$map), length.out = n),
        crop = factor(names(crop_labs)[crop_labs == cr], levels = levels(dat_r3$crop))
      )
    }

    pr <- predict(model, newdata = nd, se.fit = TRUE)

    nd %>%
      mutate(
        region = region_name,
        crop_name = cr,
        fit = plogis(pr$fit),
        lwr = plogis(pr$fit - 1.96 * pr$se.fit),
        upr = plogis(pr$fit + 1.96 * pr$se.fit)
      )
  }))
}

pred_temp <- bind_rows(lapply(levels(dat_plot$region), make_pred_grid, xvar = "mat"))
pred_prec <- bind_rows(lapply(levels(dat_plot$region), make_pred_grid, xvar = "map"))

saveRDS(list(temp = pred_temp, prec = pred_prec), file.path(gam_dir, "gam_predictions.rds"))

p_temp <- ggplot() +
  geom_point(data = dat_plot, aes(mat, prop, colour = crop_name), alpha = 0.35, size = 1) +
  geom_ribbon(data = pred_temp, aes(mat, ymin = lwr, ymax = upr, fill = crop_name), alpha = 0.15, colour = NA) +
  geom_line(data = pred_temp, aes(mat, fit, colour = crop_name), linewidth = 1.1) +
  facet_wrap(~ region, nrow = 1, scales = "free") +
  scale_colour_manual(values = crop_cols) +
  scale_fill_manual(values = crop_cols) +
  labs(title = "A. Temperature response", x = "Mean annual temperature (°C)", y = "Highly suitable area proportion") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

p_prec <- ggplot() +
  geom_point(data = dat_plot, aes(map, prop, colour = crop_name), alpha = 0.35, size = 1) +
  geom_ribbon(data = pred_prec, aes(map, ymin = lwr, ymax = upr, fill = crop_name), alpha = 0.15, colour = NA) +
  geom_line(data = pred_prec, aes(map, fit, colour = crop_name), linewidth = 1.1) +
  facet_wrap(~ region, nrow = 1, scales = "free") +
  scale_colour_manual(values = crop_cols) +
  scale_fill_manual(values = crop_cols) +
  labs(title = "B. Precipitation response", x = "Mean annual precipitation (mm)", y = "Highly suitable area proportion") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

fig_gam <- p_temp / p_prec +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Fig_GAM_climate_response.png"), fig_gam, width = 13, height = 8, dpi = 600)
ggsave(file.path(fig_dir, "Fig_GAM_climate_response.pdf"), fig_gam, width = 13, height = 8)

# ============================================================
# 6. Supplementary ECDF and diagnostics
# ============================================================

fig_ecdf <- ggplot(dat_plot, aes(prop, colour = crop_name)) +
  stat_ecdf(linewidth = 1.1) +
  facet_wrap(~ region, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = crop_cols) +
  labs(
    x = "Highly suitable area proportion",
    y = "Empirical cumulative probability"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(supp_fig_dir, "SM_GAM_ECDF.png"), fig_ecdf, width = 9, height = 6.5, dpi = 600)
ggsave(file.path(supp_fig_dir, "SM_GAM_ECDF.pdf"), fig_ecdf, width = 9, height = 6.5)

diag_dir <- file.path(supp_fig_dir, "gam_diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

for (r in names(r3_models)) {
  png(file.path(diag_dir, paste0("SM_GAM_diagnostics_", r, ".png")), width = 8, height = 8, units = "in", res = 600)
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0), mar = c(4.2, 4.2, 2.5, 1.2))
  gam.check(r3_models[[r]])
  mtext(region_full_labs[[r]], outer = TRUE, side = 3, line = 1, font = 2, cex = 1.3)
  dev.off()
}

cat("\nCorrected GAM analysis complete.\n")

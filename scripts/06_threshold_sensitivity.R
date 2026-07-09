# ============================================================
# 06_threshold_sensitivity.R
# Threshold-dependent climate sensitivity analysis
# ============================================================

# This script tests whether inferred climate sensitivity changes
# when suitable area is defined using alternative EcoCrop suitability
# thresholds. The main analysis uses thresholds 0.6, 0.7, and 0.8.

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

sens_dir <- file.path(project_dir, "data", "processed", "threshold_sensitivity")
fig_dir <- file.path(project_dir, "outputs", "figures")
supp_fig_dir <- file.path(project_dir, "outputs", "figures", "supplement")
tab_dir <- file.path(project_dir, "outputs", "tables")

dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(threshold_file)) stop("Missing input: ", threshold_file)
if (!file.exists(climate_file)) stop("Missing input: ", climate_file)

# ---------- settings ----------
region_order <- c("North", "West", "East", "Central", "South")
thresholds_main <- c("0.6", "0.7", "0.8")
thresholds_driver <- c("0.4", "0.6", "0.7", "0.8")

crop_labs <- c(
  SG = "Sorghum",
  PM = "Pearl millet",
  FM = "Finger millet"
)

region_labs <- c(
  North = "Northern Africa",
  West = "Western Africa",
  East = "Eastern Africa",
  Central = "Central Africa",
  South = "Southern Africa"
)

threshold_cols <- c(
  "Not statistically significant" = "grey70",
  "Suitability ≥ 0.6" = "#D55E00",
  "Suitability ≥ 0.7" = "#1B9E77",
  "Suitability ≥ 0.8" = "black"
)

# ============================================================
# 1. Read and prepare corrected threshold data
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

dat_base <- threshold_dat %>%
  left_join(climate_dat, by = c("crop", "region", "time_num")) %>%
  filter(!is.na(mat), !is.na(map))

# ============================================================
# 2. Main threshold-sensitivity dataset
# ============================================================

dat_r4 <- dat_base %>%
  select(
    crop, region, time_id, time_num, time_bp,
    mat, map,
    prop_0.6, prop_0.7, prop_0.8
  ) %>%
  pivot_longer(
    cols = starts_with("prop_"),
    names_to = "threshold",
    values_to = "prop"
  ) %>%
  mutate(
    threshold = gsub("prop_", "", threshold),
    threshold_label = paste0("Suitability ≥ ", threshold),
    prop_clamped = pmin(pmax(prop, 0.001), 0.999),
    logit_prop = qlogis(prop_clamped),
    crop = factor(as.character(crop), levels = c("SG", "PM", "FM")),
    region = factor(as.character(region), levels = region_order)
  ) %>%
  filter(threshold %in% thresholds_main, !is.na(prop), !is.na(logit_prop))

write.xlsx(dat_r4, file.path(sens_dir, "threshold_gam_input.xlsx"), overwrite = TRUE)

# ============================================================
# 3. Fit one GAM per threshold and region
# ============================================================

fit_model <- function(d) {
  gam(
    logit_prop ~ crop +
      s(mat, by = crop, k = 5) +
      s(map, by = crop, k = 5),
    data = d,
    method = "REML"
  )
}

model_groups <- dat_r4 %>%
  group_by(threshold, region) %>%
  group_split()

model_names <- dat_r4 %>%
  distinct(threshold, region) %>%
  arrange(threshold, region) %>%
  transmute(name = paste(threshold, region, sep = "_")) %>%
  pull(name)

r4_models <- setNames(map(model_groups, fit_model), model_names)
saveRDS(r4_models, file.path(sens_dir, "threshold_gam_models.rds"))

# ============================================================
# 4. Extract threshold-dependent effect sizes
# ============================================================

get_smooth_p <- function(model, crop_code, var_name) {

  sm <- summary(model)$s.table %>%
    as.data.frame() %>%
    rownames_to_column("smooth")

  smooth_name <- paste0("s(", var_name, "):crop", crop_code)

  sm %>%
    filter(smooth == smooth_name) %>%
    pull(`p-value`)
}

sig_code <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

extract_effects <- function(thr, reg, model, data) {

  d <- data %>% filter(threshold == thr, region == reg)

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
    region = as.character(reg),
    crop = as.character(base_dat$crop),
    threshold = thr,
    baseline_prop = as.numeric(p_base),
    temp_change = as.numeric(100 * (p_temp - p_base)),
    prec_change = as.numeric(100 * (p_prec - p_base)),
    temp_p = map_dbl(as.character(base_dat$crop), ~ get_smooth_p(model, .x, "mat")),
    prec_p = map_dbl(as.character(base_dat$crop), ~ get_smooth_p(model, .x, "map"))
  )
}

threshold_effects <- imap_dfr(r4_models, function(model, nm) {
  parts <- strsplit(nm, "_")[[1]]
  extract_effects(parts[1], parts[2], model, dat_r4)
}) %>%
  mutate(
    crop_name = recode(crop, !!!crop_labs),
    region_name = recode(region, !!!region_labs),
    threshold_label = paste0("Suitability ≥ ", threshold),
    temp_sig = sig_code(temp_p),
    prec_sig = sig_code(prec_p)
  )

threshold_table <- threshold_effects %>%
  transmute(
    Region = region_name,
    Crop = crop_name,
    Threshold = threshold,
    `Baseline suitability` = round(baseline_prop, 3),
    `Δ +1°C MAT (pp)` = round(temp_change, 2),
    `Δ +100 mm MAP (pp)` = round(prec_change, 2),
    `MAT p-value` = signif(temp_p, 3),
    `MAP p-value` = signif(prec_p, 3),
    `MAT significance` = temp_sig,
    `MAP significance` = prec_sig
  ) %>%
  arrange(Region, Crop, Threshold)

write.xlsx(threshold_table, file.path(tab_dir, "threshold_sensitivity.xlsx"), overwrite = TRUE)
write.xlsx(threshold_effects, file.path(tab_dir, "threshold_driver_effects.xlsx"), overwrite = TRUE)

# ============================================================
# 5. Main threshold-sensitivity figure
# ============================================================

region_order_short <- c("Southern", "Central", "Eastern", "Western", "Northern")
crop_order <- c("Sorghum", "Pearl millet", "Finger millet")

legend_order <- c(
  "Not statistically significant",
  "Suitability ≥ 0.6",
  "Suitability ≥ 0.7",
  "Suitability ≥ 0.8"
)

threshold_offsets <- c(
  "Suitability ≥ 0.6" = -0.18,
  "Suitability ≥ 0.7" =  0.00,
  "Suitability ≥ 0.8" =  0.18
)

region_y <- c(
  "Northern" = 1,
  "Western"  = 2,
  "Eastern"  = 3,
  "Central"  = 4,
  "Southern" = 5
)

plot_r4 <- threshold_effects %>%
  mutate(
    Region = recode(
      region_name,
      "Southern Africa" = "Southern",
      "Central Africa"  = "Central",
      "Eastern Africa"  = "Eastern",
      "Western Africa"  = "Western",
      "Northern Africa" = "Northern"
    ),
    Crop = factor(crop_name, levels = crop_order),
    Threshold = factor(
      threshold_label,
      levels = c("Suitability ≥ 0.6", "Suitability ≥ 0.7", "Suitability ≥ 0.8")
    ),
    y_base = region_y[Region],
    y_pos = y_base + threshold_offsets[as.character(Threshold)],
    temp_group = ifelse(temp_sig == "ns", "Not statistically significant", as.character(Threshold)),
    prec_group = ifelse(prec_sig == "ns", "Not statistically significant", as.character(Threshold)),
    temp_group = factor(temp_group, levels = legend_order),
    prec_group = factor(prec_group, levels = legend_order)
  )

base_theme_r4 <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_blank(),
    axis.text = element_text(colour = "grey25"),
    panel.grid.major.y = element_line(colour = "grey88", linetype = "dotted"),
    panel.grid.major.x = element_line(colour = "grey88"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

p_temp <- ggplot(plot_r4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = temp_change, y = y_pos, yend = y_pos, colour = temp_group), linewidth = 0.75) +
  geom_point(aes(x = temp_change, y = y_pos, colour = temp_group), size = 3) +
  facet_wrap(~ Crop, ncol = 1) +
  scale_y_continuous(breaks = region_y[region_order_short], labels = region_order_short, limits = c(0.5, 5.5)) +
  scale_colour_manual(values = threshold_cols, breaks = legend_order, drop = FALSE) +
  labs(title = "A. Temperature effect (+ 1°C)", x = "Change in suitable area (%)") +
  base_theme_r4 +
  theme(legend.position = "none", panel.border = element_rect(colour = "#D55E00", fill = NA, linewidth = 1.1))

p_prec <- ggplot(plot_r4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = prec_change, y = y_pos, yend = y_pos, colour = prec_group), linewidth = 0.75) +
  geom_point(aes(x = prec_change, y = y_pos, colour = prec_group), size = 3) +
  facet_wrap(~ Crop, ncol = 1) +
  scale_y_continuous(breaks = region_y[region_order_short], labels = region_order_short, limits = c(0.5, 5.5)) +
  scale_colour_manual(values = threshold_cols, breaks = legend_order, drop = FALSE) +
  labs(title = "B. Precipitation effect (+ 100 mm)", x = "Change in suitable area (%)") +
  base_theme_r4 +
  theme(panel.border = element_rect(colour = "#2166AC", fill = NA, linewidth = 1.1))

fig_threshold <- p_temp + p_prec +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Fig_threshold_sensitivity.png"), fig_threshold, width = 12.5, height = 7.5, dpi = 600)
ggsave(file.path(fig_dir, "Fig_threshold_sensitivity.pdf"), fig_threshold, width = 12.5, height = 7.5)

# ============================================================
# 6. Supplementary driver-dominance heatmap
# ============================================================

dat_driver <- dat_base %>%
  select(crop, region, time_num, mat, map, prop_0.4, prop_0.6, prop_0.7, prop_0.8) %>%
  pivot_longer(
    cols = starts_with("prop_"),
    names_to = "threshold",
    values_to = "prop"
  ) %>%
  mutate(
    threshold = gsub("prop_", "", threshold),
    prop_clamped = pmin(pmax(prop, 0.001), 0.999),
    logit_prop = qlogis(prop_clamped),
    crop = factor(as.character(crop), levels = c("SG", "PM", "FM")),
    region = factor(as.character(region), levels = region_order)
  ) %>%
  filter(threshold %in% thresholds_driver, !is.na(prop), !is.na(logit_prop))

driver_groups <- dat_driver %>%
  group_by(threshold, region) %>%
  group_split()

driver_names <- dat_driver %>%
  distinct(threshold, region) %>%
  arrange(threshold, region) %>%
  transmute(name = paste(threshold, region, sep = "_")) %>%
  pull(name)

driver_models <- setNames(map(driver_groups, fit_model), driver_names)

extract_driver_effects <- function(thr, reg, model, data) {

  d <- data %>% filter(threshold == thr, region == reg)

  base_dat <- d %>%
    group_by(crop) %>%
    summarise(
      mat = median(mat, na.rm = TRUE),
      map = median(map, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(crop = factor(crop, levels = levels(data$crop)))

  p_base <- plogis(predict(model, newdata = base_dat))
  p_temp <- plogis(predict(model, newdata = base_dat %>% mutate(mat = mat + 1)))
  p_prec <- plogis(predict(model, newdata = base_dat %>% mutate(map = map + 100)))

  tibble(
    region = as.character(reg),
    crop = as.character(base_dat$crop),
    threshold = thr,
    baseline = as.numeric(p_base),
    temp_change = as.numeric(100 * (p_temp - p_base)),
    prec_change = as.numeric(100 * (p_prec - p_base))
  )
}

driver_effects <- imap_dfr(driver_models, function(model, nm) {
  parts <- strsplit(nm, "_")[[1]]
  extract_driver_effects(parts[1], parts[2], model, dat_driver)
})

write.xlsx(driver_effects, file.path(tab_dir, "threshold_driver_dominance.xlsx"), overwrite = TRUE)

driver_plot <- driver_effects %>%
  mutate(
    crop = recode(crop, !!!crop_labs),
    region = recode(region,
                    North = "Northern",
                    West = "Western",
                    East = "Eastern",
                    Central = "Central",
                    South = "Southern"),
    driver = case_when(
      abs(temp_change) < 1 & abs(prec_change) < 1 ~ "Weak / insensitive",
      abs(temp_change) >= abs(prec_change) ~ "Temperature-sensitive",
      abs(prec_change) > abs(temp_change) ~ "Precipitation-sensitive"
    ),
    effect_sign = case_when(
      driver == "Weak / insensitive" ~ "0",
      driver == "Temperature-sensitive" & temp_change > 0 ~ "+",
      driver == "Temperature-sensitive" & temp_change < 0 ~ "−",
      driver == "Precipitation-sensitive" & prec_change > 0 ~ "+",
      driver == "Precipitation-sensitive" & prec_change < 0 ~ "−"
    ),
    region = factor(region, levels = c("Southern", "Central", "Eastern", "Western", "Northern")),
    threshold = factor(threshold, levels = thresholds_driver),
    crop = factor(crop, levels = crop_order),
    driver = factor(driver, levels = c("Temperature-sensitive", "Precipitation-sensitive", "Weak / insensitive"))
  )

fig_driver <- ggplot(driver_plot, aes(x = threshold, y = region, fill = driver)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = effect_sign), size = 4.2, fontface = "bold") +
  facet_wrap(~ crop, nrow = 1) +
  scale_fill_manual(values = c(
    "Temperature-sensitive" = "#D55E00",
    "Precipitation-sensitive" = "#2166AC",
    "Weak / insensitive" = "grey80"
  )) +
  labs(x = "Suitability threshold", y = NULL, fill = "Dominant sensitivity") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom", panel.grid = element_blank())

ggsave(file.path(supp_fig_dir, "SM_driver_dominance_heatmap.png"), fig_driver, width = 10.5, height = 6.2, dpi = 600)
ggsave(file.path(supp_fig_dir, "SM_driver_dominance_heatmap.pdf"), fig_driver, width = 10.5, height = 6.2)

cat("\nThreshold-sensitivity analysis complete.\n")

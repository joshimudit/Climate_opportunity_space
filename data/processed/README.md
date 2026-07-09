# Processed data

This folder contains processed tables and generated data products used by later scripts.

## Files

| File/folder | Description | Produced/used by |
|---|---|---|
| `ecocrop_climate.xlsx` | Climate summary table generated during EcoCrop processing. | Produced by `02_ecocrop_model_runs.R`; used for documentation and some figure workflows. |
| `ecocrop_suitability.xlsx` | EcoCrop suitability summary table. | Produced by `02_ecocrop_model_runs.R`; optionally used by `08_opportunity_space_figures.R`. |
| `threshold_metrics.xlsx` | Corrected area-weighted suitable-area metrics at multiple suitability thresholds. | Produced by `04_threshold_metrics.R`; used by `05_gam_analysis.R` and `06_threshold_sensitivity.R`. |
| `climate_for_gam.xlsx` | Climate variables matched to crop, region, and timestep for corrected GAM analyses. | Used by `05_gam_analysis.R` and `06_threshold_sensitivity.R`. |
| `ecocrop_rasters/` | Continuous EcoCrop suitability rasters for sorghum, pearl millet, and finger millet. | Produced by `02_ecocrop_model_runs.R`; used by scripts `04`, `07`, and `08`. |
| `parameter_sensitivity/` | Processed Monte Carlo parameter-sensitivity outputs. | Produced by `03_parameter_sensitivity.R`. |

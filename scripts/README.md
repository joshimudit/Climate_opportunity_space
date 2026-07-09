# Scripts

This folder contains the R scripts used to generate the climate opportunity-space analyses, EcoCrop outputs, sensitivity analyses, figures, tables, and archaeobotanical evaluation.

Run scripts from the repository root so that relative paths such as `data/processed/` and `outputs/figures/` work correctly.

## Script order

| Script | Purpose | Main inputs | Main outputs |
|---|---|---|---|
| `01_climate_data_prep.R` | Prepares Africa-region climate inputs from CHELSA-TraCE21k climate data. | External CHELSA-TraCE21k rasters; Africa region boundaries | Processed climate rasters and summary tables |
| `02_ecocrop_model_runs.R` | Runs EcoCrop suitability models for sorghum, pearl millet, and finger millet. | Processed monthly temperature and precipitation rasters | Crop suitability rasters and EcoCrop summary tables |
| `03_parameter_sensitivity.R` | Tests sensitivity of EcoCrop outputs to crop-parameter uncertainty. | Climate rasters and crop parameters | Parameter-sensitivity tables and supplementary figures |
| `04_threshold_metrics.R` | Recalculates area-weighted suitable-area metrics from continuous EcoCrop rasters. | Continuous EcoCrop suitability rasters; Africa region boundaries | `data/processed/threshold_metrics.xlsx` |
| `05_gam_analysis.R` | Runs the corrected climate-driver GAM analysis using `prop_0.8`. | `threshold_metrics.xlsx`; `climate_for_gam.xlsx` | Corrected GAM figure and tables |
| `06_threshold_sensitivity.R` | Tests whether climate sensitivity changes with suitability cutoff. | `threshold_metrics.xlsx`; `climate_for_gam.xlsx` | Threshold-sensitivity figure and table |
| `07_arch_evaluation.R` | Compares modelled suitability with archaeobotanical crop evidence from southern Africa. | Restricted archaeobotanical data; EcoCrop rasters; southern Africa boundary | Archaeobotanical evaluation tables and figures |
| `08_opportunity_space_figures.R` | Creates climate opportunity-space maps and summary figures. | EcoCrop suitability rasters; Africa region boundaries | Main and supplementary opportunity-space figures/tables |

## Notes

- Scripts are intended to document the final analysis workflow, not every exploratory test.
- Some scripts require large generated raster files in `data/processed/ecocrop_rasters/`.
- `07_arch_evaluation.R` requires a restricted archaeobotanical input file that is not included in the public repository.
- CHELSA-TraCE21k input climate data are not included here and must be downloaded separately from the CHELSA data portal.

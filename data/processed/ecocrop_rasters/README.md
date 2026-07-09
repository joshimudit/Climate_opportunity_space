# EcoCrop suitability rasters

This folder contains continuous EcoCrop suitability rasters generated for each crop and timestep.

These rasters are a key generated data product of the workflow.

## Folder structure

```text
ecocrop_rasters/
├── SG/
├── PM/
└── FM/
```

## Crop codes

| Code | Crop |
|---|---|
| `SG` | Sorghum |
| `PM` | Pearl millet |
| `FM` | Finger millet |

## File naming

Raster files follow this pattern:

```text
<CROP>_suit_t<TIME>.tif
```

Examples:

```text
SG_suit_t0.tif
PM_suit_t60.tif
FM_suit_t120.tif
```

Here, `t0` is the present-day timestep and larger timestep values represent older Holocene conditions in 100-year increments.

## Notes

- These rasters are used by:
  - `04_threshold_metrics.R`
  - `07_arch_evaluation.R`
  - `08_opportunity_space_figures.R`

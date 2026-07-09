# Region boundaries

This folder contains vector boundary files used by the analysis scripts.

## Files

| File set | Used for | Notes |
|---|---|---|
| `africa_UN_regions.*` | Africa regional summaries, climate processing, threshold metrics, and opportunity-space figures. | Shapefile components must stay together: `.shp`, `.dbf`, `.shx`, `.prj`, and any sidecar files. |
| `Southern Africa Merged.*` | Southern Africa archaeobotanical evaluation. | Used by `07_arch_evaluation.R` to sample background points and map archaeobotanical sites. |

## Important

For shapefiles, all sidecar files must have exactly the same base name. For example:

```text
Southern Africa Merged.shp
Southern Africa Merged.dbf
Southern Africa Merged.shx
Southern Africa Merged.prj
```

If one component has a typo or different spelling, R may fail to read the shapefile correctly.

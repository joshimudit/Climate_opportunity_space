# Parameter sensitivity outputs

This folder contains processed outputs from the EcoCrop parameter-sensitivity analysis.

## Files

| File/folder | Description |
|---|---|
| `mc_long.rds` | Long-format Monte Carlo output used for plotting and summary. |
| `summary.rds` | Summary statistics from the Monte Carlo sensitivity analysis. |
| `effects.rds` | Estimated parameter effects on suitability outputs. |
| `claims.rds` | Derived claim/interpretation table from the sensitivity analysis. |
| `raw/` | Raw Monte Carlo objects such as baseline outputs and sampled parameter sets. |

## Notes

The `raw/` subfolder contain large files. Hence are not uploaded, while re-running the script they should be automatically generated in the local folder.

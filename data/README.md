# Data

This folder contains small input data, processed tables, and generated data products used by the analysis scripts.

Large external raw climate data are not included. Restricted archaeobotanical site-level data are also not included.

## Folder contents

| Folder | Contents | Public in repository? |
|---|---|---|
| `raw/` | Small spatial boundary files used by the scripts. | Yes |
| `processed/` | Processed climate/EcoCrop tables and generated model outputs. | Yes, except files marked as large or sensitive |
| `restricted/` | Restricted archaeobotanical input data needed for `07_arch_evaluation.R`. | No |

## External climate data

The original CHELSA-TraCE21k climate rasters are not stored in this repository because they are large external data. They should be downloaded from the CHELSA data portal:

- CHELSA-TraCE21k model page: <https://www.chelsa-climate.org/models/chelsa-trace21k>
- CHELSA-TraCE21k-centennial dataset page: <https://www.chelsa-climate.org/datasets/chelsa-trace21k-centennial>

The scripts expect processed climate rasters/tables derived from CHELSA-TraCE21k, not the full raw external archive.

## Restricted archaeobotanical data

The raw archaeobotanical workbook used in `07_arch_evaluation.R` is not included because the dataset has not yet been openly published/licensed. It contains site-level crop occurrence information, locations, and chronology.

The data can be made available upon reasonable request to the corresponding author of the paper. Exact rerunning of the archaeobotanical evaluation requires access to this restricted input file.

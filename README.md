# Breast GeoMx immune landscapes

This project examines breast tumor regions classified as Desert, Excluded or Inflamed. The dataset contains 155 areas of illumination (AOIs) from 8 patients, profiled in tumor, fibroblast and immune compartments.

| Compartment | Desert | Excluded | Inflamed |
| --- | ---: | ---: | ---: |
| Tumor | 11 | 42 | 34 |
| Fibroblast | 4 | 14 | 10 |
| Immune | 0 | 22 | 18 |

All analysis scripts are in [`code/`](code), numbered in execution order. Start from the supplied processed matrices and metadata. Run `code/00_install_environment.R` to install and verify the locked environment, then start the analyses with script 01. See [RUNNING.md](RUNNING.md) for environment setup, the nine required input files and the complete run sequence.

The workflow produces compartment-specific expression figures, selected pathway figures, and T-cell and CAF composition panels.

Large original GeoMx and single-cell objects are not required for the downstream workflow and are not included. Processed inputs and reference files are supplied separately. Analysis definitions are recorded in the script comments; validation and interpretation limits are described in [VALIDATION.md](VALIDATION.md).

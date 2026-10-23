# Breast GeoMx immune landscapes

This project examines breast tumor regions classified as Desert, Excluded or Inflamed. The dataset contains 155 areas of illumination (AOIs) from 8 patients, profiled in tumor, fibroblast and immune compartments.

| Compartment | Desert | Excluded | Inflamed |
| --- | ---: | ---: | ---: |
| Tumor | 11 | 42 | 34 |
| Fibroblast | 4 | 14 | 10 |
| Immune | 0 | 22 | 18 |

All analysis scripts are in [`code/`](code), numbered in execution order. The five processed data and reference files are included at the repository root. Run `code/00_install_environment.R` to install and verify the locked environment, then start the analyses with script 01. See [RUNNING.md](RUNNING.md) for environment setup, input locations and the four external GMT files and the complete run sequence.

The workflow produces compartment-specific expression figures, selected pathway figures, and T-cell and CAF composition panels.

## Input files

Download or clone the complete repository and open `GeoMx_revised.Rproj`. Keep these included files at the project root, alongside `renv.lock`:

| Included file | Use |
| --- | --- |
| `counts.csv` | Compartment expression analyses and deconvolution |
| `meta_cleaned_v5.xlsx` | AOI, patient, compartment and subtype annotations |
| `q3norm_nonlog.csv.gz` | Linear Q3 expression for deconvolution; read directly by script 11 |
| `ref_scAtlasVAE.rds` | Compact T-cell reference matrix |
| `CAF_reference.rds` | Compact CAF reference matrix |

The Q3 matrix is losslessly compressed; no manual extraction is needed. Script 11 also accepts the original `q3norm_nonlog.csv`.

These files are sufficient for scripts 01–03, 05–06 and 11 after environment installation. GSEA additionally requires the four MSigDB 2025.1 GMT files listed in [gene_sets/README.md](gene_sets/README.md), placed in `gene_sets/`. GMT files are not bundled. File checksums are recorded in [environment/input_checksums.csv](environment/input_checksums.csv).

Large original GeoMx and single-cell objects are not required for this downstream workflow and are not included. Analysis definitions are recorded in the script comments; validation and interpretation limits are described in [VALIDATION.md](VALIDATION.md).

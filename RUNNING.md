# Running the analysis

This repository contains code, environment records, and five processed data/reference files. Only the four MSigDB GMT files must be obtained separately for GSEA.

## 1. Prepare the project

Clone or download the complete repository into a writable directory. Open `GeoMx_revised.Rproj` in RStudio. Keep the working directory at the repository root, where `GeoMx_revised.Rproj` and `renv.lock` are located; do not change into `code/`. A private repository requires access through the collaborator's GitHub account.

On Windows, install R 4.5.1 and [Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html) in their default locations, then restart RStudio and select R 4.5.1. Use R 4.5.1 for the recorded environment. The lockfile specifies Bioconductor 3.21 and the package versions. RStudio does not install or select the required R version automatically.

In a fresh R session, run:

```r
source("code/00_install_environment.R")
```

Script 00 installs `renv` 1.1.5 into the project, restores all 184 packages from `renv.lock`, selects the project library, verifies every recorded package version and loads the analysis dependencies. It requires no data files and does not modify the analysis definitions or lockfile. Continue only after `Environment ready` is printed. A local version table and session record are saved in `renv/setup-report/`.

Run the same command at the start of each later R session. Already matching packages are reused. If installation fails, correct the reported network or system dependency problem and rerun 00. Restart R first if it reports that a different package version is already loaded. Do not run `renv::snapshot()` or update individual packages for this reproduction.

| Component | Required version |
| --- | --- |
| R | 4.5.1 |
| Bioconductor | 3.21 |
| renv | 1.1.5 |
| edgeR | 4.6.3 |
| limma | 3.64.1 |
| fgsea | 1.34.0 |
| clusterProfiler | 4.16.0 |
| SpatialDecon | 1.18.0 |
| ggplot2 | 3.5.2 |

The complete package list is in `renv.lock`; `environment/package_versions.csv` lists all 184 locked package versions. The installer checks R and Cairo support before package installation. R itself and operating-system dependencies are not installed by this script. Install R 4.5.1 separately and select it in RStudio; see the [CRAN R archive](https://cran.r-project.org/).

Installation requires access to CRAN and Bioconductor and write access to the project library/cache. Archived package versions may require source compilation: Windows needs the matching Rtools45 toolchain, macOS needs command-line tools and any required Fortran toolchain, and Linux needs compilers and package development libraries. A missing system library must be installed through the operating system. Package version locking does not freeze operating-system libraries or guarantee bit-identical numerical results across platforms. The Windows R 4.5.1/Rtools45 restore, fresh-session verification and synthetic-data core-script test have passed on a Windows Server 2022 runner. See `VALIDATION.md` for the execution record and scope.

## 2. Input locations

| Location relative to the project root | Required for |
| --- | --- |
| `counts.csv` | Expression analyses and deconvolution |
| `meta_cleaned_v5.xlsx` | AOI, patient, slide, compartment and subtype annotations |
| `q3norm_nonlog.csv.gz` | Deconvolution; losslessly compressed linear Q3 expression |
| `ref_scAtlasVAE.rds` | Deconvolution; small T-cell reference matrix |
| `CAF_reference.rds` | Deconvolution; small CAF reference matrix |
| `gene_sets/h.all.v2025.1.Hs.symbols.gmt` | GSEA; Hallmark collection |
| `gene_sets/c2.all.v2025.1.Hs.symbols.gmt` | GSEA; C2 collection |
| `gene_sets/c5.all.v2025.1.Hs.symbols.gmt` | GSEA; C5 collection |
| `gene_sets/c7.all.v2025.1.Hs.symbols.gmt` | GSEA; C7 collection |

The five data/reference files are included at the repository root and must remain there. Script 11 reads `q3norm_nonlog.csv.gz` directly, or uses the original `q3norm_nonlog.csv` when present; no manual extraction is needed. Obtain the exact MSigDB 2025.1 human-symbol GMT files through the normal MSigDB access process; details are in `gene_sets/README.md`. `environment/input_checksums.csv` records the nine expected file identities. Do not substitute a newer GMT release or a log-transformed Q3 matrix.

The two reference RDS files are small matrices, not the original single-cell objects. `geomx_obj_processed.RDS`, the original large single-cell objects and CellPhoneDB tables are not required. The workflow starts after acquisition and probe QC and does not reconstruct those upstream steps.

## 3. Run the analyses

Script 00 prepares the software environment; the remaining numbered files are analysis scripts. There is no separate preprocessing or QC step before script 01. Each expression script reads the counts and metadata itself; input alignment checks remain in the analysis code. Some scripts clear the R environment, so save unrelated objects first.

Run these commands in order from the project root:

```r
source("code/01_tumor_TMM_voom_heatmap_PCA.R")
source("code/02_fibroblast_TMM_voom_heatmap_PCA.R")
source("code/03_immune_TMM_voom_heatmap_PCA.R")
source("code/05_tumor_gene_overlap.R")
source("code/06_prepare_GSEA_rankings.R")
source("code/07_run_GSEA.R")
source("code/08_tumor_GSEA_heatmap.R")
source("code/09_fibroblast_GSEA_heatmap.R")
source("code/10_selected_GSEA_NES_bars.R")
source("code/11_Tcell_CAF_Q3_deconvolution.R")
```

| Scripts | Inputs and dependencies | Main outputs |
| --- | --- | --- |
| 00 | R 4.5.1 and `renv.lock`; no data | Verified project environment |
| 01-03 | Counts and metadata | Three DE tables, three marker heatmaps and three PCA plots |
| 05 | Tumor DE table from 01 | Tumor gene-overlap schematic |
| 06 | DE tables from 01-03 | Thirteen gene-ranking tables in `enrichment_inputs/` |
| 07 | Rankings from 06 and four GMT files | Fifty-two GSEA calculations, combined tables and checkpoints |
| 08-10 | GSEA outputs from 07 | Tumor/fibroblast pathway heatmaps and selected NES bars |
| 11 | Counts, metadata, linear Q3 and two small references | T-cell and CAF composition panels |

For expression figures alone, scripts 01-03 and 05 need only `counts.csv` and `meta_cleaned_v5.xlsx`. Script 11 can run independently of DE and GSEA once its five inputs are present. TMM factors and voom weights are calculated within scripts 01-03; no separate TMM matrix is needed.

## 4. Outputs and reruns

Most PDFs are written to the project root. Selected NES bars are written to `quick_plots/`. The deconvolution outputs are `Tcell_IE_II_scaled_abundance_bubble.pdf` and `CAF_ID_IE_II_scaled_abundance_bubble.pdf`.

Keep the intermediate CSV tables: later scripts read them directly. A rerun replaces the script's named outputs. Do not copy unrelated historical output folders into a fresh checkout.

GSEA runs four workers using the Windows socket backend or the Unix multicore backend. Its checkpoints are tied to input hashes and package versions. Re-run script 07 to resume compatible checkpoints; to start fresh, move `gsea_results/` outside the project first. Completion is recorded in `gsea_results/RUN_STATUS.txt` as `COMPLETED: 13 ranking definitions x 4 GMT collections.` Review the per-collection warning files. Numerical Monte Carlo results need not be bit-identical across operating systems or parallel backends.

## Validation boundary

The supplied-data analysis has prior execution records described in `VALIDATION.md`. GitHub Actions checks R syntax and, when the installer, lockfile or analysis scripts change, restores the full environment on Windows and repeats verification in a fresh R session. It also runs the released core scripts on synthetic inputs from a path containing Chinese characters and a space, and exercises the Windows socket GSEA backend. The environment job does not use study data; the separate supplied-data checks are recorded in `VALIDATION.md`. These checks do not establish that a particular collaborator's workstation has working repository access, input files or system dependencies.

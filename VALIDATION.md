# Validation

The 2026-10-21 workflow was executed using the supplied processed inputs under R 4.5.1. It completed the three compartment analyses, 52 fresh GSEA calculations, T-cell and CAF deconvolution. GSEA warning logs were retained; the 18 pathways displayed in the two selected heatmaps had complete estimates. This validation does not reconstruct upstream acquisition QC or demonstrate a clean Windows environment restore.

The 2026-10-22 reorganization places all analysis scripts in `code/`. The three-compartment calculation retains its existing statistical definitions. The 2026-10-22 compatibility and equivalence checks are recorded in the development repository.

## Interpretation

- DE models use `Slide_number` as the repeated-AOI block. Slide and patient identifiers are one-to-one in the supplied cohort. This equivalence must be checked before applying the scripts to another dataset.
- Heatmap markers and overlap sets use nominal P < 0.05 as exploratory selection criteria. Full BH FDR values remain in the differential-expression tables.
- The processed matrix contains 14,833 biological genes and `NegProbe-WTX`. The existing DE/ranking implementation retains all supplied features. The negative probe is a technical control and is not interpreted as a biological marker; changing the tested feature universe requires a documented rerun.
- Deconvolution proportions are relative to the respective reference matrix. T-cell and CAF panels have separate denominators and are not absolute tissue cell counts.

The final handoff check on 2026-10-22 verified the active script paths, their input/output dependencies, direct package coverage in the lockfile and the supplied input hashes. The separate input-check script is not required and is omitted. The syntax workflow discovers the current R scripts without a fixed historical file count. These changes do not modify the analysis scripts or their numerical results. The collaborator's own workstation remains untested.

## Environment installer

On 2026-10-22 (America/Phoenix), the installer was executed with an official R 4.5.1 distribution on Ubuntu 24.04.3, starting with an empty project library and package cache. All 184 locked versions were verified, including seven matching recommended packages supplied with R. All 19 analysis dependencies listed in script 00 loaded successfully. A second, fresh R session repeated the restore and verification successfully. The package report and session record are retained in the development repository under `review/environment_restore_20261022/`.

The local test compiled packages from their recorded source versions. System development libraries were provided in an isolated local directory because the execution container could not use the normal system-package installation path. This verifies the R environment in that Linux container; it is not a Windows or macOS installation test.

The ten analysis scripts and all 184 locked package versions are unchanged. CRAN repository labels were aligned with the repository declared in the lockfile. The installer verifies matching R-supplied recommended packages as well as packages installed in the project library.

## Supplied-data check in the restored environment

A fresh R session completed scripts 01, 02, 03, 05, 06 and 11 using the nine checksum-matched inputs. It generated nine PDFs and thirteen ranking tables. All nine PDFs were opened and rendered for visual inspection. The three DE tables each contain 14,834 features; their numerical values agree with the preceding validation run to within 1e-15.

One tumor Excluded ranking was tested against Hallmark 2025.1 using the settings in script 07. All 50 pathways returned finite NES, P, adjusted P and q values. fgsea reported possible overestimation of some P-values and values below its default 1e-10 estimation floor. This check did not repeat all 52 GSEA calculations or scripts 08-10; the earlier full-workflow validation remains recorded above.

## Windows target platform

The [Windows verification run](https://github.com/YaominLi/GeoMx_revised-public/actions/runs/39882078663) completed successfully on 2026-10-22 (America/Phoenix) using Windows Server 2022 x64, R 4.5.1 UCRT and Rtools45. The clean installation verified all 184 locked package versions and loaded all 19 analysis dependencies. Verification also passed in fresh R sessions. BiocManager 1.30.26 is bootstrapped before project activation so that the installed and loaded versions agree.

The Windows job executed the unchanged scripts 01, 02, 03, 05, 06 and 11 on synthetic counts, metadata and reference matrices from a working path containing Chinese characters and a space. It generated nine nonempty PDFs and thirteen ranking tables, and completed a two-set GSEA using `SnowParam(4, type = "SOCK")`, the backend selected by script 07 on Windows. The test did not upload study data. The Windows package-version table and session records are retained in the development repository under `review/windows_validation_20261022/`.

This establishes the tested environment and core-script execution on a Windows runner. It does not represent a full study-cohort rerun on Windows, execution of scripts 08-10, all 52 study GSEA calculations, or validation of the collaborator's own workstation. The Linux supplied-data checks above provide separate evidence for the numerical outputs.

## Bundled processed inputs

The five processed data/reference files are included in the repository. The Q3 matrix is stored as `q3norm_nonlog.csv.gz`; its decompressed bytes match the original CSV checksum. Under R 4.5.1, `read.csv()` returned identical objects from the compressed and original files (14,834 features by 155 AOIs). Script 11 accepts both forms and preserves the same calculations. All R scripts passed syntax parsing after the input-path update. This distribution check did not repeat the full analysis or the Windows environment test.

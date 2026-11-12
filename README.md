# Histopathology-anchored single-cell RNA sequencing of symptomatic and asymptomatic human carotid atherosclerotic plaques

Analysis code for the carotid endarterectomy (CEA) single-cell RNA sequencing (scRNA-seq) manuscript.

This repository contains the R code used to analyse scRNA-seq data from human carotid atherosclerotic plaques, stratified by blinded histopathological grading combined with clinical presentation. It does not contain the sequencing data itself: the data are deposited in the NCBI Gene Expression Omnibus (GEO), separate from this code repository.

> [!NOTE]
> This README describes the nine scripts in this repository (cleaned UTF-8 versions uploaded 2026-11-12; the original GBK uploads were verified 2026-11-11). Two placeholders remain on purpose while the manuscript is unpublished: the GEO Series accession (assigned after GEO finishes validating the submission) and the Zenodo DOI (created when the code is archived at publication time). Fill those in when they exist; everything else in this file is complete. See "Known limitations" before publishing: the scripts depend on intermediate `.rds` files that are not in this repository.

## Study design

Single-cell RNA sequencing was performed on CEA specimens from 12 patients. Each specimen was split longitudinally: one half was used for scRNA-seq and the other half was fixed in 10% neutral-buffered formalin, paraffin-embedded, sectioned at 4 μm, and assessed with haematoxylin and eosin (H&E) and Sirius Red histopathology. Two pathologists blinded to symptom status classified the plaques with a modified American Heart Association (AHA) framework. Only clinico-pathologically concordant specimens were included.

| Group | Samples | Histopathology |
|---|---|---|
| Asymptomatic | ASYM1–ASYM6 | AHA type IV–V, non-complicated plaques |
| Symptomatic | SYM1–SYM6 | AHA type VI, complicated plaques with rupture or healed-rupture features |

No other sample is part of this dataset. In particular, sample N1 is not included in this study, in the GEO submission, or in this code.

## What the code does

The repository holds one script per figure. Each script starts from an intermediate Seurat or CellChat object (see "Required intermediate files") and produces the panels of one manuscript figure.

| Script | What it does |
|---|---|
| `Figure 2.R` | Pathology-stratified single-cell atlas: quality control, clustering and cell-type annotation of the integrated atlas object, comparison of complicated versus non-complicated lesions, and integration with the published Alsaigh 2022 dataset (GSE159677). Produces `integrated_obj_with_Alsaigh.rds`. |
| `Figure 3 and Supplement Figure 3.R` | Vascular smooth muscle cell-to-reparative fibroblast continuum: subsetting of SMC/fibroblast states, monocle3 pseudotime trajectory, contractile versus reparative module scores, and transcription-factor progression. |
| `Figure 4.R` | Rare cycling proliferating reparative fibroblasts (PRFs): pseudotime analysis of the fibroblast/PRF subset and its abundance in complicated versus non-complicated lesions. |
| `Figure 5.R` | TREM2+ mature lipid-associated macrophages (LAMs) along the monocyte-to-macrophage continuum: donor-level pseudobulk comparison and trajectory validation. Sources the shared helper file `00_common.R` (not in this repository; see "Known limitations"). |
| `Figure 6.R` | Three-tiered cell–cell communication architecture (CellChat): comparison of interaction number and strength between asymptomatic and symptomatic lesions, pathway ranking, and centrality analysis. |
| `Figure 7.R` | Clinical translation of the TREM2+ LAM signature: single-sample GSEA scoring of the signature in independent bulk transcriptomic cohorts (GSE100927, GSE28829, GSE43292) with ROC/AUC analysis. Downloads GSE28829 itself through GEOquery. |
| `Supplemental Figure 2.R` | Atlas quality panels: top-1000 variable-gene expression across cell types (raincloud, half-violin and sina plots) for the integrated atlas. Uses the `combined131` object from the R session, so run it in the same session as `Figure 2.R` (or load `combined131.rds` first). |
| `Supplemental Figure 4.R` | Macrophage donor panels, trajectory validation and GO enrichment (clusterProfiler) for foamy versus mature macrophages. Sources the shared helper file `00_common.R` (not in this repository; see "Known limitations"). |
| `Supplemental Figure 5.R` | CellChat centrality and temporal/donor validation of the communication architecture. |

## Data availability

The data generated in this study are deposited in GEO:

- GEO Series accession: **\[TO COMPLETE: add the GSE accession once GEO assigns it\]**
- Status at the time of writing: the submission has been sent and the metadata accepted; GEO is validating the raw files, after which the Series will be queued for processing and a GSE accession assigned. The record is private until 30 June 2027.
- Samples: 12 (ASYM1–ASYM6, SYM1–SYM6)

Public datasets reanalysed alongside the newly generated data:

- GSE159677 (Alsaigh 2022; downloaded and integrated by `Figure 2.R`)
- GSE100927 (bulk cohort; enters `Figure 7.R` as a pre-computed score file, see below)
- GSE28829 (bulk cohort; downloaded by `Figure 7.R` through GEOquery)
- GSE43292 (bulk cohort; enters `Figure 7.R` as pre-computed score and paired-sample files, see below)

## Input data

### Newly generated data (from GEO)

| Input | Files | Notes |
|---|---|---|
| Raw reads | 24 FASTQ files | Paired-end, 150 bp; two files per sample |
| Count matrices | 36 files | `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` for each of the 12 samples |
| Sample information | 1 table | Sample ID (ASYM1–SYM6, SYM1–SYM6), group, and AHA type for each sample |

Sequencing and primary processing, for reference: fresh plaque tissue was dissociated into a single-cell suspension, libraries were prepared with the SeekOne MM Single Cell 3' kit (SeekGene, cat. no. K00104) on a SeekOne MM chip, and sequencing was performed on an Illumina NovaSeq 6000 (paired-end, 150 bp). Reads were processed with SeekOne tools version 1.0.0, chemistry seekone v1, against the GRCh38 reference transcriptome.

### Required intermediate files (not in this repository)

The scripts do not start from the GEO count matrices. Each one starts from an intermediate object that was produced during the analysis and is not stored here. Before the code can be run, these files must be supplied (ideally by adding the upstream script that builds them from the GEO matrices; see "Known limitations").

| File | Needed by |
|---|---|
| `combined131.rds` — the integrated Seurat object used for the atlas (working name from the original analysis) | `Figure 2.R`, `Figure 3 and Supplement Figure 3.R`, `Supplemental Figure 2.R` (in-session) |
| `smc.rds` — SMC/fibroblast subset | `Figure 3 and Supplement Figure 3.R` |
| `fib_prf_v3_pseudotime.rds`, `fib_prf_v3_cds.rds` — fibroblast/PRF trajectory objects | `Figure 4.R`, `Supplemental Figure 5.R` |
| `macro.rds` — macrophage subset | `Figure 5.R`, `Supplemental Figure 4.R` |
| `macro_with_pseudotime_3grp.rds` — macrophage trajectory object | `Supplemental Figure 4.R`, `Supplemental Figure 5.R` |
| `cellchat/cellchat_asym.rds`, `cellchat/cellchat_sym.rds` — CellChat objects | `Figure 6.R`, `Supplemental Figure 5.R` |
| `Fig6_GSE100927_scores.csv`, `Fig6_GSE43292_scores.csv`, `GSE43292_paired_info.csv` — pre-computed bulk-cohort scores | `Figure 7.R` |
| `integrated_obj_with_Alsaigh.rds` — produced by `Figure 2.R` | `Supplemental Figure 2.R` |
| `00_common.R` — shared settings/helpers sourced by `Figure 5.R` and `Supplemental Figure 4.R` | `Figure 5.R`, `Supplemental Figure 4.R` |

## Repository structure

```text
CEA-scRNA-seq/
├── README.md
├── LICENSE
├── session_info.txt
├── Figure 2.R
├── Figure 3 and Supplement Figure 3.R
├── Figure 4.R
├── Figure 5.R
├── Figure 6.R
├── Figure 7.R
├── Supplemental Figure 2.R
├── Supplemental Figure 4.R
└── Supplemental Figure 5.R
```

> [!NOTE]
> The file was renamed on upload (2026-11-12) from `Figure 3 and Supplement Figure 3 .R` — the original name had a stray space before `.R` — to `Figure 3 and Supplement Figure 3.R`.

There is no `data/` folder in this repository: the count matrices come from GEO, and the intermediate objects listed above are not stored here.

## Requirements

- R: 4.3.0 (2023-04-21), Windows 11 x64 (verified from `sessionInfo()` run on 2026-11-12)
- Reference genome: GRCh38 (annotation version is not stated in the primary-processing report)

### R packages

The table lists every package the nine scripts load (`library()`/`require()` or `::`). Versions filled in below come from `sessionInfo()` run on the analysis machine on 2026-11-12 (R 4.3.0, after running Figure 2, so it covers the packages that session loaded). Two rows (`pheatmap`, `BiocManager`) record that the version was not captured on the analysis machine; run `packageVersion()` for them there and fill the versions in. `UCell` records that the package is not installed on that machine. Do not type versions from memory. Alternatively commit an `renv.lock` produced by `renv::snapshot()`.

| Package | Version | Used in |
|---|---|---|
| BiocManager | Not recorded — run `packageVersion("BiocManager")` on the analysis machine | Figure 7 |
| CellChat | 1.6.1 | Figure 6, Supplemental Figure 5 |
| clusterProfiler | 4.8.3 | Supplemental Figure 4 |
| dplyr | 1.1.2 | All nine scripts |
| edgeR | 3.42.4 | Figure 4, Figure 5, Supplemental Figure 4 |
| fgsea | 1.26.0 | Figure 5 |
| forcats | 1.0.1 | Figure 3 |
| GEOquery | 2.68.0 | Figure 2, Figure 7 |
| ggh4x | 0.3.1 | Figure 7 |
| ggplot2 | 4.0.3 | All nine scripts |
| ggbeeswarm | 0.7.3 | Supplemental Figure 2 |
| ggdist | 3.3.3 | Supplemental Figure 2 |
| gghalves | 0.1.4 | Supplemental Figure 2 |
| ggrastr | 1.0.2 | Figure 2, Figure 3, Figure 4, Supplemental Figure 2, Supplemental Figure 4 |
| ggrepel | 0.9.8 | Figure 2, Figure 6, Supplemental Figure 4, Supplemental Figure 5 |
| ggridges | 0.5.7 | Supplemental Figure 5 |
| ggsignif | 0.6.4 | Figure 4, Figure 6, Supplemental Figure 2 |
| ggpubr | 0.6.3 | Figure 4 |
| GSVA | 1.48.3 | Figure 7 |
| harmony | 2.0.3 | Figure 2 |
| igraph | 2.3.1 | Figure 3, Supplemental Figure 4 |
| Matrix | 1.6-5 | Figure 2, Figure 4, Figure 5, Supplemental Figure 4, Supplemental Figure 5 |
| mgcv | 1.8-42 | Figure 3, Supplemental Figure 5 |
| monocle3 | 1.3.1 | Figure 3, Figure 4, Supplemental Figure 4 |
| org.Hs.eg.db | 3.17.0 | Supplemental Figure 4 |
| patchwork | 1.3.2 | All nine scripts |
| pheatmap | Not recorded — run `packageVersion("pheatmap")` on the analysis machine | Figure 2, Supplemental Figure 2 |
| pROC | 1.19.0.1 | Figure 7 |
| purrr | 1.0.1 | Supplemental Figure 5 |
| RColorBrewer | 1.1-3 | Figure 2, Supplemental Figure 2 |
| readr | 2.2.0 | Figure 6, Supplemental Figure 5 |
| scales | 1.4.0 | Figure 3, Figure 4, Figure 6 |
| Seurat | 5.5.0 | Figure 2, Figure 3, Figure 4, Figure 5, Supplemental Figure 2, Supplemental Figure 4, Supplemental Figure 5 |
| SeuratWrappers | 0.4.0 | Figure 4 |
| SingleR | 2.2.0 | Figure 2 |
| tidydr | 0.0.6 | Figure 2, Figure 3 |
| tidyr | 1.3.2 | Figure 2, Figure 3, Figure 4, Figure 5, Figure 6, Supplemental Figure 4, Supplemental Figure 5 |
| tidyverse | 2.0.0 | Figure 2, Figure 3, Figure 5, Figure 6, Figure 7, Supplemental Figure 4, Supplemental Figure 5 |
| Biobase | 2.60.0 | Figure 7 |
| SeuratObject | 5.4.0 | Figure 5, Supplemental Figure 4 |
| tibble | 3.2.1 | Figure 4, Figure 7, Supplemental Figure 4, Supplemental Figure 5 |
| rstudioapi | 0.18.0 | Figure 5, Supplemental Figure 4 |
| UCell | Not installed on the analysis machine (2026-11-12: `packageVersion("UCell")` errors with "there is no package called 'UCell'"); install with `BiocManager::install("UCell")` before running Supplemental Figure 5, then record the installed version here | Supplemental Figure 5 |
| viridis | 0.6.5 | Figure 3 |

Install the packages the code needs with:

```r
# CRAN packages (BiocManager itself is on CRAN; Matrix ships with R):
install.packages(c("BiocManager", "CellChat", "dplyr", "forcats", "ggplot2", "ggbeeswarm", "ggdist", "gghalves", "ggh4x", "ggrastr", "ggrepel", "ggridges", "ggsignif", "ggpubr", "igraph", "mgcv", "patchwork", "pheatmap", "pROC", "purrr", "RColorBrewer", "readr", "rstudioapi", "scales", "Seurat", "SeuratObject", "tibble", "tidydr", "tidyr", "tidyverse", "viridis"))
# Bioconductor packages:
BiocManager::install(c("Biobase", "clusterProfiler", "edgeR", "fgsea", "GEOquery", "GSVA", "harmony", "monocle3", "org.Hs.eg.db", "SingleR", "UCell"))
# SeuratWrappers is not on CRAN or Bioconductor; install it from GitHub:
# remotes::install_github("satijalab/seurat-wrappers")
```

## How to run

The scripts are per-figure scripts, not a single linear pipeline: each one starts from an intermediate object listed under "Required intermediate files". Run them in this order only once those objects are in place (each script takes its project root from the working directory; see step 2).

1. Clone this repository and open it as your working directory.
2. Place the intermediate `.rds` files (see "Required intermediate files") in one folder, open R with that folder as the working directory, and source the scripts from there: each script takes its project root from the working directory (`cea.project_root`, set at the top of every script), so all input and output paths resolve under it.
3. Install the R packages listed above.
4. Run the scripts in this order, each in a fresh R session unless noted otherwise:
   1. `Figure 2.R` — builds `integrated_obj_with_Alsaigh.rds` (also downloads GSE159677 itself).
   2. `Supplemental Figure 2.R` — in the same R session as `Figure 2.R`, because it uses the `combined131` object from that session.
   3. `Figure 3 and Supplement Figure 3.R`
   4. `Figure 4.R`
   5. `Figure 5.R` — requires `00_common.R` (not in this repository).
   6. `Figure 6.R`
   7. `Figure 7.R`
   8. `Supplemental Figure 4.R` — requires `00_common.R` (not in this repository).
   9. `Supplemental Figure 5.R`
5. Check the output folders for the figures and tables. Each script writes its figures (PDF/PNG) next to its outputs; see the `ggsave()` calls at the end of each section.

## Reproducibility

- Keep the input file names exactly as they are in GEO. The scripts identify samples by the ASYM and SYM sample IDs in those names.
- Record the environment that produced the published results. Either commit an `renv.lock` file created with `renv::snapshot()`, or add a `session_info.txt` file created with the command below, and refer to it here.

```r
writeLines(capture.output(sessionInfo()), "session_info.txt")
```

- This repository uses `session_info.txt` (committed here, captured on the analysis machine on 2026-11-12 after running `Figure 2.R`). It does not use `renv`.

## Known limitations

Read this before publishing or making the repository public. These points were found by checking the nine uploaded scripts on 2026-11-11.

1. **Absolute paths (fixed in the uploaded files).** The original uploads used absolute Windows paths rooted at `G:/13sample/` (75 occurrences; `setwd()` in 7 of the 9 scripts). The files now in this repository replace them with a single project-root setting at the top of each script (`options(cea.project_root = normalizePath(".", mustWork = FALSE))`); all input and output paths are built from it, so run the scripts with the working directory set to the folder holding the intermediate files. No `G:/` path remains in the uploaded scripts (verified 2026-11-12).
2. **Missing upstream step.** No script in this repository builds the starting object `combined131.rds` from the GEO count matrices, and no script produces `smc.rds`, `macro.rds`, the `fib_prf_v3_*` and `macro_with_pseudotime_3grp` trajectory objects, the `cellchat_*` objects, or the Figure 7 score CSVs. Add the upstream quality-control/integration script that creates `combined131.rds` (and the scripts that derive the other objects), or state in this README where a reader can obtain them. Until then the repository reproduces the figures only from those intermediate objects, not from the GEO data alone.
3. **Missing helper file.** `Figure 5.R` and `Supplemental Figure 4.R` both `source("00_common.R")` and refer to `01_donor_panels.R` / `03_trajectory_validation.R`; none of these files is in the repository. Add `00_common.R` (or inline the helpers it provides) before publication.
4. **File encoding (fixed in the uploaded files).** The original uploads were GBK-encoded with CRLF line endings, so their Chinese comments displayed as garbled text on GitHub. The files now in this repository are UTF-8 with LF line endings (converted 2026-11-12, verified by strict re-decoding with no characters lost).
5. **Two package versions are not yet recorded, and one package is not installed.** `pheatmap` and `BiocManager` versions are not recorded (the `packageVersion()` loop stopped at the UCell error before reaching them). `UCell` is not installed on the analysis machine — Supplemental Figure 5 cannot run there until it is installed (`BiocManager::install("UCell")`); the user decided on 2026-11-12 not to install it for now, so the UCell row in the package table records that instead of a version.
6. **Unresolved path conflict for the fibroblast/PRF trajectory files.** `Figure 4.R` and `Supplemental Figure 5.R` refer to `fib_prf_v3_pseudotime.rds` under two different folders (`Fig3/Fig3_v3/...` in some places, `Fig3_v3/...` in others). Only one of them exists on the analysis machine; confirm which one and make the paths consistent before running those scripts.
7. **Sample count inside `combined131.rds` is not verified.** The working folder was called `13sample` and the object `combined131`, which suggests the object may hold 13 samples, while this study has 12 (N1 excluded). The code text does not state the sample count, and the object itself is not in this repository, so confirm that `combined131.rds` contains only the 12 study samples before relying on it.

## Data and code privacy

This repository contains analysis code only. The statements below were checked against the nine uploaded scripts on 2026-11-11, and the cleaned files now in this repository were verified to differ from them only in path lines, the added header, and encoding:

- It contains no passwords, access tokens, or other credentials (no such strings were found in any of the nine scripts).
- It contains no patient-identifiable information in the code text: no patient names, no sequencing IDs, and no ASYM/SYM sample IDs appear in the scripts; samples are referred to only by group (ASYM/SYM) and cell-type labels. (This check covers the code text only, not the contents of the `.rds` objects, which are not in this repository.)
- It contains no sequencing data, count matrices, or analysis result objects. Those belong in GEO, not in this repository.
- The scripts use a single configurable project root (`cea.project_root`, set at the top of each script) instead of absolute paths (see "Known limitations" point 1). The original uploads contained absolute `G:/13sample/` paths; none remain in the files in this repository (verified 2026-11-12).

Before making this repository public, fix points 1–4 under "Known limitations" and re-check the files for `.Rhistory`, `.RData`, and other session files, which can carry the same problems.

## Licence

This code is released under the MIT Licence. See `LICENSE` in this repository.

## Citation

**\[TO COMPLETE: after the code is archived on Zenodo, add the version DOI here and delete this note.\]**

Please cite the archived release of this code:

> \[Authors\]. Analysis code for "Histopathology-anchored single-cell RNA sequencing of symptomatic and asymptomatic human carotid atherosclerotic plaques". Zenodo. DOI: \[TO COMPLETE: 10.5281/zenodo.XXXXXXX\]

Code availability statement for the manuscript:

> Analysis code is available at https://github.com/gaoopoo/CEA-scRNA-seq and archived on Zenodo (DOI: \[TO COMPLETE: 10.5281/zenodo.XXXXXXX\]).

# GSAuto: A Faster, Cross-Platform R Package for Infinium Methylation BeadChip Data Analyses

## Overview

**GSAuto** is an open-source R package designed to process raw IDAT files from Illumina methylation beadchip arrays (27k, 450k, 850k/EPIC v1.0, 935k/EPIC v2.0, and MSA).  
It provides fast and scalable extraction of probe intensities, beta values, detection p-values, and quality metrics, optimized for large-scale cohort studies.  

Compared to Illumina GenomeStudio, GSAuto is cross-platform (Windows, macOS, Linux), automated, lightweight, and reproducible.

## Features

- Supports all major Illumina Infinium methylation platforms: **27k, 450k, 850k, 935k, and MSA**  
- Directly processes `.idat` files to generate:
  - **Beta values** (DNA methylation levels)  
  - **Detection P-values** (probe reliability)  
  - **Detection rates** (per-sample QC)  
  - **Intensity matrices** (methylated and unmethylated probe signals)  
- Built-in manifests and negative control probe sets for all platforms  
- Cross-platform, scriptable, and suitable for high-throughput pipelines  

## Installation

### Option 1: Install from GitHub (recommended)

```r
# install devtools if not already installed
install.packages("devtools")

# install GSAuto from GitHub
devtools::install_github("zzfcharlie/GSAuto")

# load package
library(GSAuto)
```

### Option 2: Source scripts directly (basic use)

```r
source("R/GSAuto_core.R")
```

## Prerequisites

Before use, please ensure the following R packages are installed:

```r
install.packages(c("illuminaio", "dplyr", "data.table"))
```

## Usage Examples

### Basic IDAT processing

```r
# Run GSAuto on IDAT files
res <- gsauto_run_from_idat(idat_dir = "extdata/GEO10223_idat", platform = "450k")

# Extract beta values
beta_mat <- res$beta_matrix

# Extract detection P-values
pval_mat <- res$pval_matrix
```

### Generate QC report

```r
gsauto_qc_report(output_dir = "extdata/GEO10223_output")
```

### Differential methylation analysis

```r
dmp_res <- gsauto_dmp(beta_matrix, group_labels)
```

## Output

The main outputs include:  
- `beta_matrix`: Matrix of beta values per sample  
- `pval_matrix`: Detection P-value matrix  
- `detect_re`: Sample-level QC summary  
- `M_matrix`: Methylated intensities  
- `U_matrix`: Unmethylated intensities  
- `Intensity_matrix`: Total probe intensities  

## File Structure

The repository includes platform-specific manifests and control probe sets:

```
data/
  ├── 27k/
  │   ├── manifest.rds
  │   ├── negControl.rds
  ├── 450k/
  │   ├── manifest.rds
  │   ├── negControl.rds
  ├── 850k/
  │   ├── manifest.rds
  │   ├── negControl.rds
  ├── 935k/
  │   ├── manifest.rds
  │   ├── negControl.rds
  └── MSA/
      ├── manifest.rds
      ├── negControl.rds
```

Example IDATs and outputs are in `extdata/`.

## Citation

If you use **GSAuto** in your research, please cite:

> Weng H, Zheng Z, Tang K, et al.  
> GSAuto: A Faster, Cross-Platform R tool for Infinium methylation beadchips data analyses.  
> *Bioinformatics Advances*. 2025.  

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

# ===============================
# Test script for GSAuto package
# ===============================

# Load all GSAuto functions (assumed to be in the same folder or sourced)
source("R/GSAuto_core.R")
source("R/GSAuto_qc.R")
source("R/GSAuto_run_from_idat.R")
source("R/GSAuto_dmp.R")

# Define test input directory
idat_dir <- "./extdata/GEO10223_idat"
platform <- "MSA"  # Change to "935k" or others as needed

# Run full pipeline
results <- gsauto_run_from_idat(
  idat_dir = idat_dir,
  platform = platform,
  run_filter = TRUE,
  run_qc = TRUE,
  sample_cutoff = 0.9,
  probe_cutoff = 0.9,
  manifest_path = file.path("data", platform, "manifest.rds")
)

# Check filtered results
print(names(results))
print(results$out_dir)

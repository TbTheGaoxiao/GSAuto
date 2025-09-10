# ===================== End-to-end runner =====================

#' Full pipeline: IDAT → matrices → (filter) → (QC report)
#'
#' @param idat_dir Directory containing *_Red.idat and *_Grn.idat files
#' @param platform One of "27k", "450k", "850k", "935k", or "MSA"
#' @param run_filter Whether to perform filtering (default TRUE)
#' @param run_qc Whether to generate a QC report (default TRUE)
#' @param sample_cutoff Detection rate threshold per sample (default 0.9)
#' @param probe_cutoff Detection rate threshold per probe (default 0.9)
#' @param manifest_path Optional path to manifest.rds (for probe type summary in log)
#'
#' @return A list with `core`, `filter`, and `out_dir`; filter is NULL if filtering is skipped
#' @export
gsauto_run_from_idat <- function(idat_dir,
                                 platform      = "450k",
                                 run_filter    = TRUE,
                                 run_qc        = TRUE,
                                 sample_cutoff = 0.9,
                                 probe_cutoff  = 0.9,
                                 manifest_path = NULL) {
  
  # Clean output directory name
  clean_dir  <- sub("/+$", "", idat_dir)
  base_name  <- basename(clean_dir)
  output_dir <- file.path(dirname(clean_dir), sub("_idat$", "_output", base_name))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Step 1: Extract matrices from IDAT files
  res_core <- GSAuto(idat_dir, platform = platform)
  
  data.table::fwrite(cbind(Probe_ID = rownames(res_core$beta_matrix),
                           res_core$beta_matrix),
                     file.path(output_dir, "beta_mat.csv"))
  data.table::fwrite(cbind(Probe_ID = rownames(res_core$pval_matrix),
                           res_core$pval_matrix),
                     file.path(output_dir, "detP_mat.csv"))
  data.table::fwrite(res_core$detect_re,
                     file.path(output_dir, "detection_rate.csv"))
  data.table::fwrite(
    data.table::as.data.table(res_core$Intensity_matrix, keep.rownames = "Probe_ID"),
    file.path(output_dir, "Intensity_matrix.csv")
  )
  
  message("✓ Matrices written to ", output_dir)
  
  # Step 2: Apply filtering (optional)
  res_filter <- NULL
  if (run_filter) {
    res_filter <- gsauto_run_filter(output_dir,
                                    sample_cutoff = sample_cutoff,
                                    probe_cutoff  = probe_cutoff,
                                    manifest_path = manifest_path,
                                    write_log     = TRUE)
  }
  
  # Step 3: Generate QC report (optional)
  if (run_qc) {
    gsauto_qc_report(output_dir, open_browser = interactive())
  }
  
  # Return structured results
  invisible(list(
    core   = res_core,
    filter = res_filter,
    out_dir = output_dir
  ))
}
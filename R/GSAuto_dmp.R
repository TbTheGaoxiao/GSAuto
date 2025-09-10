# ===================== 6. Differential Methylation Analysis =====================

#' Perform differential methylation analysis (DMP) between case and control groups
#'
#' @param case_csv CSV file with beta values for case samples (first column = Probe_ID)
#' @param ctrl_csv CSV file with beta values for control samples (same format)
#' @param delta_cut Minimum deltaBeta threshold (default 0, i.e., no filter)
#' @param fdr_cut Adjusted p-value cutoff (default 0.05)
#' @param out_dir Output directory (default = same as case_csv)
#' @param verbose Whether to print status messages (default TRUE)
#' @param strip_suffix Whether to strip suffix from Probe_IDs (default TRUE)
#' @export
gsauto_dmp <- function(case_csv,
                       ctrl_csv,
                       delta_cut = 0,
                       fdr_cut   = 0.05,
                       out_dir   = NULL,
                       verbose   = TRUE,
                       strip_suffix = TRUE) {
  
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")
  if (!requireNamespace("limma", quietly = TRUE))
    stop("Package 'limma' is required.")
  
  if (is.null(out_dir)) out_dir <- dirname(path.expand(case_csv))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_csv <- file.path(out_dir, "dmp_result.csv")
  
  # Load and average replicates
  read_beta <- function(file) {
    dt <- data.table::fread(file)
    ids <- if (strip_suffix) sub("_.*$", "", dt[[1]]) else dt[[1]]
    dt[, Probe_ID := ids]
    dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = Probe_ID]
    data.frame(dt, check.names = FALSE)
  }
  
  dt_case <- read_beta(case_csv)
  dt_ctrl <- read_beta(ctrl_csv)
  
  common <- intersect(dt_case$Probe_ID, dt_ctrl$Probe_ID)
  if (length(common) == 0)
    stop("No common Probe_IDs found after merging case and control tables.")
  
  dt_case <- dt_case[match(common, dt_case$Probe_ID), ]
  dt_ctrl <- dt_ctrl[match(common, dt_ctrl$Probe_ID), ]
  
  beta_mat <- cbind(dt_case[, -1], dt_ctrl[, -1])
  rownames(beta_mat) <- common
  beta_mat <- as.matrix(beta_mat)
  
  group <- factor(c(rep("case", ncol(dt_case) - 1),
                    rep("ctrl", ncol(dt_ctrl) - 1)))
  
  # Convert beta to M-values for analysis
  M <- log2(beta_mat / (1 - beta_mat) + 1e-3)
  fit <- limma::eBayes(limma::lmFit(M, model.matrix(~ group)))
  tt  <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  
  Ave_case <- rowMeans(beta_mat[, group == "case", drop = FALSE], na.rm = TRUE)
  Ave_ctrl <- rowMeans(beta_mat[, group == "ctrl", drop = FALSE], na.rm = TRUE)
  
  res <- data.frame(
    Probe_ID     = rownames(tt),
    logFC        = tt$logFC,
    AveBeta_case = Ave_case,
    AveBeta_ctrl = Ave_ctrl,
    deltaBeta    = Ave_case - Ave_ctrl,
    P.Value      = tt$P.Value,
    adj.P.Val    = tt$adj.P.Val,
    check.names  = FALSE
  )
  
  res <- res[abs(res$deltaBeta) >= delta_cut & res$adj.P.Val <= fdr_cut, ]
  res <- res[order(res$adj.P.Val), ]
  
  data.table::fwrite(res, out_csv)
  if (verbose) {
    message("✓ DMP result saved to ", out_csv, " (", nrow(res), " probes)")
  }
  
  invisible(res)
}
# ===================== 1. Core filtering =====================

#' Filter low-quality probes and samples
#'
#' @param beta_mat  Beta value matrix or data.table (first column is Probe_ID)
#' @param detP_mat  Detection P-value matrix (same shape as beta_mat)
#' @param sample_cutoff Minimum proportion of detected probes per sample (default 0.9)
#' @param probe_cutoff  Minimum proportion of samples passed per probe (default 0.9)
#' @return A list containing filtered matrices and removal logs
#' @export
gsauto_filter <- function(beta_mat, detP_mat,
                          sample_cutoff = 0.9,
                          probe_cutoff  = 0.9) {
  
  if (!all(dim(beta_mat) == dim(detP_mat)))
    stop("beta_mat and detP_mat must have the same dimensions")
  
  probe_ids <- beta_mat[[1]]
  beta_mat  <- as.data.frame(beta_mat[, -1])
  detP_mat  <- as.data.frame(detP_mat[, -1])
  rownames(beta_mat) <- rownames(detP_mat) <- probe_ids
  
  sample_pass <- colMeans(detP_mat <= 0.05, na.rm = TRUE)
  keep_samples <- sample_pass >= sample_cutoff
  removed_samples <- names(sample_pass)[!keep_samples]
  
  probe_pass <- rowMeans(detP_mat[, keep_samples, drop = FALSE] <= 0.05, na.rm = TRUE)
  keep_probes <- probe_pass >= probe_cutoff
  removed_probes <- rownames(detP_mat)[!keep_probes]
  
  filt_beta <- beta_mat[keep_probes, keep_samples, drop = FALSE]
  filt_detP <- detP_mat[keep_probes, keep_samples, drop = FALSE]
  
  filt_beta <- cbind(Probe_ID = rownames(filt_beta), filt_beta)
  filt_detP <- cbind(Probe_ID = rownames(filt_detP), filt_detP)
  
  list(
    filtered_beta = filt_beta,
    filtered_detP = filt_detP,
    removed_samples = removed_samples,
    removed_probes  = removed_probes
  )
}

# ===================== 2. Filtering and export =====================

#' Perform filtering and export results to disk
#'
#' @param output_dir Directory containing beta_mat.csv and detP_mat.csv
#' @param sample_cutoff Sample detection rate threshold (default 0.9)
#' @param probe_cutoff Probe pass rate threshold (default 0.9)
#' @param write_log Whether to write filtering log (default TRUE)
#' @param manifest_path Optional path to manifest.rds for probe type summary
#' @export
gsauto_run_filter <- function(output_dir,
                              sample_cutoff = 0.9,
                              probe_cutoff  = 0.9,
                              write_log     = TRUE,
                              manifest_path = NULL) {
  
  requireNamespace("data.table", quietly = TRUE)
  
  beta_path <- file.path(output_dir, "beta_mat.csv")
  detp_path <- file.path(output_dir, "detP_mat.csv")
  if (!file.exists(beta_path) || !file.exists(detp_path))
    stop("Missing beta_mat.csv or detP_mat.csv in ", output_dir)
  
  beta_dt <- data.table::fread(beta_path)
  detp_dt <- data.table::fread(detp_path)
  
  res <- gsauto_filter(beta_dt, detp_dt,
                       sample_cutoff = sample_cutoff,
                       probe_cutoff  = probe_cutoff)
  
  data.table::fwrite(res$filtered_beta, file.path(output_dir, "beta_mat_filtered.csv"))
  data.table::fwrite(res$filtered_detP, file.path(output_dir, "detP_mat_filtered.csv"))
  
  if (write_log) {
    log_f <- file.path(output_dir, "filtered_log.txt")
    cat("Removed samples:", length(res$removed_samples),  "\n",
        res$removed_samples, sep = "\n", file = log_f)
    cat("\nRemoved probes:",  length(res$removed_probes), "\n",
        head(res$removed_probes, 5), sep = "\n", file = log_f, append = TRUE)
    
    if (!is.null(manifest_path) && file.exists(manifest_path)) {
      man <- readRDS(manifest_path)
      man_sub <- man[match(res$removed_probes, man$Probe_ID), ]
      cat("\n\nRemoved probes by category\n", file = log_f, append = TRUE)
      cat_tbl <- table(
        ifelse(is.na(man_sub$col), "Type-II",
               paste0("Type-I (", man_sub$col, ")"))
      )
      utils::write.table(cat_tbl, file = log_f,
                         append = TRUE, col.names = FALSE,
                         quote = FALSE)
    }
  }
  
  message("✓ Filtered matrices saved to: ", output_dir)
  invisible(res)
}

# ===================== 3. QC Report =====================

#' Generate an HTML QC report with summary statistics and plots
#'
#' @param output_dir Directory containing beta_mat.csv and detP_mat.csv
#' @param report_name Name of the output HTML file (default "GSAuto_QC_Report.html")
#' @param open_browser Whether to open the report in browser (default interactive())
#' @export

gsauto_qc_report <- function(output_dir,
                             report_name  = "GSAuto_QC_Report.html",
                             open_browser = interactive()) {
  
  required <- c("data.table", "ggplot2", "pheatmap", "rmarkdown", "knitr", "grid")
  missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  
  beta_path <- if (file.exists(file.path(output_dir, "beta_mat_filtered.csv")))
    file.path(output_dir, "beta_mat_filtered.csv") else
      file.path(output_dir, "beta_mat.csv")
  
  detp_path <- if (file.exists(file.path(output_dir, "detP_mat_filtered.csv")))
    file.path(output_dir, "detP_mat_filtered.csv") else
      file.path(output_dir, "detP_mat.csv")
  
  beta_dt <- data.table::fread(beta_path)
  detp_dt <- data.table::fread(detp_path)
  
  sample_ids <- colnames(beta_dt)[-1]
  beta_mat   <- as.matrix(beta_dt[, -1])
  detp_mat   <- as.matrix(detp_dt[, -1])
  
  det_rate   <- colMeans(detp_mat <= 0.05, na.rm = TRUE)
  sample_cor <- if (ncol(beta_mat) >= 2) stats::cor(beta_mat, use="pairwise.complete.obs") else NULL
  
  intensity_path <- file.path(output_dir, "Intensity_matrix.csv")
  has_intensity  <- file.exists(intensity_path)
  if (has_intensity) {
    int_dt   <- data.table::fread(intensity_path)
    int_mat  <- as.matrix(int_dt[, -1])
    mean_int <- colMeans(int_mat, na.rm = TRUE)
  }
  
  log_path <- normalizePath(file.path(output_dir, "filtered_log.txt"),
                            winslash = "/", mustWork = FALSE)
  
  # Rmd 内容
  rmd_lines <- c(
    "---",
    "title: \"GSAuto QC Report\"",
    "output:",
    "  html_document:",
    "    theme: flatly",
    "    toc: false",
    "---",
    "<div class='main-box'>",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo=FALSE, message=FALSE, warning=FALSE,",
    "  fig.align='center', fig.width=8, fig.height=4, out.width='90%', dpi=120)",
    "library(ggplot2); library(pheatmap); library(knitr); library(grid)",
    "```",
    "",
    "# Summary",
    "```{r}",
    "summary_df <- data.frame(Metric = c('Samples', 'Probes', 'Mean detection rate', 'Min detection rate'),",
    "  Value = c(ncol(beta_mat), nrow(beta_mat), round(mean(det_rate), 3), round(min(det_rate), 3)))",
    "kable(summary_df, row.names = FALSE)",
    "```",
    "",
    "# Detection rate per sample",
    "```{r results='asis'}",
    "dr <- data.frame(Sample = names(det_rate), Rate = det_rate)",
    "if (all(dr$Rate >= 0.99)) {",
    "  knitr::asis_output('*All samples have detection rate ≥ 0.99. Skipping plot.*')",
    "} else {",
    "  print(ggplot(dr, aes(Sample, Rate)) +",
    "    geom_col(fill = '#3182bd') +",
    "    geom_hline(yintercept = 0.9, colour = 'red', linetype = 2) +",
    "    coord_flip() +",
    "    scale_y_continuous(limits = c(0.85, 1)) +",
    "    labs(x = NULL, y = 'Detection rate') +",
    "    theme_bw() + theme(axis.text.y = element_text(size = 7)))",
    "}",
    "```",
    "",
    "# Total intensity per sample",
    if (has_intensity) {
      c(
        "```{r}",
        "ti <- data.frame(Sample = sample_ids, Intensity = log2(mean_int))",
        "print(ggplot(ti, aes(Sample, Intensity)) +",
        "  geom_col(fill = '#a6cee3') +",
        "  geom_hline(yintercept = median(ti$Intensity) - 1, linetype = 2, colour = 'red') +",
        "  coord_flip() + labs(x = NULL, y = 'log2(mean M+U)') +",
        "  theme_bw() + theme(axis.text.y = element_text(size = 7)))",
        "```"
      )
    } else {
      c(
        "```{r results='asis'}",
        "knitr::asis_output('*Intensity_matrix.csv not found — section skipped.*')",
        "```"
      )
    },
    "",
    "# Beta value distribution",
    "```{r}",
    "beta_df <- as.data.frame(beta_mat); colnames(beta_df) <- sample_ids",
    "if (ncol(beta_df) > 40) beta_df <- beta_df[, 1:40, drop = FALSE]",
    "boxplot(beta_df, outline = FALSE, las = 2, ylab = 'Beta', main = 'Beta value distribution')",
    "```",
    "",
    "# PCA of samples (top 10k probes)",
    "```{r}",
    "top_idx <- order(apply(beta_mat, 1, var, na.rm = TRUE), decreasing = TRUE)[1:min(10000, nrow(beta_mat))]",
    "pca_res <- prcomp(t(beta_mat[top_idx, ]), scale. = TRUE)",
    "pca_df <- data.frame(pca_res$x[, 1:2], Sample = sample_ids)",
    "ggplot(pca_df, aes(PC1, PC2, label = Sample)) +",
    "  geom_point(size = 2.5, colour = '#1f78b4') +",
    "  geom_text(vjust = 1.2, hjust = 0.5, size = 3) +",
    "  theme_bw() + labs(x = 'PC1', y = 'PC2')",
    "```",
    "",
    "# Inter-sample Pearson correlation",
    "```{r, fig.width=8, fig.height=6}",
    "if (!is.null(sample_cor)) {",
    "  ph <- pheatmap::pheatmap(sample_cor,",
    "                           fontsize_col = 7,",
    "                           fontsize_row = 7,",
    "                           main = 'Inter-sample Pearson correlation',",
    "                           silent = TRUE)",
    "  grid::grid.newpage()",
    "  grid::grid.draw(ph$gtable)",
    "} else {",
    "  knitr::asis_output('*Only one sample available. Skipping heatmap.*')",
    "}",
    "```",
    "",
    "# Filtering log summary",
    "```{r results='asis'}",
    paste0("log_f <- \"", log_path, "\""),
    "if (!file.exists(log_f)) {",
    "  knitr::asis_output('*No filtering log found.*')",
    "} else {",
    "  txt <- readLines(log_f)",
    "  cat('<pre>', paste(txt, collapse='\\n'), '</pre>')",
    "}",
    "```",
    "</div>"
  )
  
  rmd <- tempfile(fileext = ".Rmd")
  writeLines(trimws(unlist(rmd_lines), "left"), rmd)
  
  out_path <- file.path(normalizePath(output_dir), report_name)
  rmarkdown::render(rmd, output_dir = dirname(out_path),
                    output_file = basename(out_path), quiet = TRUE)
  
  message("✓ QC report written to ", out_path)
  if (open_browser) utils::browseURL(out_path)
  invisible(out_path)
}

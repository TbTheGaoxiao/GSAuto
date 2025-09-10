# ===================== IDAT to matrix extraction =====================

#' Extract methylation signal matrices from IDAT files
#'
#' @param idat_dir Directory containing *.idat files (paired _Red.idat and _Grn.idat)
#' @param platform Platform identifier: "27k", "450k", "850k", "935k", "MSA"
#' @return A list with beta matrix, detection p-values, M/U signals, total intensity, detection rate and QC flag
#' @export
GSAuto_core <- function(idat_dir, platform = "450k") {
  
  suppressMessages(library(illuminaio))
  
  supported <- c("27k", "450k", "850k", "935k", "MSA")
  stopifnot(platform %in% supported)
  
  base <- file.path("data", platform)
  control_idx <- readRDS(file.path(base, "negControl.rds"))
  manifest    <- readRDS(file.path(base, "manifest.rds"))
  
  probe_ids <- manifest$Probe_ID
  n_probes  <- length(probe_ids)
  
  # Get sample stems (e.g., GSM1234_R01C01)
  idat_files <- list.files(idat_dir, pattern = "\\.idat$", full.names = FALSE)
  stems <- unique(sub("(_Grn|_Red)\\.idat$", "", idat_files, perl = TRUE))
  if (length(stems) == 0)
    stop("No .idat files found in ", idat_dir)
  
  message("Loading ", length(stems), " samples from ", idat_dir)
  
  n_samples <- length(stems)
  new_matrix <- function() {
    m <- matrix(NA_real_, nrow = n_probes, ncol = n_samples)
    rownames(m) <- probe_ids
    m
  }
  
  beta <- new_matrix()
  detP <- new_matrix()
  Mmat <- new_matrix()
  Umat <- new_matrix()
  Imat <- new_matrix()
  det_rate <- qc_flag <- numeric(n_samples)
  
  # Helper to validate paired IDATs
  pick <- function(stem, suffix) {
    f <- file.path(idat_dir, paste0(stem, suffix, ".idat"))
    if (!file.exists(f)) stop("Missing file: ", f)
    f
  }
  
  for (i in seq_along(stems)) {
    stem <- stems[i]
    idaG <- readIDAT(pick(stem, "_Grn"))
    idaR <- readIDAT(pick(stem, "_Red"))
    
    ida <- cbind(G = idaG$Quants[, "Mean"],
                 R = idaR$Quants[, "Mean"])
    rownames(ida) <- rownames(idaG$Quants)
    
    # Extract Type-I red/green and Type-II probes
    get_signals <- function(df, Mcol, Ucol) {
      M_idx <- match(df$M, rownames(ida))
      U_idx <- match(df$U, rownames(ida))
      keep <- !is.na(M_idx) & !is.na(U_idx) & !is.na(df$Probe_ID)
      data.frame(M = ida[M_idx[keep], Mcol],
                 U = ida[U_idx[keep], Ucol],
                 row.names = df$Probe_ID[keep])
    }
    
    IR <- get_signals(manifest[manifest$col == "R", ], "R", "R")
    IG <- get_signals(manifest[manifest$col == "G", ], "G", "G")
    
    II_ord <- manifest[is.na(manifest$col), ]
    II_idx <- match(II_ord$U, rownames(ida))
    II <- data.frame(M = ida[II_idx, "G"], U = ida[II_idx, "R"],
                     row.names = II_ord$Probe_ID)
    
    sset <- rbind(IR, IG, II)
    row_idx <- match(rownames(sset), probe_ids)
    
    # Write raw signals
    Mmat[row_idx, i] <- sset$M
    Umat[row_idx, i] <- sset$U
    Imat[row_idx, i] <- rowSums(sset)
    
    # Detection P-value calculation (GenomeStudio-like empirical method)
    neg <- ida[control_idx, ]
    ref_R <- sort(2 * neg[, "R"])
    ref_G <- sort(2 * neg[, "G"])
    ref_S <- sort(rowSums(neg))
    
    detP[row_idx, i] <- c(
      1 - findInterval(rowSums(IR), ref_R, left.open = TRUE) / length(ref_R),
      1 - findInterval(rowSums(IG), ref_G, left.open = TRUE) / length(ref_G),
      1 - findInterval(rowSums(II), ref_S, left.open = TRUE) / length(ref_S)
    )
    
    # Beta value calculation
    beta[, i] <- pmax(Mmat[, i], 0) / (abs(Mmat[, i]) + abs(Umat[, i]) + 100)
    
    det_rate[i] <- 1 - mean(detP[, i] > 0.05, na.rm = TRUE)
    qc_flag[i]  <- ifelse(det_rate[i] >= 0.96, "yes", "no")
  }
  
  colnames(beta) <- colnames(detP) <- colnames(Mmat) <-
    colnames(Umat) <- colnames(Imat) <- stems
  
  list(
    beta_matrix      = beta,
    pval_matrix      = detP,
    Intensity_matrix = Imat,
    detect_re        = data.frame(detection_rate = det_rate, qc = qc_flag, row.names = stems),
    M_matrix         = Mmat,
    U_matrix         = Umat
  )
}

#' Backward-compatible alias
#' @export
GSAuto <- function(idat_dir, platform = "450k") {
  GSAuto_core(idat_dir, platform)
}
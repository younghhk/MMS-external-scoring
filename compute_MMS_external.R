compute_MMS_external <- function(ext_df, IDATA_df, case_col = "case") {
  
  # IDATA_df must include:
  #   term    : metabolite ID used in IDATA
  #   beta    : coefficient from IDATA model
  #   matched : 1 if metabolite intended for external MMS, 0 otherwise
  
  stopifnot(all(c("term", "beta", "matched") %in% names(IDATA_df)))
  stopifnot(case_col %in% names(ext_df))
  
  # 1) keep intended metabolites and drop intercept
  IDATA_use <- subset(IDATA_df, term != "(Intercept)" & matched == 1)
  
  # 2) use only metabolites present in external dataset
  mets <- intersect(IDATA_use$term, names(ext_df))
  if (length(mets) == 0) {
    stop("No matched metabolites found in external dataset.")
  }
  
  # align coefficients to metabolite order
  IDATA_use <- IDATA_use[match(mets, IDATA_use$term), ]
  betas <- IDATA_use$beta
  
  # 3) extract metabolite matrix
  X <- as.matrix(ext_df[, mets, drop = FALSE])
  
  # 4) standardize using controls
  ctrl_idx <- which(ext_df[[case_col]] == 0)
  if (length(ctrl_idx) < 2) stop("Not enough controls to compute SD.")
  
  mu  <- colMeans(X[ctrl_idx, , drop = FALSE], na.rm = TRUE)
  sdv <- apply(X[ctrl_idx, , drop = FALSE], 2, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- NA_real_
  
  Z <- sweep(X, 2, mu, "-")
  Z <- sweep(Z, 2, sdv, "/")
  
  # 5) compute MMS
  MMS <- as.vector(Z %*% betas)
  MMS[!complete.cases(Z)] <- NA_real_
  
  message("MMS computed using ", length(mets),
          " metabolites; dropped ",
          sum(IDATA_df$matched == 1) - length(mets))
  
  return(MMS)
}

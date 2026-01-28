compute_MMS_external <- function(
  ext_df,
  IDATA_df,
  case_col = "case",
  use_controls = TRUE,
  use_idata_scaling = FALSE,
  na_handling = c("complete", "mean_impute"),
  # ONLY use Metabolon IDs for matching
  idata_id_cols = c("COMP_ID_IDATA", "CHEM_ID_IDATA"),
  idata_id_labels = c("COMP_ID", "CHEM_ID")
) {

  na_handling <- match.arg(na_handling)

  stopifnot(is.data.frame(ext_df), is.data.frame(IDATA_df))
  stopifnot(all(c("term", "beta") %in% names(IDATA_df)))

  # -------------------------------------------------
  # Remove intercept row (term used ONLY for this)
  # -------------------------------------------------
  IDATA_df <- IDATA_df[IDATA_df$term != "(Intercept)", , drop = FALSE]

  has_case_col <- case_col %in% names(ext_df)

  if (use_controls && !has_case_col && !use_idata_scaling) {
    warning(
      "case_col not found in ext_df. ",
      "Standardizing metabolites using all samples instead of controls."
    )
  }

  if (use_idata_scaling &&
      !all(c("mean_idata", "sd_idata") %in% names(IDATA_df))) {
    stop("use_idata_scaling = TRUE requires IDATA_df to have mean_idata and sd_idata.")
  }

  # -------------------------------------------------
  # Matching metabolites (NO term matching)
  # -------------------------------------------------
  ext_colnames <- as.character(names(ext_df))

  mapping <- data.frame(
    IDATA_row = seq_len(nrow(IDATA_df)),
    term = as.character(IDATA_df$term),
    matched_ext_col = NA_character_,
    matched_id_type = NA_character_,
    stringsAsFactors = FALSE
  )

  for (k in seq_along(idata_id_cols)) {
    idcol <- idata_id_cols[k]
    if (!idcol %in% names(IDATA_df)) next

    id_values <- as.character(IDATA_df[[idcol]])
    valid_idx <- which(!is.na(id_values) & id_values != "")

    for (i in valid_idx) {
      if (!is.na(mapping$matched_ext_col[i])) next
      if (id_values[i] %in% ext_colnames) {
        mapping$matched_ext_col[i] <- id_values[i]
        mapping$matched_id_type[i] <- idata_id_labels[k]
      }
    }
  }

  matched_rows <- which(!is.na(mapping$matched_ext_col))
  if (length(matched_rows) == 0) {
    stop(
      "No overlapping metabolites found between IDATA and ext_df. ",
      "Matching is performed using COMP_ID_IDATA and CHEM_ID_IDATA only."
    )
  }

  used_mapping <- mapping[matched_rows, , drop = FALSE]
  used_terms <- used_mapping$term
  used_ext_cols <- used_mapping$matched_ext_col

  betas <- IDATA_df$beta[used_mapping$IDATA_row]
  names(betas) <- used_terms

  X <- as.matrix(ext_df[, used_ext_cols, drop = FALSE])
  colnames(X) <- used_terms

  # -------------------------------------------------
  # Determine mean / SD (external dataset only)
  # -------------------------------------------------
  if (use_idata_scaling) {

    mu  <- IDATA_df$mean_idata[used_mapping$IDATA_row]
    sdv <- IDATA_df$sd_idata[used_mapping$IDATA_row]

  } else if (use_controls && has_case_col) {

    ctrl_idx <- which(as.numeric(ext_df[[case_col]]) == 0)
    if (length(ctrl_idx) < 2) {
      warning("Fewer than 2 controls found. Using all samples for standardization.")
      mu  <- colMeans(X, na.rm = TRUE)
      sdv <- apply(X, 2, sd, na.rm = TRUE)
    } else {
      mu  <- colMeans(X[ctrl_idx, , drop = FALSE], na.rm = TRUE)
      sdv <- apply(X[ctrl_idx, , drop = FALSE], 2, sd, na.rm = TRUE)
    }

  } else {

    mu  <- colMeans(X, na.rm = TRUE)
    sdv <- apply(X, 2, sd, na.rm = TRUE)

  }

  sdv[sdv == 0 | is.na(sdv)] <- NA_real_

  # -------------------------------------------------
  # Missing value handling
  # -------------------------------------------------
  if (na_handling == "mean_impute") {
    for (j in seq_len(ncol(X))) {
      miss <- is.na(X[, j])
      if (any(miss)) X[miss, j] <- mu[j]
    }
  }

  Z <- sweep(X, 2, mu, "-")
  Z <- sweep(Z, 2, sdv, "/")

  score <- as.vector(Z %*% betas[colnames(Z)])

  if (na_handling == "complete") {
    score[!complete.cases(Z)] <- NA_real_
  }

  details <- list(
    n_idata_terms = nrow(IDATA_df),
    n_matched = length(used_terms),
    mapping = used_mapping[, c("term", "matched_ext_col", "matched_id_type")],
    dropped_terms = IDATA_df$term[!(IDATA_df$term %in% used_terms)],
    scaling = if (use_idata_scaling) {
      "IDATA"
    } else if (use_controls && has_case_col) {
      "external_controls"
    } else {
      "external_all"
    }
  )

  return(list(
    score = score,
    used_metabolites = used_terms,
    dropped_metabolites = details$dropped_terms,
    details = details
  ))
}


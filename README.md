
# Computing MMS in External Datasets

This repository provides a standardized approach for computing **Metabolite-based Scores (MMS)** in external datasets using coefficients derived from the **IDATA** study.

In external datasets, MMS are computed using the subset of IDATA metabolites that are available on the external metabolomics platform. As a result, MMS values represent a **relative metabolite pattern within each study** and are intended for **within-study association analyses**, rather than direct numerical comparison across cohorts.

---

## Repository Contents

```
├── compute_MMS_external.R
├── IDATA_MMS_sodium_serum.xlsx
├── IDATA_MMS_sodium_urine.xlsx
├── IDATA_MMS_sodium_FMV.xlsx
├── IDATA_MMS_...
└── README.md
```

---

## Files

### `compute_MMS_external.R`

Contains a single reusable function:

```r
compute_MMS_external(ext_df, IDATA_df)
```

This function:

* Matches metabolites between IDATA and the external dataset using Metabolon identifiers
* Standardizes metabolites in the external dataset
* Computes the MMS as a weighted sum of standardized metabolites

The same function is used for all MMS; only the MMS-specific IDATA coefficient file (one per exposure and biospecimen) changes.

---

### MMS Coefficient Files (`IDATA_MMS_*.xlsx`)

One coefficient file is provided per MMS (e.g., sodium–serum, sodium–urine).

Each file must include the following columns:

| Column name     | Description                                                  |
| --------------- | ------------------------------------------------------------ |
| `term`          | Metabolite name used internally in IDATA                     |
| `beta`          | Regression coefficient from the IDATA MMS model              |
| `COMP_ID_IDATA` | Metabolon Compound ID (recommended for cross-study matching) |
| `CHEM_ID_IDATA` | Metabolon Chemical ID (alternative identifier)               |

Additional metadata columns (e.g., biochemical name, pathway) may be present and are ignored by the code.

> Rows with `term == "(Intercept)"` are automatically excluded.

---

## Matching Metabolites Across Datasets

Metabolon assigns multiple identifiers to metabolites. In IDATA, metabolites are stored with:

* a study-specific name (`term`)
* Metabolon identifiers (`COMP_ID_IDATA`, `CHEM_ID_IDATA`)

External datasets typically use **COMP ID** or **CHEM ID** as column names.

The function attempts to match metabolites using the following priority:

1. **COMP ID** (`COMP_ID_IDATA`)
2. **CHEM ID** (`CHEM_ID_IDATA`)
3. **term** (fallback only)

Only metabolites that can be successfully matched are used in the MMS. Unmatched metabolites are dropped automatically and reported in the output.

**Recommendation:**
External datasets should provide metabolite columns named by **COMP ID** (preferred) or **CHEM ID**.

---

## Requirements for the External Dataset

The external dataset (`ext_df`) must:

* Have **one row per participant**
* Contain metabolite columns named using **COMP ID**, **CHEM ID**, or (less ideally) `term`

A **case/control indicator is optional** and is only required if control-based standardization is desired.

---

## How MMS Is Computed

### 1. Metabolite selection

* Metabolites are matched automatically using Metabolon identifiers
* Only metabolites present in both IDATA and the external dataset are used

### 2. Standardization (External Dataset Only)

The IDATA models were developed using standardized metabolites, and the files provided here contain coefficients only.
To apply these coefficients, metabolites must be standardized within the external dataset.

By default, standardization is performed as follows:

- If the external dataset includes a case/control column:

  - Metabolites are standardized using controls only

- If the external dataset does not include a case/control column:

  - Metabolites are standardized using all samples

- If `use_controls = FALSE`:

  - Metabolites are standardized using all samples, regardless of whether a case/control column exists

### 3. Score calculation

MMS is computed as a weighted sum of standardized metabolites:

```
MMS_i = sum_j ( beta_j * Z_ij )
```

where `Z_ij` is the standardized value of metabolite `j` for participant `i`.

### 4. Missing data

* By default, if a participant is missing any metabolite used in the MMS, their MMS is set to `NA`
* Optionally, missing metabolite values can be mean-imputed prior to standardization

Because metabolite availability may differ across platforms, the MMS in each external dataset may be based on a **subset of the original IDATA metabolites**. The same IDATA-derived coefficients are applied to the available metabolites only.

---

## Example Usage

### External dataset **with** case/control

```r
library(readxl)

IDATA_df <- read_excel("IDATA_MMS_sodium_serum.xlsx")
ext_df   <- read.csv("external_case_control.csv")

res <- compute_MMS_external(
  ext_df   = ext_df,
  IDATA_df = IDATA_df,
  case_col = "case"   # used to standardize using controls
)

ext_df$MMS_sodium_serum <- res$score
```

### External dataset **without** case/control

```r
IDATA_df <- read_excel("IDATA_MMS_sodium_serum.xlsx")
ext_df   <- read.csv("external_cohort.csv")

# Automatically standardizes using all samples
res <- compute_MMS_external(
  ext_df   = ext_df,
  IDATA_df = IDATA_df
)

ext_df$MMS_sodium_serum <- res$score
```

### Force standardization using all samples

```r
res <- compute_MMS_external(
  ext_df        = ext_df,
  IDATA_df      = IDATA_df,
  use_controls  = FALSE
)
```

---

## Downstream Analysis Example (Case–Control)

```r
library(survival)

clogit(
  case ~ MMS_sodium_serum + covariates + strata(match_id),
  data = ext_df
)
```

---

## Notes and Assumptions

* MMS are computed using **IDATA-derived coefficients** and the subset of metabolites available in the external dataset
* Metabolite matching is performed programmatically using Metabolon identifiers
* MMS values are **study-specific** and should not be interpreted as numerically equivalent across cohorts
* MMS are intended to be used as **relative exposure variables** in within-study analyses
* Intercepts from IDATA models are intentionally excluded

---

## Contact

For questions or clarifications, please contact the study team.



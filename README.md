
# Computing MMS in External Nested Case–Control Studies

This repository provides a standardized approach for computing **Metabolite-based Scores (MMS)** in external nested case–control datasets using coefficients derived from the **IDATA** study.

In external datasets, MMS are computed using the subset of IDATA-selected metabolites that are available in the external platform. As a result, MMS values represent a relative metabolite pattern within each study and are intended for within-study association analyses rather than direct numerical comparison across cohorts.

---

## Repository Contents

```
├── compute_MMS_external.R
├── MMS_sodium_serum_stomachcancer.xlsx
├── MMS_sodium_urine_stomachcancer.xlsx
├── MMS_...
└── README.md
```

---

## Files

### `compute_MMS_external.R`

Contains a single reusable function:

```r
compute_MMS_external(ext_df, IDATA_df, case_col = "case")
```

This function:

* Matches metabolites between IDATA and the external dataset
* Standardizes metabolites using controls
* Computes the MMS as a weighted sum of metabolites

The same function is used for **all MMS**; only the coefficient file changes.

---

### MMS Coefficient Files (`MMS_*.xlsx`)

One coefficient file is provided per MMS (e.g., sodium–serum, sodium–urine).

Each file contains coefficients estimated in IDATA and must include the following columns:

| Column name | Description                                                                  |
| ----------- | ---------------------------------------------------------------------------- |
| `term`      | Metabolite ID / column name (as used in IDATA and expected in external data) |
| `beta`      | Regression coefficient from the IDATA MMS model                              |
| `matched`   | `1` = metabolite intended for external MMS; `0` = not intended / unavailable |

Additional metadata columns (e.g., biochemical name, pathway) may be present and are ignored by the code.

> Rows with `term == "(Intercept)"` are automatically excluded.

---

## Requirements for the External Dataset

The external nested case–control dataset (`ext_df`) must:

* Have **one row per participant**
* Contain metabolite columns named according to `term`
* Include a **case/control indicator**:

  * `case == 1` for cases
  * `case == 0` for controls
* (For downstream analysis) include a matching identifier (e.g., `match_id`) for conditional logistic regression

---

## How MMS Is Computed

For a given MMS:

1. **Metabolite selection**

   * Uses metabolites with `matched == 1`
   * Further restricts to metabolites that exist as columns in the external dataset

2. **Standardization**

   * For each metabolite, mean and SD are computed **using controls only**
   * Both cases and controls are standardized using these control-based parameters

3. **Score calculation**
   
MMS is computed as a weighted sum of standardized metabolites:

`MMS_i = sum_j ( beta_j * Z_ij )`

where `Z_ij` is the standardized value of metabolite `j` for participant `i`.

4. **Missing data**

   * If a participant is missing any metabolite used in the MMS, their MMS is set to `NA`

5. Because metabolite availability may differ across platforms, the MMS in each external dataset may be based on a subset of the original IDATA metabolites. The same IDATA-derived coefficients are applied to the available metabolites only.
---

## Example Usage

```r
library(readxl)

# Load external dataset
# ext_df <- read.csv("external_metabolites.csv")

# Load MMS coefficients from IDATA
IDATA_df <- read_excel("MMS_sodium_serum.xlsx")

# Compute MMS
ext_df$MMS_sodium_serum <-
  compute_MMS_external(
    ext_df   = ext_df,
    IDATA_df = IDATA_df,
    case_col = "case"
  )
```

Repeat the above steps for each MMS coefficient file.

---

## Downstream Analysis Example

The MMS can be used as an exposure in conditional logistic regression:

```r
library(survival)

clogit(
  case ~ MMS_sodium_serum + covariates + strata(match_id),
  data = ext_df
)
```

---

## Notes and Assumptions


- MMS are computed using IDATA-derived coefficients and the subset of metabolites available in the external dataset.
- As metabolite overlap may differ across studies, MMS values are **study-specific** and should not be interpreted as numerically equivalent across cohorts.
- MMS are intended to be used as **relative exposure variables** in within-study association analyses (e.g., conditional logistic regression).
- Intercepts from IDATA models are intentionally excluded.

---

## Contact

For questions or clarifications, please contact the study team.


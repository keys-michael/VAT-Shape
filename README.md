# VAT-Shape

Statistical analysis code for:

**Keys M, Cockburn R, Chung E, Gately R, Preston J, Lawson M, Wood S, 
Dissabandara L, Gopalan V, Rhee H.** Visceral adipose tissue delays kidney 
function recovery after live donor nephrectomy. [in press].

---

# Overview
This repository contains the R analysis script for a retrospective cohort study examining the association of pre-operative visceral adipose tissue (VAT) and remnant kidney cortex proportion with kidney function recovery after live donor nephrectomy.
Kidney function recovery was modelled as the ratio of post-operative eGFR to donor-specific mean pre-operative eGFR (CKD-EPI 2021), using generalised additive mixed models with tensor product interactions via the mgcv package.

**Study period:** 2020–2023, single tertiary centre, Brisbane, Australia  
**Follow-up:** Up to December 2025 (maximum 5.7 years)  
**Analytic sample:** 79 donors, 519 pre-operative creatinine observations, 521 post-operative creatinine observations

---

## Script outputs

The script reproduces all results reported in the manuscript:

- Descriptive statistics (Table 1) and missingness assessment
- Primary model and sensitivity analyses 1–6 (pre-specified)
- Sensitivity analysis 7: sex-specific recovery trajectories (post hoc, added during peer review)
- Model diagnostics (Figure S1), including a heavy-tailed (scaled t) robustness check
- Observed versus predicted calibration (Table S1)
- Figures 2–4: predicted trajectories by VAT and by remnant kidney cortex proportion,
  and the high- versus low-VAT difference curve
- Table 2 and recovery slopes by VAT level

A fixed random seed (`SEED`) is set so that `gam.check()` permutation results are reproducible.

---

## Requirements

- R (developed with version 4.x)
- Required packages:
  - dplyr
  - tidyr
  - readxl
  - stringr
  - ggplot2
  - mgcv
  - misty
  - moments
  - patchwork
  - gridExtra

Install all packages with:
```r
install.packages(c("dplyr", "tidyr", "readxl", "stringr", "ggplot2", 
                   "mgcv", "misty", "moments", "patchwork", "gridExtra"))
```

---

## Required Setup

**REQUIRED**: Before running the script, set 'GRAVES_DONOR_DOD' in the Static Variables section to the date of donation of the donor who developed Graves' disease during follow-up. This value has been replaced with a placeholder ('YYYY-MM-DD') prior to publication to prevent identification of the donor.
```r
GRAVES_DONOR_DOD <- as.Date("YYYY-MM-DD")  # Replace with actual date
```

## Data Availability

Patient data are not publicly available due to confidentiality requirements. 
Data may be shared upon reasonable request directed to the corresponding author.

The script expects a single Excel workbook (`Renal Function.xlsx`) with the 
following sheets:
- `Post Cr` — post-operative serum creatinine collections
- `Pre Cr` — pre-operative serum creatinine collections
- `L3` — CT-derived body composition variables including VAT
- `Data` — demographic, clinical, and kidney volumetry variables

---

## Author

**Michael Keys**  
Princess Alexandra Hospital, Brisbane, Australia  
ORCID: 0009-0002-5881-2637  
Email: michael.keys@health.qld.gov.au  

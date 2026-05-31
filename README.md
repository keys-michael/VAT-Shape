# VAT-Shape

Statistical analysis code for:

**Keys M, Cockburn R, Chung E, Gately R, Preston J, Lawson M, Wood S, 
Dissabandara L, Gopalan V, Rhee H.** Visceral adipose tissue delays kidney 
function recovery after live donor nephrectomy. [in press].

---

# Overview
This repository contains the R analysis script for a retrospective cohort study examining the influence of pre-operative visceral adipose tissue (VAT) on kidney function recovery after live donor nephrectomy.
Kidney function recovery was modelled as a ratio of post-operative eGFR to pre-operative mean eGFR using generalised additive mixed models with tensor product interactions via the mgcv package.

**Study period:** 2020–2023, single tertiary centre, Brisbane, Australia  
**Follow-up:** Up to December 2025 (maximum 5.7 years)  
**Analytic sample:** 79 donors, 519 pre-operative creatinine observations, 521 post-operative creatinine observations

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

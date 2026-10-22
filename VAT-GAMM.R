library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(ggplot2)
library(mgcv)
library(misty)
library(moments)
library(patchwork)
library(gridExtra)

# ---------- Static Variables --------------------------------------------------
DAYS_MIN_CUTOFF = 28

# Timepoints for prediction table (months 1-36)
TIMEPOINTS <- c(28, 183, 365, 548, 731, 913, 1096) # 1,6,12,18 months where 1 month = 365.25/12 and then rounded

# Extended timepoints including months 48 and 60
TIMEPOINTS_FULL <- c(TIMEPOINTS, 1461, 1826)
TIMEPOINTS_FULL_LABELS <- c("1", "6", "12", "18", "24", "30", "36", "48", "60")

# Time bands for observation density and observed vs predicted analysis
TIME_BAND_BREAKS <- c(DAYS_MIN_CUTOFF - 1, 182, 365, 548, 730, 913, 1096, Inf)
TIME_BAND_LABELS <- c("1-6m", "6-12m", "12-18m", "18-24m", "24-30m", "30-36m", "36m+")

# Donor with Graves' disease diagnosis during follow-up
# Post-diagnosis observations excluded as hyperthyroidism independently
# elevates eGFR via haemodynamic and creatinine production mechanisms
GRAVES_DONOR_DOD <- as.Date("YYYY-MM-DD") # Replace with actual date before running

# Seed for gam check random permutation tests
SEED <- 2026


# ---------- Helper functions --------------------------------------------------
# CKD-EPI 2021 eGFR from serum creatinine (µmol/L), sex (M/F), and age (years)
calc_egfr <- function(cr_umol, sex, age) {
  cr_mgdl <- cr_umol * 0.0113
  k <- ifelse(sex == "F", 0.7, 0.9)
  a <- ifelse(sex == "F", -0.241, -0.302)
  sex_coefficient <- ifelse(sex == "F", 1.012, 1.0)
  142 * (pmin(cr_mgdl / k, 1) ^ a) *
    (pmax(cr_mgdl / k, 1) ^ -1.2) *
    (0.9938 ^ age) * sex_coefficient
}

# Kidney function recovery: post-operative eGFR as a proportion of donor-specific pre-operative eGFR
ratio_egfr <- function(pre_egfr, post_egfr) post_egfr / pre_egfr

# gam.check() with a fixed seed, so k-index permuation p-values are reproducible per call
check_gam <- function(model, seed = SEED) {
  set.seed(seed)
  gam.check(model)
}

# One table 2 row: VAT x time interaction (edf, F, p), and remKCV beta (95% CI, p) from a fitted model
extract_row <- function(m, model_name) {
  s <- summary(m)
  
  # remKCV: parametric coefficient, SE, 95% CI (Wald, normal approximation)
  b  <- s$p.table["remKCV_s", "Estimate"]
  se <- s$p.table["remKCV_s", "Std. Error"]
  z  <- qnorm(0.975)
  
  # VAT x time interaction: the ti() smooth term (absent in SA6)
  ti_idx <- grep("^ti\\(", rownames(s$s.table))
  has_ti <- length(ti_idx) == 1
  
  tibble(
    Model     = model_name,
    edf       = if (has_ti) s$s.table[ti_idx, "edf"]     else NA_real_,
    F_stat    = if (has_ti) s$s.table[ti_idx, "F"]       else NA_real_,
    p_ti      = if (has_ti) s$s.table[ti_idx, "p-value"] else NA_real_,
    remKCV_b  = b,
    remKCV_lo = b - z * se,
    remKCV_hi = b + z * se,
    p_remKCV  = s$p.table["remKCV_s", "Pr(>|t|)"]
  )
}

# Mean oberved vs predicted recovery ratio by VAT group and time band (donor weighted by default)
# by_donor = TRUE averages within each donor first, so every donor contributes
# equally regardless of how many creatinines they had (matches the model)
calibration_table <- function(data, pred_col, by_donor = TRUE) {
  d <- data %>% mutate(pred = .data[[pred_col]])

  if (by_donor) {
    d <- d %>%
      group_by(VAT_group, time_band, ID) %>%
      summarise(egfr_ratio = mean(egfr_ratio), pred = mean(pred), .groups = "drop")
  }
  d %>%
    group_by(VAT_group, time_band) %>%
    summarise(
      n_donors  = n_distinct(ID),
      observed  = mean(egfr_ratio),
      predicted = mean(pred),
      .groups = "drop"
    ) %>%
    mutate(difference = observed - predicted)
}

# Recovery slope (change in recovery ratio per year) design vector at a given VAT level.
# Trajectories are linear (edf ~ 1), so the slope is the same over any interval.
slope_X <- function(model, vat, t0 = 28, t1 = 393) {
  nd <- data.frame(Days_Post = c(t0, t1), VAT_s = vat, remKCV_s = 0, ID = model$model$ID[1])
  X  <- predict(model, newdata = nd, type = "lpmatrix", exclude = "s(ID)")
  (X[2, ] - X[1, ]) / ((t1 - t0) / 365.25)
}

# Estimate, 95% CI and p-value for a slope (or difference in slopes) design vector
slope_ci <- function(model, Xs) {
  est <- sum(Xs * coef(model))
  se  <- sqrt(drop(t(Xs) %*% vcov(model) %*% Xs))
  round(c(slope = est, lower = est - 1.96 * se, upper = est + 1.96 * se,
          p = 2 * pnorm(-abs(est / se))), 4)
}

# ---------- Data Loading ------------------------------------------------------

# Post-operative creatinine data — base dataframe, all others joined to this
df <- read_excel("Renal Function.xlsx", sheet = "Post Cr") %>%
  mutate(ID = row_number()) %>%
  relocate(ID)

# Pre-operative creatinine data - used to calculate donor-specific mean baseline eGFR
pre_op_cr <- read_excel("Renal Function.xlsx", sheet = "Pre Cr") %>%
  mutate(ID = row_number()) %>%
  relocate(ID)

pre_op_cr_tidy <- pre_op_cr %>%
  rename_with(~ str_replace(., "^Pre-op Cr #(\\d+)\\s*$", "Cr #\\1")) %>%
  rename_with(~ str_replace(., "^Days Pre #(\\d+)\\s*$", "Days_Pre #\\1")) %>%
  rename_with(~ str_replace(., "^Pre-op Cr #(\\d+) Location\\s*$", "Location #\\1")) %>%
  rename_with(~ str_replace(., "^Pre-op Cr #(\\d+) Fasting\\?\\s*$", "Fasting #\\1")) %>%
  pivot_longer(
    cols = matches("^(Cr|Days_Pre|Location|Fasting) #\\d+$"),
    names_to = c(".value", "Collect_Count"),
    names_sep = " #"
  ) %>%
  filter(!(is.na(Cr) | Cr == 0)) %>%
  select(ID, Collect_Count, Days_Pre, Cr, Gender, Age) %>%
  mutate(
    currentAge = Age - (Days_Pre/365),
    eGFR = calc_egfr(Cr, Gender, currentAge)
  )

pre_op_cr_averages <- pre_op_cr_tidy %>%
  group_by(ID) %>%
  summarise(mean_pre_egfr = mean(eGFR, na.rm = TRUE))

# Visceral adipose tissue data, from the L3 slice
VAT_export <- read_excel("Renal Function.xlsx", sheet = "L3") %>%
  mutate(ID = row_number()) %>%
  select(ID, `VAT`)

# Load variables from the data sheet
data_sheet <- read_excel("Renal Function.xlsx", sheet = "Data") %>%
  mutate(ID = row_number()) %>%
  select(ID, `Pre-op KCV SRF`, `BMI`, `Donated Kidney Laterality`,
         `Height`, `HTN`, `Dyslipidaemia`, `Date of Donation`,
         `Diabetes`, Hyperuricaemia,
         `DL (R) Kidney Cortex Volume`, `DL (R) Kidney Medulla Volume`,
         `DL (L) Kidney Cortex Volume`, `DL (L) Kidney Medulla Volume`)

# Kidney volumetrics from deep learning segmentation (Korfiatis model)
# Remnant (preserved) kidney is the side NOT donated
# PKV/TKV include cortex + medulla (parenchymal volume); PKCV/TKCV cortex only
volumetrics <- data_sheet %>%
  transmute(
    ID,
    laterality = toupper(trimws(`Donated Kidney Laterality`)),
    R_cortex   = as.numeric(`DL (R) Kidney Cortex Volume`),
    R_medulla  = as.numeric(`DL (R) Kidney Medulla Volume`),
    L_cortex   = as.numeric(`DL (L) Kidney Cortex Volume`),
    L_medulla  = as.numeric(`DL (L) Kidney Medulla Volume`),
    TKCV = R_cortex + L_cortex,
    TKV  = R_cortex + R_medulla + L_cortex + L_medulla,
    PKCV = case_when(laterality == "LEFT"  ~ R_cortex,
                     laterality == "RIGHT" ~ L_cortex),
    PKV  = case_when(laterality == "LEFT"  ~ R_cortex + R_medulla,
                     laterality == "RIGHT" ~ L_cortex + L_medulla),
    remKCV_check = PKCV / TKCV
  )

# Sanity check: R-derived PKCV/TKCV should match Excel's Pre-op KCV SRF
kcv_mismatch <- volumetrics %>%
  left_join(data_sheet %>% select(ID, `Pre-op KCV SRF`), by = "ID") %>%
  mutate(excel_kcv = as.numeric(`Pre-op KCV SRF`)) %>%
  filter(!is.na(remKCV_check), !is.na(excel_kcv),
         abs(remKCV_check - excel_kcv) > 1e-6)

if (nrow(kcv_mismatch) == 0) {
  cat("Sanity check passed: R-derived PKCV/TKCV matches Excel remKCV.\n")
} else {
  print(kcv_mismatch)
  stop("STOP: ", nrow(kcv_mismatch), " remKCV mismatches found!")
}

volumetrics <- volumetrics %>% select(ID, PKV, TKV, PKCV, TKCV)

rem_KCV_ratio   <- data_sheet %>% select(ID, `Pre-op KCV SRF`)
BMI             <- data_sheet %>% select(ID, `BMI`)
Laterality      <- data_sheet %>% select(ID, `Donated Kidney Laterality`)
Height          <- data_sheet %>% select(ID, `Height`)
Hypertension    <- data_sheet %>% select(ID, `HTN`)
Dyslipidaemia   <- data_sheet %>% select(ID, `Dyslipidaemia`)


# Load Date of Donation from both sheets as a sanity check
# Date of Donation is unique per donor across the entire dataset
# and should match exactly after row-number-based joining
donation_date_data <- data_sheet %>% select(ID, `Date of Donation`)

donation_date_postCr <- read_excel("Renal Function.xlsx", sheet = "Post Cr") %>%
  mutate(ID = row_number()) %>%
  select(ID, `Date of Donation`)

donation_date_preCr <- read_excel("Renal Function.xlsx", sheet = "Pre Cr") %>%
  mutate(ID = row_number()) %>%
  select(ID, `Date of Donation`)

donation_date_L3 <- read_excel("Renal Function.xlsx", sheet = "L3") %>%
  mutate(ID = row_number()) %>%
  select(ID, `Date of Donation`)

date_check_all <- donation_date_data %>%
  rename(date_Data = `Date of Donation`) %>%
  left_join(donation_date_postCr %>% rename(date_PostCr = `Date of Donation`), by = "ID") %>%
  left_join(donation_date_preCr %>% rename(date_PreCr = `Date of Donation`), by = "ID") %>%
  left_join(donation_date_L3 %>% rename(date_L3 = `Date of Donation`), by = "ID")

date_mismatches <- date_check_all %>%
  filter(
    date_Data != date_PostCr |
    date_Data != date_PreCr |
    date_Data != date_L3
  )

if (nrow(date_mismatches) == 0) {
  cat("Sanity check passed: Date of donation matches across all four sheets.\n")
} else {
  stop("STOP: ", nrow(date_mismatches), " mismatches found!\n")
  print(date_mismatches)
}

# Join the above 
df <- df %>%
  left_join(VAT_export, by = "ID") %>%
  left_join(rem_KCV_ratio, by = "ID") %>%
  left_join(pre_op_cr_averages, by = "ID") %>%
  left_join(BMI, by = "ID") %>%
  left_join(Laterality, by = "ID") %>%
  left_join(Height, by = "ID") %>%
  left_join(Hypertension, by = "ID") %>%
  left_join(Dyslipidaemia, by = "ID")

df_renamed <- df %>%
  rename_with(~ str_replace(., "^Post-op Cr #(\\d+)\\s*$", "Cr #\\1")) %>%
  rename_with(~ str_replace(., "^Days Post #(\\d+)\\s*$", "Days_Post #\\1")) %>%
  rename_with(~ str_replace(., "^Post-op Cr #(\\d+) Location\\s*$", "Location #\\1")) %>%
  rename_with(~ str_replace(., "^Post-op Cr #(\\d+) Fasting\\?\\s*$", "Fasting #\\1")) %>%
  rename_with(~ str_replace(., "VAT", "VAT")) %>%
  rename_with(~ str_replace(., "Donated Kidney Laterality", "Laterality")) %>%
  mutate(VAT = as.numeric(VAT)) %>%
  mutate(BMI = as.numeric(BMI)) %>%
  mutate(Height = as.numeric(Height)) %>%
  rename_with(~ str_replace(., "Pre-op KCV SRF", "remKCV_ratio")) %>%
  mutate(remKCV_ratio = as.numeric(remKCV_ratio))

mcar_data <- df_renamed %>%
  select(VAT, BMI, remKCV_ratio, mean_pre_egfr)
na.test(mcar_data)

df_tidy <- df_renamed %>%
  pivot_longer(
    cols = matches("^(Cr|Days_Post|Location|Fasting) #\\d+"),
    names_to = c(".value", "Collect_Count"),
    names_sep = " #"
  ) %>%
  filter(!(is.na(Cr) | Cr == 0 |
           is.na(VAT) | VAT == 0 |
           is.na(remKCV_ratio) | remKCV_ratio == 0 |
           is.na(BMI) | BMI == 0 |
           is.na(Height) | Height == 0 | 
           Days_Post < DAYS_MIN_CUTOFF))

unique_ids <- unique(df_tidy$ID)
df_renamed_backup <- df_renamed[df_renamed$ID %in% unique_ids, ]

# Used for descriptive analysis later
comorbidities <- data_sheet %>%
  select(ID, `Diabetes`, `HTN`, `Dyslipidaemia`, `Hyperuricaemia`) %>%
  filter(ID %in% unique_ids)

# ---------- Analysis Setup ----------------------------------------------------

# Analysis settings
day_cap <- max(df_tidy$Days_Post, na.rm = TRUE)

donor_ref <- df_renamed_backup %>%
  summarise(
    VAT_mean  = mean(VAT),                  VAT_sd  = sd(VAT),
    KCV_mean  = mean(remKCV_ratio),         KCV_sd  = sd(remKCV_ratio),
    BMI_mean  = mean(BMI),                  BMI_sd  = sd(BMI),
    VATH_mean = mean(VAT / Height^2),       VATH_sd = sd(VAT / Height^2),
    Age_mean  = mean(Age),                  Age_sd  = sd(Age)
  )

# Create analysis variables - egfr_ratio, scaled covariates, ordered gender factor
df_tidy <- df_tidy %>%
  select(ID, Collect_Count, Days_Post, Cr, Gender, Age, mean_pre_egfr, VAT, 
         remKCV_ratio, BMI, Laterality, Height, Dyslipidaemia, HTN) %>%
  mutate(
    currentAge = Age + (Days_Post/365),
    post_egfr = calc_egfr(Cr, Gender, currentAge),     # Calculate post-op eGFR with the current age
    egfr_ratio = ratio_egfr(mean_pre_egfr, post_egfr),
    Gender_o = ordered(Gender, levels = c("F", "M"))  # Required for SA7
  ) %>%
  filter(Days_Post <= day_cap) %>%
  mutate(
    VAT_s       = (VAT - donor_ref$VAT_mean) / donor_ref$VAT_sd,
    remKCV_s    = (remKCV_ratio - donor_ref$KCV_mean) / donor_ref$KCV_sd,
    BMI_s       = (BMI - donor_ref$BMI_mean) / donor_ref$BMI_sd,
    VATHeight_s = (VAT / Height^2 - donor_ref$VATH_mean) / donor_ref$VATH_sd,
    Age_s       = (Age - donor_ref$Age_mean) / donor_ref$Age_sd
  ) %>%
  select(ID, Collect_Count, Days_Post, egfr_ratio, post_egfr, Gender_o, VAT_s, 
         remKCV_s, BMI_s, Laterality, VATHeight_s, Dyslipidaemia, HTN, Gender, Age_s)

VAT_mean <- donor_ref$VAT_mean
VAT_sd   <- donor_ref$VAT_sd

df_tidy <- df_tidy %>% mutate(ID = as.factor(ID))

# Identify the Graves' donor ID dynamically from their date of donation
graves_id <- df_renamed_backup %>%
  filter(as.Date(`Date of Donation`) == GRAVES_DONOR_DOD) %>%
  pull(ID)

# Confirm exactly one donor was matched as the Graves donor
stopifnot(length(graves_id) == 1)

# Exclude post-Graves'-diagnosis observations for the donor that developed graves
# One participant developed Graves' disease during follow-up; hyperthyroidism
# independently elevates eGFR via haemodynamic and creatinine production
# mechanisms unrelated to nephrectomy recovery. Diagnosis after day 694
# observations excluded for all results after this day
df_tidy_full <- df_tidy  # save before applying graves ID exclusion
df_tidy <- df_tidy %>%
  filter(!(ID == graves_id & Days_Post >= 694))


# ---------- Descriptive Statistics --------------------------------------------

# Descriptive statistics for analytic sample
descriptive_data <- df_renamed_backup %>%
  select(ID, Gender, Age, Laterality, VAT, BMI, remKCV_ratio, mean_pre_egfr) %>%
  left_join(volumetrics, by = "ID") %>%
  mutate(
    VAT = as.numeric(VAT),
    BMI = as.numeric(BMI),
    remKCV_ratio = as.numeric(remKCV_ratio)
  )

# N, sex split, donated laterality and comorbidities
cat("N donors:", nrow(descriptive_data), "\n")
cat("Male n (%):", sum(descriptive_data$Gender == "M"),
    "(", round(mean(descriptive_data$Gender == "M") * 100, 1), "%)\n")
cat("Left donation (%):", sum(descriptive_data$Laterality == "LEFT"),
    "(", round(mean(descriptive_data$Laterality == "LEFT") * 100, 1), "%)\n")
cat("Diabetes (%):", sum(comorbidities$Diabetes),
    "(", round(mean(comorbidities$Diabetes) * 100,1), "%)\n")
cat("Hypertension (%):", sum(comorbidities$HTN),
    "(", round(mean(comorbidities$HTN) * 100,1), "%)\n")
cat("Dyslipidaemia (%):", sum(comorbidities$Dyslipidaemia),
    "(", round(mean(comorbidities$Dyslipidaemia) * 100,1), "%)\n")
cat("Hyperuricaemia (%):", sum(comorbidities$Hyperuricaemia),
    "(", round(mean(comorbidities$Hyperuricaemia) * 100,1), "%)\n")


# Continuous variables - mean and SD
descriptive_data %>%
  select(Age, VAT, BMI, remKCV_ratio, mean_pre_egfr, PKV, TKV, PKCV, TKCV) %>%
  summarise(across(everything(), list(
    mean     = ~ round(mean(.x, na.rm = TRUE), 4),
    sd       = ~ round(sd(.x, na.rm = TRUE), 4),
    median   = ~ round(median(.x, na.rm = TRUE), 4),
    iqr      = ~ round(IQR(.x, na.rm = TRUE), 4),
    skewness = ~ round(moments::skewness(.x, na.rm = TRUE), 4),
    kurtosis = ~ round(moments::kurtosis(.x, na.rm = TRUE), 4)
  ))) %>%
  pivot_longer(everything(),
               names_to = c("variable", "stat"),
               names_sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  print()

# Pre-operative creatinine observations used to calculate the mean pre-op eGFR
pre_op_counts <- pre_op_cr_tidy %>%
  filter(ID %in% unique_ids) %>% #restrict to analytic sample only
  group_by(ID) %>%
  summarise(n_pre_cr = n())

cat("Pre-op creatinines per donor - median (IQR):",
    median(pre_op_counts$n_pre_cr), "(",
    quantile(pre_op_counts$n_pre_cr, 0.25), "-",
    quantile(pre_op_counts$n_pre_cr, 0.75), ")\n")

cat("Total pre-op creatinine measurements:", sum(pre_op_counts$n_pre_cr), "\n")

# Distribution breakdown of pre-op creatinines
table(pre_op_counts$n_pre_cr)

pre_op_days <- pre_op_cr_tidy %>%
  filter(ID %in% unique_ids) %>% # restrict to analytic sample only
  group_by(ID) %>%
  summarise(days_span = max(Days_Pre) - min(Days_Pre))

cat("Span of pre-op creatinine collection per donor - median (IQR) in years:",
    round(median(pre_op_days$days_span) / 365.25, 1), "(",
    round(quantile(pre_op_days$days_span, 0.25) / 365.25, 1), "-",
    round(quantile(pre_op_days$days_span, 0.75)/365.25, 1), ")\n")

# Vat +-1 SD in raw units (for interpreting trajectory plots)
cat("\nVAT mean:", round(mean(descriptive_data$VAT, na.rm=TRUE), 1), "cm²\n")
cat("VAT SD:", round(sd(descriptive_data$VAT, na.rm=TRUE), 1), "cm²\n")
cat("VAT -1SD:", round(mean(descriptive_data$VAT, na.rm=TRUE) - sd(descriptive_data$VAT, na.rm=TRUE), 1), "cm²\n")
cat("VAT +1SD:", round(mean(descriptive_data$VAT, na.rm=TRUE) + sd(descriptive_data$VAT, na.rm=TRUE), 1), "cm²\n")


post_op_counts <- as.numeric(table(df_tidy$ID))
cat("Total post-operative observations (primary model):", nrow(df_tidy), "\n")
cat("Observations per donor - median (IQR):",
    median(post_op_counts), "(",
    quantile(post_op_counts, 0.25), "-",
    quantile(post_op_counts, 0.75), ")\n")

post_op_days <- df_tidy %>%
  group_by(ID) %>%
  summarise(max_days = max(Days_Post))

cat("Follow-up duration per donor - median (IQR):",
    round(median(post_op_days$max_days) / 365.25, 1), "years (",
    round(quantile(post_op_days$max_days, 0.25) / 365.25, 1), "-",
    round(quantile(post_op_days$max_days, 0.75) / 365.25, 1), ")\n")
cat("Maximum follow-up:",
    round(max(post_op_days$max_days) / 365.25, 1), "years\n")

# High and Low-VAT subgroup descriptive reporting
df_renamed_backup %>%
  mutate(VAT_group = case_when(
    VAT > VAT_mean + VAT_sd ~ "High VAT",
    VAT < VAT_mean - VAT_sd ~ "Low VAT",
    TRUE ~ "Middle"
  )) %>%
  filter(VAT_group != "Middle") %>%
  group_by(VAT_group) %>%
  summarise(
    n = n(),
    n_male = sum(Gender == "M"),
    pct_male = round(100 * mean(Gender == "M"),1),
    n_HTN = sum(HTN == 1, na.rm = TRUE),
    n_Dys = sum(Dyslipidaemia == 1, na.rm = TRUE),
    mean_age = mean(Age, na.rm = TRUE),
    sd_age = sd(Age, na.rm = TRUE),
    mean_egfr = mean(mean_pre_egfr, na.rm = TRUE),
    sd_pre_egfr = sd(mean_pre_egfr, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))

# Count observations by time band for each VAT group
df_tidy %>%
  mutate(
    VAT_group = case_when(
      VAT_s > 0.5 ~ "High VAT",
      VAT_s < -0.5 ~ "Low VAT",
      TRUE ~ "Middle"
    ),
    time_band = cut(Days_Post,
                    breaks = TIME_BAND_BREAKS,
                    labels = TIME_BAND_LABELS,
                    right = TRUE)
  ) %>%
  filter(VAT_group != "Middle") %>%
  group_by(VAT_group, time_band) %>%
  summarise(n_obs = n(), n_donors = n_distinct(ID), .groups = "drop") %>%
  pivot_wider(names_from = VAT_group, 
              values_from = c(n_obs, n_donors))


# ---------- Determine Missingness Reasons -------------------------------------
# Full 93-donor set before filtering
all_ids <- df_renamed$ID

# Identify excluded donors
excluded_ids <- setdiff(all_ids, unique_ids)

cat("Total donors screened:", length(all_ids), "\n")
cat("Analytic sample:", length(unique_ids), "\n")
cat("Excluded:", length(excluded_ids), "\n\n")

# For each excluded donor, identify which variable(s) caused exclusion
df_renamed_wide <- df_renamed %>%
  filter(ID %in% excluded_ids)

# Check each exclusion criterion
cat("Excluded donors breakdown:\n")
df_renamed_wide %>%
  rowwise() %>%
  mutate(
    missing_VAT = is.na(VAT) | VAT == 0,
    missing_remKCV = is.na(remKCV_ratio) | remKCV_ratio == 0,
    missing_BMI = is.na(BMI) | BMI == 0,
    missing_Height = is.na(Height) | Height == 0,
    missing_post_cr = all(is.na(c_across(matches("^Cr #\\d+$"))) | 
                            c_across(matches("^Cr #\\d+$")) == 0),
    only_early_cr = {
      days <- c_across(matches("^Days_Post #\\d+$"))
      crs  <- c_across(matches("^Cr #\\d+$"))
      valid_cr <- !is.na(crs) & crs != 0
      has_any_cr <- any(valid_cr)
      all_within <- all(days[valid_cr] < DAYS_MIN_CUTOFF, na.rm = TRUE)
      has_any_cr & all_within
    }
  ) %>%
  ungroup() %>%
  select(ID, missing_VAT, missing_remKCV, missing_BMI, 
         missing_Height, missing_post_cr, only_early_cr) %>%
  print()


# ---------- Fit Models --------------------------------------------------------

# Fit primary GAMM model:
# - remKCV_s as a fixed linear population-level effect
# - s(Days_Post) and s(VAT_s) are fixed non-linear population-level effects
# - ti(Days_Post, VAT_s) are fixed population-level smooths that show how VAT modifies recovery curve over time
# - s(ID) is a non-population level (i.e., vary per donor), random effect. Accounts for repeated measurements from the same person
m_gam <- gam(
  egfr_ratio ~ remKCV_s +
    s(Days_Post, k = 10) +
    s(VAT_s, k = 5) +
    ti(Days_Post, VAT_s, k = c(10, 5)) +
    s(ID, bs = "re"),
  data = df_tidy,
  method = "REML"
)
check_gam(m_gam)
summary(m_gam)
concurvity(m_gam, full = FALSE)


# Sensitivity Analysis 1: BMI as a smooth covariate
m_gam_BMI_sens <- gam(
  egfr_ratio ~ remKCV_s +
    s(Days_Post, k = 10) +
    s(VAT_s, k = 5) +
    s(BMI_s, k = 5) +
    ti(Days_Post, VAT_s, k = c(10, 5)) +
    s(ID, bs = "re"),
  data = df_tidy,
  method = "REML"
)
check_gam(m_gam_BMI_sens)
summary(m_gam_BMI_sens)
concurvity(m_gam_BMI_sens, full = FALSE)

# Sensitivity Analysis 2: Replace VAT with (VAT / Height^2)
m_gam_height_sens <- gam(
  egfr_ratio ~ remKCV_s +
    s(Days_Post, k = 10) +
    s(VATHeight_s, k = 5) +
    ti(Days_Post, VATHeight_s, k = c(10, 5)) +
    s(ID, bs = "re"),
  data = df_tidy,
  method = "REML"
)
check_gam(m_gam_height_sens)
summary(m_gam_height_sens)
concurvity(m_gam_height_sens, full = FALSE)

# Sensitivity Analysis 3: All post-op Cr's (including post-Graves' diagnosis)
m_gam_graves_sens <- gam(
  egfr_ratio ~ remKCV_s +
    s(Days_Post, k = 10) +
    s(VAT_s, k = 5) +
    ti(Days_Post, VAT_s, k = c(10, 5)) +
    s(ID, bs = "re"),
  data = df_tidy_full,
  method = "REML"
)
check_gam(m_gam_graves_sens)
summary(m_gam_graves_sens)
concurvity(m_gam_graves_sens, full = FALSE)

# Sensitivity Analysis 4: Comorbidities (with >5% prevalence) as fixed linear covariates
m_gam_comorbidity_sens <- gam(
  egfr_ratio ~ remKCV_s + HTN + Dyslipidaemia +
    s(Days_Post, k = 10) +
    s(VAT_s, k = 5) +
    ti(Days_Post, VAT_s, k = c(10, 5)) +
    s(ID, bs = "re"),
  data = df_tidy,
  method = "REML"
)
check_gam(m_gam_comorbidity_sens)
summary(m_gam_comorbidity_sens)
concurvity(m_gam_comorbidity_sens, full = FALSE)

# Sensitivity analysis 5: Add Gender and Age as fixed linear covariates
m_gam_demographic_fixed_sens <- gam(
  egfr_ratio ~ remKCV_s + Gender + Age_s +
    s(Days_Post, k = 10) +
    s(VAT_s, k = 5) +
    ti(Days_Post, VAT_s, k = c(10, 5)) +
    s(ID, bs = "re"),
  data = df_tidy,
  method = "REML"
)
check_gam(m_gam_demographic_fixed_sens)
summary(m_gam_demographic_fixed_sens)
concurvity(m_gam_demographic_fixed_sens, full = FALSE)

# Sensitivity analysis 6: VAT, Gender, Age and remKCV as fixed linear covariates
m_gam_fixed <- gam(
  egfr_ratio ~ VAT_s + Age_s + Gender + remKCV_s +
    s(Days_Post, k = 10) +
    s(ID, bs = "re"),
  data = df_tidy,
  method = "REML",
)
check_gam(m_gam_fixed)
summary(m_gam_fixed)
concurvity(m_gam_fixed, full = FALSE)

# Sensitivity analysis 7: sex-specific recovery trajectories
# POST-HOC, added during peer review. High-VAT donors were predominantly male.
# This tests whether the VAT x time interaction persists when men and women are
# permitted separate recovery curves (not just a level shift, as in SA5)
m_gam_sex_time_sens <- gam(
  egfr_ratio ~ remKCV_s + Gender_o +
    s(Days_Post, k = 10) +
    s(Days_Post, by = Gender_o, k = 10) +   # male-minus-female difference in trajectory
    s(VAT_s, k = 5) +
    ti(Days_Post, VAT_s, k = c(10, 5)) +
    s(ID, bs = "re"),
  data = df_tidy,
  method = "REML"
)
check_gam(m_gam_sex_time_sens)
summary(m_gam_sex_time_sens)
concurvity(m_gam_sex_time_sens, full = FALSE)


# ---------- Model Diagnostics (including supplementary figure) ----------------
# Supplementary Figure: residual diagnostics for the primary model
# Panels drawn individually (rather than via gam.check) to allow full titles
png("Supplementary_Figure_diagnostics.png",
    width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2))

dev_resid <- residuals(m_gam, type = "deviance")

# A: Normal Q-Q plot (same method gam.check uses, with red reference line)
qq.gam(m_gam, rep = 0, rl.col = "red",
       main = "A. Normal Q-Q plot of residuals",
       pch = 20, cex = 0.6)

# B: Residuals vs linear predictor
plot(m_gam$linear.predictors, dev_resid,
     main = "B. Residuals vs linear predictor",
     xlab = "Linear predictor", ylab = "Deviance residuals")
abline(h = 0, lty = 2)

# C: Histogram of residuals
hist(dev_resid,
     main = "C. Histogram of residuals",
     xlab = "Deviance residuals", ylab = "Frequency")

# D: Observed vs fitted
plot(fitted(m_gam), m_gam$y,
     main = "D. Observed vs fitted values",
     xlab = "Fitted values", ylab = "Observed kidney function recovery ratio")
abline(0, 1, lty = 2)

dev.off()


# Primary model  robustness check: heavy-tailed scaled-t family
m_gam_scat <- gam(
  egfr_ratio ~ remKCV_s +
    s(Days_Post, k = 10) +
    s(VAT_s, k = 5) +
    ti(Days_Post, VAT_s, k = c(10, 5)) +
    s(ID, bs = "re"),
  data = df_tidy,
  family = scat(),
  method = "REML"
)
summary(m_gam_scat)

# ---------- Observed vs Predicted: Supplementary Table S1 --------------------
# Compares model predictions with observed data by VAT subgroup (same +/-1 SD
# cut-offs used for the subgroup baseline eGFRs). Predictions are population-
# level (random effect excluded), made at each observation's actual covariate
# values - a like-for-like comparison. High and low VAT are tabulated in full;
# intermediate VAT is summarised in the Table S1 footnote.

# Donor-level VAT subgroup assignment
vat_groups <- df_renamed_backup %>%
  transmute(
    ID_chr = as.character(ID),
    VAT_group = factor(
      case_when(
        VAT > VAT_mean + VAT_sd ~ "High VAT",
        VAT < VAT_mean - VAT_sd ~ "Low VAT",
        TRUE                    ~ "Intermediate VAT"
      ),
      levels = c("High VAT", "Low VAT", "Intermediate VAT")
    )
  )

# Observation-level data with primary-model predictions, all donors
obs_pred <- df_tidy %>%
  mutate(
    pred_primary = as.numeric(predict(m_gam, newdata = df_tidy, exclude = "s(ID)")),
    ID_chr       = as.character(ID),
    time_band    = cut(Days_Post, breaks = TIME_BAND_BREAKS, labels = TIME_BAND_LABELS)
  ) %>%
  left_join(vat_groups, by = "ID_chr")

# High and low VAT only, for the subgroup comparisons
obs_pred_extremes <- obs_pred %>% filter(VAT_group != "Intermediate VAT")


# a) Observed eGFR in the first 12 months, high vs low VAT donors

# Observation level
obs_pred_extremes %>%
  filter(Days_Post <= 365) %>%
  group_by(VAT_group) %>%
  summarise(
    n_donors            = n_distinct(ID),
    n_obs               = n(),
    mean_egfr_obs       = mean(post_egfr),
    median_egfr_obs     = median(post_egfr),
    pct_obs_below_60    = 100 * mean(post_egfr < 60),
    donors_any_below_60 = n_distinct(ID[post_egfr < 60]),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1))) %>%
  print(width = Inf)

# Donor level (each donor's mean first-year eGFR, weighted equally)
obs_pred_extremes %>%
  filter(Days_Post <= 365) %>%
  group_by(VAT_group, ID) %>%
  summarise(donor_mean = mean(post_egfr), .groups = "drop") %>%
  group_by(VAT_group) %>%
  summarise(
    n_donors              = n(),
    mean_of_donor_means   = round(mean(donor_mean), 1),
    median_of_donor_means = round(median(donor_mean), 1),
    donors_mean_below_60  = sum(donor_mean < 60),
    .groups = "drop"
  ) %>%
  print(width = Inf)

# b) Calibration: observed vs predicted recovery ratio by time band (Table S1)

calib_primary <- calibration_table(obs_pred, "pred_primary", by_donor = TRUE)

# Table S1 body: high and low VAT, all time bands
calib_primary %>%
  filter(VAT_group != "Intermediate VAT") %>%
  mutate(across(where(is.double), ~ round(.x, 3))) %>%
  arrange(VAT_group, time_band) %>%
  print(n = Inf, width = Inf)

# Table S1 footnote: mean and maximum absolute difference, all three groups
calib_primary %>%
  group_by(VAT_group) %>%
  summarise(
    n_donors           = n_distinct(vat_groups$ID_chr[vat_groups$VAT_group == first(VAT_group)]),
    mean_difference    = round(mean(difference), 3),
    max_abs_difference = round(max(abs(difference)), 3),
    .groups = "drop"
  )


# ---------- High-low Difference Curve - Figure 4 ------------------------------

# Draws the high-low difference curve for the PRIMARY model only

# 1) Create a sequence of days and two grids for the two VAT_s values
days_seq <- seq(min(df_tidy$Days_Post, na.rm=TRUE),
                max(df_tidy$Days_Post, na.rm=TRUE))

make_new <- function(MV) {
  data.frame(
    Days_Post   = days_seq,
    VAT_s = MV,
    remKCV_s    = 0,
    ID = df_tidy$ID[1] # dummy value, will be excluded
  )
}
nd_low <- make_new(-1)
nd_high <- make_new( 1)

# 2) Get design (lpmatrix) for both conditions on the link scale (identity link here)
X_low  <- predict(m_gam, newdata = nd_low,  type = "lpmatrix", exclude = "s(ID)")
X_high <- predict(m_gam, newdata = nd_high, type = "lpmatrix", exclude = "s(ID)")

# 3) Difference, standard errors, and CIs using the full covariance
beta <- coef(m_gam)
V    <- vcov(m_gam)
Xdiff <- X_high - X_low

fit_diff <- as.vector(Xdiff %*% beta)
se_diff  <- sqrt(pmax(0, diag(Xdiff %*% V %*% t(Xdiff))))
lwr <- fit_diff - 1.96 * se_diff
upr <- fit_diff + 1.96 * se_diff
sig <- (lwr > 0) | (upr < 0)

diff_df <- tibble(
  Days_Post = days_seq,
  diff      = fit_diff,
  lower     = lwr,
  upper     = upr,
  sig       = sig
)

# 4) Plot: shaded where significantly different
ggplot(diff_df, aes(x = Days_Post, y = diff)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sig), alpha = 0.25) +
  geom_line(linewidth = 1.1) +
  scale_x_continuous(
    breaks = TIMEPOINTS_FULL,
    labels = TIMEPOINTS_FULL_LABELS,
    name = "Months post-nephrectomy",
    minor_breaks = NULL
  ) +
  scale_fill_manual(values = c(`TRUE` = "darkgreen", `FALSE` = "grey80"),
                    breaks = c(FALSE, TRUE),
                    labels = c("Includes 0", "Excludes 0"),
                    name = "Pointwise 95% CI") +
  labs(y = "Difference in fitted kidney function recovery ratio\nHigh (+1 SD) – Low (-1 SD) VAT",
       title = "High vs low VAT: fitted difference over time") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
  )
ggsave("VAT_difference_curve.png", width = 8, height = 5, dpi = 300, bg="white")

# First and last day where CI excludes zero
sig_day_first <- min(diff_df$Days_Post[diff_df$sig])
sig_day_last  <- max(diff_df$Days_Post[diff_df$sig])
cat("Significance window:", sig_day_first, "to", sig_day_last, "days\n")

cross_day <- diff_df$Days_Post[which(diff_df$diff > 0)[1]]
cat("Trajectories cross at day", cross_day, "(", round(cross_day / 30.4375, 1), "months )\n")


# ---------- Predicted Trajectories by VAT - Figure 2, part 1 ------------------

# Draws the predicted trajectories for the PRIMARY model only

# Prediction grid for VAT_s = -1, 0, +1 SD at average remKCV_s
adiposity_levels <- c(-1, 0, 1)
pred_grid <- expand.grid(
  Days_Post = seq(min(df_tidy$Days_Post, na.rm = TRUE),
                    max(df_tidy$Days_Post, na.rm = TRUE)),
  VAT_s     = adiposity_levels,
  remKCV_s  = 0,
  ID        = df_tidy$ID[1]
)

# Legend labels with raw VAT values (cohort donor-level mean and SD)
vat_raw    <- donor_ref$VAT_mean + adiposity_levels * donor_ref$VAT_sd
vat_labels <- c(sprintf("-1 SD (%.0f cm²)", vat_raw[1]),
                sprintf("Mean (%.0f cm²)",  vat_raw[2]),
                sprintf("+1 SD (%.0f cm²)", vat_raw[3]))

# Make sure types match and rows are finite
pred_grid <- pred_grid %>% mutate(
  Days_Post   = as.numeric(Days_Post),
  VAT_s = as.numeric(VAT_s),
  remKCV_s    = as.numeric(remKCV_s)
)
stopifnot(all(is.finite(pred_grid$Days_Post)),
          all(is.finite(pred_grid$VAT_s)),
          all(is.finite(pred_grid$remKCV_s)))

pred_grid$ID <- df_tidy$ID[1]  # Any valid ID value, doesn't matter which
# Predict from the GAMM 'gam' part (population-level smooth)
p <- predict(m_gam, newdata = pred_grid, se.fit = TRUE, exclude = "s(ID)")

pred_grid$pred <- as.numeric(p$fit)
pred_grid$se   <- as.numeric(p$se.fit)
pred_grid <- pred_grid %>% mutate(
  group      = factor(VAT_s, levels = adiposity_levels,
                      labels = vat_labels),
  lower      = pred - 1.96 * se,
  upper      = pred + 1.96 * se
)

# Sanity check - ensure we actually have finite predictions
if (!any(is.finite(pred_grid$pred))) {
  stop("All predictions are NA/Inf. Check variable names, types, and k values in ti().")
}

# Plot: predicted lines + ribbons for VAT_s = -1,0,1 + non-significant period shaded
p_pred_traj <- ggplot(pred_grid, aes(x = Days_Post, y = pred, color = group)) +
  annotate("rect", xmin = sig_day_last, xmax = max(pred_grid$Days_Post), ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  annotate("text",
           x = (sig_day_last + max(pred_grid$Days_Post)) / 2, y = Inf,
           label = "Difference not statistically significant",
           vjust = 1.5, size = 3.5, colour = "grey30") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(
    breaks = TIMEPOINTS_FULL,
    labels = TIMEPOINTS_FULL_LABELS,
    name = "Months post-nephrectomy",
    expand = expansion(mult = c(0.01, 0)) 
  ) +
  labs(
    y = "Kidney function recovery\n(post-/pre-operative eGFR ratio)",
    title = "Predicted recovery trajectories by VAT level",
    color = NULL,
    fill  = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
print(p_pred_traj)
ggsave("VAT_trajectory_plot.png", plot = p_pred_traj,
       width = 8, height = 5, dpi = 300, bg="white")


# ---------- Predicted eGFR Ratios by VAT at timepoints - Figure 2 in total ----

# From the primary model

pred_table <- expand.grid(
  Days_Post = TIMEPOINTS,
  VAT_s = c(-1,1),
  remKCV_s = 0,
  ID = df_tidy$ID[1]
)

p_table <- predict(m_gam, newdata = pred_table, se.fit = TRUE, exclude = "s(ID)")
pred_table$pred <- as.numeric(p_table$fit)
pred_table$se <- as.numeric(p_table$se.fit)

baseline_egfr <- df_renamed_backup %>%
  mutate(VAT_group = case_when(
    VAT > VAT_mean + VAT_sd ~ 1,   # High VAT = +1SD, matches VAT_s = 1
    VAT < VAT_mean - VAT_sd ~ -1,  # Low VAT = -1SD, matches VAT_s = -1
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(VAT_group)) %>%
  group_by(VAT_group) %>%
  summarise(baseline = mean(mean_pre_egfr, na.rm = TRUE)) %>%
  rename(VAT_s = VAT_group)

# Join the baselines and calculate absolute eGFR using group-specific means
pred_table <- pred_table %>%
  left_join(baseline_egfr, by = "VAT_s") %>%
  mutate(
    eGFR_predicted = pred * baseline,
    lower_CI = (pred - 1.96*se) * baseline,
    upper_CI = (pred + 1.96*se) * baseline,
    VAT_group = ifelse(VAT_s == -1, "LowVAT", "HighVAT"),
    Month = Days_Post / 30.4375
  )

# Pivot wider and add diff eGFR for each month between the groups
pred_table_wide <- pred_table %>%
  select(VAT_group, Month, Days_Post, pred, se, eGFR_predicted, lower_CI, upper_CI) %>%
  mutate(
    lower_ratio = pred - 1.96*se,
    upper_ratio = pred + 1.96*se
  ) %>%
  select(-se, -Days_Post) %>% # drop se and Days_Post before pivoting
  pivot_wider(
    names_from = VAT_group,
    values_from = c(pred, lower_ratio, upper_ratio, eGFR_predicted, lower_CI, upper_CI)
  ) %>%
  mutate(
    diff_eGFR = `eGFR_predicted_LowVAT` - `eGFR_predicted_HighVAT`
  ) %>%
  select(
    Month,
    `pred_LowVAT`, `lower_ratio_LowVAT`, `upper_ratio_LowVAT`,
    `eGFR_predicted_LowVAT`, `lower_CI_LowVAT`, `upper_CI_LowVAT`,
    `pred_HighVAT`, `lower_ratio_HighVAT`, `upper_ratio_HighVAT`,
    `eGFR_predicted_HighVAT`, `lower_CI_HighVAT`, `upper_CI_HighVAT`,
    diff_eGFR
  )

pred_table_wide %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  transmute(
    Month = Month,
    `Low VAT ratio (95% CI)` = paste0(pred_LowVAT, " (", lower_ratio_LowVAT, "–", upper_ratio_LowVAT, ")"),
    `Low VAT eGFR (95% CI)` = paste0(eGFR_predicted_LowVAT, " (", lower_CI_LowVAT, "–", upper_CI_LowVAT, ")"),
    `High VAT ratio (95% CI)` = paste0(pred_HighVAT, " (", lower_ratio_HighVAT, "–", upper_ratio_HighVAT, ")"),
    `High VAT eGFR (95% CI)` = paste0(eGFR_predicted_HighVAT, " (", lower_CI_HighVAT, "–", upper_CI_HighVAT, ")"),
    `Difference in eGFR` = round(diff_eGFR, 1)
  ) %>%
  print()

print(pred_table_wide %>% 
        mutate(across(where(is.numeric), ~ round(.x, 3))), 
      width = Inf)

# Build a small table from the computed differences above
fig_table <- pred_table_wide %>%
  mutate(Month = round(Month)) %>%
  mutate(
    Month_label = as.character(Month),
    eGFR_diff = round(diff_eGFR, 1)
  ) %>%
  select(Month_label, eGFR_diff) %>%
  rename(
    `\nMonth` = Month_label,
    `eGFR Difference\n(Low minus High VAT)` = eGFR_diff
  )

italic_rows <- which(round(pred_table_wide$Month) >= 24)

tbl_grob <- tableGrob(
  fig_table,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 12,
    core = list(
      fg_params = list(
        fontface = ifelse(seq_len(nrow(fig_table)) %in% italic_rows,
                          "italic", "plain")
      )
    )
  )
)

# Combine the above table to the predicted trajectories plot
p_combined <- p_pred_traj + wrap_elements(tbl_grob) +
  plot_layout(widths = c(3,1)) 
# +
# plot_annotation(
#    caption = "Italic rows indicate timepoints beyond the period of statistical significance (~21 months)",
#    theme = theme(plot.caption = element_text(size=9, hjust = 1, face = "italic"))
#  )
print(p_combined)
ggsave("VAT_trajectory_with_table.png", plot = p_combined,
       width = 11, height = 5, dpi = 300, bg = "white")


# ---------- Predicted Trajectories by remKCV - Figure 3 -----------------------
# Primary model only; VAT held at its mean, remKCV at -1, 0, +1 SD

# 1) Prediction grid
remkcv_levels <- c(-1, 0, 1)
remkcv_pct    <- 100 * (donor_ref$KCV_mean + remkcv_levels * donor_ref$KCV_sd)
remkcv_labels <- c(sprintf("-1 SD (%.1f%%)", remkcv_pct[1]),
                   sprintf("Mean (%.1f%%)",  remkcv_pct[2]),
                   sprintf("+1 SD (%.1f%%)", remkcv_pct[3]))

pred_grid_kcv <- expand.grid(
  Days_Post = seq(min(df_tidy$Days_Post), max(df_tidy$Days_Post)),
  remKCV_s  = remkcv_levels,
  VAT_s     = 0
)
pred_grid_kcv$ID <- df_tidy$ID[1]  # dummy value, random effect excluded below

p_kcv <- predict(m_gam, newdata = pred_grid_kcv, se.fit = TRUE, exclude = "s(ID)")

pred_grid_kcv <- pred_grid_kcv %>%
  mutate(
    pred  = as.numeric(p_kcv$fit),
    se    = as.numeric(p_kcv$se.fit),
    lower = pred - 1.96 * se,
    upper = pred + 1.96 * se,
    group = factor(remKCV_s, levels = remkcv_levels, labels = remkcv_labels)
  )

# 2) Plot: no raw points, so the y-axis zooms to the trajectories (reviewer request)
p_pred_kcv <- ggplot(pred_grid_kcv, aes(x = Days_Post, y = pred, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_viridis_d(end = 0.85, direction = -1) +
  scale_fill_viridis_d(end = 0.85, direction = -1) +
  scale_x_continuous(breaks = TIMEPOINTS_FULL, labels = TIMEPOINTS_FULL_LABELS,
                     name = "Months post-nephrectomy") +
  labs(
    y = "Kidney function recovery\n(post-/pre-operative eGFR ratio)",
    title = "Predicted recovery trajectories by remnant kidney cortex proportion",
    colour = NULL,
    fill   = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())

print(p_pred_kcv)
ggsave("Figure_remKCV_trajectory.png", plot = p_pred_kcv,
       width = 8, height = 5, dpi = 300, bg = "white")

# 3) Numbers for the results text: constant +1 vs -1 SD gap
kcv_gap_ratio <- 2 * coef(m_gam)["remKCV_s"]
cat("+1 vs -1 SD remKCV difference in recovery ratio:", round(kcv_gap_ratio, 3), "\n")
cat("Approximate eGFR equivalent at cohort mean baseline:",
    round(kcv_gap_ratio * mean(descriptive_data$mean_pre_egfr), 1), "mL/min/1.73m²\n")


# ---------- Table 2: VAT interaction + remKCV across models ------------------

models <- list(
  "Primary model"             = m_gam,
  "BMI adjustment"            = m_gam_BMI_sens,
  "VAT body-size indexing"    = m_gam_height_sens,
  "Post-Graves' inclusion"    = m_gam_graves_sens,
  "Comorbidities adjustment"  = m_gam_comorbidity_sens,
  "Age and sex adjustment"    = m_gam_demographic_fixed_sens,
  "Fixed linear VAT (SA6)"    = m_gam_fixed,
  "Sex-specific trajectories" = m_gam_sex_time_sens
)

table2_raw <- bind_rows(lapply(names(models), function(n) extract_row(models[[n]], n)))

# Full precision, for checking against the manuscript
print(table2_raw, width = Inf)

# Formatted, ready to paste
table2_raw %>%
  transmute(
    Model,
    edf  = ifelse(is.na(edf), "n/a", sprintf("%.3f", edf)),
    F    = ifelse(is.na(F_stat), "n/a", sprintf("%.1f", F_stat)),
    P_ti = case_when(is.na(p_ti) ~ "n/a", p_ti < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p_ti)),
    `remKCV β (95% CI)` = sprintf("%.3f (%.3f to %.3f)", remKCV_b, remKCV_lo, remKCV_hi),
    P_remKCV = ifelse(p_remKCV < 0.001, "<0.001", sprintf("%.3f", p_remKCV))
  ) %>%
  print(width = Inf)

# ---------- Recovery slopes by VAT level (response to reviewer Comment 17) ----

X_low  <- slope_X(m_gam, -1)
X_high <- slope_X(m_gam,  1)

slope_ci(m_gam, X_low)            # low VAT (-1 SD) slope per year
slope_ci(m_gam, X_high)           # high VAT (+1 SD) slope per year
slope_ci(m_gam, X_high - X_low)   # difference in slopes

# *** END ***

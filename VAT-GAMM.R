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

# Donor with Graves' disease diagnosis during follow-up
# Post-diagnosis observations excluded as hyperthyroidism independently
# elevates eGFR via haemodynamic and creatinine production mechanisms
GRAVES_DONOR_DOD <- as.Date("YYYY-MM-DD") # Replace with actual date before running


# ---------- Helper functions --------------------------------------------------
calc_egfr <- function(cr_umol, sex, age) {
  cr_mgdl <- cr_umol * 0.0113
  k <- ifelse(sex == "F", 0.7, 0.9)
  a <- ifelse(sex == "F", -0.241, -0.302)
  sex_coefficient <- ifelse(sex == "F", 1.012, 1.0)
  142 * (pmin(cr_mgdl / k, 1) ^ a) *
    (pmax(cr_mgdl / k, 1) ^ -1.2) *
    (0.9938 ^ age) * sex_coefficient
}
ratio_egfr <- function(pre_egfr, post_egfr) post_egfr / pre_egfr

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
         `Diabetes`, Hyperuricaemia)

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

# Create analysis variables - egfr_ratio, scaled covariates
df_tidy <- df_tidy %>%
  select(ID, Collect_Count, Days_Post, Cr, Gender, Age, mean_pre_egfr, VAT, remKCV_ratio, BMI, Laterality, Height, Dyslipidaemia, HTN) %>%
  mutate(
    currentAge = Age + (Days_Post/365),
    post_egfr = calc_egfr(Cr, Gender, currentAge),     # Calculate post-op eGFR with the current age
    egfr_ratio = ratio_egfr(mean_pre_egfr, post_egfr)
  ) %>%
  filter(Days_Post <= day_cap) %>%
  mutate(
    VAT_s = as.numeric(scale(VAT, center = TRUE, scale = TRUE)),
    remKCV_s    = as.numeric(scale(remKCV_ratio, center = TRUE, scale = TRUE)),
    BMI_s = as.numeric(scale(BMI, center = TRUE, scale = TRUE)),
    VATHeight_s = as.numeric(scale(VAT/Height^2, center = TRUE, scale = TRUE)),
    Dyslipid_s = as.numeric(scale(Dyslipidaemia, center = TRUE, scale = TRUE)),
    HTN_s = as.numeric(scale(HTN, center = TRUE, scale = TRUE)),
    Age_s = as.numeric(scale(Age, center = TRUE, scale = TRUE))
  ) %>%
  select(ID, Collect_Count, Days_Post, egfr_ratio, VAT_s, remKCV_s, BMI_s, Laterality, VATHeight_s, Dyslipidaemia, HTN, Gender, Age_s)

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

# ---------- Fit Models --------------------------------------------------------

# Fit primary GAMM model:
# - remKCV_s as a fixed linear population-level effect
# - s(Days_Post) and s(VAT_s) are fixed non-linear population-level effects
# - ti(Days_Post, VAT_s) is are fixed population-level smooths that show how VAT modifies recovery curve over time
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
gam.check(m_gam)
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
gam.check(m_gam_BMI_sens)
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
gam.check(m_gam_height_sens)
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
gam.check(m_gam_graves_sens)
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
gam.check(m_gam_comorbidity_sens)
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
gam.check(m_gam_demographic_fixed_sens)
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
gam.check(m_gam_fixed)
summary(m_gam_fixed)
concurvity(m_gam_fixed, full = FALSE)

# ---------- Predicted Trajectories --------------------------------------------

# Draws the predicted trajectories for the PRIMARY model only

# Prediction grid for VAT_s = -1, 0, +1 SD at average remKCV_s
adiposity_levels <- c(-1, 0, 1)
pred_grid <- expand.grid(
  Days_Post   = seq(min(df_tidy$Days_Post, na.rm = TRUE),
                    max(df_tidy$Days_Post, na.rm = TRUE)),
  VAT_s = adiposity_levels,
  remKCV_s    = mean(df_tidy$remKCV_s, na.rm = TRUE),
  ID = df_tidy$ID[1]
)

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
                      labels = c("-1SD", "0SD", "+1SD")),
  lower      = pred - 1.96 * se,
  upper      = pred + 1.96 * se
)

# Sanity check - ensure we actually have finite predictions
if (!any(is.finite(pred_grid$pred))) {
  stop("All predictions are NA/Inf. Check variable names, types, and k values in t2().")
}

# Plot: predicted lines + ribbons for VAT_s = -1,0,1; raw datapoints behind
p_pred_traj <- ggplot(pred_grid, aes(x = Days_Post, y = pred, color = group)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.18, color = NA) +
  geom_point(data = df_tidy,
             aes(x = Days_Post, y = egfr_ratio),
             inherit.aes = FALSE, alpha = 0.15, size = 1) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(
    breaks = TIMEPOINTS_FULL,
    labels = TIMEPOINTS_FULL_LABELS,
    name = "Months post-nephrectomy"
  ) +
  labs(
    y = "Kidney Function Recovery",
    title = "Predicted recovery trajectories by VAT level",
    color = "VAT level",
    fill  = "VAT level"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
print(p_pred_traj)
ggsave("VAT_trajectory_plot.png", plot = p_pred_traj,
       width = 8, height = 5, dpi = 300, bg="white")


# ---------- High-low Difference Curve -----------------------------------------

# Draws the high-low difference curve for the PRIMARY model only

# 1) Create a sequence of days and two grids for the two VAT_s values
days_seq <- seq(min(df_tidy$Days_Post, na.rm=TRUE),
                max(df_tidy$Days_Post, na.rm=TRUE))

make_new <- function(MV) {
  data.frame(
    Days_Post   = days_seq,
    VAT_s = MV,
    remKCV_s    = mean(df_tidy$remKCV_s, na.rm = TRUE),
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
  scale_fill_manual(values = c(`TRUE` = "tomato", `FALSE` = "grey80"),
                    breaks = c(FALSE, TRUE),
                    labels = c("Includes 0", "Excludes 0"),
                    name = "Pointwise 95% CI") +
  labs(y = "High (+1 SD) – Low (-1 SD) VAT\nDifference in fitted kidney function recovery",
       title = "High vs low VAT: fitted difference over time") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
  )
ggsave("VAT_difference_curve.png", width = 8, height = 5, dpi = 300, bg="white")

# First and last day where CI excludes zero
sig_days <- diff_df %>% filter(sig == TRUE)
cat("Significance window:", min(sig_days$Days_Post), "to", max(sig_days$Days_Post), "days\n")


# ---------- Predicted eGFR Ratios at timepoints -------------------------------

# From the primary model

pred_table <- expand.grid(
  Days_Post = TIMEPOINTS,
  VAT_s = c(-1,1),
  remKCV_s = mean(df_tidy$remKCV_s, na.rm=TRUE),
  ID = df_tidy$ID[1]
)

p_table <- predict(m_gam, newdata = pred_table, se.fit = TRUE, exclude = "s(ID)")
pred_table$pred <- as.numeric(p_table$fit)
pred_table$se <- as.numeric(p_table$se.fit)

# Calculate group-specific baseline eGFR means from the data
VAT_mean <- mean(df_renamed_backup$VAT, na.rm = TRUE)
VAT_sd <- sd(df_renamed_backup$VAT, na.rm = TRUE)

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
#  plot_annotation(
#    caption = "Italic rows indicate timepoints beyond the period of statistical significance (~21 months)",
#    theme = theme(plot.caption = element_text(size=9, hjust = 1, face = "italic"))
#  )
print(p_combined)
ggsave("VAT_trajectory_with_table.png", plot = p_combined,
       width = 11, height = 5, dpi = 300, bg = "white")


# ---------- Descriptive Statistics --------------------------------------------

# Descriptive statistics for analytic sample
descriptive_data <- df_renamed_backup %>%
  select(ID, Gender, Age, Laterality, VAT, BMI, remKCV_ratio, mean_pre_egfr) %>%
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
  select(Age, VAT, BMI, remKCV_ratio, mean_pre_egfr) %>%
  summarise(across(everything(), list(
    mean     = ~ round(mean(.x, na.rm = TRUE), 4),
    sd       = ~ round(sd(.x, na.rm = TRUE), 4),
    median   = ~ round(median(.x, na.rm = TRUE), 4),
    iqr      = ~ round(IQR(.x, na.rm = TRUE), 4),
    skewness = ~ round(moments::skewness(.x, na.rm = TRUE), 4),
    kurtosis = ~ round(moments::kurtosis(.x, na.rm = TRUE), 4)
  ))) %>%
  pivot_longer(everything(),
               names_to = c("variable",  "stat"),
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
cat("\nVAT mean:", round(mean(descriptive_data$VAT, na.rm=TRUE), 1), "cm2\n")
cat("VAT SD:", round(sd(descriptive_data$VAT, na.rm=TRUE), 1), "cm2\n")
cat("VAT -1SD:", round(mean(descriptive_data$VAT, na.rm=TRUE) - sd(descriptive_data$VAT, na.rm=TRUE), 1), "cm2\n")
cat("VAT +1SD:", round(mean(descriptive_data$VAT, na.rm=TRUE) + sd(descriptive_data$VAT, na.rm=TRUE), 1), "cm2\n")


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

# *** END ***
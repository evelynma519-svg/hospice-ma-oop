# =============================================================================
# 00_setup.R
# Packages, paths, parameters, and helpers. Sourced by run_all.R.
# Run everything from the repository root.
# =============================================================================

PKGS <- c("haven", "data.table", "sandwich", "lmtest", "ggplot2")
missing_pkgs <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Run once:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))", call. = FALSE)
}
for (p in PKGS) suppressPackageStartupMessages(library(p, character.only = TRUE))

if (!file.exists("README.md") || !dir.exists("R")) {
  stop("Please run the scripts from the repository root (the folder with README.md).",
       call. = FALSE)
}

DIR_RAW     <- file.path("data", "raw")
DIR_DERIVED <- file.path("data", "derived")
DIR_TABLES  <- file.path("output", "tables")
DIR_FIGS    <- file.path("output", "figures")
for (d in c(DIR_DERIVED, DIR_TABLES, DIR_FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

ANALYSIS_FILE <- file.path(DIR_DERIVED, "analysis.rds")          # controls per ORIG_CONTROL_RULE
ANALYSIS_ALIGNED <- file.path(DIR_DERIVED, "analysis_aligned.rds") # controls from last core wave

# ---- HRS Exit interview files -----------------------------------------------
# Section N (health services and insurance). Wave letter Q = 2018, R = 2020.
# Each Exit wave is linked to the deceased respondent's LAST core interview:
# 2018 Exit -> 2016 core (RAND wave 13); 2020 Exit -> 2018 core (wave 14).
EXIT_WAVES <- data.table(
  exit_year = c(2018L, 2020L),
  file      = c("X18N_R", "X20N_R"),
  prefix    = c("XQN", "XRN"),
  wave      = c(14L, 15L),          # wave label used in the original code
  core_wave = c(13L, 14L))          # RAND wave of the last core interview

EXIT_ITEMS <- c(ma_hmo       = "009",   # Medicare HMO / Medicare Advantage
                oop_nursing  = "119",   # out-of-pocket nursing home costs
                life_support = "306",   # life support in last two years
                oop_hospice  = "328",   # out-of-pocket hospice costs
                hospice_days = "437")   # days in hospice

# Values to treat as missing (don't know / refused). Check these against the
# HRS codebook for each item before relying on the results.
MISSING_CODES <- list(
  ma_hmo       = c(8, 9),
  life_support = c(8, 9),
  oop_hospice  = c(99997, 99998, 99999),
  oop_nursing  = c(999997, 999998, 999999),
  hospice_days = c(997, 998, 999))

RAND_FILE_PATTERN <- "^randhrs1992.*\\.dta$"     # RAND HRS Longitudinal File
RAND_IMP_PATTERN  <- "^randhrsimp.*\\.dta$"      # RAND HRS Detailed Imputations File

# RAND end-of-life hospice OOP cost (Exit interviews from 2014 on = wave 12+):
#   REMAHSPC amount, REOAHSPC any OOP expense (0/1), REFAHSPC imputation flag:
#   1 continuous value, 2 complete bracket, 3 incomplete bracket,
#   5 no value/bracket (fully imputed), 6 no expense.
# (REMCHSPC/REFCHSPC cover facility hospice for Exit waves 3-11 only.)
RAND_FLAG_LABELS <- c("1" = "Reported amount", "2" = "Complete bracket",
                      "3" = "Incomplete bracket", "5" = "Fully imputed",
                      "6" = "No expense")

# Which core wave the ORIGINAL analysis took its controls from is not recorded
# (the cleaning was done interactively). 01_build.R reports the sample size
# under each rule; set this to the rule that reproduces the paper's N.
#   "aligned": each person's last core wave (13 for 2018 Exit, 14 for 2020)
#   "13" / "14": the same wave for everyone
# Income in the original: the do-file keep list names HwIFTOT (RAND's
# imputation FLAG for income), but the published income coefficient
# (0.006 per dollar in Table 3) shows that dollar income HwITOT was used.
# 02_replicate.R reports the baseline model under both so this can be checked.
ORIG_INCOME <- "itot"     # "itot" = HwITOT (dollars); "iftot" = HwIFTOT (flag)

ORIG_CONTROL_RULE <- "13"   # reproduces the paper's N most closely (291 vs. 294)

# ---- Helpers ----------------------------------------------------------------

# Find a file anywhere under data/raw (case-insensitive), or NULL.
find_raw <- function(pattern) {
  hits <- list.files(DIR_RAW, pattern = pattern, recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) NULL else hits[1]
}

# Stata "reg ..., robust" = HC1.
vcov_robust <- function(fit) sandwich::vcovHC(fit, type = "HC1")

# Stata "glm ..., vce(robust)" = sandwich HC0 scaled by n / (n - 1).
vcov_glm_robust <- function(fit) {
  n <- nobs(fit)
  sandwich::vcovHC(fit, type = "HC0") * n / (n - 1)
}

coef_row <- function(fit, V, term, label) {
  est <- unname(coef(fit)[term]); se <- sqrt(V[term, term])
  data.table(model = label, term = term, estimate = est, se = se,
             p = 2 * pt(-abs(est / se), fit$df.residual), n = nobs(fit))
}

save_csv <- function(x, name) fwrite(x, file.path(DIR_TABLES, paste0(name, ".csv")))

save_fig <- function(plot, name, width = 7, height = 4.5) {
  ggsave(file.path(DIR_FIGS, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300)
}

# =============================================================================
# 02_replicate.R
# Replication of the original Stata analysis (final_do_file.do) in R, using
# the original-style controls (spouse demographics, as in the original keep
# list). Specifications follow the do-file line by line:
#   * OLS with robust (HC1) SEs; race entered as a single linear variable
#   * gamma GLM with log link and robust SEs (quasi-likelihood, so zeros are
#     allowed, as in Stata's glm)
#   * log(OOP) model on positive costs; interactions; wave fixed effects
# No methodological changes here; see 03_improved.R.
# =============================================================================

d <- copy(readRDS(ANALYSIS_FILE))
# Income version used in the original (see ORIG_INCOME in 00_setup.R).
d[, income_flag := income_orig]
if (ORIG_INCOME == "iftot") d[, income := income_orig]
d[, income_orig := NULL]
setnames(d, c("gender_orig", "race_orig", "eduyear_orig"), c("gender", "race", "eduyear"))
d[, log_oophosp := fifelse(oophosp > 0, log(oophosp), NA_real_)]

# ---- Table 2: summary statistics --------------------------------------------
sum_vars <- c("oophosp", "maenroll", "gender", "race", "income", "eduyear",
              "oopmdre", "oopmdsp", "renrsnit")
summ <- rbindlist(lapply(sum_vars, function(v) {
  x <- d[[v]]
  data.table(variable = v, mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
             min = suppressWarnings(min(x, na.rm = TRUE)),
             max = suppressWarnings(max(x, na.rm = TRUE)), n = sum(!is.na(x)))
}))
save_csv(summ, "table2_summary_statistics")

# ---- Models -----------------------------------------------------------------
gamma_log <- quasi(link = "log", variance = "mu^2")   # = Stata glm, family(gamma) link(log)
ctrl <- "gender + race + eduyear + income"

specs <- list(
  "T3 OLS no controls"      = list("ols", "oophosp ~ maenroll"),
  "T3 GLM no controls"      = list("glm", "oophosp ~ maenroll"),
  "T3 OLS baseline"         = list("ols", paste("oophosp ~ maenroll +", ctrl)),
  "T3a GLM baseline"        = list("glm", paste("oophosp ~ maenroll +", ctrl)),
  "T3b OLS extra controls"  = list("ols", paste("oophosp ~ maenroll +", ctrl,
                                               "+ oopnurs + lifesupt + renrsnit + hospday")),
  "T3c log OOP"             = list("ols", paste("log_oophosp ~ maenroll +", ctrl)),
  "T4 MA x income"          = list("ols", "oophosp ~ maenroll * income + gender + race + eduyear"),
  "T4 MA x life support"    = list("ols", paste("oophosp ~ maenroll * lifesupt +", ctrl)),
  "T4 MA x OOP nursing"     = list("ols", paste("oophosp ~ maenroll * oopnurs +", ctrl)),
  "T4 MA x gender"          = list("ols", "oophosp ~ maenroll * gender + race + eduyear + income + oopnurs + renrsnit"),
  "T4 MA x hospice days"    = list("ols", paste("oophosp ~ maenroll * hospday +", ctrl, "+ oopnurs + renrsnit")),
  "T5 wave FE"              = list("ols", paste("oophosp ~ maenroll +", ctrl,
                                               "+ oopmdre + oopmdsp + renrsnit + factor(wave)")),
  "T5 wave FE, MA x gender" = list("ols", "oophosp ~ maenroll * gender + race + eduyear + income + oopnurs + renrsnit + factor(wave)"))

fit_spec <- function(label, spec) {
  f <- as.formula(spec[[2]])
  vars <- all.vars(f)
  empty <- vars[vapply(vars, function(v) all(is.na(d[[v]])), logical(1))]
  if (length(empty)) {
    message("Skipping '", label, "': no data for ", paste(empty, collapse = ", "))
    return(NULL)
  }
  df <- na.omit(d[, ..vars])
  if (spec[[1]] == "ols") {
    fit <- lm(f, data = df); V <- vcov_robust(fit)
  } else {
    # Start from the sample mean: with zero costs, R's default start (mu = y)
    # implies log(0). Stata's glm starts from the mean as well.
    y0 <- df[[all.vars(f)[1]]]
    fit <- glm(f, data = df, family = gamma_log, mustart = rep(mean(y0), length(y0)),
               control = glm.control(maxit = 100))
    V <- vcov_glm_robust(fit)
  }
  if (label == "T3 OLS no controls") V <- vcov(fit)   # original used default SEs here
  terms <- intersect(c("maenroll", grep("^maenroll:", names(coef(fit)), value = TRUE),
                       "gender", "race", "eduyear", "income"),
                     names(coef(fit)))
  rbindlist(lapply(terms, function(t) coef_row(fit, V, t, label)))
}
results <- rbindlist(Map(fit_spec, names(specs), specs))
save_csv(results, "replication_models")

# Percentage interpretation for the log and GLM models (exp(b) - 1).
pct <- results[model %in% c("T3a GLM baseline", "T3c log OOP", "T3 GLM no controls") & term == "maenroll",
               .(model, estimate, pct_change = 100 * (exp(estimate) - 1))]
save_csv(pct, "replication_pct_interpretation")

# ---- Which income variable reproduces the paper? ----------------------------
income_check <- rbindlist(lapply(c("income", "income_flag"), function(v) {
  f <- as.formula(paste("oophosp ~ maenroll + gender + race + eduyear +", v))
  df <- na.omit(d[, all.vars(f), with = FALSE])
  fit <- lm(f, data = df)
  data.table(income_variable = fifelse(v == "income", "HwITOT (dollars)", "HwIFTOT (flag)"),
             maenroll = coef(fit)[["maenroll"]], income_coef = coef(fit)[[v]], n = nobs(fit),
             paper_maenroll = 157.189, paper_income_coef = 0.006, paper_n = 294)
}))
save_csv(income_check, "replication_income_variable_check")
print(income_check, digits = 4)

# ---- Compare with the paper (config/paper_values.csv) -----------------------
paper <- fread(file.path("config", "paper_values.csv"))
if (!"term" %in% names(paper)) paper[, term := "maenroll"]   # older one-row-per-model format
paper <- paper[!is.na(estimate)]
check <- merge(paper, results[, .(model, term, rep_estimate = estimate, rep_se = se, rep_n = n)],
               by = c("model", "term"), all.x = TRUE, sort = FALSE)
check[, `:=`(diff = rep_estimate - estimate,
             same_sign = sign(rep_estimate) == sign(estimate),
             same_significance_5pct = (abs(rep_estimate / rep_se) > 1.96) == (abs(estimate / se) > 1.96))]
save_csv(check, "replication_check")
print(check, digits = 4)
print(results, digits = 3)
print(pct, digits = 3)

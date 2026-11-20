# =============================================================================
# 03_improved.R
# Improved analysis of the association between Medicare Advantage (MA)
# enrollment and out-of-pocket (OOP) hospice costs at the end of life.
#
# Changes relative to the original:
#   1. Outcome: RAND's harmonized, imputed hospice OOP cost (REMAHSPC, RAND HRS
#      Detailed Imputations File), which uses bracket answers when the next of
#      kin could not give an exact amount. The raw Exit item is kept as a
#      robustness check.
#   2. Controls are the deceased respondent's OWN sex, race/ethnicity and
#      education (not the spouse's), and household income (HwITOT, not its
#      imputation flag), from their last core interview.
#   3. Race/ethnicity enters as categories, not as one linear variable.
#   4. Within each block, all models use the same estimation sample.
#   5. Two-part model for a cost outcome with many zeros: (a) probability of
#      any OOP hospice cost (logit), (b) amount among those with a positive
#      cost (gamma GLM, log link), combined into the average difference in
#      expected OOP cost, with bootstrap standard errors.
#   7. Results for all decedents and for hospice users only.
#   8. Minimum detectable effects (80% power, 5% two-sided test) to show how
#      large a difference the data could have detected.
#   9. Robustness: keep decedents with missing income (income set to 0 plus a
#      missing-income indicator) instead of dropping them.
# MA enrollment is self-selected, so all estimates are associations.
# =============================================================================

d <- readRDS(ANALYSIS_ALIGNED)
d[, `:=`(log_income = log1p(pmax(income, 0)),
         female     = 1L - gender_resp,
         late_exit  = as.integer(wave == 15L),
         hospice_user = fifelse(!is.na(hospice_stay), as.integer(hospice_stay == 1),
                                as.integer(hospday > 0)),
         imputed    = as.integer(oop_hospice_flag %in% c(2, 3, 5)),
         income_missing = as.integer(is.na(income)))]
d[, log_income_filled := fifelse(is.na(log_income), 0, log_income)]

# Race/ethnicity as dummies (White = reference). Dummies rather than a factor
# so that predictions still work when a group has no positive costs in a
# (bootstrap) sample. Groups with fewer than 10 decedents with a positive cost
# (on either outcome) are merged into "Other" so the amount model is estimable.
pos_counts <- d[(oophosp > 0 | oop_hospice_rand > 0) & !is.na(race_resp), .N, by = race_resp]
small <- setdiff(pos_counts[N < 10, race_resp], c("White", "Other"))
if (length(small)) {
  message("Merging into 'Other' (fewer than 10 positive-cost decedents): ",
          paste(small, collapse = ", "))
  d[race_resp %in% small, race_resp := "Other"]
}
# Dummy names use their own prefix so they cannot be confused with other
# race variables in the data (e.g. race_orig, the SPOUSE's race).
for (g in setdiff(unique(na.omit(d$race_resp)), "White")) {
  d[, (paste0("rdum_", tolower(g))) := as.integer(race_resp == g)]
}
race_dummies <- grep("^rdum_", names(d), value = TRUE)
message("Race/ethnicity dummies: ", paste(race_dummies, collapse = ", "))
CTRL_SETS <- list(
  "Complete cases" = c("female", race_dummies, "eduyear_resp", "log_income", "late_exit"),
  "Missing-income indicator" = c("female", race_dummies, "eduyear_resp",
                                 "log_income_filled", "income_missing", "late_exit"))
ctrl_vars <- CTRL_SETS[["Complete cases"]]
ctrl <- paste(ctrl_vars, collapse = " + ")

gamma_log <- quasi(link = "log", variance = "mu^2")

# Two-part model on outcome column y; returns the average difference in
# expected OOP cost between MA = 1 and MA = 0, holding everything else fixed.
two_part <- function(df) {
  p1 <- glm(as.formula(paste("any_oop ~ maenroll +", ctrl)), data = df, family = binomial())
  pos <- df[y > 0]
  p2 <- glm(as.formula(paste("y ~ maenroll +", ctrl)), data = pos,
            family = gamma_log, mustart = rep(mean(pos$y), nrow(pos)),
            control = glm.control(maxit = 100))
  pr1 <- function(ma) suppressWarnings(predict(p1, copy(df)[, maenroll := ma], type = "response"))
  pr2 <- function(ma) suppressWarnings(predict(p2, copy(df)[, maenroll := ma], type = "response"))
  dollar <- mean(pr1(1) * pr2(1) - pr1(0) * pr2(0))
  # With few positive costs the amount model can diverge in a resample;
  # such draws are flagged and dropped from the bootstrap.
  ok <- p2$converged && is.finite(dollar) && abs(dollar) < 10 * max(df$y)
  list(p1 = p1, p2 = p2, ok = ok,
       ame_prob   = mean(pr1(1) - pr1(0)),
       ame_dollar = if (ok) dollar else NA_real_)
}
MIN_POSITIVE_FOR_AMOUNT <- 30   # below this, the amount model is not reported

OUTCOMES <- c("RAND imputed (REMAHSPC)" = "oop_hospice_rand",
              "Raw Exit item"           = "oophosp")
B <- 500
MDE_FACTOR <- qnorm(0.975) + qnorm(0.80)   # = 2.80: 80% power, 5% two-sided
SPECS <- list(list(outcome = "RAND imputed (REMAHSPC)", controls = "Complete cases"),
              list(outcome = "Raw Exit item",           controls = "Complete cases"),
              list(outcome = "RAND imputed (REMAHSPC)", controls = "Missing-income indicator"))
out <- list()
attrition <- list()
for (sp in SPECS) {
  o <- sp$outcome
  ctrl_vars <- CTRL_SETS[[sp$controls]]
  ctrl <- paste(ctrl_vars, collapse = " + ")
  base <- d[!is.na(get(OUTCOMES[[o]]))]
  base[, y := get(OUTCOMES[[o]])]
  # Sample attrition: how many observations each requirement removes, in order.
  steps <- c("y", "maenroll", "female", "race_resp", "eduyear_resp",
             if (sp$controls == "Complete cases") "log_income")
  keep <- rep(TRUE, nrow(base)); n_prev <- nrow(base)
  att <- data.table(outcome = o, step = "Exit respondents with outcome", n = n_prev)
  for (v in steps[-1]) {
    keep <- keep & !is.na(base[[v]])
    att <- rbind(att, data.table(outcome = o, step = paste("non-missing", v), n = sum(keep)))
  }
  base <- na.omit(base, cols = c("y", "maenroll", ctrl_vars))
  att <- rbind(att, data.table(outcome = o, step = "estimation sample (all controls)", n = nrow(base)))
  att[, controls := sp$controls]
  attrition[[paste(o, sp$controls)]] <- att
  base[, any_oop := as.integer(y > 0)]
  blocks <- list("All decedents" = base, "Hospice users" = base[hospice_user == 1])
  for (s in names(blocks)) {
    df <- blocks[[s]]
    if (nrow(df) < 30 || sum(df$y > 0) < 5 || uniqueN(df$maenroll) < 2) {
      message("Skipping '", o, " / ", s, "': too few observations (n = ", nrow(df),
              ", positive = ", sum(df$y > 0), ").")
      next
    }
    ols <- lm(as.formula(paste("y ~ maenroll +", ctrl)), data = df)
    tp  <- two_part(df)
    amount_ok <- sum(df$y > 0) >= MIN_POSITIVE_FOR_AMOUNT && tp$ok
    boot <- replicate(B, {
      bd <- df[sample.int(nrow(df), replace = TRUE)]
      r <- tryCatch(two_part(bd), error = function(e) NULL)
      if (is.null(r)) c(NA, NA) else c(r$ame_prob, r$ame_dollar)
    })
    if (!amount_ok) boot[2, ] <- NA
    V2 <- vcov_glm_robust(tp$p2)
    q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
    out[[length(out) + 1]] <- data.table(
      outcome = o, controls = sp$controls, sample = s,
      n = nrow(df), n_positive = sum(df$y > 0), share_ma = mean(df$maenroll),
      share_imputed = if (OUTCOMES[[o]] == "oop_hospice_rand") mean(df$imputed) else NA_real_,
      mean_oop = mean(df$y),
      ols_dollars = coef(ols)[["maenroll"]], ols_se = sqrt(vcov_robust(ols)["maenroll", "maenroll"]),
      ols_ci_low = coef(ols)[["maenroll"]] - 1.96 * sqrt(vcov_robust(ols)["maenroll", "maenroll"]),
      ols_ci_high = coef(ols)[["maenroll"]] + 1.96 * sqrt(vcov_robust(ols)["maenroll", "maenroll"]),
      part1_prob_diff = tp$ame_prob,
      part1_ci_low = q(boot[1, ], 0.025), part1_ci_high = q(boot[1, ], 0.975),
      part2_pct_diff = if (amount_ok) 100 * (exp(coef(tp$p2)[["maenroll"]]) - 1) else NA_real_,
      part2_se = if (amount_ok) sqrt(V2["maenroll", "maenroll"]) else NA_real_,
      twopart_dollars = if (amount_ok) tp$ame_dollar else NA_real_,
      twopart_ci_low = q(boot[2, ], 0.025), twopart_ci_high = q(boot[2, ], 0.975),
      bootstrap_reps_used = sum(!is.na(boot[2, ])),
      share_positive = mean(df$y > 0),
      mde_ols_dollars = MDE_FACTOR * sqrt(vcov_robust(ols)["maenroll", "maenroll"]),
      mde_part1_prob  = MDE_FACTOR * sd(boot[1, ], na.rm = TRUE))
  }
}
improved <- rbindlist(out)
save_csv(improved, "improved_two_part_model")
save_csv(rbindlist(attrition), "improved_sample_attrition")
print(rbindlist(attrition))
print(improved, digits = 3)
if (any(is.na(improved$twopart_dollars))) {
  message("Amount model (part 2) not reported where fewer than ", MIN_POSITIVE_FOR_AMOUNT,
          " decedents have a positive cost.")
}

# ---- Figure: distribution of positive OOP hospice costs by MA status --------
pos <- d[!is.na(maenroll) & oop_hospice_rand > 0]
if (nrow(pos) > 0) {
  p <- ggplot(pos, aes(oop_hospice_rand,
                       fill = factor(maenroll, levels = 0:1,
                                     labels = c("Traditional Medicare", "Medicare Advantage")))) +
    geom_histogram(bins = 25, alpha = 0.6, position = "identity") +
    scale_x_log10(labels = scales::label_dollar()) +
    labs(title = "Out-of-pocket hospice costs among decedents with any cost",
         subtitle = "RAND imputed measure (REMAHSPC), HRS Exit interviews 2018 and 2020",
         x = "OOP hospice cost (log scale)", y = "Decedents", fill = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  save_fig(p, "fig_oop_distribution_by_ma")
}

# Medicare Advantage and Out-of-Pocket Hospice Costs at the End of Life: Replication and Extension

A replication, in R, of my econometrics research project on whether Medicare Advantage (MA) enrollment is associated with higher out-of-pocket (OOP) hospice costs in the last year of life, using the Health and Retirement Study (HRS) Exit interviews. The original analysis was written in Stata as an individual project (Tufts University, Econometrics). This repository rebuilds the dataset from the raw files, reproduces the original specifications, and extends them.

---

## Research question

Hospice care is carved out of MA: once an MA enrollee elects hospice, traditional Medicare pays for it. Differences in OOP hospice costs between MA and traditional Medicare decedents may still arise through care coordination, provider networks, or selection into MA. The analysis compares OOP hospice costs reported by next of kin in the HRS Exit interviews for decedents enrolled in MA versus traditional Medicare. **MA enrollment is self-selected, so all estimates are associations, not causal effects.**

## What this repository does

| Step | Script | Notes |
|---|---|---|
| Rebuild the dataset | `R/01_build.R` | The original cleaning was run interactively and only partly saved. This script rebuilds every step from the raw Exit files and RAND HRS, and reports the sample size under each possible construction (`build_sample_sizes_by_rule.csv`). |
| Original specifications | `R/02_replicate.R` | Every model in the original do-file: OLS, gamma GLM, log model, interactions, wave fixed effects. |
| Improved analysis | `R/03_improved.R` | Corrected controls, a common estimation sample, a two-part model, and results for hospice users. |

### Changes in the improved analysis

1. **A better-measured outcome.** The main outcome is RAND's harmonized end-of-life hospice OOP cost, `REMAHSPC`, from the RAND HRS Detailed Imputations File. When the next of kin could not give an exact amount, RAND uses their answers to the follow-up dollar brackets to impute one, and a flag (`REFAHSPC`) records which values are reported, bracketed, or fully imputed. The raw Exit item used in the original analysis is kept as a robustness check. `build_raw_vs_rand_outcome.csv` and `build_rand_hospice_flags.csv` compare the two.
2. **Own demographics, not the spouse's.** The original keep list took gender, race, and education from RAND's `s*` variables, which describe the **spouse**. That also drops every unmarried decedent. The improved analysis uses the decedent's own `rasex`, `raracem`, `rahispan`, and `raedyrs`.
3. **Income in logs.** Household income (`HwITOT`) is right-skewed, so it enters as log(1 + income) rather than in dollars.
4. **Controls from the right interview.** Each decedent is linked to their last core interview: 2016 (wave 13) for the 2018 Exit, 2018 (wave 14) for the 2020 Exit.
5. **Race/ethnicity as categories** (White, Black, Other, Hispanic) instead of a single linear variable.
6. **One estimation sample across models**, so differences between columns are not driven by different samples.
7. **Two-part model.** Many decedents have zero OOP hospice cost. Part 1 models whether any cost occurs (logit). Part 2 models the amount among those with a positive cost (gamma GLM, log link). The two parts are combined into the average difference in expected OOP cost between MA and traditional Medicare, with bootstrap standard errors.
8. **Hospice users.** Results are also shown for decedents who received any hospice care (RAND `REHPCSTY`).
9. **Interpretation of log coefficients** as percentage differences, exp(b) − 1.
10. **Minimum detectable effects** (80% power, 5% two-sided test), to show how large a difference the data could have detected.
11. **Robustness to missing income:** decedents without income data are kept, with income set to zero and a missing-income indicator added, instead of being dropped.

## Repository structure

```
hospice-ma-oop/
├── README.md
├── run_all.R              # runs everything in order
├── config/
│   └── paper_values.csv   # the original paper's estimates, for comparison
├── R/
│   ├── 00_setup.R         # packages, paths, parameters, helpers
│   ├── 01_build.R         # raw HRS files -> analysis data
│   ├── 02_replicate.R     # original specifications
│   └── 03_improved.R      # improved analysis
├── data/
│   ├── raw/               # HRS files go here (never tracked)
│   └── derived/           # created by the scripts (never tracked)
└── output/
    ├── tables/
    └── figures/
```

## Data

HRS data are free with registration at https://hrsdata.isr.umich.edu. The HRS conditions of use do not allow redistribution, so no data are included here.

Place these files anywhere under `data/raw/`. Subfolders are fine, because the scripts search recursively:

| File | Source |
|---|---|
| `X18N_R` (2018 Exit, Section N) | HRS 2018 Exit. Either `X18N_R.dta`, or `X18N_R.dct` (from the "sta" download) plus `X18N_R.da` (from the "da" download) |
| `X20N_R` (2020 Exit, Section N) | HRS 2020 Exit, same options |
| `randhrs1992_2022v2.dta` | RAND HRS Longitudinal File 2022 (V2), Stata version |
| `randhrsimp1992_2022v2.dta` | RAND HRS Detailed Imputations File 2022 (V2), Stata version |

The original project used RAND HRS 1992–2020 (V2) and did not use the Detailed Imputations File. This repository uses the 2022 (V2) release, which contains the same waves but may include minor data revisions.

## How to run

Requires R (≥ 4.1). Install the packages once:

```r
install.packages(c("haven", "data.table", "sandwich", "lmtest", "ggplot2"))
```

Then, from the repository root:

```bash
Rscript run_all.R
```

**Matching the original sample.** Which core wave the original controls came from was not recorded. `output/tables/build_sample_sizes_by_rule.csv` shows the baseline sample size under each rule; wave 13 for all decedents reproduces the paper, so `ORIG_CONTROL_RULE` in `R/00_setup.R` is set to `"13"`. The paper's estimates are stored in `config/paper_values.csv`, and `replication_check.csv` compares them coefficient by coefficient.

## Results

### Replication

`R/02_replicate.R` re-runs every specification in the original do-file. The original cleaning was done interactively and only partly saved, so the construction had to be reconstructed. Taking gender, race, and education from the **spouse's** records in wave 13 for all decedents, as the original keep list does, reproduces the original sample sizes: 31 and 23 observations in Tables 3b and 3c, exactly as published, and 291 in Table 3 (294 published).

| MA coefficient | Paper | Replication |
|---|---|---|
| Table 3, OLS | 157.2 (480.9), N = 294 | 128.3 (484.3), N = 291 |
| Table 3a, gamma GLM | 0.779 (0.733), N = 294 | 0.326 (0.721), N = 291 |
| Table 3b, OLS with extra controls | −58.9 (314.9), N = 31 | −30.9 (324.5), N = 31 |
| Table 3c, log OOP | 2.109** (0.733), N = 23 | 1.961** (0.804), N = 23 |

Robust standard errors in parentheses; ** p < 0.05. The MA coefficient has the same sign and the same significance at the 5% level as in the paper in all four models. The remaining differences in magnitude are consistent with data revisions between RAND HRS 2020 (V2), used originally, and 2022 (V2). One difference is unresolved: the race coefficient has the opposite sign in all four models, which suggests the original recoded race in a step that was not saved. The full comparison is in `output/tables/replication_check.csv`.

Three features of the original construction explain its small samples and fragile estimates:

- **Spouse instead of own demographics.** Using the spouse's records drops every decedent without a surviving spouse with core data. Under the same wave rule, the decedent's own demographics leave 668 observations instead of 291.
- **Raw Exit item asked only of hospice users.** The raw OOP hospice question is asked only when hospice care was received, so the original sample is implicitly limited to hospice users. RAND's measure records zero for decedents without hospice care and fills in bracketed or missing amounts.
- **The one significant result rests on 23 people.** The log-OOP estimate (Table 3c) uses only decedents with positive costs and complete spouse data. The improved specifications below, on samples 20 to 70 times larger, find no comparable difference.

### Improved analysis

Main specification: RAND imputed hospice OOP cost, decedent's own demographics, all decedents with complete data (**N = 1,177**; `output/tables/improved_two_part_model.csv`).

| MA vs. traditional Medicare | Estimate | 95% CI | Minimum detectable effect |
|---|---|---|---|
| Probability of any OOP hospice cost (sample mean 7.8%) | +0.2 pp | −2.9 to +3.2 pp | 4.3 pp |
| Expected OOP hospice cost, two-part model (sample mean $576) | −$199 | −$705 to +$278 | |
| Expected OOP hospice cost, OLS | −$81 | −$606 to +$444 | $751 |

- **OOP hospice costs are rare.** About 90% of decedents in the 2018 and 2020 Exit interviews had no OOP hospice expense (`output/tables/build_rand_hospice_flags.csv`). This is consistent with the Medicare hospice benefit's minimal cost sharing, and it leaves little variation for MA status to explain.
- **No detectable difference by MA status.** Every estimate is small relative to its uncertainty.
- **What the data can and cannot rule out.** For the probability of any OOP cost, the confidence interval excludes an increase of more than about 3 percentage points. For dollar amounts, the minimum detectable effect ($751) exceeds the average cost itself ($576), so modest dollar differences cannot be distinguished from zero. The appropriate conclusion is "no evidence of a difference, with limited power for dollar amounts", not "no relationship".

**Robustness**

| Specification | N | Any OOP cost (pp) | Two-part ($) | OLS ($) |
|---|---|---|---|---|
| Main (complete cases) | 1,177 | +0.2 [−2.9, +3.2] | −199 [−705, +278] | −81 [−606, +444] |
| Hospice users only | 632 | +0.9 [−4.0, +6.0] | −337 [−1,393, +599] | −128 [−1,081, +825] |
| Keep missing income (indicator) | 1,554 | +2.5 [−0.0, +5.2] | +46 [−346, +421] | +89 [−345, +523] |
| Keep missing income, hospice users | 834 | +4.9 [−0.3, +9.8] | +99 [−652, +771] | +176 [−615, +967] |
| Raw Exit item (hospice users) | 579 | −1.7 [−6.0, +2.3] | +107 [−545, +801] | −11 [−512, +490] |

95% confidence intervals in brackets: bootstrap percentile intervals for the probability and two-part estimates, robust (HC1) intervals for OLS.

The dollar estimates change sign across specifications and are never distinguishable from zero. The one pattern worth noting is in the probability of any OOP cost. When decedents with missing income are kept, MA enrollees are 2.5 to 4.9 percentage points more likely to have any OOP hospice cost, and the confidence intervals end just below zero. This is not robust: it is close to zero in the main sample and negative with the raw Exit item. It is best read as a question for larger data, not as a finding.

### Limitations

- MA status at death is missing for about 30% of Exit respondents (690 of 2,246), mostly next-of-kin "don't know" answers; these decedents are excluded.
- Household income from the last core interview is missing for 377 of the remaining decedents; the missing-income specification shows how much this matters.
- MA enrollment is self-selected; estimates are associations.
- Only about 90 to 130 decedents have positive costs, so the amount model (part 2) is imprecise. Some bootstrap resamples do not converge; these are dropped and counted in `bootstrap_reps_used` (about 430 of 500 kept).

## Notes

- **Missing-value codes** for the Exit items (don't know/refused) are listed in `R/00_setup.R` and should be checked against the HRS codebook.
- **Why Exit interviews.** Hospice costs at the end of life are only observed in the Exit interview, answered by a proxy (usually a family member) after the respondent's death. Controls come from the respondent's last core interview, linked through RAND HRS. MA status at death is not harmonized by RAND, so it comes from the raw Exit item (`XQN009` / `XRN009`).
- **Hospice cost variables by Exit wave.** For Exit interviews from 2014 onward, RAND stores hospice costs in any location in `REMAHSPC` / `REFAHSPC`. `REMCHSPC` covers hospice in a facility for earlier Exit waves only.
- **Income variable in the original.** The original do-file's keep list names `HwIFTOT`, which in RAND HRS is the imputation flag for household income. The published income coefficient (0.006 per dollar) shows that dollar income `HwITOT` was actually used. The replication uses `HwITOT`; `replication_income_variable_check.csv` compares both.
- **Variable names** follow RAND HRS 2022 (V2): sex is `RASEX` / `SwSEX` (older releases used `*GENDER`; the code accepts either). `RENRSNIT` is the number of nights in a nursing home reported in the Exit interview.
- The original also merged each Exit section with its preload (`PR`) file. That merge only added identifiers already present in Section N, so it is not repeated here.

## Contact

Evelyn Ma · evelynma519@gmail.com · [LinkedIn](https://www.linkedin.com/in/-evelynma-)

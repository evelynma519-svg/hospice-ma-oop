# =============================================================================
# 01_build.R
# Build the analysis file from the HRS Exit interviews (2018, 2020) and the
# RAND HRS longitudinal file.
# Output: data/derived/analysis.rds and output/tables/build_*.csv
#
# The original cleaning was run interactively in Stata and only partly saved
# (data_cleaning_do.do). This script rebuilds every step and keeps two versions
# of the demographic controls side by side:
#   *_orig : SPOUSE variables (s*gender, s*racem, s*edyrs), as in the original
#            keep list. In RAND HRS, the s-prefix refers to the spouse.
#   *_resp : the deceased respondent's own variables (ragender, raracem,
#            rahispan, raedyrs).
# =============================================================================

# ---- Read one Exit section: .dta if available, else .dct + .da --------------
read_dct_da <- function(dct_path, da_path, keep) {
  dct <- readLines(dct_path, warn = FALSE)
  rx  <- "_column\\(\\s*(\\d+)\\s*\\)\\s+(\\w+)\\s+(\\w+)\\s+%(\\d+)"
  m   <- regmatches(dct, regexec(rx, dct))
  m   <- do.call(rbind, m[lengths(m) == 5])
  spec <- data.table(start = as.integer(m[, 2]), type = m[, 3],
                     name = toupper(m[, 4]), width = as.integer(m[, 5]))
  spec <- spec[name %in% keep]
  lines <- readLines(da_path, warn = FALSE)
  out <- as.data.table(lapply(seq_len(nrow(spec)), function(i) {
    v <- trimws(substr(lines, spec$start[i], spec$start[i] + spec$width[i] - 1))
    v[v == ""] <- NA
    if (startsWith(spec$type[i], "str")) v else as.numeric(v)
  }))
  setnames(out, spec$name)
  out
}

read_exit <- function(file, prefix) {
  keep <- c("HHID", "PN", paste0(prefix, EXIT_ITEMS))
  dta <- find_raw(paste0("^", file, "\\.dta$"))
  if (!is.null(dta)) {
    x <- as.data.table(haven::zap_labels(haven::read_dta(dta)))
    setnames(x, toupper(names(x)))
  } else {
    dct <- find_raw(paste0("^", file, "\\.dct$"))
    da  <- find_raw(paste0("^", file, "\\.da$"))
    if (is.null(dct) || is.null(da)) {
      stop("Could not find ", file, ".dta, or ", file, ".dct + ", file, ".da, under data/raw/. ",
           "The .da data files come in the HRS 'da' download (e.g. x18da.zip); ",
           "the .dct dictionaries in the 'sta' download.", call. = FALSE)
    }
    x <- read_dct_da(dct, da, keep)
  }
  absent <- setdiff(keep, names(x))
  if (length(absent)) stop(file, " is missing: ", paste(absent, collapse = ", "), call. = FALSE)
  x <- x[, ..keep]
  setnames(x, paste0(prefix, EXIT_ITEMS), names(EXIT_ITEMS))
  x
}

exit <- rbindlist(lapply(seq_len(nrow(EXIT_WAVES)), function(i) {
  w <- EXIT_WAVES[i]
  x <- read_exit(w$file, w$prefix)
  x[, `:=`(exit_year = w$exit_year, wave = w$wave, core_wave = w$core_wave)]
  x
}))
exit[, hhidpn := as.numeric(HHID) * 1000 + as.numeric(PN)]
if (anyDuplicated(exit$hhidpn)) stop("Duplicate hhidpn in Exit files.", call. = FALSE)

for (v in names(MISSING_CODES)) exit[get(v) %in% MISSING_CODES[[v]], (v) := NA]
exit[, `:=`(maenroll = fcase(ma_hmo == 1, 1L, ma_hmo == 5, 0L),
            lifesupt = fcase(life_support == 1, 1L, life_support == 5, 0L),
            oophosp  = oop_hospice,
            oopnurs  = oop_nursing,
            hospday  = hospice_days)]

# ---- RAND HRS: read only the needed columns ---------------------------------
rand_path <- find_raw(RAND_FILE_PATTERN)
if (is.null(rand_path)) stop("RAND HRS .dta file not found under data/raw/.", call. = FALSE)
rand_names <- names(haven::read_dta(rand_path, n_max = 0))

pick <- function(candidates) {
  hit <- intersect(candidates, rand_names)
  if (length(hit) == 0) NA_character_ else hit[1]
}
# Variable names follow RAND HRS 2022 (V2): sex is RASEX / SwSEX (older
# releases used *GENDER). HwITOT is total household income; HwIFTOT is only its
# imputation FLAG (0/1/2). The original code used HwIFTOT as "income", so it is
# kept as income_orig for the exact replication and replaced in 03_improved.R.
wave_vars <- function(w) c(
  sp_gender   = pick(paste0("s", w, c("sex", "gender"))),
  sp_race     = pick(paste0("s", w, "racem")),
  sp_edyrs    = pick(paste0("s", w, "edyrs")),
  income_orig = pick(paste0("h", w, "iftot")),
  income      = pick(paste0("h", w, "itot")),
  oopmd_r     = pick(paste0("r", w, "oopmd")),
  oopmd_s     = pick(paste0("s", w, "oopmd")))
fixed_vars <- c(resp_gender = pick(c("rasex", "ragender")), resp_race = pick("raracem"),
                resp_hisp = pick("rahispan"), resp_edyrs = pick("raedyrs"),
                renrsnit = pick("renrsnit"),
                rand_exit_wave = pick("rexitwv"),       # wave of the Exit interview
                hospice_stay   = pick("rehpcsty"),      # any hospice care before death
                hospice_nights = pick("rehpcnit"))
all_wave_vars <- c(wave_vars(13), wave_vars(14))
needed <- c("hhidpn", fixed_vars, all_wave_vars)
not_found <- c(names(fixed_vars)[is.na(fixed_vars)],
               paste0(names(all_wave_vars), " (wave ", rep(13:14, each = 7), ")")[is.na(all_wave_vars)])
if (length(not_found)) warning("Not found in RAND file (left missing): ",
                               paste(not_found, collapse = ", "), call. = FALSE)

message("Reading ", sum(!is.na(needed)), " columns from ", rand_path, " ...")
cols <- unname(needed[!is.na(needed)])
rand <- as.data.table(haven::zap_labels(haven::read_dta(rand_path, col_select = tidyselect::all_of(cols))))
rand <- rand[hhidpn %in% exit$hhidpn]

# Controls from a given core wave, renamed to generic names.
controls_from_wave <- function(w) {
  vars <- wave_vars(w)
  out <- data.table(hhidpn = rand$hhidpn)
  for (nm in names(vars)) out[, (nm) := if (is.na(vars[[nm]])) NA_real_ else rand[[vars[[nm]]]]]
  out
}
attach_controls <- function(rule) {
  if (rule == "aligned") {
    rbindlist(lapply(13:14, function(w)
      merge(exit[core_wave == w], controls_from_wave(w), by = "hhidpn", all.x = TRUE)))
  } else {
    merge(exit, controls_from_wave(as.integer(rule)), by = "hhidpn", all.x = TRUE)
  }
}

fixed <- data.table(hhidpn = rand$hhidpn)
for (nm in names(fixed_vars)) fixed[, (nm) := if (is.na(fixed_vars[[nm]])) NA_real_ else rand[[fixed_vars[[nm]]]]]

finalize <- function(dt) {
  dt <- merge(dt, fixed, by = "hhidpn", all.x = TRUE)
  dt[, `:=`(
    # Original-style controls (spouse variables). gender: 1 = male, 0 = female.
    gender_orig  = fcase(sp_gender == 1, 1L, sp_gender == 2, 0L),
    race_orig    = sp_race,
    eduyear_orig = sp_edyrs,
    # Respondent's own controls.
    gender_resp  = fcase(resp_gender == 1, 1L, resp_gender == 2, 0L),
    race_resp    = fcase(resp_hisp == 1, "Hispanic",
                         resp_race == 1, "White", resp_race == 2, "Black",
                         resp_race == 3, "Other"),
    eduyear_resp = resp_edyrs,
    oopmdre      = oopmd_r,
    oopmdsp      = oopmd_s)]
  dt
}

# ---- RAND Detailed Imputations: imputed hospice OOP cost --------------------
imp_path <- find_raw(RAND_IMP_PATTERN)
if (is.null(imp_path)) {
  warning("RAND HRS Detailed Imputations file not found: imputed hospice costs ",
          "will be missing and 03_improved.R falls back to the raw Exit item.", call. = FALSE)
  fixed[, `:=`(oop_hospice_rand = NA_real_, oop_hospice_flag = NA_real_)]
} else {
  imp_names <- names(haven::read_dta(imp_path, n_max = 0))
  imp_cols  <- intersect(c("hhidpn", "remahspc", "refahspc"), imp_names)
  if (!all(c("remahspc", "refahspc") %in% imp_cols)) {
    stop("REMAHSPC / REFAHSPC not found in ", imp_path, call. = FALSE)
  }
  message("Reading ", paste(imp_cols, collapse = ", "), " from ", imp_path, " ...")
  imp <- as.data.table(haven::zap_labels(
    haven::read_dta(imp_path, col_select = tidyselect::all_of(imp_cols))))
  imp <- imp[hhidpn %in% exit$hhidpn]
  setnames(imp, c("remahspc", "refahspc"), c("oop_hospice_rand", "oop_hospice_flag"))
  fixed <- merge(fixed, imp, by = "hhidpn", all.x = TRUE)
}

# ---- Sample-size diagnostics under each control-wave rule -------------------
base_orig <- c("oophosp", "maenroll", "gender_orig", "race_orig", "eduyear_orig",
               if (ORIG_INCOME == "itot") "income" else "income_orig")
base_resp <- c("oophosp", "maenroll", "gender_resp", "race_resp", "eduyear_resp", "income")
diag <- rbindlist(lapply(c("aligned", "13", "14"), function(rule) {
  dt <- finalize(attach_controls(rule))
  data.table(control_wave_rule = rule,
             exit_respondents = nrow(dt),
             with_oophosp_and_ma = sum(complete.cases(dt[, .(oophosp, maenroll)])),
             baseline_N_spouse_controls = sum(complete.cases(dt[, ..base_orig])),
             baseline_N_own_controls    = sum(complete.cases(dt[, ..base_resp])),
             has_spouse_data = sum(!is.na(dt$sp_gender)))
}))
save_csv(diag, "build_sample_sizes_by_rule")
print(diag)

analysis <- finalize(attach_controls(ORIG_CONTROL_RULE))
saveRDS(analysis, ANALYSIS_FILE)
saveRDS(finalize(attach_controls("aligned")), ANALYSIS_ALIGNED)

# ---- Raw Exit item vs. RAND imputed hospice OOP cost ------------------------
al <- readRDS(ANALYSIS_ALIGNED)
exit_wave_check <- al[, .(n = .N, rand_exit_wave_matches = sum(rand_exit_wave == wave, na.rm = TRUE)),
                      by = exit_year]
save_csv(exit_wave_check, "build_exit_wave_check")
flag_tab <- al[, .N, by = .(exit_year, flag = RAND_FLAG_LABELS[as.character(oop_hospice_flag)])]
flag_tab[is.na(flag), flag := "Missing in RAND"]
save_csv(dcast(flag_tab, flag ~ exit_year, value.var = "N", fill = 0), "build_rand_hospice_flags")
both <- al[!is.na(oophosp) & !is.na(oop_hospice_rand)]
compare <- data.table(
  raw_item_nonmissing   = sum(!is.na(al$oophosp)),
  rand_nonmissing       = sum(!is.na(al$oop_hospice_rand)),
  both_nonmissing       = nrow(both),
  identical_when_both   = sum(abs(both$oophosp - both$oop_hospice_rand) < 1),
  rand_with_ma_nonmiss  = sum(!is.na(al$oop_hospice_rand) & !is.na(al$maenroll)),
  raw_with_ma_nonmiss   = sum(!is.na(al$oophosp) & !is.na(al$maenroll)))
save_csv(compare, "build_raw_vs_rand_outcome")
print(exit_wave_check); print(dcast(flag_tab, flag ~ exit_year, value.var = "N", fill = 0)); print(compare)

# Item-level missingness, for the README / paper.
miss <- analysis[, lapply(.SD, function(v) sum(!is.na(v))),
                 .SDcols = c("maenroll", "oophosp", "oopnurs", "lifesupt", "hospday",
                             "gender_orig", "gender_resp", "race_resp", "eduyear_resp",
                             "income_orig", "income", "oopmdre", "oopmdsp", "renrsnit",
                             "oop_hospice_rand", "hospice_stay")]
save_csv(melt(miss, measure.vars = names(miss), variable.name = "variable",
              value.name = "non_missing"), "build_nonmissing_counts")
message("Saved ", ANALYSIS_FILE, " (", nrow(analysis), " Exit respondents; control rule = ",
        ORIG_CONTROL_RULE, ")")

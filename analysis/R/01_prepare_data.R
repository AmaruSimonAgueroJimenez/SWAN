# =============================================================================
# 01 - Read the SWAN visit 10 public-use file (ICPSR 32961) and build the
#      analysis data set. The raw file is never modified.
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
suppressPackageStartupMessages(library(readstata13))

# readstata13 warns about three symptom-count variables and one raw impedance
# variable whose value labels cannot be attached; none of them is used here.
raw <- withCallingHandlers(
  read.dta13(CFG$data_file),
  warning = function(w) {
    if (grepl("NUMHOTF10|IMPERAW10", conditionMessage(w))) invokeRestart("muffleWarning")
  }
)
stopifnot(nrow(raw) == 2245L, !anyDuplicated(raw$SWANID))

yes <- function(x) as.character(x) == "Yes"

# ---- recode -----------------------------------------------------------------
status_map <- c(
  "Pre-menopausal"                     = "Premenopause",
  "Early perimenopause"                = "Early perimenopause",
  "Late perimenopause"                 = "Late perimenopause",
  "Post-menopausal"                    = "Natural post-menopause",
  "Hysterectomy/both ovaries removed"  = "Surgical menopause (BSO)",
  "Unknown due to hormones (HT) use"   = "Undetermined (hormone therapy)",
  "Unknown due to hysterectomy"        = "Undetermined (hysterectomy)"
)
activity_map <- c(
  "No, not limited at all" = "Not limited",
  "Yes, limited a little"  = "Limited a little",
  "Yes, limited a lot"     = "Limited a lot"
)
race_map <- c("1" = "Black", "2" = "Chinese", "3" = "Japanese", "4" = "White")

status_chr   <- as.character(raw$STATUS10)
activity_chr <- as.character(raw$V_ACTI10)
exercise_chr <- as.character(raw$EXERCIS10)
stopifnot(all(na.omit(status_chr) %in% names(status_map)))
stopifnot(all(na.omit(activity_chr) %in% c(names(activity_map), "Missing")))
stopifnot(all(na.omit(exercise_chr) %in% c("No", "Yes")))
stopifnot(all(as.character(raw$RACE) %in% names(race_map)))

# exogenous sex-hormone use reported since the previous visit
hormone_vars <- c("ESTROG110", "ESTROG210", "ESTRNJ110", "ESTRNJ210", "COMBIN110",
                  "COMBIN210", "PROGES110", "PROGES210", "BCP110", "BCP210")
hormone_yes <- sapply(hormone_vars, function(v) yes(raw[[v]]))
hormone_any <- apply(hormone_yes, 1, function(r) any(r, na.rm = TRUE))
hormone_all_missing <- apply(hormone_yes, 1, function(r) all(is.na(r)))
hormone_any[hormone_all_missing] <- NA

dat <- data.frame(
  id            = raw$SWANID,
  shbg          = raw$SHBG10,       # nmol/L
  smm           = raw$SKELMM10,     # kg, Janssen (2000) bioimpedance equation
  testosterone  = raw$T10,          # ng/dL
  fsh           = raw$FSH10,        # mIU/mL
  estradiol     = raw$E2AVE10,      # pg/mL (mean of duplicate assays)
  tsh           = raw$TSH10,        # uIU/mL
  dheas         = raw$DHAS10,       # ug/dL
  age           = raw$AGE10,        # years
  pbf           = raw$PBFBIA10,     # % body fat (Sowers bioimpedance equation)
  ffm           = raw$FFMNHAN10,    # kg, fat-free mass (RJL/NHANES III equation)
  height_cm     = raw$HEIGHT10,
  weight_kg     = raw$WEIGHT10,
  waist_cm      = raw$WAIST10,
  resistance    = raw$CONDRAW10,    # ohms
  status        = factor(unname(status_map[status_chr]), levels = unname(status_map)),
  activity_limit = factor(unname(activity_map[activity_chr]), levels = unname(activity_map)),
  exercise      = factor(exercise_chr, levels = c("No", "Yes")),
  race          = factor(unname(race_map[as.character(raw$RACE)]), levels = c("White", "Black", "Chinese", "Japanese")),
  hormone_use   = hormone_any,
  e2_flag       = yes(raw$FLGCV10),
  stringsAsFactors = FALSE
)

# ---- derived variables ------------------------------------------------------
dat$bmi <- dat$weight_kg / (dat$height_cm / 100)^2             # kg/m2
# BMI supplied by SWAN must agree with the value computed from height and weight
both <- !is.na(raw$BMI10) & !is.na(dat$bmi)
stopifnot(max(abs(raw$BMI10[both] - dat$bmi[both])) < 1e-4)

dat$fat_mass_kg    <- dat$weight_kg * dat$pbf / 100            # kg
dat$smi            <- dat$smm / (dat$height_cm / 100)^2        # kg/m2
dat$smm_pct_weight <- 100 * dat$smm / dat$weight_kg            # % of body weight

# The SWAN skeletal muscle mass variable must be reproduced by the Janssen
# equation for women: SMM = 0.401 * height^2 / R - 0.071 * age + 5.102
janssen <- 0.401 * dat$height_cm^2 / dat$resistance - 0.071 * dat$age + 5.102
ok <- !is.na(janssen) & !is.na(dat$smm)
stopifnot(cor(janssen[ok], dat$smm[ok]) > 0.9999, abs(mean(janssen[ok] - dat$smm[ok])) < 0.01)

# reference categories
dat$status         <- relevel(dat$status, ref = CFG$ref_status)
dat$activity_limit <- relevel(dat$activity_limit, ref = CFG$ref_activity)
dat$exercise       <- relevel(dat$exercise, ref = CFG$ref_exercise)
dat$race           <- relevel(dat$race, ref = CFG$ref_race)

# ---- analytic sample --------------------------------------------------------
model_vars <- MODEL_VARIABLES
body_vars  <- c("weight_kg")          # weight, for the alternative measures of body size

flow <- data.frame(
  step = c("Assessed at visit 10",
           "SHBG measured",
           "Skeletal muscle mass estimated (bioimpedance)",
           "SHBG and skeletal muscle mass available",
           "Complete data on all model covariates (analytic sample)",
           "Analytic sample with weight (alternative measures of body size)"),
  n = c(nrow(dat),
        sum(!is.na(dat$shbg)),
        sum(!is.na(dat$smm)),
        sum(!is.na(dat$shbg) & !is.na(dat$smm)),
        sum(complete.cases(dat[, model_vars])),
        sum(complete.cases(dat[, c(model_vars, body_vars)])))
)
dat$analytic <- complete.cases(dat[, model_vars])
# percent body fat is needed only for the model of the first submission; requiring it must not
# exclude anyone who could be analysed with the primary model
stopifnot(identical(complete.cases(dat[, setdiff(model_vars, "pbf")]), dat$analytic))
analytic <- droplevels(dat[dat$analytic, ])

# every sensitivity model must be estimable on the full analytic sample
stopifnot(nrow(analytic) == CFG$expected_primary$n,
          flow$n[5] == flow$n[6])
stopifnot(all(analytic$smm > 0), all(analytic$shbg > 0))

dir.create(file.path(CFG$out_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CFG$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CFG$out_dir, "data"),    recursive = TRUE, showWarnings = FALSE)
saveRDS(dat,      file.path(CFG$out_dir, "data", "visit10_all.rds"))
saveRDS(analytic, file.path(CFG$out_dir, "data", "analytic.rds"))
write_table(flow, "flow_participants", CFG$out_dir)
print(flow)

# ---- missing data by variable (STROBE item 14b) -------------------------------
both <- !is.na(dat$shbg) & !is.na(dat$smm)
missing_tab <- data.frame(
  variable = model_vars,
  missing_visit10 = vapply(model_vars, function(v) sum(is.na(dat[[v]])), integer(1)),
  missing_with_shbg_and_smm = vapply(model_vars, function(v) sum(is.na(dat[[v]][both])), integer(1)),
  row.names = NULL, stringsAsFactors = FALSE
)
stopifnot(sum(both) - sum(complete.cases(dat[both, model_vars])) == flow$n[4] - flow$n[5])
write_table(missing_tab, "tableS_missing_by_variable", CFG$out_dir)
print(missing_tab)

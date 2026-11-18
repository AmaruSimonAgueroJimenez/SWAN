# =============================================================================
# 03 - Primary model: skeletal muscle mass ~ log2(SHBG) + covariates
#      Gamma family, log link. Covariates: age, height, waist circumference,
#      race/ethnicity, menopausal status, log2 of testosterone, estradiol, FSH,
#      DHEA-S and TSH, limitation of vigorous activities and exercise.
#      The model of the first submission (SKELMMM10) is also fitted here and kept
#      as a sensitivity analysis.
#   table2_primary_model.csv     ratios of means, 95% CI and P values
#   tableS_submitted_model.csv   the same for the model of the first submission
#   tableS_primary_joint_tests.csv  joint F tests of the terms of the primary model
#   tableS_predicted_smm.csv     adjusted mean skeletal muscle mass at selected SHBG values
#   tableS_diagnostics.csv       dispersion, collinearity, influence, link test
#   tableS_vif.csv               adjusted generalised variance inflation factors
#   tableS_diagnostics_initial.csv  the same diagnostics for the initial model
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
suppressPackageStartupMessages(library(car))
analytic <- readRDS(file.path(CFG$out_dir, "data", "analytic.rds"))
cl <- CFG$conf_level

f_primary <- make_formula("smm", PRIMARY_EXPOSURE, PRIMARY_COVARIATES)
m_primary <- fit_gamma_log(f_primary, analytic)

# ---- regression tests -----------------------------------------------------------
exp_ <- CFG$expected_primary
stopifnot(nobs(m_primary) == exp_$n,
          abs(AIC(m_primary) - exp_$aic) < 1e-3,
          abs(coef(m_primary)[["log2(shbg)"]] - exp_$beta_log2_shbg) < 1e-6)
message("Primary model reproduces the expected n, AIC and SHBG coefficient.")

f_submitted <- make_formula("smm", SUBMITTED_EXPOSURE, SUBMITTED_COVARIATES)
m_submitted <- fit_gamma_log(f_submitted, analytic)
exs <- CFG$expected_submitted
stopifnot(nobs(m_submitted) == exs$n,
          abs(AIC(m_submitted) - exs$aic) < 1e-3,
          abs(coef(m_submitted)[["shbg"]] - exs$beta_shbg) < 1e-6)
message("The model of the first submission is reproduced (n, AIC and SHBG coefficient).")

# ---- tables of ratios of means ----------------------------------------------------
cat_labels <- c(
  "raceBlack" = "  Black", "raceChinese" = "  Chinese", "raceJapanese" = "  Japanese",
  "statusPremenopause" = "  Premenopause",
  "statusEarly perimenopause" = "  Early perimenopause",
  "statusLate perimenopause" = "  Late perimenopause",
  "statusSurgical menopause (BSO)" = "  Surgical menopause (BSO)",
  "statusUndetermined (hormone therapy)" = "  Undetermined (hormone therapy)",
  "statusUndetermined (hysterectomy)" = "  Undetermined (hysterectomy)",
  "activity_limitLimited a little" = "  Limited a little",
  "activity_limitLimited a lot" = "  Limited a lot",
  "exerciseYes" = "  Yes")
labels_primary <- c(
  "log2(shbg)" = "SHBG, per doubling",
  "log2(testosterone)" = "Testosterone, per doubling",
  "log2(estradiol)" = "Estradiol, per doubling",
  "log2(fsh)" = "FSH, per doubling",
  "log2(dheas)" = "DHEA-S, per doubling",
  "log2(tsh)" = "TSH, per doubling",
  age = "Age, per 1 year",
  height_cm = "Height, per 5 cm",
  waist_cm = "Waist circumference, per 5 cm",
  cat_labels)
labels_submitted <- c(
  shbg = "SHBG, per 10 nmol/L",
  testosterone = "Testosterone, per 10 ng/dL",
  estradiol = "Estradiol, per 10 pg/mL",
  fsh = "FSH, per 10 mIU/mL",
  dheas = "DHEA-S, per 50 ug/dL",
  tsh = "TSH, per 1 uIU/mL",
  age = "Age, per 1 year",
  pbf = "Body fat, per 5 percentage points",
  cat_labels[!grepl("^race", names(cat_labels))])

# reporting unit of a term: 1 for log2 terms (a doubling) and for indicator variables
unit_of <- function(term) if (term %in% names(CFG$unit)) CFG$unit[[term]] else 1
ratio_table <- function(model, labels) {
  stopifnot(setequal(names(labels), setdiff(names(coef(model)), "(Intercept)")))
  tab <- do.call(rbind, lapply(names(labels), function(tm) ratio_for_term(model, tm, unit_of(tm), cl)))
  tab$label <- unname(labels[tab$term]); tab$pct_difference <- pct(tab$ratio)
  tab[, c("label", "term", "unit", "beta", "se", "ratio", "lcl", "ucl", "pct_difference", "p", "n")]
}
tab2 <- ratio_table(m_primary, labels_primary)
write_table(tab2, "table2_primary_model", CFG$out_dir)
tab_sub <- ratio_table(m_submitted, labels_submitted)
write_table(tab_sub, "tableS_submitted_model", CFG$out_dir)

# joint (multi-degree-of-freedom) tests of the terms of the primary model
joint <- car::Anova(m_primary, type = 2, test.statistic = "F")
joint_tab <- data.frame(term = rownames(joint), F = joint$F, df = joint$Df, df_residual = m_primary$df.residual, p = joint$`Pr(>F)`)
write_table(joint_tab, "tableS_primary_joint_tests", CFG$out_dir)

# ---- adjusted means at selected SHBG concentrations --------------------------
# Marginal standardisation: every participant is assigned the same SHBG value,
# predictions are averaged over the observed covariate distribution.
q <- quantile(analytic$shbg, c(0.10, 0.25, 0.50, 0.75, 0.90), names = FALSE)
pred_tab <- cbind(shbg_percentile = c(10, 25, 50, 75, 90), marginal_curve(m_primary, analytic, q, cl))
names(pred_tab) <- c("shbg_percentile", "shbg_nmol_L", "adjusted_mean_smm_kg", "lcl", "ucl")
pred_tab$difference_vs_p25_kg <- pred_tab$adjusted_mean_smm_kg - pred_tab$adjusted_mean_smm_kg[2]
write_table(pred_tab, "tableS_predicted_smm", CFG$out_dir)

# average change in skeletal muscle mass (kg) for a doubling of SHBG
ame_doubling <- mean(predict(m_primary, transform(analytic, shbg = 2 * shbg), type = "response") -
                     predict(m_primary, analytic, type = "response"))

# ---- diagnostics --------------------------------------------------------------
# The same quantities are computed for the primary and for the initial model.
diagnostics_of <- function(model, formula, data) {
  gv <- car::vif(model)                            # generalised VIF
  vif_adj <- gv[, "GVIF^(1/(2*Df))"]^2             # comparable to an ordinary VIF
  cd <- cooks.distance(model)
  # Pregibon link test: the squared linear predictor should add nothing
  eta <- predict(model, type = "link")
  m_link <- fit_gamma_log(update(formula, . ~ . + I(eta^2)), cbind(data, eta = eta))
  fitted_ <- fitted(model)
  list(vif_adj = vif_adj, cooks = cd,
       table = data.frame(
         quantity = c("Dispersion parameter", "Residual deviance", "Residual df",
                      "Deviance explained (1 - deviance/null deviance)",
                      "RMSE, kg", "Correlation observed vs fitted",
                      "Largest adjusted GVIF", "Term with largest adjusted GVIF",
                      "Observations with Cook's distance > 4/n", "Largest Cook's distance", "Largest leverage",
                      "Link test P (squared linear predictor)"),
         value = c(signif(summary(model)$dispersion, 4), round(deviance(model), 3), model$df.residual,
                   round(1 - deviance(model) / model$null.deviance, 4),
                   round(sqrt(mean((data$smm - fitted_)^2)), 3), round(cor(data$smm, fitted_), 3),
                   round(max(vif_adj), 2), names(which.max(vif_adj)),
                   sum(cd > 4 / nobs(model)), signif(max(cd), 3), signif(max(hatvalues(model)), 3),
                   signif(summary(m_link)$coefficients["I(eta^2)", "Pr(>|t|)"], 3))))
}
phi <- summary(m_primary)$dispersion
dg_primary <- diagnostics_of(m_primary, f_primary, analytic)
dg_initial <- diagnostics_of(m_submitted, f_submitted, analytic)
vif_adj <- dg_primary$vif_adj
cd <- dg_primary$cooks
diag_tab <- rbind(dg_primary$table,
                  data.frame(quantity = "Average change in SMM for a doubling of SHBG, kg", value = round(ame_doubling, 3)))
write_table(diag_tab, "tableS_diagnostics", CFG$out_dir)
write_table(dg_initial$table, "tableS_diagnostics_initial", CFG$out_dir)
write_table(data.frame(term = names(vif_adj), adjusted_gvif = round(vif_adj, 2)), "tableS_vif", CFG$out_dir)

saveRDS(list(model = m_primary, formula = f_primary, table2 = tab2, pred = pred_tab, ame_doubling = ame_doubling,
             joint = joint_tab, phi = phi, cooks = cd,
             submitted = m_submitted, formula_submitted = f_submitted, table_submitted = tab_sub,
             cooks_submitted = cooks.distance(m_submitted)),
        file.path(CFG$out_dir, "data", "primary_model.rds"))
print(tab2[, c("label", "ratio", "lcl", "ucl", "p")], digits = 4, row.names = FALSE)
print(joint_tab, digits = 4, row.names = FALSE)
print(pred_tab, digits = 4, row.names = FALSE)
print(diag_tab, row.names = FALSE)

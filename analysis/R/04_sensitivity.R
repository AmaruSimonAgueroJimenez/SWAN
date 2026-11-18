# =============================================================================
# 04 - Sensitivity analyses for the association between SHBG and skeletal
#      muscle mass. Every estimate is the ratio of means per doubling of SHBG
#      (coefficient of log2(SHBG)) from a gamma log-link model unless stated.
#   tableS_sensitivity.csv              all specifications (groups A to F)
#   tableS_nonlinearity.csv             shape of the association (log2, spline, linear term)
#   tableS_adjusted_means_by_model.csv  adjusted mean muscle mass across SHBG (curves of Fig 1)
#   tableS_full_regression_output.csv   coefficients of the primary and of the submitted model
#   tableS_fit_measures.csv             fit measures of the same two models
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
suppressPackageStartupMessages({library(sandwich); library(splines)})
analytic <- readRDS(file.path(CFG$out_dir, "data", "analytic.rds"))
pm       <- readRDS(file.path(CFG$out_dir, "data", "primary_model.rds"))
cl <- CFG$conf_level
EXPO <- PRIMARY_EXPOSURE                       # "log2(shbg)": ratios are per doubling of SHBG

est <- function(model, group, label, note = "", vcov_matrix = NULL) {
  r <- ratio_for_term(model, EXPO, 1, cl, vcov_matrix)
  cbind(data.frame(group = group, specification = label, note = note, stringsAsFactors = FALSE),
        r[, c("ratio", "lcl", "ucl", "p", "n")])
}
cov_p <- PRIMARY_COVARIATES
g <- function(outcome, covs, data = analytic) fit_gamma_log(make_formula(outcome, EXPO, covs), data)

res  <- list()
mods <- list()      # fitted models that are reused below

# ---- A. sequential adjustment ---------------------------------------------------
gA <- "A. Sequential adjustment"
res$A1 <- est(g("smm", character(0)), gA, "Unadjusted")
res$A2 <- est(g("smm", "age"), gA, "Adjusted for age")
res$A3 <- est(g("smm", c("age", "race", "height_cm")), gA, "Adjusted for age, race/ethnicity and height")
res$A4 <- est(g("smm", c("age", "race", "height_cm", "waist_cm")), gA, "Adjusted for age, race/ethnicity, height and waist circumference")
res$A5 <- est(pm$model, gA, "Primary model", "further adjusted for menopausal status, hormones, activity limitation and exercise")
mods$A5 <- pm$model

# ---- B. one covariate removed from the primary model ------------------------------
# Shows which covariates of the primary model the SHBG estimate depends on.
gB <- "B. One covariate removed from the primary model"
hormone_terms <- paste0("log2(", HORMONES, ")")
drop_sets <- list(
  B1 = list(vars = "waist_cm",       label = "waist circumference"),
  B2 = list(vars = "height_cm",      label = "height"),
  B3 = list(vars = "race",           label = "race/ethnicity"),
  B4 = list(vars = "age",            label = "age"),
  B5 = list(vars = "status",         label = "menopausal status"),
  B6 = list(vars = hormone_terms,    label = "the five hormone covariates"),
  B7 = list(vars = "activity_limit", label = "activity limitation"),
  B8 = list(vars = "exercise",       label = "exercise"))
stopifnot(setequal(unlist(lapply(drop_sets, `[[`, "vars")), cov_p))      # every covariate is removed once
for (id in names(drop_sets)) {
  mods[[id]] <- g("smm", setdiff(cov_p, drop_sets[[id]]$vars))
  res[[id]]  <- est(mods[[id]], gB, paste("Primary model without", drop_sets[[id]]$label))
}

# ---- C. other measures of body size ------------------------------------------------
gC <- "C. Other measures of body size"
no_size <- setdiff(cov_p, c("height_cm", "waist_cm"))
res$C1 <- est(g("smm", c(no_size, "bmi")), gC, "Body mass index in place of height and waist circumference")
res$C2 <- est(g("smm", c(no_size, "height_cm", "weight_kg")), gC, "Weight in place of waist circumference")
res$C3 <- est(g("smm", c(no_size, "height_cm", "fat_mass_kg")), gC, "Fat mass in place of waist circumference",
              "fat mass from the same impedance measurement as the outcome")
res$C4 <- est(g("smm", c(cov_p, "pbf")), gC, "Primary model + percent body fat")
res$C5 <- est(g("smm", c(cov_p, "bmi")), gC, "Primary model + body mass index")
res$C6 <- est(g("smm", c(cov_p, "fat_mass_kg")), gC, "Primary model + fat mass")
mods$C7 <- g("smm", SUBMITTED_COVARIATES)
res$C7 <- est(mods$C7, gC, "Covariates of the initial model",
              "hormones untransformed and percent body fat; no race/ethnicity, height or waist circumference")

# ---- D. alternative scaling of the outcome -------------------------------------------
gD <- "D. Alternative scaling of muscle mass"
res$D1 <- est(g("smi", setdiff(cov_p, "height_cm")), gD, "Skeletal muscle mass index (kg/m2)", "height is not a covariate")
res$D2 <- est(g("smm_pct_weight", cov_p), gD, "Skeletal muscle mass as % of body weight")

# ---- E. model specification ------------------------------------------------------------
gE <- "E. Model specification"
res$E1 <- est(pm$model, gE, "Primary model, robust (HC3) standard errors", vcov_matrix = sandwich::vcovHC(pm$model, type = "HC3"))
m_gauss <- glm(pm$formula, data = analytic, family = gaussian(link = "log"), start = coef(pm$model))
stopifnot(m_gauss$converged)
res$E2 <- est(m_gauss, gE, "Gaussian family, log link")
res$E3 <- est(lm(update(pm$formula, log(smm) ~ .), data = analytic), gE, "Linear regression of log skeletal muscle mass", "ratio of geometric means")
res$E4 <- est(g("smm", c(HORMONES, setdiff(cov_p, hormone_terms))), gE, "Hormone covariates untransformed")

# ---- F. sample restrictions ---------------------------------------------------------
gF <- "F. Sample restrictions"
keep_cook <- pm$cooks <= 4 / nrow(analytic)
res$F1 <- est(fit_gamma_log(pm$formula, analytic[keep_cook, ]), gF, "Excluding observations with Cook's distance > 4/n")
stopifnot(!anyNA(analytic$hormone_use))
no_ht <- !analytic$hormone_use & analytic$status != "Undetermined (hormone therapy)"
res$F2 <- est(fit_gamma_log(pm$formula, droplevels(analytic[no_ht, ])), gF, "Excluding users of exogenous sex hormones")
natural <- analytic$status %in% c("Premenopause", "Early perimenopause", "Late perimenopause", "Natural post-menopause")
res$F3 <- est(fit_gamma_log(pm$formula, droplevels(analytic[natural, ])), gF, "Natural menopausal stages only")
post <- analytic$status == "Natural post-menopause"
res$F4 <- est(fit_gamma_log(update(pm$formula, . ~ . - status), analytic[post, ]), gF, "Natural post-menopause only")
top_t <- which.max(analytic$testosterone)
stopifnot(analytic$testosterone[top_t] > 1500, which.max(pm$cooks_submitted) == top_t)   # the influential point of the first submission
res$F5 <- est(fit_gamma_log(pm$formula, analytic[-top_t, ]), gF, "Excluding the participant with the extreme testosterone value")

sens <- do.call(rbind, res)
sens$id <- names(res)
sens$pct_difference <- pct(sens$ratio)
stopifnot(all(sens$n[!grepl("^F", sens$id)] == nrow(analytic)))          # same women in every unrestricted model
write_table(sens[, c("id", "group", "specification", "note", "n", "ratio", "lcl", "ucl", "pct_difference", "p")],
            "tableS_sensitivity", CFG$out_dir)

# ---- shape of the association ----------------------------------------------------------
# Primary model: natural cubic spline of log2(SHBG) against the linear term in log2(SHBG),
# and SHBG as an untransformed linear term (reported per 10 nmol/L).
u10 <- CFG$unit[["shbg"]]
m_spline <- fit_gamma_log(update(pm$formula, . ~ . - log2(shbg) + ns(log2(shbg), df = 3)), analytic)
nl <- anova(pm$model, m_spline, test = "F")
m_linear <- fit_gamma_log(update(pm$formula, . ~ . - log2(shbg) + shbg), analytic)
r_lin <- ratio_for_term(m_linear, "shbg", u10, cl)
# Covariates of the model submitted originally: spline of SHBG against the linear term
m_spline_sub <- fit_gamma_log(update(pm$formula_submitted, . ~ . - shbg + ns(shbg, df = 3)), analytic)
nl_sub <- anova(pm$submitted, m_spline_sub, test = "F")
r_sub  <- ratio_for_term(pm$submitted, "shbg", u10, cl)
nonlin <- data.frame(
  model = c("Primary model", "Covariates of the initial model"),
  term_tested = c("log2(SHBG)", "SHBG, linear"),
  aic_log2 = c(AIC(pm$model), AIC(mods$C7)),
  aic_linear = c(AIC(m_linear), AIC(pm$submitted)),
  aic_spline_3df = c(AIC(m_spline), AIC(m_spline_sub)),
  F_nonlinearity = c(nl$F[2], nl_sub$F[2]), df1 = c(nl$Df[2], nl_sub$Df[2]), df2 = c(nl$`Resid. Df`[2], nl_sub$`Resid. Df`[2]),
  p_nonlinearity = c(nl$`Pr(>F)`[2], nl_sub$`Pr(>F)`[2]),
  ratio_per_10 = c(r_lin$ratio, r_sub$ratio), lcl_per_10 = c(r_lin$lcl, r_sub$lcl), ucl_per_10 = c(r_lin$ucl, r_sub$ucl),
  p_per_10 = c(r_lin$p, r_sub$p),
  stringsAsFactors = FALSE)
write_table(nonlin, "tableS_nonlinearity", CFG$out_dir)

# ---- adjusted mean skeletal muscle mass across SHBG (curves of Fig 1) ------------------
# Same standardisation for every model: all women are assigned the same SHBG and the
# predictions are averaged over the observed covariates of the analytic sample.
shbg_range <- quantile(analytic$shbg, c(0.025, 0.975), names = FALSE)
shbg_grid  <- seq(shbg_range[1], shbg_range[2], length.out = 60)
curve_models <- list(A5 = mods$A5, B1 = mods$B1, C7 = mods$C7)
curves <- do.call(rbind, lapply(names(curve_models), function(id)
  cbind(id = id, marginal_curve(curve_models[[id]], analytic, shbg_grid, cl), stringsAsFactors = FALSE)))
write_table(curves, "tableS_adjusted_means_by_model", CFG$out_dir)

# ---- full regression output (log scale) and fit measures of the two main models --------
write_table(rbind(full_output(pm$model, "Primary model"),
                  full_output(pm$submitted, "Initial model")), "tableS_full_regression_output", CFG$out_dir)
write_table(rbind(fit_measures(pm$model, "Primary model"),
                  fit_measures(pm$submitted, "Initial model")), "tableS_fit_measures", CFG$out_dir)

saveRDS(list(sens = sens, curves = curves, shbg_range = shbg_range, nonlin = nonlin,
             m_spline = m_spline, m_linear = m_linear, m_spline_sub = m_spline_sub),
        file.path(CFG$out_dir, "data", "sensitivity.rds"))
options(width = 200)
print(sens[, c("id", "specification", "n", "pct_difference", "ratio", "lcl", "ucl", "p")], digits = 4, row.names = FALSE)
print(nonlin, digits = 4)

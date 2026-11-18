# =============================================================================
# PROPOSAL of 17 November 2026, ADOPTED the same day as the primary model of the
# manuscript (see analysis/R/functions.R and analysis/R/03_primary_model.R).
# Kept as the record of how the model was arrived at; run_all.R does not call it.
# In this script "submitted" is the model of the first submission, which the
# pipeline now calls the initial model.
#
# Rationale, fixed before the models below were fitted:
#   1. SHBG and the hormone covariates on the log2 scale (right-skewed
#      concentrations; one extreme testosterone value has leverage 0.77 on the
#      original scale).
#   2. Race/ethnicity: related to SHBG and to body composition; a design variable
#      of SWAN.
#   3. Height: determinant of absolute muscle mass and a term of the impedance
#      equation that defines the outcome; not affected by SHBG.
#   4. Waist circumference: central adiposity, a common cause of low SHBG and of
#      larger absolute muscle mass, measured independently of the impedance device.
#   5. Percent body fat is dropped: it comes from the same impedance measurement
#      as the outcome and, together with any weight-based measure, it determines
#      fat-free mass, of which skeletal muscle is the main part.
#
# Run from the repository root, after analysis/run_all.R:
#     Rscript analysis/proposal/primary_model_proposal.R
# Output: analysis/proposal/output/*.csv
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
suppressPackageStartupMessages(library(splines))
analytic <- readRDS(file.path(CFG$out_dir, "data", "analytic.rds"))
pm       <- readRDS(file.path(CFG$out_dir, "data", "primary_model.rds"))
out_dir  <- file.path("analysis", "proposal", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cl <- CFG$conf_level
n  <- nrow(analytic)
stopifnot(n == CFG$expected_primary$n, !anyNA(analytic[, c("height_cm", "weight_kg", "waist_cm", "race", "bmi", "fat_mass_kg")]))

hormones <- c("testosterone", "fsh", "estradiol", "tsh", "dheas")
stopifnot(all(sapply(analytic[hormones], function(x) all(x > 0))), all(analytic$shbg > 0))
log2h <- paste0("log2(", hormones, ")")
other <- c("age", "status", "activity_limit", "exercise")

# ---- one row per model: SHBG as a linear term and as log2(SHBG) ----------------
describe <- function(label, covariates, data = analytic) {
  m_lin <- fit_gamma_log(make_formula("smm", "shbg", covariates), data)
  m_log <- fit_gamma_log(make_formula("smm", "log2(shbg)", covariates), data)
  r10 <- ratio_for_term(m_lin, "shbg", CFG$unit[["shbg"]], cl)
  rdb <- ratio_for_term(m_log, "log2(shbg)", 1, cl)
  data.frame(
    model = label, n = nobs(m_log), coefficients = length(coef(m_log)),
    pct_per_10 = 100 * (r10$ratio - 1), pct_per_10_lcl = 100 * (r10$lcl - 1), pct_per_10_ucl = 100 * (r10$ucl - 1), p_per_10 = r10$p,
    pct_per_doubling = 100 * (rdb$ratio - 1), pct_doubling_lcl = 100 * (rdb$lcl - 1), pct_doubling_ucl = 100 * (rdb$ucl - 1), p_doubling = rdb$p,
    aic_linear = AIC(m_lin), aic_log2 = AIC(m_log),
    dispersion_log2 = summary(m_log)$dispersion,
    deviance_explained_log2 = 1 - deviance(m_log) / m_log$null.deviance,
    stringsAsFactors = FALSE)
}

submitted <- c(hormones, "age", "pbf", other[-1])
stopifnot(setequal(submitted, SUBMITTED_COVARIATES))
proposed  <- c(log2h, other, "race", "height_cm", "waist_cm")
stopifnot(setequal(proposed, PRIMARY_COVARIATES))      # the proposal is the primary model of the pipeline

path <- rbind(
  describe("0. Submitted model",                    submitted),
  describe("1. Hormone covariates on the log2 scale", c(log2h, "pbf", other)),
  describe("2. + race/ethnicity",                   c(log2h, "pbf", other, "race")),
  describe("3. + height",                           c(log2h, "pbf", other, "race", "height_cm")),
  describe("4. + waist circumference",              c(log2h, "pbf", other, "race", "height_cm", "waist_cm")),
  describe("5. - percent body fat (proposed)",      proposed)
)
# the first row must reproduce the submitted model
stopifnot(abs(path$aic_linear[1] - CFG$expected_submitted$aic) < 1e-3)
# the last row must reproduce the primary model of the pipeline
stopifnot(abs(path$aic_log2[nrow(path)] - CFG$expected_primary$aic) < 1e-3)
write_table <- function(df, name) { p <- file.path(out_dir, paste0(name, ".csv")); write.csv(df, p, row.names = FALSE, fileEncoding = "UTF-8"); message("  wrote ", p) }
write_table(path, "path_to_proposed_model")

base <- c(log2h, other, "race")
alternatives <- rbind(
  describe("Proposed: height + waist circumference",            proposed),
  describe("Body mass index instead of height and waist",       c(base, "bmi")),
  describe("Height + weight instead of height and waist",       c(base, "height_cm", "weight_kg")),
  describe("Height + fat mass instead of height and waist",     c(base, "height_cm", "fat_mass_kg")),
  describe("Height only (no adiposity measure)",                c(base, "height_cm")),
  describe("Proposed + percent body fat",                       c(proposed, "pbf")),
  describe("Proposed + fat mass",                               c(proposed, "fat_mass_kg")),
  describe("Proposed + body mass index",                        c(proposed, "bmi")),
  describe("Proposed, hormone covariates untransformed",        c(hormones, other, "race", "height_cm", "waist_cm")),
  describe("Proposed, without hormone covariates",              c(other, "race", "height_cm", "waist_cm")),
  describe("Proposed, log(height) and log(waist)",              c(log2h, other, "race", "log(height_cm)", "log(waist_cm)")),
  describe("Proposed, excluding the extreme testosterone value", proposed, analytic[-which.max(analytic$testosterone), ]),
  describe("Proposed, natural post-menopause only",             setdiff(proposed, "status"), droplevels(analytic[analytic$status == "Natural post-menopause", ]))
)
write_table(alternatives, "alternatives_to_proposed_model")

# ---- the proposed model in detail -------------------------------------------------
f_prop <- make_formula("smm", "log2(shbg)", proposed)
m_prop <- fit_gamma_log(f_prop, analytic)
co <- summary(m_prop)$coefficients
crit <- qt(1 - (1 - cl) / 2, m_prop$df.residual)
coef_tab <- data.frame(term = rownames(co), beta = co[, 1], se = co[, 2], t = co[, 3], p = co[, 4],
                       ratio = exp(co[, 1]), lcl = exp(co[, 1] - crit * co[, 2]), ucl = exp(co[, 1] + crit * co[, 2]), row.names = NULL)
write_table(coef_tab, "proposed_model_coefficients")

# shape of the association: spline of log2(SHBG) against the linear term in log2(SHBG)
m_spl <- fit_gamma_log(update(f_prop, . ~ . - log2(shbg) + ns(log2(shbg), df = 3)), analytic)
nl <- anova(m_prop, m_spl, test = "F")
# diagnostics
cd <- cooks.distance(m_prop); gv <- car::vif(m_prop); vif_adj <- gv[, "GVIF^(1/(2*Df))"]^2
eta <- predict(m_prop, type = "link")
m_link <- fit_gamma_log(update(f_prop, . ~ . + I(eta^2)), cbind(analytic, eta = eta))
park_d <- data.frame(r2 = (analytic$smm - fitted(m_prop))^2, log_mu = log(fitted(m_prop)))
m_park <- glm(r2 ~ log_mu, data = park_d, family = Gamma(link = "log")); stopifnot(m_park$converged)
se_park <- sqrt(sandwich::vcovHC(m_park, type = "HC3")["log_mu", "log_mu"])
# other families with the proposed linear predictor
m_ig <- glm(f_prop, data = analytic, family = inverse.gaussian(link = "log"), start = coef(m_prop)); stopifnot(m_ig$converged)
m_ln <- lm(update(f_prop, log(smm) ~ .), data = analytic)
dbl <- function(m) 100 * (exp(coef(m)[["log2(shbg)"]]) - 1)
# average difference in kg for a doubling of SHBG
ame_dbl <- mean(predict(m_prop, transform(analytic, shbg = 2 * shbg), type = "response") - fitted(m_prop))

diag <- data.frame(
  quantity = c("n", "coefficients", "residual df", "dispersion", "deviance explained", "AIC",
               "non-linearity of log2(SHBG): F", "non-linearity: df1", "non-linearity: df2", "non-linearity: P",
               "largest Cook's distance", "largest leverage", "observations with Cook's distance > 4/n",
               "largest adjusted GVIF", "term with largest adjusted GVIF", "link test P",
               "Park test lambda", "Park test lambda, lower 95% limit", "Park test lambda, upper 95% limit",
               "inverse Gaussian: % per doubling", "inverse Gaussian: AIC", "log-normal: % per doubling", "log-normal: AIC (scale of muscle mass)",
               "average difference in kg per doubling of SHBG", "Spearman SHBG - waist", "Spearman SHBG - height", "Spearman SMM - waist", "Spearman SMM - height",
               "Pearson SMM - fat-free mass from weight and percent body fat"),
  value = c(nobs(m_prop), length(coef(m_prop)), m_prop$df.residual, summary(m_prop)$dispersion, 1 - deviance(m_prop) / m_prop$null.deviance, AIC(m_prop),
            nl$F[2], nl$Df[2], nl$`Resid. Df`[2], nl$`Pr(>F)`[2],
            max(cd), max(hatvalues(m_prop)), sum(cd > 4 / n),
            max(vif_adj), names(which.max(vif_adj)), summary(m_link)$coefficients["I(eta^2)", "Pr(>|t|)"],
            coef(m_park)[["log_mu"]], coef(m_park)[["log_mu"]] - qnorm(0.975) * se_park, coef(m_park)[["log_mu"]] + qnorm(0.975) * se_park,
            dbl(m_ig), AIC(m_ig), dbl(m_ln), AIC(m_ln) + 2 * sum(log(analytic$smm)),
            ame_dbl, cor(analytic$shbg, analytic$waist_cm, method = "spearman"), cor(analytic$shbg, analytic$height_cm, method = "spearman"),
            cor(analytic$smm, analytic$waist_cm, method = "spearman"), cor(analytic$smm, analytic$height_cm, method = "spearman"),
            cor(analytic$smm, analytic$weight_kg * (1 - analytic$pbf / 100))),
  stringsAsFactors = FALSE)
write_table(diag, "proposed_model_diagnostics")

options(width = 200)
print(path[, c("model", "coefficients", "pct_per_10", "pct_per_10_lcl", "pct_per_10_ucl", "p_per_10", "pct_per_doubling", "pct_doubling_lcl", "pct_doubling_ucl", "p_doubling", "aic_linear", "aic_log2", "deviance_explained_log2")], digits = 3, row.names = FALSE)
print(alternatives[, c("model", "n", "pct_per_doubling", "pct_doubling_lcl", "pct_doubling_ucl", "p_doubling", "pct_per_10", "p_per_10", "aic_log2", "deviance_explained_log2")], digits = 3, row.names = FALSE)
print(coef_tab, digits = 3, row.names = FALSE)
print(diag, row.names = FALSE)

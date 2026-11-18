# =============================================================================
# 06 - Complementary model with SHBG as the dependent variable (gamma, log link).
#      (a) covariates of the primary model; (b) covariates of the model of the
#      first submission (SHBGM18), kept as a regression test and for comparison.
#   tableS_shbg_outcome.csv               ratios of mean SHBG, both sets of covariates
#   tableS_shbg_outcome_joint_tests.csv   joint F tests of the terms of model (a)
#   tableS_shbg_outcome_full_output.csv   coefficients, standard errors and t statistics of model (a)
#   tableS_shbg_outcome_fit.csv           fit measures of model (a)
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
suppressPackageStartupMessages(library(car))
analytic <- readRDS(file.path(CFG$out_dir, "data", "analytic.rds"))
cl <- CFG$conf_level

m_a <- fit_gamma_log(make_formula("shbg", "smm", PRIMARY_COVARIATES), analytic)
m_b <- fit_gamma_log(make_formula("shbg", "smm", SUBMITTED_COVARIATES), analytic)

# regression test against Table 2 of the first submission (SHBGM18)
stopifnot(abs(coef(m_b)["smm"] - (-0.032982)) < 1e-5, abs(AIC(m_b) - 15164.54) < 1e-2)

unit_of <- function(term) if (term %in% names(CFG$unit)) CFG$unit[[term]] else 1
tidy_all <- function(m, label) {
  tm <- setdiff(names(coef(m)), "(Intercept)")
  out <- do.call(rbind, lapply(tm, function(t) ratio_for_term(m, t, unit_of(t), cl)))
  cbind(model = label, out)
}
tab <- rbind(tidy_all(m_a, "Covariates of the primary model"), tidy_all(m_b, "Covariates of the initial model"))
tab$pct_difference <- pct(tab$ratio)
write_table(tab, "tableS_shbg_outcome", CFG$out_dir)

# Model (a) is Table 3 of the article: joint tests, full regression output and fit measures, as for the primary model
COMPLEMENTARY <- "Complementary model"
joint <- car::Anova(m_a, type = 2, test.statistic = "F")
joint_tab <- data.frame(term = rownames(joint), F = joint$F, df = joint$Df, df_residual = m_a$df.residual, p = joint$`Pr(>F)`)
joint_tab <- joint_tab[joint_tab$term != "Residuals", ]
stopifnot(setequal(joint_tab$term, c("smm", PRIMARY_COVARIATES)), !anyNA(joint_tab$F))
write_table(joint_tab, "tableS_shbg_outcome_joint_tests", CFG$out_dir)
write_table(full_output(m_a, COMPLEMENTARY), "tableS_shbg_outcome_full_output", CFG$out_dir)
fit_a <- fit_measures(m_a, COMPLEMENTARY)
write_table(fit_a, "tableS_shbg_outcome_fit", CFG$out_dir)

saveRDS(list(m_a = m_a, m_b = m_b, tab = tab, joint = joint_tab, fit = fit_a), file.path(CFG$out_dir, "data", "shbg_outcome.rds"))
print(joint_tab, digits = 4, row.names = FALSE)
print(tab[, c("model", "term", "unit", "ratio", "lcl", "ucl", "p")], digits = 4, row.names = FALSE)

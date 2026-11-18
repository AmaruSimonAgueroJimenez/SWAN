# =============================================================================
# 05 - Alternative error distributions and link functions for skeletal muscle
#      mass, all fitted to the same analytic sample with the primary covariates.
#      The SHBG association is expressed on a common scale: the average change
#      in skeletal muscle mass (kg) for a doubling of SHBG.
#   tableS_model_comparison.csv
#      Two further checks of the choice of a gamma model:
#   tableS_alternative_families.csv  inverse Gaussian and log-normal models
#   tableS_variance_function.csv     modified Park test of the variance function
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
analytic <- readRDS(file.path(CFG$out_dir, "data", "analytic.rds"))
pm       <- readRDS(file.path(CFG$out_dir, "data", "primary_model.rds"))
f <- pm$formula

# starting values: inverse-link fits need a feasible start (mu > 0)
start_inverse <- c(1 / mean(analytic$smm), rep(0, length(coef(pm$model)) - 1))
specs <- list(
  list(label = "Gamma, log link (primary)", family = Gamma(link = "log"),          start = NULL),
  list(label = "Gamma, identity link",      family = Gamma(link = "identity"),     start = NULL),
  list(label = "Gamma, inverse link",       family = Gamma(link = "inverse"),      start = start_inverse),
  list(label = "Gaussian, identity link",   family = gaussian(link = "identity"),  start = NULL),
  list(label = "Gaussian, log link",        family = gaussian(link = "log"),       start = coef(pm$model)),
  list(label = "Gaussian, inverse link",    family = gaussian(link = "inverse"),   start = start_inverse)
)

# average change in skeletal muscle mass (kg) when SHBG is multiplied by `factor`
ame_shbg <- function(model, data, factor) {
  mean(predict(model, transform(data, shbg = shbg * factor), type = "response") -
       predict(model, data, type = "response"))
}
EXPO <- PRIMARY_EXPOSURE

rows <- lapply(specs, function(s) {
  m <- glm(f, data = analytic, family = s$family, start = s$start)
  stopifnot(m$converged, nobs(m) == nrow(analytic))
  co <- summary(m)$coefficients[EXPO, ]
  data.frame(
    model = s$label, n = nobs(m),
    aic = AIC(m), bic = BIC(m),
    rmse_kg = sqrt(mean((analytic$smm - fitted(m))^2)),
    shbg_coefficient = co[["Estimate"]], shbg_se = co[["Std. Error"]], shbg_p = co[[4]],
    change_smm_kg_per_doubling = ame_shbg(m, analytic, 2),
    stringsAsFactors = FALSE)
})
cmp <- do.call(rbind, rows)
cmp$delta_aic <- cmp$aic - min(cmp$aic)
cmp <- cmp[order(cmp$aic), ]
write_table(cmp, "tableS_model_comparison", CFG$out_dir)
saveRDS(cmp, file.path(CFG$out_dir, "data", "model_comparison.rds"))
print(cmp, digits = 5, row.names = FALSE)

# ---- further families with a log-scale mean model ------------------------------
# Same sample and linear predictor. AIC is always on the scale of skeletal muscle
# mass (kg): for the log-normal model the Jacobian term 2 * sum(log(y)) is added to
# the AIC of the regression of log(y).
m_gamma  <- pm$model
m_invg   <- glm(f, data = analytic, family = inverse.gaussian(link = "log"), start = coef(pm$model))
m_lognor <- lm(update(f, log(smm) ~ .), data = analytic)
stopifnot(m_invg$converged, nobs(m_invg) == nrow(analytic), nobs(m_lognor) == nrow(analytic))

ci_pct <- function(model) {            # percent difference per doubling of SHBG, Wald t interval
  b <- coef(model)[[EXPO]]; se <- sqrt(vcov(model)[EXPO, EXPO])
  crit <- qt(1 - (1 - CFG$conf_level) / 2, df.residual(model))
  100 * (exp(c(b, b - crit * se, b + crit * se)) - 1)
}
# average change in kg; the log-normal mean is exp(xb + sigma^2 / 2)
ame_lognormal <- function(model, data, factor) {
  s2 <- summary(model)$sigma^2
  mean(exp(predict(model, transform(data, shbg = shbg * factor)) + s2 / 2) -
       exp(predict(model, data) + s2 / 2))
}
alt <- data.frame(
  model = c("Gamma, log link (primary)", "Inverse Gaussian, log link", "Log-normal (linear regression of log muscle mass)"),
  variance_function = c("proportional to mean^2", "proportional to mean^3", "proportional to mean^2 (constant on the log scale)"),
  estimand = c("ratio of arithmetic means", "ratio of arithmetic means", "ratio of geometric means"),
  n = c(nobs(m_gamma), nobs(m_invg), nobs(m_lognor)),
  aic = c(AIC(m_gamma), AIC(m_invg), AIC(m_lognor) + 2 * sum(log(analytic$smm))),
  rbind(ci_pct(m_gamma), ci_pct(m_invg), ci_pct(m_lognor)),
  change_smm_kg_per_doubling = c(ame_shbg(m_gamma, analytic, 2), ame_shbg(m_invg, analytic, 2),
                                 ame_lognormal(m_lognor, analytic, 2)),
  stringsAsFactors = FALSE
)
names(alt)[6:8] <- c("pct_difference", "pct_lcl", "pct_ucl")
alt$delta_aic_vs_primary <- alt$aic - alt$aic[1]
stopifnot(abs(alt$aic[1] - CFG$expected_primary$aic) < 1e-3)
write_table(alt, "tableS_alternative_families", CFG$out_dir)
print(alt, digits = 5, row.names = FALSE)

# ---- variance function: modified Park test (Manning and Mullahy, 2001) ---------
# Regress the squared raw residuals of the primary model on the log of the fitted
# mean with a gamma log-link model. The slope lambda estimates the power of the
# variance function, Var(Y) proportional to mean^lambda: 0 Gaussian, 1 Poisson-like,
# 2 gamma, 3 inverse Gaussian. Robust (HC3) standard error.
park_dat <- data.frame(r2 = (analytic$smm - fitted(m_gamma))^2, log_mu = log(fitted(m_gamma)))
m_park   <- glm(r2 ~ log_mu, data = park_dat, family = Gamma(link = "log"))
stopifnot(m_park$converged)
se_park  <- sqrt(sandwich::vcovHC(m_park, type = "HC3")["log_mu", "log_mu"])
lambda   <- coef(m_park)[["log_mu"]]
zcrit    <- qnorm(1 - (1 - CFG$conf_level) / 2)
park <- data.frame(
  quantity = c("lambda (power of the variance function)", "lower 95% limit", "upper 95% limit",
               "P value, lambda = 2 (gamma)", "P value, lambda = 0 (Gaussian)", "P value, lambda = 3 (inverse Gaussian)"),
  value = c(lambda, lambda - zcrit * se_park, lambda + zcrit * se_park,
            2 * pnorm(-abs((lambda - 2) / se_park)), 2 * pnorm(-abs(lambda / se_park)),
            2 * pnorm(-abs((lambda - 3) / se_park))),
  stringsAsFactors = FALSE
)
write_table(park, "tableS_variance_function", CFG$out_dir)
print(park, digits = 4, row.names = FALSE)
saveRDS(list(alt = alt, park = park), file.path(CFG$out_dir, "data", "model_families.rds"))

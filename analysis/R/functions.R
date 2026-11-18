# =============================================================================
# Helper functions shared by every analysis script
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ---- model fitting ----------------------------------------------------------
# All ratio-scale models go through this single function so that the family and
# link are guaranteed to be identical across the primary and sensitivity models.
fit_gamma_log <- function(formula, data) {
  m <- glm(formula, data = data, family = Gamma(link = "log"))
  stopifnot(m$converged)
  m
}

# Hormone covariates (all strictly positive concentrations).
HORMONES <- c("testosterone", "fsh", "estradiol", "tsh", "dheas")

# Primary model (revised in November 2026). SHBG and the hormone covariates enter
# on the log2 scale; body size is represented by height and waist circumference.
PRIMARY_EXPOSURE   <- "log2(shbg)"
PRIMARY_COVARIATES <- c(paste0("log2(", HORMONES, ")"), "age", "height_cm", "waist_cm", "race",
                        "status", "activity_limit", "exercise")

# Initial model: the model of the first submission (SKELMMM10; only the reference level of
# menopausal status differs, which does not change the fit). Kept as a sensitivity analysis
# and called "initial model" in every table and figure.
SUBMITTED_EXPOSURE   <- "shbg"
SUBMITTED_COVARIATES <- c(HORMONES, "age", "pbf", "status", "activity_limit", "exercise")

# Variables that must be complete for a woman to be analysed: those of the primary
# model and of the submitted model, so that every model uses the same sample.
MODEL_VARIABLES <- c("smm", "shbg", HORMONES, "age", "height_cm", "waist_cm", "race",
                     "status", "activity_limit", "exercise", "pbf")

make_formula <- function(outcome, exposure, covariates) {
  reformulate(c(exposure, covariates), response = outcome)
}

# ---- extraction -------------------------------------------------------------
# Ratio of means (exp(beta * unit)) with a Wald interval based on the
# t distribution with the residual degrees of freedom, i.e. the same reference
# distribution that summary.glm() uses for families with estimated dispersion.
ratio_for_term <- function(model, term, unit, conf_level, vcov_matrix = NULL) {
  stopifnot(term %in% names(coef(model)))
  beta <- unname(coef(model)[term])
  V    <- if (is.null(vcov_matrix)) vcov(model) else vcov_matrix
  se   <- sqrt(V[term, term])
  df   <- model$df.residual
  crit <- qt(1 - (1 - conf_level) / 2, df)
  tval <- beta / se
  data.frame(
    term     = term,
    unit     = unit,
    beta     = beta,
    se       = se,
    ratio    = exp(unit * beta),
    lcl      = exp(unit * (beta - crit * se)),
    ucl      = exp(unit * (beta + crit * se)),
    p        = 2 * pt(abs(tval), df, lower.tail = FALSE),
    n        = nobs(model),
    stringsAsFactors = FALSE
  )
}

# Marginally standardised mean of a log-link model across values of SHBG, with
# delta-method confidence limits. Every participant is assigned SHBG = x and the
# predictions are averaged over the observed covariates:
#   m(x) = mean_i exp(X_i(x) b);   dm/db = mean_i exp(X_i(x) b) X_i(x)
marginal_curve <- function(model, data, grid, conf_level) {
  stopifnot(family(model)$link == "log", "shbg" %in% all.vars(formula(model)))
  tt <- delete.response(terms(model)); V <- vcov(model); b <- coef(model)
  crit <- qt(1 - (1 - conf_level) / 2, model$df.residual)
  do.call(rbind, lapply(grid, function(x) {
    nd <- data; nd$shbg <- x
    X  <- model.matrix(tt, nd, xlev = model$xlevels)
    mu <- as.vector(exp(X %*% b))
    g  <- colMeans(mu * X)
    se <- sqrt(as.numeric(t(g) %*% V %*% g))
    data.frame(shbg = x, mean = mean(mu), lcl = mean(mu) - crit * se, ucl = mean(mu) + crit * se)
  }))
}

# Full regression output on the scale of the linear predictor
full_output <- function(model, label) {
  co <- summary(model)$coefficients
  data.frame(model = label, term = rownames(co), beta = co[, "Estimate"], se = co[, "Std. Error"],
             t = co[, "t value"], p = co[, "Pr(>|t|)"], row.names = NULL, stringsAsFactors = FALSE)
}
fit_measures <- function(model, label) {
  data.frame(model = label, n = nobs(model), residual_df = model$df.residual,
             dispersion = summary(model)$dispersion, null_deviance = model$null.deviance,
             residual_deviance = deviance(model), log_likelihood = as.numeric(logLik(model)),
             aic = AIC(model), bic = BIC(model), stringsAsFactors = FALSE)
}

# percent difference implied by a ratio
pct <- function(ratio) 100 * (ratio - 1)

# ---- formatting -------------------------------------------------------------
fmt_num <- function(x, digits) formatC(x, format = "f", digits = digits, big.mark = ",")

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}

fmt_ci <- function(est, lcl, ucl, digits) {
  sprintf("%s (%s to %s)", fmt_num(est, digits), fmt_num(lcl, digits), fmt_num(ucl, digits))
}

mean_sd <- function(x, digits = 1) {
  stopifnot(!anyNA(x))
  sprintf("%s (%s)", fmt_num(mean(x), digits), fmt_num(sd(x), digits))
}

median_iqr <- function(x, digits = 1) {
  stopifnot(!anyNA(x))
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
  sprintf("%s (%s to %s)", fmt_num(q[2], digits), fmt_num(q[1], digits), fmt_num(q[3], digits))
}

n_pct <- function(x, level, digits = 1) {
  stopifnot(!anyNA(x))
  k <- sum(x == level)
  sprintf("%s (%s)", fmt_num(k, 0), fmt_num(100 * k / length(x), digits))
}

# ---- output -----------------------------------------------------------------
write_table <- function(df, name, out_dir) {
  path <- file.path(out_dir, "tables", paste0(name, ".csv"))
  write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8")
  message("  wrote ", path)
  invisible(path)
}

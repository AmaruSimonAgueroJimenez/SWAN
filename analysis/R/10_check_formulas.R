# =============================================================================
# 10 - Numerical check of the formulas written in Text A of S1 File.
#      Each equation is computed by hand from the design matrix, the outcome and
#      the fitted values of the primary model and compared with what R returns.
#      The pipeline stops if any of them disagrees.
#   tableS_formula_checks.csv
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
suppressPackageStartupMessages({library(sandwich); library(splines)})
analytic <- readRDS(file.path(CFG$out_dir, "data", "analytic.rds"))
pm <- readRDS(file.path(CFG$out_dir, "data", "primary_model.rds"))
m  <- pm$model
X  <- model.matrix(m); y <- analytic$smm; mu <- unname(fitted(m)); n <- nrow(X); k <- ncol(X)
stopifnot(n == CFG$expected_primary$n, length(y) == n, all(weights(m, "working") == 1))   # working weights are 1 (gamma, log link)

checks <- list()
check <- function(equation, quantity, by_hand, from_r, tol) {
  d <- max(abs(by_hand - from_r))
  checks[[length(checks) + 1]] <<- data.frame(equation = equation, quantity = quantity, max_abs_difference = d,
                                              tolerance = tol, agrees = d < tol, stringsAsFactors = FALSE)
}
XtXi <- solve(crossprod(X))
r    <- (y - mu) / mu                                   # Pearson residual of the gamma model
phi  <- sum(r^2) / (n - k)
V    <- phi * XtXi
h    <- rowSums((X %*% XtXi) * X)                       # diagonal of X (X'X)^-1 X'

# S5: the score equations hold at the estimate (up to the convergence tolerance of glm)
check("S5", "score equations, relative to sum(|r x_j|)", drop(crossprod(X, r)) / colSums(abs(r * X)), 0, 1e-4)
check("S6", "dispersion from the Pearson statistic", phi, summary(m)$dispersion, 1e-12)
check("S6", "covariance matrix of the coefficients", V, vcov(m), 1e-10)

# S9: confidence interval of the ratio of means per doubling of SHBG
b1 <- coef(m)[["log2(shbg)"]]; se1 <- sqrt(V["log2(shbg)", "log2(shbg)"]); tq <- qt(0.975, n - k)
t2 <- pm$table2[pm$table2$term == "log2(shbg)", ]
check("S7, S9", "ratio of means per doubling and its 95% CI", exp(c(b1, b1 - tq * se1, b1 + tq * se1)), c(t2$ratio, t2$lcl, t2$ucl), 1e-10)

# S10: joint F test of race/ethnicity
m0 <- fit_gamma_log(update(pm$formula, . ~ . - race), analytic)
check("S10", "joint F test, race/ethnicity", ((deviance(m0) - deviance(m)) / 3) / phi, pm$joint$F[pm$joint$term == "race"], 1e-8)

# S11, S12: standardized mean at the median SHBG and its delta-method interval
s  <- median(analytic$shbg); Xs <- X; Xs[, "log2(shbg)"] <- log2(s)
mus <- exp(drop(Xs %*% coef(m))); g <- colMeans(mus * Xs); se_s <- sqrt(drop(t(g) %*% V %*% g))
mc <- marginal_curve(m, analytic, s, 0.95)
check("S11, S12", "standardized mean at the median SHBG and its 95% CI",
      c(mean(mus), mean(mus) - tq * se_s, mean(mus) + tq * se_s), c(mc$mean, mc$lcl, mc$ucl), 1e-8)

# S13: average difference in kg for a doubling of SHBG
check("S13", "average difference in kg per doubling", (exp(b1) - 1) * mean(mu), pm$ame_doubling, 1e-10)

# S16: HC3 covariance matrix
HC3 <- XtXi %*% crossprod(X * (r / (1 - h))) %*% XtXi
check("S16", "HC3 covariance matrix", HC3, sandwich::vcovHC(m, type = "HC3"), 1e-10)

# S19: F test of the spline against the single term
ms <- fit_gamma_log(update(pm$formula, . ~ . - log2(shbg) + ns(log2(shbg), df = 3)), analytic)
nl <- read.csv(file.path(CFG$out_dir, "tables", "tableS_nonlinearity.csv"))
check("S19", "F test of the spline of log2(SHBG)", ((deviance(m) - deviance(ms)) / 2) / summary(ms)$dispersion, nl$F_nonlinearity[1], 1e-8)

# S22: AIC from the log-likelihood of S4. R evaluates it at the dispersion deviance / n.
ph <- deviance(m) / n
ll <- sum(-lgamma(1 / ph) - (1 / ph) * log(ph * mu) + (1 / ph - 1) * log(y) - y / (ph * mu))
check("S4, S22", "AIC from the gamma log-likelihood", -2 * ll + 2 * (k + 1), AIC(m), 1e-6)

# S25, S26, S28: leverage, deviance residual, Cook's distance
check("S25", "leverage", h, hatvalues(m), 1e-10)
check("S26", "squared deviance residual", 2 * ((y - mu) / mu - log(y / mu)), residuals(m, "deviance")^2, 1e-10)
check("S28", "Cook's distance", r^2 * h / (k * phi * (1 - h)^2), cooks.distance(m), 1e-10)

tab <- do.call(rbind, checks)
write_table(tab, "tableS_formula_checks", CFG$out_dir)
print(tab, digits = 3, row.names = FALSE)
if (!all(tab$agrees)) stop("A formula of Text A in S1 File does not reproduce the value computed by R: ",
                           paste(tab$equation[!tab$agrees], collapse = ", "))

# =============================================================================
# 08 - Audit of the statistics quoted in the first submission
#      (Manuscript SHBG.BSD.2sep26.Vcorrected.docx) against the data.
#      Each row gives the value printed in the manuscript, the value recomputed
#      here and whether they agree.
#   audit_first_submission.csv
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
all10    <- readRDS(file.path(CFG$out_dir, "data", "visit10_all.rds"))
analytic <- readRDS(file.path(CFG$out_dir, "data", "analytic.rds"))
pm       <- readRDS(file.path(CFG$out_dir, "data", "primary_model.rds"))
suppressPackageStartupMessages(library(readstata13))
raw <- suppressWarnings(read.dta13(CFG$data_file))

rows <- list()
chk <- function(item, reported, recomputed, tol, comment = "") {
  agree <- if (is.numeric(reported) && is.numeric(recomputed)) abs(reported - recomputed) <= tol else NA
  rows[[length(rows) + 1]] <<- data.frame(item = item, reported = as.character(reported),
    recomputed = as.character(signif(recomputed, 6)), agrees = agree, comment = comment, stringsAsFactors = FALSE)
}

# ---- descriptive statistics (full visit 10 file, as in the first submission) ----
d <- function(v, f) f(raw[[v]], na.rm = TRUE)
chk("Age, mean",              55.97, d("AGE10", mean), 0.005)
chk("Age, SD",                2.69,  d("AGE10", sd), 0.005)
chk("BMI, mean",              29.02, d("BMI10", mean), 0.005)
chk("BMI, SD",                7.21,  d("BMI10", sd), 0.005)
chk("Fat-free mass, mean",    45.48, d("FFMNHAN10", mean), 0.005)
chk("Percent body fat, mean", 36.07, d("PBFBIA10", mean), 0.005)
chk("Skeletal muscle mass, mean", 20.3, d("SKELMM10", mean), 0.05)
chk("Skeletal muscle mass, SD",   3.23, d("SKELMM10", sd), 0.005)
chk("SHBG, mean",             45.64, d("SHBG10", mean), 0.005)
chk("SHBG, SD",               27.44, d("SHBG10", sd), 0.005)
chk("DHEA-S, mean",           122.14, d("DHAS10", mean), 0.005)
chk("DHEA-S, SD (text)",      76.29, d("DHAS10", sd), 0.005, "Text says 76.29, Additional file 2 says 76.4")
chk("Estradiol, SD",          40.33, d("E2AVE10", sd), 0.005)
chk("FSH, mean",              110,   d("FSH10", mean), 0.5)
chk("TSH, SD",                2.75,  d("TSH10", sd), 0.005)
chk("Testosterone, SD",       41.8,  d("T10", sd), 0.05)
chk("Post-menopausal, %",     75.4,  100 * sum(as.character(raw$STATUS10) == "Post-menopausal", na.rm = TRUE) / nrow(raw), 0.05)

# ---- Spearman correlations quoted in the text (pairwise complete, full file) -----
sp <- function(a, b) cor(raw[[a]], raw[[b]], use = "pairwise.complete.obs", method = "spearman")
chk("Spearman testosterone - fat-free mass", 0.37, sp("T10", "FFMNHAN10"), 0.005,
    "0.37 is the testosterone - DHEA-S correlation")
chk("Spearman testosterone - DHEA-S (not quoted)", NA, sp("T10", "DHAS10"), 0)
chk("Spearman testosterone - hip BMD", 0.21, sp("T10", "HPBMDT10"), 0.005)
chk("Spearman testosterone - SHBG", -0.14, sp("T10", "SHBG10"), 0.005)
chk("Spearman estradiol - percent body fat", 0.24, sp("E2AVE10", "PBFBIA10"), 0.005)
chk("Spearman estradiol - hip BMD", 0.32, sp("E2AVE10", "HPBMDT10"), 0.005)
chk("Spearman estradiol - skeletal muscle mass", 0.26, sp("E2AVE10", "SKELMM10"), 0.005)

# ---- ANOVA by menopausal status (all available observations) --------------------
st <- droplevels(raw$STATUS10)
chk("ANOVA skeletal muscle mass by status, P", 1.7e-7, anova(lm(raw$SKELMM10 ~ st))$`Pr(>F)`[1], 0.05e-7)
chk("ANOVA SHBG by status, P", 0.073, anova(lm(raw$SHBG10 ~ st))$`Pr(>F)`[1], 0.0005)

# pairwise t tests on skeletal muscle mass: unadjusted versus Bonferroni
pw_un  <- pairwise.t.test(raw$SKELMM10, st, p.adjust.method = "none", pool.sd = FALSE)$p.value
pw_bon <- pairwise.t.test(raw$SKELMM10, st, p.adjust.method = "bonferroni", pool.sd = FALSE)$p.value
chk("Pairwise SMM comparisons with unadjusted P < 0.05 (of 21)", NA, sum(pw_un < 0.05, na.rm = TRUE), 0,
    "The figure displayed unadjusted P values although the Methods state Bonferroni correction")
chk("Pairwise SMM comparisons with Bonferroni P < 0.05 (of 21)", NA, sum(pw_bon < 0.05, na.rm = TRUE), 0)

# ---- model of the first submission (SKELMMM10) ------------------------------------
m <- pm$submitted
chk("SKELMMM10 observations", 1693, nobs(m), 0)
chk("SKELMMM10 AIC", 8383.171, AIC(m), 0.001)
chk("SKELMMM10 log likelihood", -4173.585, as.numeric(logLik(m)), 0.001,
    "stargazer derives the log likelihood without the dispersion parameter; logLik() differs by 1")
chk("SKELMMM10 SHBG coefficient", -0.001, coef(m)[["shbg"]], 0.0005)
chk("SKELMMM10 testosterone coefficient", 0.0003, coef(m)[["testosterone"]], 0.00005)
chk("SKELMMM10 RMSE", 2.956, sqrt(mean(residuals(m, "response")^2)), 0.0005)
chk("SKELMMM10 Sigma", 0.145, sqrt(summary(m)$dispersion), 0.0005)
chk("Largest Cook's distance, SKELMMM10", NA, max(pm$cooks_submitted), 0,
    "Text states that no observation was unusually influential")
chk("Largest adjusted GVIF, SKELMMM10", NA, max(car::vif(m)[, 3]^2), 0,
    "Text states that TSH exceeded the VIF threshold of 5")

audit <- do.call(rbind, rows)
write_table(audit, "audit_first_submission", CFG$out_dir)
print(audit, row.names = FALSE)

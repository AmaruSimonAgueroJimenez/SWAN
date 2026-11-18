# =============================================================================
# 02 - Descriptive statistics
#   table1_characteristics.csv      analytic sample, overall and by SHBG tertile
#   tableS_included_excluded.csv    analytic sample vs visit 10 participants not analysed
#   tableS_by_status.csv            SHBG and skeletal muscle mass by menopausal status
#   tableS_status_pairwise.csv      multiplicity-adjusted pairwise comparisons
#   tableS_spearman.csv             Spearman correlations (analytic sample)
#   tableS_smm_by_shbg_and_waist_tertile.csv  mean skeletal muscle mass by tertile of SHBG,
#                                   in all women and within tertiles of waist circumference (Fig 1)
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
all10    <- readRDS(file.path(CFG$out_dir, "data", "visit10_all.rds"))
analytic <- readRDS(file.path(CFG$out_dir, "data", "analytic.rds"))

# ---- Table 1 ----------------------------------------------------------------
cuts <- quantile(analytic$shbg, c(0, 1/3, 2/3, 1), names = FALSE)
analytic$shbg_tertile <- cut(analytic$shbg, cuts, include.lowest = TRUE,
                             labels = c("Lowest", "Middle", "Highest"))
stopifnot(!anyNA(analytic$shbg_tertile))
tertile_range <- tapply(analytic$shbg, analytic$shbg_tertile, function(x) sprintf("%.1f to %.1f", min(x), max(x)))

groups <- c(list(Overall = analytic), split(analytic, analytic$shbg_tertile))

cont_mean <- list(   # label = c(variable, digits)   -> mean (SD)
  "Age, years"                          = c("age", 1),
  "Height, cm"                          = c("height_cm", 1),
  "Weight, kg"                          = c("weight_kg", 1),
  "Body mass index, kg/m2"              = c("bmi", 1),
  "Waist circumference, cm"             = c("waist_cm", 1),
  "Body fat, %"                         = c("pbf", 1),
  "Skeletal muscle mass, kg"            = c("smm", 1),
  "Skeletal muscle mass index, kg/m2"   = c("smi", 2),
  "Skeletal muscle mass, % of body weight" = c("smm_pct_weight", 1)
)
cont_median <- list( # label = c(variable, digits)   -> median (IQR)
  "SHBG, nmol/L"            = c("shbg", 1),
  "Testosterone, ng/dL"     = c("testosterone", 1),
  "Estradiol, pg/mL"        = c("estradiol", 1),
  "FSH, mIU/mL"             = c("fsh", 1),
  "DHEA-S, ug/dL"           = c("dheas", 1),
  "TSH, uIU/mL"             = c("tsh", 2)
)
cat_vars <- list(
  "Race/ethnicity"                    = "race",
  "Menopausal status"                 = "status",
  "Health limits vigorous activities" = "activity_limit",
  "Exercise used for health, past 12 months" = "exercise"
)

rows <- list()
add_row <- function(label, values) rows[[length(rows) + 1]] <<- c(Characteristic = label, values)

add_row("n", sapply(groups, function(g) fmt_num(nrow(g), 0)))
add_row("SHBG range, nmol/L", c(Overall = sprintf("%.1f to %.1f", min(analytic$shbg), max(analytic$shbg)), tertile_range))
for (lab in names(cont_mean)) {
  v <- cont_mean[[lab]][1]; dg <- as.integer(cont_mean[[lab]][2])
  add_row(paste0(lab, ", mean (SD)"), sapply(groups, function(g) mean_sd(g[[v]], dg)))
}
for (lab in names(cont_median)) {
  v <- cont_median[[lab]][1]; dg <- as.integer(cont_median[[lab]][2])
  add_row(paste0(lab, ", median (IQR)"), sapply(groups, function(g) median_iqr(g[[v]], dg)))
  if (v %in% names(CFG$si)) {   # same quantity in SI units (PLOS ONE asks for SI values)
    si <- CFG$si[[v]]
    add_row(sprintf("%s, %s, median (IQR)", sub(",.*$", "", lab), si$unit),
            sapply(groups, function(g) median_iqr(g[[v]] * si$factor, si$digits)))
  }
}
for (lab in names(cat_vars)) {
  v <- cat_vars[[lab]]
  add_row(paste0(lab, ", n (%)"), setNames(rep("", length(groups)), names(groups)))
  lev_order <- if (v == "status") c("Premenopause", "Early perimenopause", "Late perimenopause",
                                    "Natural post-menopause", "Surgical menopause (BSO)",
                                    "Undetermined (hormone therapy)", "Undetermined (hysterectomy)")
               else if (v == "race") c("Black", "Chinese", "Japanese", "White")
               else levels(analytic[[v]])
  if (v == "exercise") lev_order <- "Yes"
  for (lv in lev_order) add_row(paste0("  ", lv), sapply(groups, function(g) n_pct(g[[v]], lv)))
}
table1 <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
write_table(table1, "table1_characteristics", CFG$out_dir)

# ---- included vs not included -----------------------------------------------
excluded <- all10[!all10$analytic, ]
desc_avail <- function(x, digits, type) {
  x <- x[!is.na(x)]
  val <- if (type == "mean") mean_sd(x, digits) else median_iqr(x, digits)
  sprintf("%s [n = %s]", val, fmt_num(length(x), 0))
}
vars_ie <- list(
  "Age, years, mean (SD)"                    = c("age", 1, "mean"),
  "Body mass index, kg/m2, mean (SD)"        = c("bmi", 1, "mean"),
  "Skeletal muscle mass, kg, mean (SD)"      = c("smm", 1, "mean"),
  "Body fat, %, mean (SD)"                   = c("pbf", 1, "mean"),
  "SHBG, nmol/L, median (IQR)"               = c("shbg", 1, "median"),
  "Testosterone, ng/dL, median (IQR)"        = c("testosterone", 1, "median"),
  "Estradiol, pg/mL, median (IQR)"           = c("estradiol", 1, "median"),
  "FSH, mIU/mL, median (IQR)"                = c("fsh", 1, "median")
)
ie <- lapply(names(vars_ie), function(lab) {
  s <- vars_ie[[lab]]
  data.frame(Characteristic = lab,
             Analysed = desc_avail(analytic[[s[1]]], as.integer(s[2]), s[3]),
             Not_analysed = desc_avail(excluded[[s[1]]], as.integer(s[2]), s[3]))
})
cat_ie <- function(v, label, levs) {
  out <- data.frame(Characteristic = paste0(label, ", n (%)"), Analysed = "", Not_analysed = "")
  xe <- excluded[[v]][!is.na(excluded[[v]])]
  for (lv in levs) out <- rbind(out, data.frame(
    Characteristic = paste0("  ", lv),
    Analysed = n_pct(analytic[[v]], lv),
    Not_analysed = sprintf("%s (%s)", sum(xe == lv), fmt_num(100 * sum(xe == lv) / length(xe), 1))))
  out
}
ie_tab <- rbind(
  data.frame(Characteristic = "n", Analysed = fmt_num(nrow(analytic), 0), Not_analysed = fmt_num(nrow(excluded), 0)),
  do.call(rbind, ie),
  cat_ie("race", "Race/ethnicity", c("Black", "Chinese", "Japanese", "White")),
  cat_ie("status", "Menopausal status", c("Premenopause", "Early perimenopause", "Late perimenopause",
         "Natural post-menopause", "Surgical menopause (BSO)", "Undetermined (hormone therapy)", "Undetermined (hysterectomy)"))
)
write_table(ie_tab, "tableS_included_excluded", CFG$out_dir)

# ---- SHBG and skeletal muscle mass by menopausal status ---------------------
status_levels <- c("Premenopause", "Early perimenopause", "Late perimenopause", "Natural post-menopause",
                   "Surgical menopause (BSO)", "Undetermined (hormone therapy)", "Undetermined (hysterectomy)")
by_status <- do.call(rbind, lapply(status_levels, function(lv) {
  g <- analytic[analytic$status == lv, ]
  data.frame(Menopausal_status = lv, n = nrow(g),
             SMM_mean_SD = mean_sd(g$smm, 1), SMM_median_IQR = median_iqr(g$smm, 1),
             SHBG_mean_SD = mean_sd(g$shbg, 1), SHBG_median_IQR = median_iqr(g$shbg, 1))
}))
# global tests: Welch one-way ANOVA (unequal variances) and Kruskal-Wallis
w_smm <- oneway.test(smm ~ status, analytic, var.equal = FALSE)
w_shbg <- oneway.test(shbg ~ status, analytic, var.equal = FALSE)
k_smm <- kruskal.test(smm ~ status, analytic)
k_shbg <- kruskal.test(shbg ~ status, analytic)
global <- data.frame(
  variable = c("Skeletal muscle mass", "SHBG"),
  anova_classic_p = c(anova(lm(smm ~ status, analytic))$`Pr(>F)`[1], anova(lm(shbg ~ status, analytic))$`Pr(>F)`[1]),
  welch_F         = c(w_smm$statistic[[1]], w_shbg$statistic[[1]]),
  welch_df1       = c(w_smm$parameter[[1]], w_shbg$parameter[[1]]),
  welch_df2       = c(w_smm$parameter[[2]], w_shbg$parameter[[2]]),
  anova_welch_p   = c(w_smm$p.value, w_shbg$p.value),
  kruskal_chi2    = c(k_smm$statistic[[1]], k_shbg$statistic[[1]]),
  kruskal_df      = c(k_smm$parameter[[1]], k_shbg$parameter[[1]]),
  kruskal_p       = c(k_smm$p.value, k_shbg$p.value)
)
write_table(by_status, "tableS_by_status", CFG$out_dir)
write_table(global, "tableS_by_status_global_tests", CFG$out_dir)

# all 21 pairwise Welch t-tests, unadjusted and Bonferroni/Holm adjusted
pairwise_tab <- function(var) {
  pr <- combn(status_levels, 2, simplify = FALSE)
  out <- do.call(rbind, lapply(pr, function(p) {
    a <- analytic[[var]][analytic$status == p[1]]; b <- analytic[[var]][analytic$status == p[2]]
    tt <- t.test(a, b, var.equal = FALSE)
    data.frame(variable = var, group1 = p[1], group2 = p[2], n1 = length(a), n2 = length(b),
               mean1 = mean(a), mean2 = mean(b), diff = mean(a) - mean(b), p_unadjusted = tt$p.value)
  }))
  out$p_bonferroni <- p.adjust(out$p_unadjusted, "bonferroni")
  out$p_holm       <- p.adjust(out$p_unadjusted, "holm")
  out
}
pw <- rbind(pairwise_tab("smm"), pairwise_tab("shbg"))
write_table(pw, "tableS_status_pairwise", CFG$out_dir)

# ---- Spearman correlations --------------------------------------------------
cv <- c(SHBG = "shbg", `Skeletal muscle mass` = "smm", `Body mass index` = "bmi", `Body fat %` = "pbf",
        `Waist circumference` = "waist_cm", Height = "height_cm", Age = "age", Testosterone = "testosterone",
        Estradiol = "estradiol", FSH = "fsh", `DHEA-S` = "dheas", TSH = "tsh")
rho <- cor(analytic[, cv], method = "spearman")
dimnames(rho) <- list(names(cv), names(cv))
write_table(data.frame(Variable = rownames(rho), round(rho, 2), check.names = FALSE), "tableS_spearman", CFG$out_dir)

# ---- Skeletal muscle mass by tertile of SHBG, overall and within tertiles of waist circumference ----
# Unadjusted data behind Fig 1: does the gradient across SHBG remain within strata of central adiposity?
waist_cuts <- quantile(analytic$waist_cm, c(0, 1/3, 2/3, 1), names = FALSE)
analytic$waist_tertile <- cut(analytic$waist_cm, waist_cuts, include.lowest = TRUE, labels = c("Lowest", "Middle", "Highest"))
stopifnot(!anyNA(analytic$waist_tertile))
waist_range <- tapply(analytic$waist_cm, analytic$waist_tertile, function(x) sprintf("%.1f to %.1f", min(x), max(x)))
cell_mean <- function(d, waist_label) {                 # mean with a t interval for one cell
  x <- d$smm; n <- length(x); se <- sd(x) / sqrt(n); tq <- qt(1 - (1 - CFG$conf_level) / 2, n - 1)
  data.frame(waist_tertile = waist_label, shbg_tertile = as.character(d$shbg_tertile[1]), n = n,
             mean_smm = mean(x), sd_smm = sd(x), lcl = mean(x) - tq * se, ucl = mean(x) + tq * se, stringsAsFactors = FALSE)
}
strata <- do.call(rbind, c(
  lapply(split(analytic, analytic$shbg_tertile), cell_mean, waist_label = "All women"),
  unlist(lapply(levels(analytic$waist_tertile), function(w) {
    d <- analytic[analytic$waist_tertile == w, ]
    lapply(split(d, d$shbg_tertile), cell_mean, waist_label = w)
  }), recursive = FALSE)))
rownames(strata) <- NULL
strata$pct_of_shbg_tertile <- 100 * strata$n / as.vector(table(analytic$shbg_tertile)[strata$shbg_tertile])
stopifnot(nrow(strata) == 12, all(strata$n >= 30),                                   # every cell supports a mean with a t interval
          sum(strata$n[strata$waist_tertile != "All women"]) == nrow(analytic),      # the strata partition the sample
          all(abs(tapply(strata$pct_of_shbg_tertile[strata$waist_tertile != "All women"],
                         strata$shbg_tertile[strata$waist_tertile != "All women"], sum) - 100) < 1e-8))
write_table(strata, "tableS_smm_by_shbg_and_waist_tertile", CFG$out_dir)

saveRDS(list(cuts = cuts, tertile_range = tertile_range, global = global, pairwise = pw, rho = rho,
             strata = strata, waist_cuts = waist_cuts, waist_range = waist_range),
        file.path(CFG$out_dir, "data", "descriptives.rds"))
print(table1, row.names = FALSE); print(by_status, row.names = FALSE); print(global)
print(strata, digits = 4, row.names = FALSE)
print(subset(pw, p_unadjusted < 0.05)[, c("variable","group1","group2","diff","p_unadjusted","p_bonferroni","p_holm")], row.names = FALSE)

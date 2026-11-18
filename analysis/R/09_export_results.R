# =============================================================================
# 09 - Collect every number quoted in the manuscript into one JSON file, plus
#      the session information, so that the text can be checked against the
#      analysis output.
#   results.json, session_info.txt
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
suppressPackageStartupMessages(library(jsonlite))
rd <- function(name) readRDS(file.path(CFG$out_dir, "data", paste0(name, ".rds")))
all10 <- rd("visit10_all"); analytic <- rd("analytic"); pm <- rd("primary_model")
sens <- rd("sensitivity"); desc <- rd("descriptives"); cmp <- rd("model_comparison"); so <- rd("shbg_outcome"); fam <- rd("model_families")
flow <- read.csv(file.path(CFG$out_dir, "tables", "flow_participants.csv"))

row <- function(df, key, col) { r <- df[df[[col]] == key, ]; stopifnot(nrow(r) == 1); as.list(r) }
t2 <- pm$table2; ts <- pm$table_submitted; s <- sens$sens
so_lab <- "Covariates of the primary model"; so_p <- so$tab[so$tab$model == so_lab, ]
stopifnot(nrow(so_p) > 0, "smm" %in% so_p$term)
q_shbg <- quantile(analytic$shbg, c(0.25, 0.5, 0.75), names = FALSE)
pw <- desc$pairwise

res <- list(
  flow = setNames(as.list(flow$n), c("visit10", "shbg", "smm", "both", "analytic", "analytic_body")),
  sample = list(
    n = nrow(analytic),
    age_mean = mean(analytic$age), age_sd = sd(analytic$age), age_min = min(analytic$age), age_max = max(analytic$age),
    shbg_median = q_shbg[2], shbg_q1 = q_shbg[1], shbg_q3 = q_shbg[3],
    shbg_mean = mean(analytic$shbg), shbg_sd = sd(analytic$shbg),
    smm_mean = mean(analytic$smm), smm_sd = sd(analytic$smm),
    bmi_mean = mean(analytic$bmi), bmi_sd = sd(analytic$bmi),
    pct_post = 100 * mean(analytic$status == "Natural post-menopause"),
    n_hormone_users_excluded = nrow(analytic) - s$n[s$id == "F2"],
    n_cook_excluded = nrow(analytic) - s$n[s$id == "F1"],
    spearman_shbg_waist = cor(analytic$shbg, analytic$waist_cm, method = "spearman"),
    spearman_smm_waist = cor(analytic$smm, analytic$waist_cm, method = "spearman"),
    spearman_shbg_height = cor(analytic$shbg, analytic$height_cm, method = "spearman"),
    spearman_smm_height = cor(analytic$smm, analytic$height_cm, method = "spearman"),
    tertile_waist = as.list(tapply(analytic$waist_cm, cut(analytic$shbg, desc$cuts, include.lowest = TRUE, labels = c("Lowest", "Middle", "Highest")), mean)),
    tertile_range = as.list(desc$tertile_range),
    tertile_bmi = as.list(tapply(analytic$bmi, cut(analytic$shbg, desc$cuts, include.lowest = TRUE, labels = c("Lowest", "Middle", "Highest")), mean)),
    tertile_smm = as.list(tapply(analytic$smm, cut(analytic$shbg, desc$cuts, include.lowest = TRUE, labels = c("Lowest", "Middle", "Highest")), mean)),
    tertile_smm_pct = as.list(tapply(analytic$smm_pct_weight, cut(analytic$shbg, desc$cuts, include.lowest = TRUE, labels = c("Lowest", "Middle", "Highest")), mean)),
    tertile_weight = as.list(tapply(analytic$weight_kg, cut(analytic$shbg, desc$cuts, include.lowest = TRUE, labels = c("Lowest", "Middle", "Highest")), mean)),
    spearman_shbg_smm = desc$rho["SHBG", "Skeletal muscle mass"],
    spearman_shbg_bmi = desc$rho["SHBG", "Body mass index"],
    spearman_smm_bmi  = desc$rho["Skeletal muscle mass", "Body mass index"],
    spearman_smm_pbf  = desc$rho["Skeletal muscle mass", "Body fat %"],
    spearman_shbg_pbf = desc$rho["SHBG", "Body fat %"],
    # fat-free mass implied by weight and percent body fat (both from impedance) against skeletal muscle mass
    pearson_smm_ffm = cor(analytic$smm, analytic$weight_kg - analytic$fat_mass_kg),
    # unadjusted means by tertile of SHBG, in all women and within tertiles of waist circumference (Fig 1)
    strata = lapply(seq_len(nrow(desc$strata)), function(i) as.list(desc$strata[i, ])),
    waist_cuts = round(desc$waist_cuts[2:3], 1), waist_range = as.list(desc$waist_range),
    spearman_smm_fat_mass = cor(analytic$smm, analytic$fat_mass_kg, method = "spearman"),
    spearman_shbg_fat_mass = cor(analytic$shbg, analytic$fat_mass_kg, method = "spearman")
  ),
  status = list(
    global = desc$global,
    smm_post_vs_late = row(pw[pw$variable == "smm" & pw$group1 == "Late perimenopause", ], "Natural post-menopause", "group2"),
    smm_post_vs_bso  = row(pw[pw$variable == "smm" & pw$group1 == "Natural post-menopause", ], "Surgical menopause (BSO)", "group2"),
    smm_post_vs_early = row(pw[pw$variable == "smm" & pw$group1 == "Early perimenopause", ], "Natural post-menopause", "group2"),
    n_smm_pairs_bonferroni = sum(pw$p_bonferroni[pw$variable == "smm"] < 0.05),
    n_shbg_pairs_bonferroni = sum(pw$p_bonferroni[pw$variable == "shbg"] < 0.05)
  ),
  primary = list(
    terms = lapply(split(t2[, c("label", "ratio", "lcl", "ucl", "pct_difference", "p")], t2$term), as.list),
    joint = pm$joint,
    ame_kg_per_doubling = pm$ame_doubling,
    predicted = pm$pred,
    dispersion = pm$phi, aic = AIC(pm$model), n_coefficients = length(coef(pm$model)),
    deviance_explained = 1 - deviance(pm$model) / pm$model$null.deviance,
    max_cook = max(pm$cooks), n_cook_gt_4n = sum(pm$cooks > 4 / nrow(analytic)),
    max_leverage = max(hatvalues(pm$model))
  ),
  # model of the first submission, now a sensitivity analysis
  submitted = list(
    terms = lapply(split(ts[, c("label", "ratio", "lcl", "ucl", "pct_difference", "p")], ts$term), as.list),
    aic = AIC(pm$submitted), dispersion = summary(pm$submitted)$dispersion,
    deviance_explained = 1 - deviance(pm$submitted) / pm$submitted$null.deviance,
    max_cook = max(pm$cooks_submitted), max_leverage = max(hatvalues(pm$submitted)),
    max_testosterone_ng_dl = max(analytic$testosterone),
    max_testosterone_nmol_l = max(analytic$testosterone) * CFG$si$testosterone$factor
  ),
  sensitivity = lapply(split(s[, c("specification", "n", "ratio", "lcl", "ucl", "pct_difference", "p")], s$id), as.list),
  si = CFG$si,
  alternative_families = fam$alt,
  variance_function = setNames(as.list(fam$park$value), c("lambda", "lcl", "ucl", "p_gamma", "p_gaussian", "p_inverse_gaussian")),
  nonlinearity = sens$nonlin,
  model_comparison = cmp[, c("model", "aic", "delta_aic", "rmse_kg", "change_smm_kg_per_doubling")],
  shbg_outcome = list(
    primary = lapply(split(so_p[, c("ratio", "lcl", "ucl", "p")], so_p$term), as.list),
    submitted_smm = as.list(so$tab[so$tab$model != so_lab & so$tab$term == "smm", c("ratio", "lcl", "ucl", "p")]),
    # the model with the covariates of the primary model is Table 3 of the article
    joint = so$joint,
    n_coefficients = length(coef(so$m_a)), residual_df = so$m_a$df.residual, dispersion = summary(so$m_a)$dispersion,
    deviance_explained = 1 - deviance(so$m_a) / so$m_a$null.deviance
  ),
  software = list(r = R.version.string)
)
writeLines(toJSON(res, auto_unbox = TRUE, digits = 8, pretty = TRUE), file.path(CFG$out_dir, "results.json"))
writeLines(capture.output(sessionInfo()), file.path(CFG$out_dir, "session_info.txt"))
message("  wrote results.json and session_info.txt")

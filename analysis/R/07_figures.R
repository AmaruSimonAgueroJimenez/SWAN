# =============================================================================
# 07 - Figures (PLOS ONE specifications: TIFF, LZW, 600 dpi, width <= 7.5 in,
#      height <= 8.75 in, sans-serif text 8-12 pt). PNG copies are previews.
#   Fig1  the unadjusted data: (A) mean muscle mass by tertile of SHBG, in all women and within tertiles of
#         waist circumference, (B) distribution of the tertiles of waist circumference in each tertile of SHBG
#   Fig2  the association and its dependence on the covariates: (A) adjusted mean muscle mass across SHBG,
#         (B) estimate when one covariate is removed, (C) estimate with other measures of body size
#   Fig2_overview  panel A of Fig2 alone (used in the report)
#   S1_Fig  skeletal muscle mass and SHBG by menopausal status
#   S2_Fig  diagnostics of the primary model
#   S3_Fig  shape of the SHBG association (term of the model versus spline)
#   S4_Fig  SHBG association under sequential adjustment, other scalings, specifications and samples
# =============================================================================
source("analysis/config.R"); source("analysis/R/functions.R")
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(splines)})
analytic <- readRDS(file.path(CFG$out_dir, "data", "analytic.rds"))
pm   <- readRDS(file.path(CFG$out_dir, "data", "primary_model.rds"))
sens <- readRDS(file.path(CFG$out_dir, "data", "sensitivity.rds"))
fig_dir <- file.path(CFG$out_dir, "figures")

FONT <- "Liberation Sans"            # metric-compatible with Arial
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"
COL_PRIMARY <- "#2a78d6"; COL_BODY <- "#eb6834"   # validated colour-blind-safe pair

theme_fig <- function(base = 9) {
  theme_minimal(base_size = base, base_family = FONT) +
    theme(text = element_text(colour = INK),
          axis.text = element_text(colour = INK2),
          axis.title = element_text(colour = INK),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          axis.line = element_line(colour = INK2, linewidth = 0.3),
          axis.ticks = element_line(colour = INK2, linewidth = 0.3),
          plot.title = element_text(face = "bold", size = base + 1),
          plot.tag = element_text(face = "bold", size = base + 3),
          legend.position = "none",
          plot.margin = margin(6, 8, 6, 6))
}

save_fig <- function(plot, name, width_in, height_in) {
  stopifnot(width_in <= 7.5, height_in <= 8.75)
  ggsave(file.path(fig_dir, paste0(name, ".tif")), plot, width = width_in, height = height_in,
         dpi = 600, device = "tiff", compression = "lzw", type = "cairo", bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width = width_in, height = height_in,
         dpi = 200, type = "cairo", bg = "white")
  message("  wrote ", name)
}

# ---- adjusted means across SHBG (curves computed in 04_sensitivity.R) --------------
cu   <- sens$curves
rng  <- sens$shbg_range
curve_of <- function(id) { d <- cu[cu$id == id, ]; stopifnot(nrow(d) > 0); d }
# the standardised curve must average close to the observed mean at the observed SHBG
stopifnot(abs(mean(fitted(pm$model)) - mean(analytic$smm)) < 0.05)
rug_data <- analytic[analytic$shbg >= rng[1] & analytic$shbg <= rng[2], ]
stopifnot(l10n_info()[["UTF-8"]])   # the minus sign and superscripts need a UTF-8 locale
minus <- function(x) gsub("-", "−", x, fixed = TRUE)
HEAD_SIZE <- 2.75                   # size of the bold panel headers (mm)

# Primary model against the same model without waist circumference
curve_panel <- function(title = NULL) {
  c1 <- cbind(curve_of("A5"), model = "primary")
  c2 <- cbind(curve_of("B1"), model = "nowaist")
  lab <- data.frame(model = c("primary", "nowaist"), x = max(c1$shbg),
                    y = c(tail(c1$mean, 1), tail(c2$mean, 1)),
                    label = c("Primary model", "Primary model without\nwaist circumference"))
  ggplot(rbind(c1, c2), aes(shbg, mean, colour = model, fill = model)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.18, colour = NA) +
    geom_line(aes(linetype = model), linewidth = 0.8) +
    geom_rug(data = rug_data, aes(x = shbg), inherit.aes = FALSE, sides = "b", alpha = 0.12, linewidth = 0.2,
             colour = INK2, length = unit(0.03, "npc")) +
    geom_text(data = lab, aes(x, y, label = label), colour = INK, hjust = 0, nudge_x = 1.5,
              size = 2.9, family = FONT, lineheight = 0.95, inherit.aes = FALSE) +
    scale_colour_manual(values = c(primary = COL_PRIMARY, nowaist = COL_BODY)) +
    scale_fill_manual(values = c(primary = COL_PRIMARY, nowaist = COL_BODY)) +
    scale_linetype_manual(values = c(primary = "solid", nowaist = "longdash")) +
    scale_x_continuous(breaks = seq(20, 100, 20), expand = expansion(mult = c(0.02, 0.34))) +
    labs(x = "SHBG (nmol/L)", y = "Adjusted mean skeletal\nmuscle mass (kg)", title = title) +
    theme_fig(8.5) + theme(plot.title = element_text(face = "bold", size = 7.9, colour = INK, margin = margin(0, 0, 4, 0)),
                           plot.title.position = "plot")
}

# ---- forest plots of the SHBG estimate ----------------------------------------------
# rows:   named character vector, id of the sensitivity table -> label shown in the plot
# groups: named list, header of each block -> ids it contains (in display order)
# ref:    estimate of the primary model, drawn as a dashed reference line (or NULL)
forest <- function(rows, groups, xlim, xbreaks, ref = NULL) {
  stopifnot(setequal(names(rows), unlist(groups)), !anyDuplicated(unlist(groups)), !anyDuplicated(rows))
  all <- sens$sens; stopifnot(all(names(rows) %in% all$id))
  fd <- all[match(unlist(groups), all$id), ]
  fd$pct <- pct(fd$ratio); fd$pl <- pct(fd$lcl); fd$pu <- pct(fd$ucl)
  stopifnot(min(fd$pl) >= xlim[1], max(fd$pu) <= xlim[2])           # nothing may fall outside the axis
  fd$grp   <- factor(rep(names(groups), lengths(groups)), levels = names(groups))
  fd$label <- factor(rows[fd$id], levels = rev(unname(rows[unlist(groups)])))
  fd$txt   <- minus(sprintf("%.2f (%.2f to %.2f)", fd$pct, fd$pl, fd$pu))
  fd$primary <- fd$id == "A5"
  heads <- data.frame(grp = factor(names(groups), levels = names(groups)), label = names(groups))
  facet_theme <- theme(strip.text = element_blank(), panel.spacing.y = unit(15, "pt"),
                       panel.grid.major.y = element_blank())
  p_forest <- ggplot(fd, aes(pct, label)) +
    geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3) +
    { if (!is.null(ref)) geom_vline(xintercept = ref, colour = COL_PRIMARY, linewidth = 0.3, linetype = "22") } +
    geom_errorbarh(aes(xmin = pl, xmax = pu), height = 0, linewidth = 0.5, colour = INK) +
    geom_point(aes(fill = primary), shape = 21, size = 1.9, colour = INK, stroke = 0.4) +
    geom_text(data = heads, aes(x = -Inf, y = Inf, label = label), hjust = 0, vjust = -0.6, size = HEAD_SIZE,
              fontface = "bold", family = FONT, colour = INK, inherit.aes = FALSE) +
    scale_fill_manual(values = c(`TRUE` = COL_PRIMARY, `FALSE` = "white")) +
    scale_x_continuous(breaks = xbreaks, limits = xlim,
                       labels = function(b) minus(formatC(b, format = "f", digits = 0))) +
    facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
    coord_cartesian(clip = "off") +
    labs(x = "Difference in mean skeletal muscle mass\nper doubling of SHBG (%)", y = NULL) +
    theme_fig(8.5) + facet_theme +
    theme(axis.text.y = element_text(colour = INK, size = 8), plot.margin = margin(16, 2, 6, 6))
  p_text <- ggplot(fd, aes(0, label)) +
    geom_text(aes(label = txt), hjust = 0, size = 2.7, family = FONT, colour = INK) +
    geom_text(data = data.frame(grp = factor(names(groups)[1], levels = names(groups))),
              aes(x = 0, y = Inf, label = "% difference (95% CI)"), hjust = 0, vjust = -0.6, size = HEAD_SIZE,
              fontface = "bold", family = FONT, colour = INK, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
    coord_cartesian(clip = "off") +
    theme_void(base_family = FONT) + facet_theme + theme(plot.margin = margin(16, 6, 6, 0))
  p_forest + p_text + plot_layout(widths = c(3.1, 1.25))
}
primary_pct <- pct(sens$sens$ratio[sens$sens$id == "A5"])

# ---- Fig 1: the unadjusted data, in all women and within tertiles of waist circumference ----
# Tertiles of waist circumference are an ordered scale: one hue, light to dark. Identity is carried by
# direct labels (panel A) and by the same order of the segments (panel B), not by colour alone.
desc <- readRDS(file.path(CFG$out_dir, "data", "descriptives.rds"))
st <- desc$strata
RAMP <- c(Lowest = "#6aa6e8", Middle = "#2a78d6", Highest = "#103f7d")
TERT <- c("Lowest", "Middle", "Highest")
stopifnot(setequal(unique(st$waist_tertile), c("All women", TERT)), setequal(unique(st$shbg_tertile), TERT))
wc <- formatC(desc$waist_cuts[2:3], format = "f", digits = 0)
waist_lab <- c(`All women` = "All women", Lowest = paste0("Waist up to ", wc[1], " cm"),
               Middle = paste0("Waist ", wc[1], " to ", wc[2], " cm"), Highest = paste0("Waist above ", wc[2], " cm"))
shbg_lab <- setNames(paste0(TERT, "\n", desc$tertile_range[TERT]), TERT)
st$x   <- factor(shbg_lab[st$shbg_tertile], levels = shbg_lab)
st$grp <- factor(st$waist_tertile, levels = c("All women", TERT))
end_lab <- st[st$shbg_tertile == "Highest", ]
end_lab$label <- waist_lab[as.character(end_lab$grp)]
stopifnot(min(diff(sort(end_lab$mean_smm))) > 0.4)              # the direct labels do not collide
p1a <- ggplot(st, aes(x, mean_smm, group = grp, colour = grp)) +
  geom_line(aes(linewidth = grp, linetype = grp)) +
  geom_linerange(aes(ymin = lcl, ymax = ucl), linewidth = 0.4) +
  geom_point(aes(size = grp), shape = 21, fill = "white", stroke = 0.7) +
  geom_text(data = end_lab, aes(label = label), hjust = 0, nudge_x = 0.13, size = 2.8, family = FONT, colour = INK) +
  scale_colour_manual(values = c(`All women` = INK, RAMP)) +
  scale_linewidth_manual(values = c(`All women` = 0.9, Lowest = 0.6, Middle = 0.6, Highest = 0.6)) +
  scale_linetype_manual(values = c(`All women` = "22", Lowest = "solid", Middle = "solid", Highest = "solid")) +
  scale_size_manual(values = c(`All women` = 2.2, Lowest = 1.8, Middle = 1.8, Highest = 1.8)) +
  scale_x_discrete(expand = expansion(add = c(0.3, 1.6))) +
  labs(title = "A. Skeletal muscle mass by tertile of SHBG", x = "Tertile of SHBG (nmol/L)", y = "Mean skeletal muscle mass (kg)") +
  theme_fig(8.5) + theme(panel.grid.major.x = element_blank(),
                         plot.title = element_text(face = "bold", size = 7.9, colour = INK, margin = margin(0, 0, 6, 0)),
                         plot.title.position = "plot")
comp <- st[st$waist_tertile != "All women", ]
comp$seg <- factor(comp$waist_tertile, levels = rev(TERT))      # highest tertile of waist on top
p1b <- ggplot(comp, aes(x, pct_of_shbg_tertile, fill = seg)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", pct_of_shbg_tertile), colour = seg == "Lowest"), position = position_stack(vjust = 0.5),
            size = 2.7, family = FONT) +
  scale_colour_manual(values = c(`TRUE` = INK, `FALSE` = "white")) +     # dark text on the lightest segment
  scale_fill_manual(values = RAMP) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "B. Waist circumference by tertile of SHBG", x = "Tertile of SHBG (nmol/L)", y = "Women (%)") +
  theme_fig(8.5) + theme(panel.grid.major.x = element_blank(),
                         plot.title = element_text(face = "bold", size = 7.9, colour = INK, margin = margin(0, 0, 6, 0)),
                         plot.title.position = "plot")
save_fig(p1a + p1b + plot_layout(widths = c(1.3, 1)), "Fig1", 7.2, 3.5)

# ---- Fig 2: the association and its dependence on the covariates ----
rows_cov <- c(
  A5 = "Primary model",
  B1 = "Without waist circumference", B2 = "Without height", B3 = "Without race/ethnicity", B4 = "Without age",
  B5 = "Without menopausal status", B6 = "Without the five hormones", B7 = "Without activity limitation",
  B8 = "Without exercise",
  C1 = "BMI in place of height and waist", C2 = "Weight in place of waist", C3 = "Fat mass in place of waist",
  C4 = "+ percent body fat", C5 = "+ BMI", C6 = "+ fat mass", C7 = "Covariates of the initial model")
groups_cov <- list(
  "B. One covariate removed from the primary model" = c("A5", paste0("B", 1:8)),
  "C. Body size measured in other ways"             = paste0("C", 1:7))
fig2 <- curve_panel("A. Adjusted mean skeletal muscle mass across SHBG") /
  forest(rows_cov, groups_cov, xlim = c(-4.75, 1), xbreaks = seq(-4, 1, 1), ref = primary_pct) +
  plot_layout(heights = c(1.9, 5.1))
save_fig(fig2, "Fig2", 7.2, 7.6)

# the curves alone, for the report
save_fig(curve_panel(), "Fig2_overview", 5.6, 3.3)

# ---- S4 Fig: the other specifications (sequential adjustment, scaling, specification, samples) ----
rows_spec <- c(
  A1 = "Unadjusted", A2 = "Adjusted for age", A3 = "+ race/ethnicity + height", A4 = "+ waist circumference",
  A5 = "Primary model",
  D1 = "Muscle mass index (kg/m²)", D2 = "Muscle mass as % of body weight",
  E1 = "Robust standard errors", E2 = "Gaussian family, log link", E3 = "Linear model of log muscle mass",
  E4 = "Hormone covariates untransformed",
  F1 = "Excluding Cook's distance > 4/n", F2 = "Excluding sex-hormone users", F3 = "Natural menopausal stages only",
  F4 = "Natural post-menopause only", F5 = "Excluding extreme testosterone value (n = 1)")
groups_spec <- list(
  "Sequential adjustment"                    = paste0("A", 1:5),
  "Alternative scaling of muscle mass"       = c("D1", "D2"),
  "Primary model, alternative specification" = paste0("E", 1:4),
  "Primary model, restricted samples"        = paste0("F", 1:5))
# every row of the sensitivity table appears in exactly one of the two forest plots (A5 in both)
stopifnot(setequal(union(names(rows_cov), names(rows_spec)), sens$sens$id))
s4 <- forest(rows_spec, groups_spec, xlim = c(-5.5, 1), xbreaks = seq(-5, 1, 1), ref = primary_pct)
save_fig(s4, "S4_Fig", 7.2, 5.4)

# ---- S1 Fig: by menopausal status -------------------------------------------------
lv <- c("Premenopause", "Early perimenopause", "Late perimenopause", "Natural post-menopause",
        "Surgical menopause (BSO)", "Undetermined (hormone therapy)", "Undetermined (hysterectomy)")
sd_ <- analytic; sd_$status <- factor(as.character(sd_$status), levels = lv)
nlab <- as.data.frame(table(sd_$status)); names(nlab) <- c("status", "n")
wrap <- c("Premenopause" = "Premenopause", "Early perimenopause" = "Early\nperimenopause",
          "Late perimenopause" = "Late\nperimenopause", "Natural post-menopause" = "Natural\npost-menopause",
          "Surgical menopause (BSO)" = "Surgical\nmenopause (BSO)",
          "Undetermined (hormone therapy)" = "Undetermined\n(hormone therapy)",
          "Undetermined (hysterectomy)" = "Undetermined\n(hysterectomy)")
xlabs <- setNames(sprintf("%s\n(n = %s)", wrap[as.character(nlab$status)], format(nlab$n, big.mark = ",", trim = TRUE)), nlab$status)
box <- function(var, ylab) {
  ggplot(sd_, aes(status, .data[[var]])) +
    geom_boxplot(width = 0.55, outlier.size = 0.5, outlier.alpha = 0.35, linewidth = 0.35,
                 fill = "#dbe8f8", colour = INK) +
    scale_x_discrete(labels = xlabs) + labs(x = NULL, y = ylab) + theme_fig(8) +
    theme(axis.text.x = element_text(size = 6.6, colour = INK), panel.grid.major.x = element_blank())
}
s1 <- (box("smm", "Skeletal muscle mass (kg)") + labs(tag = "A")) /
      (box("shbg", "SHBG (nmol/L)") + labs(tag = "B"))
save_fig(s1, "S1_Fig", 7.2, 6.4)

# ---- S2 Fig: diagnostics ------------------------------------------------------------
m <- pm$model; mu <- fitted(m); phi <- pm$phi
# quantile residuals of a gamma GLM: Phi^-1{F(y; shape = 1/phi, scale = mu*phi)}
qres <- qnorm(pgamma(analytic$smm, shape = 1 / phi, scale = mu * phi))
dd <- data.frame(mu = mu, dres = residuals(m, "deviance") / sqrt(phi), qres = qres,
                 cook = pm$cooks, idx = seq_along(mu), smm = analytic$smm)
pA <- ggplot(dd, aes(mu, dres)) + geom_point(size = 0.4, alpha = 0.35, colour = INK2) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = INK) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = COL_PRIMARY, linewidth = 0.7) +
  labs(x = "Fitted skeletal muscle mass (kg)", y = "Standardised deviance residual", tag = "A") + theme_fig(8)
pB <- ggplot(dd, aes(sample = qres)) + stat_qq(size = 0.4, alpha = 0.35, colour = INK2) +
  stat_qq_line(colour = COL_PRIMARY, linewidth = 0.7) +
  labs(x = "Standard normal quantile", y = "Quantile residual", tag = "B") + theme_fig(8)
pC <- ggplot(dd, aes(idx, cook)) + geom_segment(aes(xend = idx, yend = 0), linewidth = 0.25, colour = INK2) +
  geom_hline(yintercept = 4 / nrow(dd), linetype = "dashed", linewidth = 0.3, colour = COL_BODY) +
  scale_y_continuous(trans = "sqrt", breaks = c(0.002, 0.01, 0.02, 0.05, 0.1)) +
  labs(x = "Observation", y = "Cook's distance (square-root scale)", tag = "C") + theme_fig(8)
dens <- rbind(data.frame(x = analytic$smm, what = "Observed"), data.frame(x = mu, what = "Fitted"))
pD <- ggplot(dd, aes(mu, smm)) + geom_point(size = 0.4, alpha = 0.3, colour = INK2) +
  geom_abline(slope = 1, intercept = 0, colour = COL_PRIMARY, linewidth = 0.7) +
  labs(x = "Fitted skeletal muscle mass (kg)", y = "Observed skeletal muscle mass (kg)", tag = "D") + theme_fig(8)
save_fig((pA + pB) / (pC + pD), "S2_Fig", 7.2, 6.2)

# ---- S3 Fig: shape of the association (term of the model against a spline) ---------------
shape_curve <- function(m_term, m_spl, term_label, label) {
  ref <- median(analytic$shbg)
  gridx <- seq(quantile(analytic$shbg, 0.01), quantile(analytic$shbg, 0.99), length.out = 80)
  # ratio relative to the median SHBG, other covariates cancel under the log link
  base <- analytic[1, ]
  nd <- base[rep(1, length(gridx)), ]; nd$shbg <- gridx
  nr <- base; nr$shbg <- ref
  rel <- function(m) exp(predict(m, nd, type = "link") - predict(m, nr, type = "link"))
  rbind(data.frame(shbg = gridx, ratio = rel(m_term), term = term_label, model = label),
        data.frame(shbg = gridx, ratio = rel(m_spl), term = "Natural cubic spline (3 df)", model = label))
}
lab_p <- "Primary model"; lab_s <- "Covariates of the initial model"
sc <- rbind(shape_curve(pm$model, sens$m_spline, "Term of the model", lab_p),
            shape_curve(pm$submitted, sens$m_spline_sub, "Term of the model", lab_s))
sc$model <- factor(sc$model, levels = c(lab_p, lab_s),
                   labels = c("A. Primary model: log2(SHBG)", "B. Initial model: SHBG, linear"))
s3 <- ggplot(sc, aes(shbg, ratio, linetype = term, colour = term)) +
  geom_hline(yintercept = 1, colour = INK2, linewidth = 0.3) +
  geom_line(linewidth = 0.8) + facet_wrap(~model) +
  scale_colour_manual(values = c("Term of the model" = COL_PRIMARY, "Natural cubic spline (3 df)" = COL_BODY)) +
  scale_linetype_manual(values = c("Term of the model" = "solid", "Natural cubic spline (3 df)" = "longdash")) +
  labs(x = "SHBG (nmol/L)", y = "Ratio of mean skeletal muscle mass\n(reference: median SHBG)", colour = NULL, linetype = NULL) +
  theme_fig(8.5) + theme(legend.position = "bottom", strip.text = element_text(face = "bold", hjust = 0, colour = INK),
                         legend.key.width = unit(26, "pt"))
save_fig(s3, "S3_Fig", 7.2, 3.6)

# numbers used in the figure legends
saveRDS(list(fig1_range = rng, curves = cu), file.path(CFG$out_dir, "data", "figures.rds"))

# Shared helpers of the report pages. Everything shown in the report is read from
# analysis/output (results.json, tables/*.csv, figures/*.png); nothing is refitted here.
suppressPackageStartupMessages(library(jsonlite))
OUT <- file.path("..", "analysis", "output")
stopifnot(file.exists(file.path(OUT, "results.json")))   # run `Rscript analysis/run_all.R` first
R   <- fromJSON(file.path(OUT, "results.json"), simplifyVector = FALSE)
tab <- function(name) read.csv(file.path(OUT, "tables", paste0(name, ".csv")),
                               check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
fig <- function(name) { p <- file.path(OUT, "figures", paste0(name, ".png")); stopifnot(file.exists(p)); p }

MINUS <- "−"
fmt <- function(x, d = 2) {                         # fixed decimals, thousands separator, true minus sign
  x <- as.numeric(x)
  s <- formatC(abs(x), format = "f", digits = d, big.mark = ",")
  neg <- !is.na(x) & x < 0 & as.numeric(formatC(abs(x), format = "f", digits = d)) != 0
  ifelse(is.na(x), "", ifelse(neg, paste0(MINUS, s), s))
}
pct   <- function(ratio) 100 * (as.numeric(ratio) - 1)             # percent difference from a ratio of means
ci    <- function(l, u, d = 2) paste0(fmt(l, d), " to ", fmt(u, d))
pval  <- function(p) { p <- as.numeric(p); ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)) }
ptext <- function(p) { p <- as.numeric(p); ifelse(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = 3))) }
ufix  <- function(s) {                                              # units as they should be displayed
  s <- gsub("kg/m2", "kg/m²", s, fixed = TRUE); s <- gsub("ug/dL", "µg/dL", s, fixed = TRUE)
  s <- gsub("uIU/mL", "µIU/mL", s, fixed = TRUE); gsub("umol/L", "µmol/L", s, fixed = TRUE)
}
rows  <- function(x) {                                              # list of JSON records -> data frame
  keys <- unique(unlist(lapply(x, names)))                          # a field that is NA in R is absent or null in the JSON
  do.call(rbind, lapply(x, function(r) as.data.frame(
    setNames(lapply(keys, function(k) if (is.null(r[[k]])) NA else r[[k]]), keys), stringsAsFactors = FALSE)))
}
term_label <- function(term, label) {                                # readable label of a model term
  pre <- ifelse(grepl("^status", term), "Menopausal status: ", ifelse(grepl("^activity_limit", term), "Vigorous activities: ",
         ifelse(grepl("^exercise", term), "Exercise for health: ", ifelse(grepl("^race", term), "Race/ethnicity: ", ""))))
  lab <- ifelse(is.na(label), sub("^race", "", term), trimws(label))
  lab <- ifelse(term == "smm", "Skeletal muscle mass, per 1 kg", lab)
  ufix(paste0(pre, lab))
}
kt    <- function(df, ...) knitr::kable(df, row.names = FALSE, ...)
est   <- function(o, d = 2) paste0(fmt(pct(o$ratio), d), "% (95% CI ", ci(pct(o$lcl), pct(o$ucl), d), ")")   # percent difference of a model
S <- R$sample; PR <- R$primary; SUB <- R$submitted; SE <- R$sensitivity

# ---- code of the pipeline, shown folded on the pages -------------------------------------------------------------
# The pages only read the output of the pipeline. show_code() displays the scripts that computed that output,
# folded and without running them. The scripts are read from the repository when the page is rendered, so
# `Rscript analysis/run_all.R` must have been run with the scripts as they are shown.
SCRIPTS <- c(
  "analysis/run_all.R"                         = "runs the ten scripts of the pipeline in order",
  "analysis/config.R"                          = "paths, reporting units, reference categories, expected values of the models",
  "analysis/R/functions.R"                     = "definition of the primary and initial models, model fitting, ratios of means, standardized means",
  "analysis/R/01_prepare_data.R"               = "reads the raw file, derives the variables, defines the analytic sample",
  "analysis/R/02_descriptives.R"               = "characteristics of the sample, comparison by menopausal status, correlations, means by tertile",
  "analysis/R/03_primary_model.R"              = "primary model, initial model, adjusted means, diagnostics of both",
  "analysis/R/04_sensitivity.R"                = "sensitivity models, shape of the association, full regression output",
  "analysis/R/05_model_comparison.R"           = "distributions and link functions, other families, variance function",
  "analysis/R/06_shbg_as_outcome.R"            = "complementary model with SHBG as the dependent variable",
  "analysis/R/07_figures.R"                    = "figures",
  "analysis/R/08_audit_first_submission.R"     = "recomputes the numbers reported in the first submission",
  "analysis/R/09_export_results.R"             = "collects every number quoted in the text into results.json",
  "analysis/R/10_check_formulas.R"             = "computes the formulas by hand and compares them with what R returns",
  "analysis/proposal/primary_model_proposal.R" = "analysis that led to the proposal of the primary model (not part of the pipeline)",
  "docs/_common.R"                             = "helpers shared by the pages of this report"
)
show_code <- function(...) {
  files <- c(...)
  stopifnot(length(files) > 0, !anyDuplicated(files), all(files %in% names(SCRIPTS)), all(file.exists(file.path("..", files))))
  chunk <- function(f) c("```{r}", "#| eval: false", "#| echo: true",
                         sprintf("#| file: \"%s\"", file.path("..", f)),
                         sprintf("#| code-summary: \"%s · %s\"", f, SCRIPTS[[f]]), "```", "")
  cat(knitr::knit_child(text = unlist(lapply(files, chunk)), quiet = TRUE), sep = "\n")
}

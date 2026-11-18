# =============================================================================
# Reproduces every table, figure and number of the PLOS ONE submission.
# Run from the repository root:   Rscript analysis/run_all.R
# =============================================================================
stopifnot(file.exists("analysis/config.R"))          # must be run from the repository root
stopifnot(l10n_info()[["UTF-8"]])                    # figures use the minus sign and superscripts

scripts <- c("01_prepare_data.R", "02_descriptives.R", "03_primary_model.R", "04_sensitivity.R",
             "05_model_comparison.R", "06_shbg_as_outcome.R", "07_figures.R",
             "08_audit_first_submission.R", "09_export_results.R", "10_check_formulas.R")
for (s in scripts) {
  message("== ", s)
  # each script runs in its own environment so that no object leaks between steps
  sys.source(file.path("analysis", "R", s), envir = new.env())
}
message("Done. Output in analysis/output/")

# Analysis for the PLOS ONE submission

Reproduces every table, figure and number of

> Sex hormone-binding globulin and skeletal muscle mass in midlife women: a
> cross-sectional analysis of the Study of Women's Health Across the Nation (SWAN)

from the public-use file of SWAN visit 10 (ICPSR study 32961, version 2).

## How to run

1. Obtain the visit 10 dataset from ICPSR (<https://doi.org/10.3886/ICPSR32961.v2>)
   and place the Stata file at `ICPSR_32961/DS0001/32961-0001-Data.dta`.
2. From the repository root, in a UTF-8 locale:

   ```sh
   Rscript analysis/run_all.R
   ```

R packages: `readstata13`, `dplyr`, `tidyr`, `car`, `sandwich`, `ggplot2`,
`patchwork`, `jsonlite` (the versions used are listed in
`output/session_info.txt`; R 4.3.3). The run takes under a minute.

The pipeline stops with an error if a formula of Text A in S1 File does not
reproduce the value computed by R (`R/10_check_formulas.R`), or if either of two
regression tests fails:

- the primary model must reproduce n = 1,693, AIC = 7,309.813 and a coefficient
  of log2(SHBG) of -0.004054;
- the initial model must reproduce n = 1,693, AIC = 8,383.171 and a coefficient
  of SHBG of -0.001069, the values of the model reported in the first
  submission (`SKELMMM10` in `docs/exploratory.qmd`).

## Scripts

| Script | What it does | Output used in the article |
|---|---|---|
| `config.R` | Paths, reporting units, SI conversion factors, reference categories, expected values of the primary and initial models | |
| `R/functions.R` | Definition of the primary and initial models; shared helpers: model fitting, ratios of means with confidence intervals, standardized means | |
| `R/01_prepare_data.R` | Reads the raw file (never modified), recodes, derives BMI, fat mass and muscle indices, checks the Janssen equation, defines the analytic sample, counts missing values by variable | participant flow; Text A in S1 File |
| `R/02_descriptives.R` | Characteristics by SHBG tertile, included vs not included, by menopausal status with multiplicity-adjusted comparisons, Spearman correlations, mean muscle mass by tertile of SHBG within tertiles of waist circumference | Table 1; Fig 1; Tables A to C and J in S1 File |
| `R/03_primary_model.R` | Primary model (gamma, log link, log2 SHBG), initial model, joint tests, adjusted means, diagnostics of both models | Table 2; Tables E and H in S1 File |
| `R/04_sensitivity.R` | Sequential adjustment, one covariate removed at a time, other measures of body size, alternative scaling of muscle mass, model specification, sample restrictions, shape of the association, full regression output and fit measures | Fig 2; Tables D and I in S1 File |
| `R/05_model_comparison.R` | Gamma and Gaussian families with log, identity and inverse links on the same sample; inverse Gaussian and log-normal models; modified Park test of the variance function | Table F in S1 File |
| `R/06_shbg_as_outcome.R` | Complementary model with SHBG as the dependent variable, with the covariates of the primary model (ratios, joint tests, coefficients, fit measures) and of the initial model | Table 3; Tables G and I in S1 File |
| `R/07_figures.R` | Figures as 600 dpi TIFF (LZW) with PNG previews | Fig 1, Fig 2; Figs A to D in S1 File |
| `R/08_audit_first_submission.R` | Recomputes the statistics quoted in the first submission and flags disagreements | change log |
| `R/09_export_results.R` | Collects every number quoted in the text into `output/results.json`; writes the session information | manuscript text |
| `R/10_check_formulas.R` | Computes by hand the formulas written in Text A of S1 File (dispersion, covariance matrix, confidence interval, joint F test, standardized mean, HC3, spline test, AIC, leverage, residuals, Cook's distance) and stops if any differs from what R returns | Text A in S1 File |

## Output

- `output/tables/*.csv`: one file per table (full precision).
- `output/figures/*.tif`, `*.png`: `Fig1`, `Fig2` and `S1_Fig` to `S4_Fig` (Figs A
  to D in S1 File). `Fig1` shows the unadjusted data: mean muscle mass by tertile
  of SHBG, in all women and within tertiles of waist circumference, and the
  distribution of the tertiles of waist circumference in each tertile of SHBG.
  `Fig2` has three panels: adjusted means across SHBG in the primary model with
  and without waist circumference, the estimate when one covariate is removed,
  and the estimate when body size is measured in other ways. `Fig2_overview` is
  its first panel alone, used in the report. `S4_Fig` shows sequential
  adjustment, scaling, specifications and restricted samples. `Fig1` and `Fig2`
  are copied to `docs/figures/` by `docs/render_report.sh`.
- `output/results.json`: the numbers quoted in the manuscript.
- `output/session_info.txt`: R and package versions.
- `output/data/`: individual-level derived data, not under version control.

## Design notes

- **Primary model.** Gamma family, log link, skeletal muscle mass on log2(SHBG),
  adjusted for age, height, waist circumference, race/ethnicity, menopausal
  status, log2 of testosterone, estradiol, FSH, DHEA-S and TSH, activity
  limitation and exercise (22 coefficients). It is defined once, in
  `R/functions.R` (`PRIMARY_EXPOSURE`, `PRIMARY_COVARIATES`).
- **Initial model.** The model of the first submission: SHBG and the hormones on
  their original scale, age, percent body fat, menopausal status, activity
  limitation and exercise (`SUBMITTED_EXPOSURE`, `SUBMITTED_COVARIATES`). Only
  the reference category of menopausal status differs from the submitted code
  (natural post-menopause instead of surgical menopause), which leaves the fit
  and the SHBG coefficient unchanged. It is reported in full and is model C7 of
  the sensitivity analyses.
- **The primary model was revised after the initial model had been fitted**
  (17 November 2026). The reasons are in the report (*Model definition*) and in
  Text A of S1 File; the steps that led to it are in `proposal/`.
- **Reporting.** Ratios of means per doubling for SHBG and the hormones,
  exp(beta), and exp(beta x unit) for the other covariates, with Wald 95% confidence
  intervals based on the t distribution with the residual degrees of freedom.
  Units are set in `config.R`.
- **Units.** Hormones are analyzed in the units of the SWAN file. The SI
  conversion factors used for Table 1 and the text are in `config.R`.
- **STROBE.** `tableS_missing_by_variable.csv`, `tableS_full_regression_output.csv`,
  `tableS_fit_measures.csv`, `tableS_diagnostics_initial.csv` and the three
  `tableS_shbg_outcome_*.csv` files support the reporting items of the STROBE
  checklist (S2 File) and of the PLOS ONE statistical reporting guidelines.
- **Same sample for every model.** All sensitivity models use the 1,693 women
  of the analytic sample unless a restriction is the point of the analysis.
- **No correction for multiple testing.** The sensitivity analyses describe the
  robustness of a single association.

## Relation to `docs/exploratory.qmd`

`docs/exploratory.qmd` (named `docs/index.qmd` until November 2026) is the exploratory analysis behind the first submission and is
kept unchanged. Three points about it matter when comparing results:

1. The candidate models compared there were fitted to slightly different samples
   (1,693 or 1,694 women), so their AIC and BIC weights are not comparable. The
   gamma model with inverse link was not fitted; with starting values it
   converges (see `R/05_model_comparison.R`).
2. The pairwise comparisons by menopausal status were displayed with unadjusted
   P values. `R/02_descriptives.R` reports Bonferroni and Holm adjusted values.
3. The chunks labelled k-fold cross-validation fit the model on the held-out
   fold and correlate its own fitted values with the observations, so they do
   not measure out-of-sample performance. They are not used in the article.

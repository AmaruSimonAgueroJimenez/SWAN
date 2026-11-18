# SHBG and skeletal muscle mass in midlife women (SWAN)

Analysis code, report and manuscript files for the study of sex hormone-binding
globulin (SHBG) and skeletal muscle mass in midlife women, using the public-use
file of the 10th visit (2006 to 2008) of the Study of Women's Health Across the
Nation (SWAN; ICPSR study 32961).

Authors: Carla Basualto-Alarcón, Amaru Agüero, Manuel Estrada, Luis Toro.

## Start here

| I want to | Go to |
|---|---|
| Read the results with tables, figures and the model equations | The analysis report: <https://amarusimonaguerojimenez.github.io/SWAN/> (source in [`docs/`](docs/)) |
| See how the models are defined mathematically | [Model definition](https://amarusimonaguerojimenez.github.io/SWAN/models.html) in the report |
| Rerun the analysis | [`analysis/`](analysis/) and the section *Reproduce* below |
| Find the manuscript files | `manuscripts/`, one folder per version, named by date (working files, not in the public repository for now) |
| Know which script produces a given table or figure | [`analysis/README.md`](analysis/README.md) |
| See the code behind a table or figure | Folded on each page of the report (*Code*); every script of the pipeline is on the page [Reproducibility](https://amarusimonaguerojimenez.github.io/SWAN/reproducibility.html#code-of-the-pipeline) |

## Repository map

```
SWAN/
├── analysis/              Reproducible pipeline in R
│   ├── run_all.R            runs everything
│   ├── config.R             paths, units, reference categories, regression tests
│   ├── R/                   functions.R, 01_prepare_data.R ... 10_check_formulas.R
│   ├── output/              tables (CSV), figures (TIFF, PNG), results.json
│   └── proposal/            record of how the primary model was arrived at
├── docs/                  Analysis report (HTML pages served by GitHub Pages)
│   ├── index.qmd ...        overview, data, model definition, results, sensitivity,
│   │                        diagnostics, reproducibility; the code is folded on each page
│   ├── render_report.sh     renders the pages from analysis/output
│   ├── figures/             Fig 1 and Fig 2 of the manuscript (PNG and TIFF)
│   └── exploratory.qmd      exploratory analysis of 2025, kept unchanged
├── manuscripts/           Manuscript versions, one dated folder each (working files, not in the
│   └── 2026-11-17_PLOS_ONE/   public repository for now)
├── ICPSR_32961/           Place of the SWAN visit 10 dataset and its documentation
└── LICENSE
```

## Reproduce

```sh
Rscript analysis/run_all.R      # tables, figures and results.json in analysis/output/
bash docs/render_report.sh      # HTML report in docs/ (optional; needs Quarto)
```

Requirements: R 4.3 or later with the packages `readstata13`, `dplyr`, `tidyr`,
`car`, `sandwich`, `ggplot2`, `patchwork` and `jsonlite` (plus `knitr` for the
report), and the SWAN visit 10 file in `ICPSR_32961/DS0001/`.

The pipeline stops if the primary model or the initial model (the model of the
first submission) does not reproduce its expected sample size, AIC and SHBG
coefficient. The report and the manuscript read every number from
`analysis/output/`, so they cannot drift apart from the code.

## What the analysis shows

**Question.** Higher SHBG has been related to lower muscle mass in
post-menopausal women. SHBG falls with adiposity and absolute muscle mass rises
with body size, so it was unclear whether that association is independent of
body size.

**Answer.** In 1,693 women aged 51 to 63 years, each doubling of SHBG was
associated with 4.3% lower skeletal muscle mass without adjustment and with a
difference of -0.4% (95% CI -1.05 to 0.25) in the primary model, a gamma
regression with log link adjusted for age, height, waist circumference,
race/ethnicity, menopausal status, five other hormones, activity limitation and
exercise. Waist circumference was the only covariate whose removal changed the
estimate appreciably (-3.5% without it). The result was the same with body mass
index or weight as the measure of body size, but not when fat mass from the same
impedance measurement replaced waist circumference (-2.5%). Read in the other
direction, with SHBG as the dependent variable, SHBG was related to waist
circumference, estradiol and DHEA-S and not to muscle mass. SHBG appears to
reflect adiposity rather than muscle mass. Details, with confidence intervals,
are in the report.

The primary model was revised on 17 November 2026, after the results of the
initial model were known (-1.06% per 10 nmol/L, with percent body fat as the
only measure of body composition). Both models are reported in full.

## Data

Sutton-Tyrrell K, Selzer F, Sowers M, Finkelstein J, Powell L, Gold E, et al.
Study of Women's Health Across the Nation (SWAN), 2006-2008: Visit 10 Dataset.
Ann Arbor, MI: Inter-university Consortium for Political and Social Research;
2018. <https://doi.org/10.3886/ICPSR32961.v2>

The data are distributed by ICPSR under its terms of use and must be obtained
from ICPSR.

## License

See `LICENSE`.

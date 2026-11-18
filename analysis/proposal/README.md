# Proposal: revised primary model (adopted)

**Status: adopted on 17 November 2026.** The model proposed here is now the
primary model of the manuscript and of the pipeline (`analysis/R/functions.R`,
`analysis/R/03_primary_model.R`). This folder is kept as the record of how the
model was arrived at. `run_all.R` does not call the script, and the numbers of
the manuscript come from `analysis/output/`, not from here.

In this folder "submitted model" means the model of the first submission, which
the pipeline and the manuscript now call the initial model.

`primary_model_proposal.R` fits a candidate for a revised primary model on the
same 1,693 women and writes four tables to `output/`:

| File | Content |
|---|---|
| `path_to_proposed_model.csv` | The SHBG estimate after each change, from the submitted model to the proposed one |
| `alternatives_to_proposed_model.csv` | The proposed model with other measures of body size, other covariate sets and restricted samples |
| `proposed_model_coefficients.csv` | Coefficients, standard errors, t statistics and ratios of means of the proposed model |
| `proposed_model_diagnostics.csv` | Fit, shape of the association, influence, collinearity, variance function, other families |

Proposed model: gamma regression with log link of skeletal muscle mass on
log2(SHBG), adjusted for age, height, waist circumference, race/ethnicity,
menopausal status, log2 of testosterone, estradiol, FSH, DHEA-S and TSH,
limitation of vigorous activities and exercise.

Run from the repository root, after `Rscript analysis/run_all.R`:

```sh
Rscript analysis/proposal/primary_model_proposal.R
```

The script checks that its first model reproduces the initial model and that
its last model reproduces the primary model of the pipeline.

The body-size sensitivity results were already known when this proposal was
written, so the choice is post hoc. The manuscript says so in the Methods and in
the limitations.

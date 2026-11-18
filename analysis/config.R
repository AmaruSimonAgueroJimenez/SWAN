# =============================================================================
# Configuration for the SHBG / skeletal muscle mass analysis (PLOS ONE revision)
# Single source of truth for paths, units and analysis settings.
# All scripts are run from the repository root (see analysis/run_all.R).
# =============================================================================

CFG <- list(
  # ---- paths ---------------------------------------------------------------
  data_file = file.path("ICPSR_32961", "DS0001", "32961-0001-Data.dta"),
  out_dir   = file.path("analysis", "output"),

  # ---- inference -----------------------------------------------------------
  conf_level = 0.95,
  seed       = 20261117L,   # only used for the randomised quantile residuals

  # ---- reporting units -----------------------------------------------------
  # Associations are reported as ratios of means per `unit` of each covariate.
  unit = c(
    shbg         = 10,   # nmol/L
    testosterone = 10,   # ng/dL
    fsh          = 10,   # mIU/mL
    estradiol    = 10,   # pg/mL
    tsh          = 1,    # uIU/mL
    dheas        = 50,   # ug/dL
    age          = 1,    # years
    pbf          = 5,    # percentage points
    bmi          = 1,    # kg/m2
    height_cm    = 5,    # cm
    weight_kg    = 5,    # kg
    fat_mass_kg  = 5,    # kg
    waist_cm     = 5,    # cm
    smm          = 1     # kg (used when SHBG is the dependent variable)
  ),

  # ---- SI conversion (dataset unit x factor = SI unit) ------------------------
  # Factors follow from the molar masses (4 significant digits): testosterone 288.42 g/mol,
  # estradiol 272.38 g/mol, DHEA-S 368.49 g/mol. FSH (mIU/mL = IU/L) and
  # TSH (uIU/mL = mIU/L) are numerically identical in SI units.
  si = list(
    testosterone = list(factor = 0.03467, unit = "nmol/L", digits = 2L),
    estradiol    = list(factor = 3.671,  unit = "pmol/L", digits = 1L),
    dheas        = list(factor = 0.02714, unit = "umol/L", digits = 2L)
  ),

  # ---- reference categories --------------------------------------------------
  ref_status   = "Natural post-menopause",
  ref_activity = "Not limited",
  ref_exercise = "No",
  ref_race     = "White",

  # ---- regression tests -------------------------------------------------------
  # Primary model (revised in November 2026): gamma log-link model of skeletal
  # muscle mass on log2(SHBG), adjusted for age, height, waist circumference,
  # race/ethnicity, menopausal status, log2 hormones, activity limitation and
  # exercise. The pipeline stops if these values are not reproduced.
  expected_primary = list(n = 1693L, aic = 7309.813, beta_log2_shbg = -0.004054),
  # Model of the first submission (SKELMMM10 in docs/exploratory.qmd and
  # docs/resultados_modelos_SKELMMM.xlsx), now a sensitivity analysis.
  expected_submitted = list(n = 1693L, aic = 8383.171, beta_shbg = -0.001069)
)

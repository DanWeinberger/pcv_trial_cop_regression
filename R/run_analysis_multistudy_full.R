# =====================================================================
# Full multi-study COP model (cop_eiv_model_multistudy.jags): PCV7 serotypes
# AND the PCV13-ADDITIONAL serotypes (1, 3, 6A, 7F, 19A), pooled into one fit.
#
# Extends R/run_analysis_multistudy.R (PCV7-only demo) by adding the WISSPAR
# head-to-head trials as BOTH new outcome studies and new immunogenicity
# sources, so the shared serotype-level hierarchy (b1[s], x_u[s], x_i[s]) and
# the study x serotype intercept a[k,s] now span all 12 serotypes across 7
# PCV7 + 13-valent studies.
#
#   OUTCOME (ragged, study x serotype intercept a[k,s]):
#     whitney      -- Whitney IPD (NCKP surveillance, US), PCV7
#     andrews      -- Andrews 2019 PCV7 IPD (England & Wales)
#     vdl_pcv7     -- van der Linden 2016 PCV7 IPD (Germany)
#     andrews_pcv13-- Andrews 2019 PCV13-additional IPD (England & Wales)
#     vdl_pcv13    -- van der Linden 2016 PCV13-additional IPD (Germany)
#
#   IMMUNOGENICITY (pooled onto ONE shared latent per serotype/arm):
#     NCKP, Am_Indian (Navajo), South_Africa           -- PCV7 serotypes
#     NCT00366340 (Germany), NCT00688870 (Taiwan),
#     NCT00689351 (South Korea)  [WISSPAR head-to-head] -- PCV13-additional
#     serotypes (NCT00205803 excluded, as in R/config.R: WISSPAR export
#     carries no 95% CI for it, and the EIV model requires that precision).
#
# NOTE -- disconnected serotype sets: the 7 PCV7 serotypes and the 5
# PCV13-additional serotypes never co-occur in the same outcome study. This
# is no longer a "connectedness" problem the way it was under the old
# a[s] + delta[k] design (there is no reference study / fixed effect that
# needs a bridging serotype to be identified) -- every a[k,s] is estimated
# from its own study/serotype's data, shrunk toward the shared mu_a/tau_a
# hyperprior. In practice this means the PCV7 and PCV13-additional blocks'
# baselines (a[k,s]) are linked only INDIRECTLY, through that shared
# hyperprior, not through direct data linkage. The serotype slope b1[s] --
# the actual correlate-of-protection quantity -- is estimated the same way it
# already is in the single-study RE model (hierarchical shrinkage across
# serotypes, each serotype informed by only 2 arms per study); check
# diagnostics.pdf / posterior_summary.csv for wide credible intervals or poor
# mixing on a[]/b1[] before over-interpreting those specific parameters.
# =====================================================================

source(file.path("R", "config.R"))
source(file.path("R", "cop_model.R"))

outcome_sources <- list(
  list(data_file = DEFAULT_DATA_FILE,   study_id = "whitney",       label = OUTCOME_WHITNEY),
  list(data_file = ANDREWS_DATA_FILE,   study_id = "andrews",       label = OUTCOME_ANDREWS),
  list(data_file = VDL_PCV7_DATA_FILE,  study_id = "vdl_pcv7",      label = OUTCOME_VDL_PCV7),
  list(data_file = WISSPAR_DATA_FILE,   study_id = "andrews_pcv13", label = OUTCOME_WISSPAR),
  list(data_file = VDL_PCV13_DATA_FILE, study_id = "vdl_pcv13",     label = OUTCOME_VDL_PCV13)
)

immuno_sources <- list(
  list(data_file = DEFAULT_DATA_FILE, study = "NCKP",         label = "NCKP (US, Whitney/Kaiser)"),
  list(data_file = DEFAULT_DATA_FILE, study = "Am_Indian",    label = "Navajo / American Indian"),
  list(data_file = DEFAULT_DATA_FILE, study = "South_Africa", label = "South Africa"),
  list(data_file = WISSPAR_DATA_FILE, study = "NCT00366340",  label = "WISSPAR NCT00366340 (Germany)"),
  list(data_file = WISSPAR_DATA_FILE, study = "NCT00688870",  label = "WISSPAR NCT00688870 (Taiwan)"),
  list(data_file = WISSPAR_DATA_FILE, study = "NCT00689351",  label = "WISSPAR NCT00689351 (South Korea)")
)

out_dir <- file.path("results", "multistudy_full_demo")

prep <- prepare_cop_data_multistudy(outcome_sources, immuno_sources)
fit  <- fit_cop_multistudy(prep, out_dir)

slope_summary(fit$samp, prep,
             analysis_id     = "multistudy_full_demo",
             predictor_label = "Pooled NCKP + Navajo + South Africa + WISSPAR (DE/TW/KR) immunogenicity",
             outcome_label   = "Pooled Whitney + Andrews + van der Linden PCV7 and PCV13-additional IPD",
             out_dir         = out_dir)

study_summary(fit$samp, prep, out_dir = out_dir)

plot_cop_absolute_multistudy(prep, fit$samp, out_dir,
                             title_suffix = "PCV7 (NCKP/Navajo/S.Africa) + PCV13-additional (WISSPAR DE/TW/KR) immunogenicity")

cat(sprintf("\nDone (max PSRF %.3f, min ESS %d). Outputs in %s\n",
            fit$max_psrf, round(fit$min_ess), out_dir))

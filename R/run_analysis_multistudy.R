# =====================================================================
# Example driver: multi-study COP model (cop_eiv_model_multistudy.jags).
#
# Demonstrates the generalization in R/cop_model.R over the existing PCV7
# data already used by the single-study analyses in R/config.R:
#
#   OUTCOME (ragged, study-specific fixed effect delta[k]):
#     whitney  -- Whitney IPD (NCKP surveillance, US)          [reference]
#     andrews  -- Andrews 2019 PCV7 IPD (England & Wales)
#     vdl_pcv7 -- van der Linden 2016 PCV7 IPD (Germany)
#   All three report the same 7 PCV7 serotypes (4, 6B, 9V, 14, 18C, 19F,
#   23F), so this genuinely exercises the ragged/study-fixed-effect design
#   even though, here, the serotype sets happen to coincide.
#
#   IMMUNOGENICITY (pooled onto ONE shared latent per serotype/arm):
#     NCKP, Am_Indian (Navajo), South_Africa -- three DIFFERENT populations/
#     studies each estimating the (assumed shared) serotype-specific GMC
#     distribution for this vaccine. In the single-study analyses these were
#     three SEPARATE analyses (nckp/navajo/south_africa); here they are
#     combined into one pooled measurement-error predictor.
#
# This is one self-contained example, not wired into the ANALYSES/
# COMPARISONS registry in R/config.R -- that registry's shape (one outcome +
# one predictor per analysis) doesn't match the multi-study design. Treat
# this script as a template: swap outcome_sources/immuno_sources for a
# different set of studies to explore other combinations.
# =====================================================================

source(file.path("R", "config.R"))
source(file.path("R", "cop_model.R"))

outcome_sources <- list(
  list(data_file = DEFAULT_DATA_FILE,   study_id = "whitney",  label = OUTCOME_WHITNEY),
  list(data_file = ANDREWS_DATA_FILE,   study_id = "andrews",  label = OUTCOME_ANDREWS),
  list(data_file = VDL_PCV7_DATA_FILE,  study_id = "vdl_pcv7", label = OUTCOME_VDL_PCV7)
)

immuno_sources <- list(
  list(data_file = DEFAULT_DATA_FILE, study = "NCKP",         label = "NCKP (US, Whitney/Kaiser)"),
  list(data_file = DEFAULT_DATA_FILE, study = "Am_Indian",    label = "Navajo / American Indian"),
  list(data_file = DEFAULT_DATA_FILE, study = "South_Africa", label = "South Africa")
)

out_dir <- file.path("results", "multistudy_pcv7_demo")

prep <- prepare_cop_data_multistudy(outcome_sources, immuno_sources)
fit  <- fit_cop_multistudy(prep, out_dir)

slope_summary(fit$samp, prep,
             analysis_id     = "multistudy_pcv7_demo",
             predictor_label = "Pooled NCKP + Navajo + South Africa immunogenicity",
             outcome_label   = "Pooled Whitney + Andrews + van der Linden PCV7 IPD",
             out_dir         = out_dir)

study_summary(fit$samp, prep, out_dir = out_dir)

plot_cop_absolute_multistudy(prep, fit$samp, out_dir,
                             title_suffix = "Pooled NCKP + Navajo + South Africa immunogenicity")

cat(sprintf("\nDone (max PSRF %.3f, min ESS %d). Outputs in %s\n",
            fit$max_psrf, round(fit$min_ess), out_dir))

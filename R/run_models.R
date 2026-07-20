# =====================================================================
# Driver: fit one or more model specs from R/config.R's MODEL_SPECS
# registry (cop_eiv_model_multistudy.jags -- the only model this project
# carries forward).
#
# Usage (from the project root):
#   Rscript R/run_models.R                          # fit every spec in the registry
#   Rscript R/run_models.R all                      # same
#   Rscript R/run_models.R pooled                   # fit a spec by id (currently the only one)
#   Rscript R/run_models.R pooled predictor_error=sd   # population-SD EIV
#
# The optional predictor_error=<se|sd> token selects the EIV predictor error
# (see prepare_gmc_predictor() in R/cop_model.R): "se" (default) is the SE of
# the mean log-GMC from the reported 95% CI; "sd" scales it to the population
# SD of individual log-titers (SE * sqrt(N_immuno)), which requires every
# immuno source in the spec to carry an N_immuno column (currently only the
# WISSPAR sources do). A "sd" run writes to results/<id>_sd/ so it never
# overwrites the default "se" results.
#
# For each spec <id> it writes results/<id>/:
#   posterior_summary.csv  mcmc.rds  diagnostics.pdf  waic.csv
#   slope_summary.csv  study_summary.csv
#   cop_scatter_absolute_risk_gmc_multistudy.png
#
# All modelling logic lives in R/cop_model.R; the registry in R/config.R
# decides WHICH specs exist. This script just wires them together.
# =====================================================================

source(file.path("R", "config.R"))
source(file.path("R", "cop_model.R"))

args <- commandArgs(trailingOnly = TRUE)

perr_arg <- grep("^predictor_error=", args, value = TRUE)
predictor_error <- if (length(perr_arg)) sub("^predictor_error=", "", tail(perr_arg, 1)) else "se"
ids_args <- setdiff(args, perr_arg)
ids <- if (length(ids_args) == 0 || identical(ids_args, "all")) {
  names(MODEL_SPECS)
} else ids_args

invisible(lapply(ids, function(id) {
  spec    <- get_model_spec(id)
  out_dir <- spec_out_dir(id, predictor_error)

  cat("\n=====================================================================\n")
  cat(sprintf("Model spec '%s'  [predictor_error: %s]\n", id, predictor_error))
  cat(sprintf("  outcomes : %s\n", paste(spec$outcomes, collapse = ", ")))
  cat(sprintf("  immuno   : %s\n", paste(spec$immuno, collapse = ", ")))
  cat(sprintf("  output   : %s\n", out_dir))
  cat("=====================================================================\n")

  prep <- prepare_cop_data_multistudy(spec_outcome_sources(spec), spec_immuno_sources(spec),
                                      predictor_error = predictor_error)
  fit  <- fit_cop_multistudy(prep, out_dir)

  predictor_label <- paste(vapply(spec$immuno, function(x) get_immuno_source(x)$label, character(1)),
                          collapse = " + ")
  outcome_label   <- paste(vapply(spec$outcomes, function(x) get_outcome_source(x)$label, character(1)),
                          collapse = " + ")
  slope_summary(fit$samp, prep, id, predictor_label, outcome_label, out_dir = out_dir)
  study_summary(fit$samp, prep, out_dir = out_dir)
  plot_cop_absolute_multistudy(prep, fit$samp, out_dir, title_suffix = sprintf("Predictor: %s", predictor_label))

  cat(sprintf("Done '%s' (max PSRF %.3f, min ESS %d, WAIC %.2f).\n",
              id, fit$max_psrf, round(fit$min_ess), fit$waic$waic))
}))

cat("\nAll requested model specs complete.\n")

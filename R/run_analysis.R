# =====================================================================
# Driver: fit + plot + summarise one or more COP analyses.
#
# Usage (from the project root):
#   Rscript R/run_analysis.R                 # run every analysis in the registry
#   Rscript R/run_analysis.R all             # same
#   Rscript R/run_analysis.R nckp navajo     # run selected analyses by id
#   Rscript R/run_analysis.R param=ratio     # use the ratio parameterization
#
# The optional param=<centered|ratio> token selects the JAGS parameterization
# (default centered). Both are algebraically equivalent and produce canonically-
# named posteriors, so downstream plots/comparisons are identical either way.
#
# For each analysis <id> it writes results/<id>/:
#   posterior_summary.csv  mcmc.rds  diagnostics.pdf
#   cop_fit_proportion.pdf (+ .png)  cop_scatter_gmr_rr.pdf (+ .png)
#   slope_summary.csv
#
# All modelling logic lives in R/cop_model.R; the registry in R/config.R
# decides WHICH analyses exist. This script just wires them together.
# =====================================================================

source(file.path("R", "config.R"))
source(file.path("R", "cop_model.R"))

args <- commandArgs(trailingOnly = TRUE)

# Optional param=<centered|ratio> token; the rest are analysis ids.
param_arg <- grep("^param=", args, value = TRUE)
parameterization <- if (length(param_arg)) {
  sub("^param=", "", tail(param_arg, 1))
} else "centered"
ids_args <- setdiff(args, param_arg)
ids <- if (length(ids_args) == 0 || identical(ids_args, "all")) {
  names(ANALYSES)
} else ids_args

invisible(lapply(ids, function(id) {
  cfg     <- get_analysis(id)
  out_dir <- analysis_out_dir(id)

  cat("\n=====================================================================\n")
  cat(sprintf("Analysis '%s'  [parameterization: %s]\n", id, parameterization))
  cat(sprintf("  predictor: %s  (Study == '%s')\n", cfg$predictor_label, cfg$predictor_study))
  cat(sprintf("  outcome  : %s\n", cfg$outcome_label))
  cat(sprintf("  output   : %s\n", out_dir))
  cat("=====================================================================\n")

  prep <- prepare_cop_data(cfg$data_file, cfg$predictor_study)
  fit  <- fit_cop(prep, out_dir, parameterization = parameterization)
  plot_cop(prep, fit$samp, out_dir,
           title_suffix = sprintf("Predictor: %s", cfg$predictor_label))
  plot_cop_scatter(prep, fit$samp, out_dir,
                   title_suffix = sprintf("Predictor: %s", cfg$predictor_label))
  slope_summary(fit$samp, prep, id, cfg$predictor_label, cfg$outcome_label,
                out_dir = out_dir)

  cat(sprintf("Done '%s' (max PSRF %.3f, min ESS %d).\n",
              id, fit$max_psrf, round(fit$min_ess)))
}))

cat("\nAll requested analyses complete.\n")

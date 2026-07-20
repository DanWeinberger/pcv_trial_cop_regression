# =====================================================================
# Independent check of the CI -> SE -> population-SD logic in
# R/cop_model.R (ci_se(), log_se(), and the SE*sqrt(N) step used by
# prepare_gmc_predictor(predictor_error = "sd")).
#
# Data (siber_nckp_table5_table6.csv): Siber et al. 2007, Vaccine 25:3816-26,
# NCKP 7vPnC-immunized infants, n = 188, "without 22F absorption" columns
# (the single-absorption ELISA, matching the assay behind siber_table3_tidy.csv):
#   GMC/Lower_CL/Upper_CL          -> Table 5
#   pct_responders_ge_0.35         -> Table 6 (% subjects with antibody
#                                     >= 0.35 ug/mL, the same n = 188 sample)
#
# Method: for the SAME sample, estimate the population SD of log-titer two
# independent ways and compare.
#   Way 1 (the code's own logic): SE of the mean log-GMC from the reported
#     95% CI (log_se()), scaled to a population SD via SE * sqrt(N) -- this
#     is exactly the predictor_error = "sd" path in prepare_gmc_predictor().
#   Way 2 (independent): back out the population SD from GMC + % responders
#     >= 0.35 ug/mL, assuming log-normal titers:
#       pct_responders = P(titer >= 0.35) = 1 - Phi((log(0.35) - log(GMC)) / SD)
#       => SD = (log(0.35) - log(GMC)) / qnorm(1 - pct_responders)
#
# Usage: Rscript validation/validate_ci_se.R   (run from the project root)
# =====================================================================

source(file.path("R", "cop_model.R"))

d <- read.csv(file.path("validation", "siber_nckp_table5_table6.csv"),
             stringsAsFactors = FALSE)

d$se_logGMC <- log_se(d$Lower_CL, d$Upper_CL)
d$sd_code   <- d$se_logGMC * sqrt(d$N_immuno)

thresh <- log(0.35)
d$z          <- qnorm(1 - d$pct_responders_ge_0.35 / 100)
d$sd_implied <- (thresh - log(d$GMC)) / d$z

d$ratio <- d$sd_code / d$sd_implied

print(d[, c("Serotype", "GMC", "se_logGMC", "sd_code",
           "pct_responders_ge_0.35", "sd_implied", "ratio")],
     digits = 3, row.names = FALSE)

cat(sprintf("\nMedian ratio (code / independent): %.2f\n", median(d$ratio)))

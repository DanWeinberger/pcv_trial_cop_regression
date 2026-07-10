# =====================================================================
# Build data/wisspar_andrews_merged.csv
#
# Turns the WISSPAR head-to-head PCV7-vs-PCV13 immunogenicity slice into COP
# predictor sources - ONE PER TRIAL - in the project's standard merged-CSV
# schema, so each trial fits with the existing engine (R/cop_model.R) exactly
# like the SIBER predictors, and their slopes can then be contrasted with
# R/compare_slopes.R.
#
# Design:
#   * Restricted to the PCV13-ADDITIONAL serotypes (1, 3, 5, 6A, 7F, 19A),
#     where the PCV7 arm carries no antigen and is a genuine no-antigen
#     ("Unimmunized") comparator within the same trial.
#   * Arm map:  PCV13 (Pfizer) -> Immunized ,  PCV7 -> Unimmunized.
#   * NO pooling: each head-to-head trial keeps its own GMC + 95% CI and enters
#     as a separate `Study` (the trial's NCT id). The engine's own SE floor /
#     asymmetry handling deals with rounded near-zero-width CIs.
#   * Outcome = Andrews 2019 Table 3 PCV13 IPD, ">=1 dose" (England & Wales),
#     broadcast identically to every trial (Vaccinated -> Immunized counts,
#     Unvaccinated -> Unimmunized counts). This mirrors the SIBER design, where
#     the Whitney outcome is identical across the three predictor studies;
#     only the immunogenicity predictor differs by Study.
#   * A serotype is kept for a trial only if that trial measured BOTH arms AND
#     an Andrews outcome exists: outcome set is {1, 3, 6A, 7F, 19A}
#     (5 has no outcome; 6C has an outcome but no head-to-head GMC).
#
# NOTE: OPA is not usable here - the head-to-head OPA trials cover no
# PCV13-additional serotype, so these predictor sources are GMC-only. OPA stays
# in data/wisspar_pcv7_pcv13_child_postprimary.csv for descriptive use.
#
# Run from the project root (after R/build_wisspar_head2head.R):
#   Rscript R/build_wisspar_andrews_merged.R
# =====================================================================

H2H_FILE     <- file.path("data", "wisspar_pcv7_pcv13_child_postprimary.csv")
ANDREWS_FILE <- file.path("data", "andrews2019_table2_3_tidy.csv")
OUT_FILE     <- file.path("data", "wisspar_andrews_merged.csv")

PCV13_ADDITIONAL <- c("1", "3", "5", "6A", "7F", "19A")

# ---- Outcome: Andrews Table 3 PCV13, ">=1 dose" (fixed across trials) --------
a <- read.csv(ANDREWS_FILE, stringsAsFactors = FALSE)
a <- subset(a, Table == "Table3" & Vaccine == "PCV13" &
               VE_Measure == ">=1 dose" & Serotype %in% PCV13_ADDITIONAL)
arm_of <- c(Vaccinated = "Immunized", Unvaccinated = "Unimmunized")
a$Vaccine_Arm <- arm_of[a$Vaccine_Status]
stopifnot(!any(is.na(a$Vaccine_Arm)))
out_key <- paste(a$Vaccine_Arm, a$Serotype, sep = "|")
stopifnot(!any(duplicated(out_key)))
outcome_st <- unique(a$Serotype)

# ---- Predictor: per-trial WISSPAR GMC (additional serotypes, both arms) ------
h <- read.csv(H2H_FILE, stringsAsFactors = FALSE)
h <- subset(h, assay == "GMC" & serotype %in% PCV13_ADDITIONAL)
h$Vaccine_Arm <- ifelse(h$vaccine == "PCV7", "Unimmunized",
                 ifelse(h$vaccine == "PCV13 (Pfizer)", "Immunized", NA))
stopifnot(!any(is.na(h$Vaccine_Arm)))

# The EIV model needs a usable 95% CI (its measurement-error input). Some
# trials report no CI (limits recorded as 0). Drop rows without a valid CI and
# note which trials that removes.
valid_ci <- !is.na(h$lower_limit) & !is.na(h$upper_limit) &
            h$lower_limit > 0 & h$upper_limit > 0 & h$upper_limit >= h$lower_limit
no_ci_trials <- setdiff(unique(h$study_id), unique(h$study_id[valid_ci]))
h <- h[valid_ci, ]

# Attach the (broadcast) Andrews outcome by arm x serotype.
m <- match(paste(h$Vaccine_Arm, h$serotype, sep = "|"), out_key)
h <- h[!is.na(m), ]                       # keep only serotypes with an outcome
h$Cases       <- a$Cases[m[!is.na(m)]]
h$Total_Cases <- a$Total_Cases[m[!is.na(m)]]

# Keep, per trial, only serotypes measured in BOTH arms.
both <- do.call(rbind, lapply(split(h, h$study_id), function(g) {
  st <- intersect(g$serotype[g$Vaccine_Arm == "Immunized"],
                  g$serotype[g$Vaccine_Arm == "Unimmunized"])
  g[g$serotype %in% st, ]
}))

out <- data.frame(
  Study          = both$study_id,
  Vaccine_Arm    = both$Vaccine_Arm,
  Serotype       = both$serotype,
  GMC            = both$value,
  Lower_CL       = both$lower_limit,
  Upper_CL       = both$upper_limit,
  vaccine_group  = ifelse(both$Vaccine_Arm == "Immunized", "PCV13", "PCV7"),
  Vaccine_Status = ifelse(both$Vaccine_Arm == "Immunized", ">=1 dose", "PCV7 comparator"),
  Cases          = both$Cases,
  Total_Cases    = both$Total_Cases,
  stringsAsFactors = FALSE
)
so  <- suppressWarnings(as.numeric(gsub("[^0-9].*$", "", out$Serotype)))
out <- out[order(out$Study, out$Vaccine_Arm, so, out$Serotype), ]

write.csv(out, OUT_FILE, row.names = FALSE)

# ---- Report per-trial serotype coverage --------------------------------------
cov <- tapply(out$Serotype[out$Vaccine_Arm == "Immunized"],
              out$Study[out$Vaccine_Arm == "Immunized"],
              function(s) paste(sort(s), collapse = ", "))
ncov <- vapply(strsplit(cov, ", "), length, integer(1))

cat("Wrote", OUT_FILE, "with", nrow(out), "rows.\n")
cat("Outcome: Andrews 2019 Table 3 PCV13 >=1 dose; outcome serotypes:",
    paste(sort(outcome_st), collapse = ", "), "\n")
cat("Dropped serotype 5 (WISSPAR GMC but no Andrews PCV13 outcome).\n")
if (length(no_ci_trials)) {
  cat("Dropped trial(s) with no usable 95% CI:",
      paste(no_ci_trials, collapse = ", "), "\n")
}
cat("\n")
cat("Per-trial additional-serotype coverage (n serotypes usable for a COP fit):\n")
for (s in names(cov)) cat(sprintf("  %-13s (%d): %s\n", s, ncov[s], cov[s]))
cat("\nTrials with >= 3 serotypes are good COP-fit candidates.\n")

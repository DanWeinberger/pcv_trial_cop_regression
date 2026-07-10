# =====================================================================
# Build data/wisspar_vanderlinden_pcv13_merged.csv
#
# Direct analog of R/build_wisspar_andrews_merged.R, but the PCV13 IPD OUTCOME
# comes from van der Linden 2016 Table 2 (Germany, PCV13, "at least one dose")
# instead of Andrews 2019 Table 3. The predictor is unchanged: the WISSPAR
# head-to-head PCV7-vs-PCV13 GMC slice, one Study per trial.
#
# Design (mirrors the WISSPAR/Andrews build):
#   * PCV13-additional serotypes, PCV7 arm = no-antigen comparator.
#   * Arm map:  PCV13 (Pfizer) -> Immunized ,  PCV7 -> Unimmunized.
#   * Outcome broadcast identically to every trial so slopes are comparable.
#   * Outcome serotype set = {1, 3, 6A, 7F, 19A}. Serotype 5 is dropped: it has
#     a head-to-head GMC but zero IPD cases in the van der Linden PCV13 period
#     (0:0 in both arms), so it carries no outcome signal - the same serotype
#     that dropped from the WISSPAR/Andrews build (there for lack of a PCV13
#     outcome). 6C has an Andrews outcome but no head-to-head GMC and is absent
#     from van der Linden's single-serotype rows; it plays no part here.
#
# Run from the project root (after R/build_wisspar_head2head.R):
#   Rscript R/build_vanderlinden_pcv13_merged.R
# =====================================================================

H2H_FILE <- file.path("data", "wisspar_pcv7_pcv13_child_postprimary.csv")
VDL_FILE <- file.path("data", "vanderlinden_tables_tidy.csv")
OUT_FILE <- file.path("data", "wisspar_vanderlinden_pcv13_merged.csv")

PCV13_ADDITIONAL <- c("1", "3", "5", "6A", "7F", "19A")

# ---- Outcome: van der Linden Table 2 PCV13, "at least one dose" ---------------
v <- read.csv(VDL_FILE, stringsAsFactors = FALSE)
v <- subset(v, Table == "Table2" & Vaccine == "PCV13" &
               VE_Measure == "at least one dose" & Serotype %in% PCV13_ADDITIONAL)
arm_of <- c(Vaccinated = "Immunized", Unvaccinated = "Unimmunized")
v$Vaccine_Arm <- arm_of[v$Vaccine_Status]
stopifnot(!any(is.na(v$Vaccine_Arm)))

# Drop serotypes with no outcome signal (zero cases in BOTH arms, e.g. ST5): a
# zero-case arm gives an infinite log-proportion init and no COP information.
tot_by_st <- tapply(v$Cases, v$Serotype, sum)
drop_zero <- names(tot_by_st)[tot_by_st == 0]
v <- v[!v$Serotype %in% drop_zero, ]

out_key <- paste(v$Vaccine_Arm, v$Serotype, sep = "|")
stopifnot(!any(duplicated(out_key)))
outcome_st <- unique(v$Serotype)

# ---- Predictor: per-trial WISSPAR GMC (additional serotypes, both arms) -------
h <- read.csv(H2H_FILE, stringsAsFactors = FALSE)
h <- subset(h, assay == "GMC" & serotype %in% PCV13_ADDITIONAL)
h$Vaccine_Arm <- ifelse(h$vaccine == "PCV7", "Unimmunized",
                 ifelse(h$vaccine == "PCV13 (Pfizer)", "Immunized", NA))
stopifnot(!any(is.na(h$Vaccine_Arm)))

valid_ci <- !is.na(h$lower_limit) & !is.na(h$upper_limit) &
            h$lower_limit > 0 & h$upper_limit > 0 & h$upper_limit >= h$lower_limit
no_ci_trials <- setdiff(unique(h$study_id), unique(h$study_id[valid_ci]))
h <- h[valid_ci, ]

m <- match(paste(h$Vaccine_Arm, h$serotype, sep = "|"), out_key)
h <- h[!is.na(m), ]                       # keep only serotypes with an outcome
h$Cases       <- v$Cases[m[!is.na(m)]]
h$Total_Cases <- v$Total_Cases[m[!is.na(m)]]

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
  Vaccine_Status = ifelse(both$Vaccine_Arm == "Immunized", "at least one dose", "PCV7 comparator"),
  Cases          = both$Cases,
  Total_Cases    = both$Total_Cases,
  stringsAsFactors = FALSE
)
so  <- suppressWarnings(as.numeric(gsub("[^0-9].*$", "", out$Serotype)))
out <- out[order(out$Study, out$Vaccine_Arm, so, out$Serotype), ]

write.csv(out, OUT_FILE, row.names = FALSE)

cov  <- tapply(out$Serotype[out$Vaccine_Arm == "Immunized"],
               out$Study[out$Vaccine_Arm == "Immunized"],
               function(s) paste(sort(s), collapse = ", "))
ncov <- vapply(strsplit(cov, ", "), length, integer(1))

cat("Wrote", OUT_FILE, "with", nrow(out), "rows.\n")
cat("Outcome: van der Linden 2016 Table 2 PCV13 'at least one dose'; outcome serotypes:",
    paste(sort(outcome_st), collapse = ", "), "\n")
if (length(drop_zero)) {
  cat("Dropped serotype(s) with zero cases in both arms (no COP signal):",
      paste(drop_zero, collapse = ", "), "\n")
}
if (length(no_ci_trials)) {
  cat("Dropped trial(s) with no usable 95% CI:", paste(no_ci_trials, collapse = ", "), "\n")
}
cat("\nPer-trial additional-serotype coverage (n serotypes usable for a COP fit):\n")
for (s in names(cov)) cat(sprintf("  %-13s (%d): %s\n", s, ncov[s], cov[s]))

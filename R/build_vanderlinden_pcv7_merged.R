# =====================================================================
# Build data/siber_vanderlinden_pcv7_merged.csv
#
# Direct analog of R/build_andrews_merged.R, but the PCV7 IPD OUTCOME comes from
# van der Linden 2016 Table 1 (Germany, PCV7, "at least one dose") instead of
# Andrews 2019. Reuses the SIBER immunogenicity predictors (GMC + 95% CI per
# Study x arm x serotype) unchanged; only the Cases / Total_Cases are swapped.
#   Arm map:  Vaccinated -> Immunized ,  Unvaccinated -> Unimmunized.
#   Serotypes: the 7 PCV7 valencies (4, 6B, 9V, 14, 18C, 19F, 23F). 6A is in the
#   van der Linden table (cross-protection) but is not a PCV7 antigen and has no
#   SIBER predictor, so it is excluded - matching the Andrews PCV7 analysis.
#
# Run from the project root:  Rscript R/build_vanderlinden_pcv7_merged.R
# =====================================================================

siber <- read.csv(file.path("data", "siber_whitney_merged.csv"), stringsAsFactors = FALSE)
vdl   <- read.csv(file.path("data", "vanderlinden_tables_tidy.csv"), stringsAsFactors = FALSE)

# ---- van der Linden PCV7 outcome: Table 1, "at least one dose", real STs ------
v <- subset(vdl,
            Table == "Table1" &
            Vaccine == "PCV7" &
            VE_Measure == "at least one dose" &
            Serotype %in% unique(siber$Serotype))
stopifnot(nrow(v) > 0)

arm_of <- c(Vaccinated = "Immunized", Unvaccinated = "Unimmunized")
v$Vaccine_Arm <- arm_of[v$Vaccine_Status]
stopifnot(!any(is.na(v$Vaccine_Arm)))

v_key <- paste(v$Vaccine_Arm, v$Serotype, sep = "|")
stopifnot(!any(duplicated(v_key)))

# ---- Merge: keep SIBER predictor columns, swap in van der Linden counts -------
out <- siber
key <- paste(out$Vaccine_Arm, out$Serotype, sep = "|")
m   <- match(key, v_key)
if (any(is.na(m))) {
  stop("No van der Linden PCV7 outcome for: ",
       paste(unique(key[is.na(m)]), collapse = ", "))
}
out$Cases       <- v$Cases[m]
out$Total_Cases <- v$Total_Cases[m]

out_file <- file.path("data", "siber_vanderlinden_pcv7_merged.csv")
write.csv(out, out_file, row.names = FALSE)

cat("Wrote", out_file, "with", nrow(out), "rows (",
    length(unique(out$Study)), "studies x",
    length(unique(out$Serotype)), "serotypes x 2 arms).\n")
cat("van der Linden PCV7 'at least one dose' outcome (identical across predictor studies):\n")
print(unique(out[, c("Vaccine_Arm", "Serotype", "Cases", "Total_Cases")]))

# =====================================================================
# Build data/siber_andrews_merged.csv
#
# Reuses the SIBER immunogenicity predictors (GMC + 95% CI, per Study x arm x
# serotype) already assembled in siber_whitney_merged.csv, but replaces the
# OUTCOME case counts with the Andrews 2019 PCV7 trial (Table 2, ">=1 dose"
# schedule, England & Wales). Arm mapping: Vaccinated -> Immunized,
# Unvaccinated -> Unimmunized. Serotypes match the 7 PCV7 valencies exactly.
#
# Run from the project root:  Rscript R/build_andrews_merged.R
# =====================================================================

siber   <- read.csv(file.path("data", "siber_whitney_merged.csv"),
                    stringsAsFactors = FALSE)
andrews <- read.csv(file.path("data", "andrews2019_table2_3_tidy.csv"),
                    stringsAsFactors = FALSE)

# ---- Andrews outcome: PCV7 (Table 2), ">=1 dose", real serotype rows --------
a <- subset(andrews,
            Vaccine == "PCV7" &
            VE_Measure == ">=1 dose" &
            Serotype %in% unique(siber$Serotype))
stopifnot(nrow(a) > 0)

arm_of <- c(Vaccinated = "Immunized", Unvaccinated = "Unimmunized")
a$Vaccine_Arm <- arm_of[a$Vaccine_Status]
stopifnot(!any(is.na(a$Vaccine_Arm)))

# key = arm x serotype -> (Cases, Total_Cases)
a_key <- paste(a$Vaccine_Arm, a$Serotype, sep = "|")
stopifnot(!any(duplicated(a_key)))

# ---- Merge: keep SIBER predictor columns, swap in Andrews counts -------------
out <- siber
key <- paste(out$Vaccine_Arm, out$Serotype, sep = "|")
m   <- match(key, a_key)

if (any(is.na(m))) {
  stop("No Andrews outcome for: ",
       paste(unique(key[is.na(m)]), collapse = ", "))
}

out$Cases       <- a$Cases[m]
out$Total_Cases <- a$Total_Cases[m]

out_file <- file.path("data", "siber_andrews_merged.csv")
write.csv(out, out_file, row.names = FALSE)

cat("Wrote", out_file, "with", nrow(out), "rows (",
    length(unique(out$Study)), "studies x",
    length(unique(out$Serotype)), "serotypes x 2 arms).\n")
cat("Andrews PCV7 >=1 dose outcome (identical across the 3 predictor studies):\n")
print(unique(out[, c("Vaccine_Arm", "Serotype", "Cases", "Total_Cases")]))

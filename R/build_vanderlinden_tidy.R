# =====================================================================
# Build data/vanderlinden_tables_tidy.csv
#
# Extracts van der Linden et al. 2016 (PLOS ONE, doi:10.1371/journal.pone.0161257,
# "Effectiveness of Pneumococcal Conjugate Vaccines (PCV7 and PCV13) against
# Invasive Pneumococcal Disease among Children under Two Years of Age in
# Germany"):
#   Table 1 -> PCV7  VE, Germany, July 2006 - June 2010
#   Table 2 -> PCV13 VE, Germany, July 2010 - June 2015
#
# The paper is an indirect-cohort / screening-method case-control study. For each
# serotype x dose-schedule it reports, as "vaccinated:unvaccinated" pairs:
#   cases    vaccinated:unvaccinated   (IPD cases of that serotype)
#   controls vaccinated:unvaccinated   (IPD cases of non-study serotypes; the
#                                        control pool is FIXED per dose schedule)
#
# To match the Andrews 2019 tidy schema (data/andrews2019_table2_3_tidy.csv),
# which the COP engine consumes as arm-specific Cases + Total_Cases, we use the
# indirect-cohort denominator exactly as Andrews does:
#   Vaccinated   row: Cases = cases_vacc  , Total_Cases = cases_vacc  + controls_vacc
#   Unvaccinated row: Cases = cases_unvacc, Total_Cases = cases_unvacc + controls_unvacc
# so Cases/Total_Cases is the case-proportion among {this-serotype cases + controls}
# within each vaccination arm - the same quantity the Andrews rows encode.
#
# Output columns match Andrews exactly:
#   Table, Vaccine, Serotype, VE_Measure, Age_Range, Vaccine_Status, Cases, Total_Cases
#
# Run from the project root:  Rscript R/build_vanderlinden_tidy.R
# =====================================================================

OUT_FILE  <- file.path("data", "vanderlinden_tables_tidy.csv")
AGE_RANGE <- "74 to 729 days"   # "under two years of age" (paper inclusion age)

# Control pool (vaccinated:unvaccinated), fixed within each table x dose schedule.
controls <- list(
  Table1 = list("at least one dose" = c(v = 94,  u = 60),
                "post primary"       = c(v = 38,  u = 27),
                "post booster"       = c(v = 11,  u = 19)),
  Table2 = list("at least one dose" = c(v = 194, u = 43),
                "post primary"       = c(v = 74,  u = 20),
                "post booster"       = c(v = 33,  u = 16))
)

# Cases (vaccinated:unvaccinated) per row, transcribed from the PDF tables.
# Each entry: Table, Vaccine, Serotype, VE_Measure, cases_vacc, cases_unvacc.
raw <- read.csv(text = "
Table,Vaccine,Serotype,VE_Measure,cv,cu
Table1,PCV7,PCV7 serotypes + 6A,at least one dose,20,81
Table1,PCV7,PCV7 serotypes + 6A,post primary,1,49
Table1,PCV7,PCV7 serotypes + 6A,post booster,0,44
Table1,PCV7,4,at least one dose,0,1
Table1,PCV7,4,post primary,0,1
Table1,PCV7,4,post booster,0,1
Table1,PCV7,6B,at least one dose,2,20
Table1,PCV7,6B,post primary,0,15
Table1,PCV7,6B,post booster,0,14
Table1,PCV7,9V,at least one dose,0,4
Table1,PCV7,9V,post primary,0,2
Table1,PCV7,9V,post booster,0,2
Table1,PCV7,14,at least one dose,2,24
Table1,PCV7,14,post primary,0,12
Table1,PCV7,14,post booster,0,15
Table1,PCV7,18C,at least one dose,4,6
Table1,PCV7,18C,post primary,0,3
Table1,PCV7,18C,post booster,0,4
Table1,PCV7,19F,at least one dose,5,11
Table1,PCV7,19F,post primary,1,6
Table1,PCV7,19F,post booster,0,4
Table1,PCV7,23F,at least one dose,3,8
Table1,PCV7,23F,post primary,0,6
Table1,PCV7,23F,post booster,0,3
Table1,PCV7,6A,at least one dose,4,7
Table1,PCV7,6A,post primary,0,4
Table1,PCV7,6A,post booster,0,1
Table2,PCV13,PCV13 serotypes,at least one dose,25,55
Table2,PCV13,PCV13 serotypes,post primary,10,22
Table2,PCV13,PCV13 serotypes,post booster,2,13
Table2,PCV13,PCV13-non-PCV7 serotypes,at least one dose,23,43
Table2,PCV13,PCV13-non-PCV7 serotypes,post primary,10,16
Table2,PCV13,PCV13-non-PCV7 serotypes,post booster,2,12
Table2,PCV13,PCV7 serotypes in PCV13,at least one dose,2,12
Table2,PCV13,PCV7 serotypes in PCV13,post primary,0,6
Table2,PCV13,PCV7 serotypes in PCV13,post booster,0,1
Table2,PCV13,1,at least one dose,2,5
Table2,PCV13,1,post primary,1,2
Table2,PCV13,1,post booster,0,1
Table2,PCV13,3,at least one dose,6,5
Table2,PCV13,3,post primary,1,2
Table2,PCV13,3,post booster,1,2
Table2,PCV13,5,at least one dose,0,0
Table2,PCV13,5,post primary,0,0
Table2,PCV13,5,post booster,0,0
Table2,PCV13,6A,at least one dose,0,4
Table2,PCV13,6A,post primary,0,1
Table2,PCV13,6A,post booster,0,2
Table2,PCV13,7F,at least one dose,1,12
Table2,PCV13,7F,post primary,0,2
Table2,PCV13,7F,post booster,0,0
Table2,PCV13,19A,at least one dose,14,17
Table2,PCV13,19A,post primary,8,9
Table2,PCV13,19A,post booster,1,7
", stringsAsFactors = FALSE, strip.white = TRUE, colClasses = c(
  Table = "character", Vaccine = "character", Serotype = "character",
  VE_Measure = "character", cv = "integer", cu = "integer"))

# ---- Consistency check: single serotypes must sum to the aggregate rows ------
single_pcv7  <- c("4", "6B", "9V", "14", "18C", "19F", "23F", "6A")   # -> "PCV7 serotypes + 6A"
single_add   <- c("1", "3", "5", "6A", "7F", "19A")                    # -> "PCV13-non-PCV7 serotypes"
check_sum <- function(df, singles, agg, tbl) {
  for (m in unique(df$VE_Measure)) {
    s  <- df[df$Table == tbl & df$Serotype %in% singles & df$VE_Measure == m, ]
    a  <- df[df$Table == tbl & df$Serotype == agg     & df$VE_Measure == m, ]
    if (nrow(a) &&
        (sum(s$cv) != a$cv || sum(s$cu) != a$cu)) {
      stop(sprintf("%s / %s / %s: singles sum (%d:%d) != aggregate (%d:%d)",
                   tbl, agg, m, sum(s$cv), sum(s$cu), a$cv, a$cu))
    }
  }
}
check_sum(raw, single_pcv7, "PCV7 serotypes + 6A",       "Table1")
check_sum(raw, single_add,  "PCV13-non-PCV7 serotypes",  "Table2")

# ---- Expand each row into a Vaccinated + Unvaccinated pair (Andrews schema) ---
rows <- list()
for (k in seq_len(nrow(raw))) {
  r  <- raw[k, ]
  ct <- controls[[r$Table]][[r$VE_Measure]]
  rows[[length(rows) + 1]] <- data.frame(
    Table = r$Table, Vaccine = r$Vaccine, Serotype = r$Serotype,
    VE_Measure = r$VE_Measure, Age_Range = AGE_RANGE,
    Vaccine_Status = "Vaccinated", Cases = r$cv, Total_Cases = r$cv + ct[["v"]],
    stringsAsFactors = FALSE)
  rows[[length(rows) + 1]] <- data.frame(
    Table = r$Table, Vaccine = r$Vaccine, Serotype = r$Serotype,
    VE_Measure = r$VE_Measure, Age_Range = AGE_RANGE,
    Vaccine_Status = "Unvaccinated", Cases = r$cu, Total_Cases = r$cu + ct[["u"]],
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)

write.csv(out, OUT_FILE, row.names = FALSE)
cat("Wrote", OUT_FILE, "with", nrow(out), "rows (", nrow(raw), "case:control pairs x 2 arms ).\n")
cat("Internal aggregate-consistency checks passed.\n")

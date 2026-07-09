# Merge Siber correlates-of-protection GMC data with Whitney serotype-specific
# vaccine-effectiveness case data, joining on serotype and vaccine category.
#
# The two source files label vaccination status differently:
#   siber_table3_tidy.csv:        Vaccine_Arm    = "Immunized" / "Unimmunized"
#   whitney_figure_serotypes.csv: Vaccine_Status = ">=1 dose"  / "Unvaccinated"
# These are recoded to a common `vaccine_group` factor ("vaccinated" /
# "unvaccinated") before joining. Serotype codes already match between files
# (e.g. "4", "6B", "9V", "14", "18C", "19F", "23F"); siber's "Aggregated"
# pseudo-serotype has no counterpart in whitney and will not match.

library(dplyr)
library(readr)

siber <- read_csv("data/siber_table3_tidy.csv", show_col_types = FALSE) %>%
  mutate(
    vaccine_group = recode(Vaccine_Arm,
      "Immunized"   = "vaccinated",
      "Unimmunized" = "unvaccinated"
    ),
    Serotype = as.character(Serotype)
  )

whitney <- read_csv("data/whitney_figure_serotypes.csv", show_col_types = FALSE) %>%
  mutate(
    vaccine_group = recode(Vaccine_Status,
      ">=1 dose"     = "vaccinated",
      "Unvaccinated" = "unvaccinated"
    ),
    Serotype = as.character(Serotype)
  )

merged <- inner_join(
  siber, whitney,
  by = c("Serotype", "vaccine_group"),
  suffix = c("_siber", "_whitney")
)

# Rows present in one file but not matched in the other (e.g. siber's
# "Aggregated" serotype, or serotypes whitney has that siber lacks).
unmatched_siber   <- anti_join(siber,   whitney, by = c("Serotype", "vaccine_group"))
unmatched_whitney <- anti_join(whitney, siber,   by = c("Serotype", "vaccine_group"))

cat(sprintf(
  "Matched %d rows. Unmatched: %d from siber, %d from whitney.\n",
  nrow(merged), nrow(unmatched_siber), nrow(unmatched_whitney)
))

write_csv(merged, "data/siber_whitney_merged.csv")

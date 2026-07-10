# =====================================================================
# Build data/wisspar_pcv7_pcv13_child_postprimary.csv
#
# Imports GMC and OPA immunogenicity from the WISSPAR database and extracts
# the FIRST-PASS slice: head-to-head PCV7 vs PCV13 (Pfizer) trials in children
# at the post-primary timepoint.
#
# Source: WISSPAR_v2 build-time export (flat array of per-arm serotype rows),
#   https://github.com/PneumococcalCapsules/WISSPAR_v2  (data/wisspar_export.json)
#   The upstream export already applies the WISSPAR transforms
#   (Wyeth->Pfizer, IgG->GMC, PCV13->"PCV13 (Pfizer)", etc.).
#
# "Head-to-head" = a single trial (study_id) that measured BOTH a PCV7 arm and
# a PCV13 (Pfizer) arm in the same population/timepoint, so the two products are
# a naturally paired within-study comparison (cf. WISSPAR app/src/arms.ts).
#
# Run from the project root:  Rscript R/build_wisspar_head2head.R
# =====================================================================

suppressWarnings(suppressMessages(library(jsonlite)))

# ---- Config ---------------------------------------------------------
EXPORT_URL <- paste0("https://raw.githubusercontent.com/",
                     "PneumococcalCapsules/WISSPAR_v2/main/data/wisspar_export.json")
EXPORT_CACHE <- file.path("data", "wisspar_export.json")   # local raw snapshot
OUT_FILE     <- file.path("data", "wisspar_pcv7_pcv13_child_postprimary.csv")

REF_VACCINE  <- "PCV7"
COMP_VACCINE <- "PCV13 (Pfizer)"
# Post-primary timepoints in children (dose_description as coded by WISSPAR).
POSTPRIMARY  <- c("1m post primary child", "1m post 2nd primary dose child")

# ---- Fetch the raw export (cache locally; refetch only if missing) --
if (!file.exists(EXPORT_CACHE)) {
  message("Downloading WISSPAR export -> ", EXPORT_CACHE)
  dir.create(dirname(EXPORT_CACHE), showWarnings = FALSE, recursive = TRUE)
  utils::download.file(EXPORT_URL, EXPORT_CACHE, mode = "wb", quiet = TRUE)
}
d <- fromJSON(EXPORT_CACHE, flatten = TRUE)

# ---- Filter to the first-pass slice ---------------------------------
d <- subset(d,
            assay %in% c("GMC", "OPA") &
            vaccine %in% c(REF_VACCINE, COMP_VACCINE) &
            dose_description %in% POSTPRIMARY &
            !is.na(value) & value > 0)

# Keep only trials that carry BOTH products at the SAME timepoint & assay
# (true within-study head-to-head, decided per assay x dose_description).
key <- with(d, paste(study_id, assay, dose_description, sep = "|"))
has_both <- tapply(d$vaccine, key,
                   function(v) all(c(REF_VACCINE, COMP_VACCINE) %in% v))
d <- d[has_both[key], ]
if (nrow(d) == 0) stop("No head-to-head PCV7 vs PCV13 post-primary rows found.")

# ---- Collapse any duplicate reads (matches WISSPAR arm aggregation) --
# One estimate per study x assay x timepoint x schedule x vaccine x serotype;
# average duplicate value reads and drop the CI once aggregated (>1 read).
grp <- with(d, paste(study_id, assay, dose_description, schedule,
                     vaccine, serotype, sep = "|"))
agg <- do.call(rbind, lapply(split(d, grp), function(g) {
  n <- nrow(g)
  data.frame(
    study_id         = g$study_id[1],
    assay            = g$assay[1],
    dose_description = g$dose_description[1],
    schedule         = g$schedule[1],
    vaccine          = g$vaccine[1],
    serotype         = g$serotype[1],
    value            = mean(g$value),
    lower_limit      = if (n == 1) g$lower_limit[1] else NA_real_,
    upper_limit      = if (n == 1) g$upper_limit[1] else NA_real_,
    phase            = g$phase[1],
    sponsor          = g$sponsor[1],
    location         = g$location_continent[1],
    n_reads          = n,
    stringsAsFactors = FALSE
  )
}))
rownames(agg) <- NULL

# Order: assay, serotype (numeric-aware), study, vaccine
so <- suppressWarnings(as.numeric(gsub("[^0-9].*$", "", agg$serotype)))
agg <- agg[order(agg$assay, so, agg$serotype, agg$study_id, agg$vaccine), ]

# ---- Write ----------------------------------------------------------
write.csv(agg, OUT_FILE, row.names = FALSE)

cat("Wrote", OUT_FILE, "with", nrow(agg), "rows.\n\n")
cat("Head-to-head trials (PCV7 & PCV13 Pfizer, children, post-primary):\n")
tab <- as.data.frame(table(agg$study_id, agg$assay, agg$vaccine))
names(tab) <- c("study_id", "assay", "vaccine", "n_serotypes")
print(subset(tab, n_serotypes > 0), row.names = FALSE)
cat("\nSerotypes per assay:\n")
for (a in unique(agg$assay)) {
  cat(" ", a, ":", paste(sort(unique(agg$serotype[agg$assay == a])), collapse = ", "), "\n")
}

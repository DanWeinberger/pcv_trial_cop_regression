# =====================================================================
# Analysis registry for the COP regression project.
#
# This is the ONLY file you edit to add a new comparison. It declares:
#   ANALYSES    - each a single COP model fit: one outcome (case counts)
#                 regressed on one immunogenicity predictor source.
#   COMPARISONS - named sets of analyses whose global slopes (mu_b1) we
#                 contrast against each other.
#
# The run driver (R/run_analysis.R) and the slope comparison
# (R/compare_slopes.R) read this file; neither hard-codes a study name.
#
# ---- How to add a new analysis --------------------------------------
# Append a list() to ANALYSES with:
#   id              short slug -> results/<id>/ output folder
#   data_file       merged CSV containing serotype x arm rows
#   predictor_study value of the `Study` column supplying the predictor GMCs
#   predictor_label human-readable predictor name (used in plots/tables)
#   outcome_label   human-readable outcome name (metadata / titles)
#
# In the current merged file the OUTCOME (Cases / Total_Cases) is the Whitney
# IPD data and is identical across studies; changing `predictor_study` swaps
# only the immunogenicity source. When a future dataset carries a genuinely
# different outcome, add it as a new merged file + a new analysis entry.
# =====================================================================

DEFAULT_DATA_FILE     <- file.path("data", "siber_whitney_merged.csv")
ANDREWS_DATA_FILE     <- file.path("data", "siber_andrews_merged.csv")
WISSPAR_DATA_FILE     <- file.path("data", "wisspar_andrews_merged.csv")
VDL_PCV7_DATA_FILE    <- file.path("data", "siber_vanderlinden_pcv7_merged.csv")
VDL_PCV13_DATA_FILE   <- file.path("data", "wisspar_vanderlinden_pcv13_merged.csv")
OUTCOME_WHITNEY   <- "Whitney IPD case counts (NCKP surveillance)"
OUTCOME_ANDREWS   <- "Andrews 2019 PCV7 IPD case counts (England & Wales, >=1 dose)"
OUTCOME_WISSPAR   <- "Andrews 2019 PCV13 additional-serotype IPD (England & Wales, >=1 dose)"
OUTCOME_VDL_PCV7  <- "van der Linden 2016 PCV7 IPD (Germany, at least one dose)"
OUTCOME_VDL_PCV13 <- "van der Linden 2016 PCV13 additional-serotype IPD (Germany, at least one dose)"

ANALYSES <- list(
  list(
    id              = "nckp",
    data_file       = DEFAULT_DATA_FILE,
    predictor_study = "NCKP",
    predictor_label = "NCKP immunogenicity (US, Whitney/Kaiser)",
    outcome_label   = OUTCOME_WHITNEY
  ),
  list(
    id              = "navajo",
    data_file       = DEFAULT_DATA_FILE,
    predictor_study = "Am_Indian",
    predictor_label = "Navajo / American Indian immunogenicity",
    outcome_label   = OUTCOME_WHITNEY
  ),
  list(
    id              = "south_africa",
    data_file       = DEFAULT_DATA_FILE,
    predictor_study = "South_Africa",
    predictor_label = "South Africa immunogenicity",
    outcome_label   = OUTCOME_WHITNEY
  ),

  # ---- Same immunogenicity predictors, Andrews 2019 PCV7 outcome ----------
  list(
    id              = "nckp_andrews",
    data_file       = ANDREWS_DATA_FILE,
    predictor_study = "NCKP",
    predictor_label = "NCKP immunogenicity (US, Whitney/Kaiser)",
    outcome_label   = OUTCOME_ANDREWS
  ),
  list(
    id              = "navajo_andrews",
    data_file       = ANDREWS_DATA_FILE,
    predictor_study = "Am_Indian",
    predictor_label = "Navajo / American Indian immunogenicity",
    outcome_label   = OUTCOME_ANDREWS
  ),
  list(
    id              = "south_africa_andrews",
    data_file       = ANDREWS_DATA_FILE,
    predictor_study = "South_Africa",
    predictor_label = "South Africa immunogenicity",
    outcome_label   = OUTCOME_ANDREWS
  ),

  # ---- WISSPAR head-to-head PCV13-vs-PCV7 GMC, one analysis per trial --------
  # Children, post-primary, PCV13-additional serotypes (1, 3, 6A, 7F, 19A).
  # PCV13 -> Immunized, PCV7 -> Unimmunized (no-antigen comparator). Outcome is
  # the Andrews 2019 PCV13 IPD data, broadcast identically across trials, so
  # the slopes are directly comparable (only the predictor GMC differs by trial).
  # (NCT00205803 (USA) is excluded: the WISSPAR export carries no 95% CI for it,
  # and the error-in-variables model requires that measurement-error input.)
  list(
    id              = "wisspar_de",
    data_file       = WISSPAR_DATA_FILE,
    predictor_study = "NCT00366340",
    predictor_label = "WISSPAR NCT00366340 PCV13/PCV7 GMC (Germany)",
    outcome_label   = OUTCOME_WISSPAR
  ),
  list(
    id              = "wisspar_tw",
    data_file       = WISSPAR_DATA_FILE,
    predictor_study = "NCT00688870",
    predictor_label = "WISSPAR NCT00688870 PCV13/PCV7 GMC (Taiwan)",
    outcome_label   = OUTCOME_WISSPAR
  ),
  list(
    id              = "wisspar_kr",
    data_file       = WISSPAR_DATA_FILE,
    predictor_study = "NCT00689351",
    predictor_label = "WISSPAR NCT00689351 PCV13/PCV7 GMC (South Korea)",
    outcome_label   = OUTCOME_WISSPAR
  ),

  # ---- SIBER predictors, van der Linden 2016 PCV7 outcome (Germany) ----------
  # Same three immunogenicity sources as nckp/navajo/south_africa, but the PCV7
  # IPD outcome is van der Linden Table 1 ("at least one dose"), over the seven
  # PCV7 serotypes. Parallels the *_andrews PCV7 analyses.
  list(
    id              = "nckp_vdl_pcv7",
    data_file       = VDL_PCV7_DATA_FILE,
    predictor_study = "NCKP",
    predictor_label = "NCKP immunogenicity (US, Whitney/Kaiser)",
    outcome_label   = OUTCOME_VDL_PCV7
  ),
  list(
    id              = "navajo_vdl_pcv7",
    data_file       = VDL_PCV7_DATA_FILE,
    predictor_study = "Am_Indian",
    predictor_label = "Navajo / American Indian immunogenicity",
    outcome_label   = OUTCOME_VDL_PCV7
  ),
  list(
    id              = "south_africa_vdl_pcv7",
    data_file       = VDL_PCV7_DATA_FILE,
    predictor_study = "South_Africa",
    predictor_label = "South Africa immunogenicity",
    outcome_label   = OUTCOME_VDL_PCV7
  ),

  # ---- WISSPAR head-to-head GMC, van der Linden 2016 PCV13 outcome -----------
  # PCV13-additional serotypes (1, 3, 6A, 7F, 19A; 5 has zero cases in Germany).
  # PCV13 -> Immunized, PCV7 -> Unimmunized. Outcome is van der Linden Table 2
  # ("at least one dose"), broadcast identically across trials. Parallels the
  # wisspar_de/tw/kr PCV13 analyses (Andrews outcome).
  list(
    id              = "wisspar_de_vdl_pcv13",
    data_file       = VDL_PCV13_DATA_FILE,
    predictor_study = "NCT00366340",
    predictor_label = "WISSPAR NCT00366340 PCV13/PCV7 GMC (Germany)",
    outcome_label   = OUTCOME_VDL_PCV13
  ),
  list(
    id              = "wisspar_tw_vdl_pcv13",
    data_file       = VDL_PCV13_DATA_FILE,
    predictor_study = "NCT00688870",
    predictor_label = "WISSPAR NCT00688870 PCV13/PCV7 GMC (Taiwan)",
    outcome_label   = OUTCOME_VDL_PCV13
  ),
  list(
    id              = "wisspar_kr_vdl_pcv13",
    data_file       = VDL_PCV13_DATA_FILE,
    predictor_study = "NCT00689351",
    predictor_label = "WISSPAR NCT00689351 PCV13/PCV7 GMC (South Korea)",
    outcome_label   = OUTCOME_VDL_PCV13
  )
)

# Named comparisons. `analyses` lists analysis ids (in display order);
# `reference` (optional) is the id used as the baseline for pairwise slope
# differences (defaults to the first listed).
COMPARISONS <- list(
  list(
    id        = "predictor_source",
    label     = "COP slope by immunogenicity source (outcome = Whitney IPD)",
    analyses  = c("nckp", "navajo", "south_africa"),
    reference = "nckp"
  ),
  list(
    id        = "predictor_source_andrews",
    label     = "PCV7 COP slope by immunogenicity source and outcome (Andrews 2019 vs van der Linden 2016)",
    analyses  = c("nckp_andrews", "navajo_andrews", "south_africa_andrews",
                  "nckp_vdl_pcv7", "navajo_vdl_pcv7", "south_africa_vdl_pcv7"),
    reference = "nckp_andrews"
  ),
  list(
    id        = "outcome_nckp",
    label     = "COP slope by outcome for NCKP immunogenicity (Whitney vs Andrews vs van der Linden)",
    analyses  = c("nckp", "nckp_andrews", "nckp_vdl_pcv7"),
    reference = "nckp"
  ),
  list(
    id        = "outcome_navajo",
    label     = "COP slope by outcome for Navajo immunogenicity (Whitney vs Andrews vs van der Linden)",
    analyses  = c("navajo", "navajo_andrews", "navajo_vdl_pcv7"),
    reference = "navajo"
  ),
  list(
    id        = "outcome_south_africa",
    label     = "COP slope by outcome for South Africa immunogenicity (Whitney vs Andrews vs van der Linden)",
    analyses  = c("south_africa", "south_africa_andrews", "south_africa_vdl_pcv7"),
    reference = "south_africa"
  ),
  list(
    id        = "wisspar_by_trial",
    label     = "PCV13 COP slope by WISSPAR trial and outcome (Andrews 2019 vs van der Linden 2016, additional serotypes)",
    analyses  = c("wisspar_de", "wisspar_tw", "wisspar_kr",
                  "wisspar_de_vdl_pcv13", "wisspar_tw_vdl_pcv13", "wisspar_kr_vdl_pcv13"),
    reference = "wisspar_de"
  ),

  # ---- van der Linden 2016 outcome comparisons ------------------------------
  list(
    id        = "predictor_source_vdl_pcv7",
    label     = "COP slope by immunogenicity source (outcome = van der Linden 2016 PCV7)",
    analyses  = c("nckp_vdl_pcv7", "navajo_vdl_pcv7", "south_africa_vdl_pcv7"),
    reference = "nckp_vdl_pcv7"
  ),
  list(
    id        = "vdl_pcv13_by_trial",
    label     = "COP slope by WISSPAR head-to-head trial (outcome = van der Linden 2016 PCV13)",
    analyses  = c("wisspar_de_vdl_pcv13", "wisspar_tw_vdl_pcv13", "wisspar_kr_vdl_pcv13"),
    reference = "wisspar_de_vdl_pcv13"
  ),

  # ---- PCV7 outcome: Andrews (England & Wales) vs van der Linden (Germany) ----
  list(
    id        = "pcv7_outcome_nckp",
    label     = "PCV7 COP slope by outcome for NCKP immunogenicity (Andrews vs van der Linden)",
    analyses  = c("nckp_andrews", "nckp_vdl_pcv7"),
    reference = "nckp_andrews"
  ),
  list(
    id        = "pcv7_outcome_navajo",
    label     = "PCV7 COP slope by outcome for Navajo immunogenicity (Andrews vs van der Linden)",
    analyses  = c("navajo_andrews", "navajo_vdl_pcv7"),
    reference = "navajo_andrews"
  ),
  list(
    id        = "pcv7_outcome_south_africa",
    label     = "PCV7 COP slope by outcome for South Africa immunogenicity (Andrews vs van der Linden)",
    analyses  = c("south_africa_andrews", "south_africa_vdl_pcv7"),
    reference = "south_africa_andrews"
  ),

  # ---- PCV13 outcome: Andrews (England & Wales) vs van der Linden (Germany) ---
  list(
    id        = "pcv13_outcome_de",
    label     = "PCV13 COP slope by outcome for NCT00366340 GMC (Andrews vs van der Linden)",
    analyses  = c("wisspar_de", "wisspar_de_vdl_pcv13"),
    reference = "wisspar_de"
  ),
  list(
    id        = "pcv13_outcome_tw",
    label     = "PCV13 COP slope by outcome for NCT00688870 GMC (Andrews vs van der Linden)",
    analyses  = c("wisspar_tw", "wisspar_tw_vdl_pcv13"),
    reference = "wisspar_tw"
  ),
  list(
    id        = "pcv13_outcome_kr",
    label     = "PCV13 COP slope by outcome for NCT00689351 GMC (Andrews vs van der Linden)",
    analyses  = c("wisspar_kr", "wisspar_kr_vdl_pcv13"),
    reference = "wisspar_kr"
  )
)

# ---- Registry helpers (no need to edit below) -----------------------
names(ANALYSES)    <- vapply(ANALYSES,    function(a) a$id, character(1))
names(COMPARISONS) <- vapply(COMPARISONS, function(c) c$id, character(1))

get_analysis <- function(id) {
  if (!id %in% names(ANALYSES)) {
    stop("Unknown analysis id '", id, "'. Known: ", paste(names(ANALYSES), collapse = ", "))
  }
  ANALYSES[[id]]
}

get_comparison <- function(id) {
  if (!id %in% names(COMPARISONS)) {
    stop("Unknown comparison id '", id, "'. Known: ", paste(names(COMPARISONS), collapse = ", "))
  }
  COMPARISONS[[id]]
}

# "centered" and "ratio" are algebraically-equivalent reparameterizations (see
# R/cop_model.R PARAMETERIZATIONS), so they intentionally share results/<id>/
# unsuffixed. Any OTHER parameterization (e.g. "RE", a genuinely different
# model) gets its own results/<id>_<parameterization>/ so it never clobbers
# the centered/ratio fit for the same analysis.
analysis_out_dir <- function(id, parameterization = "centered") {
  base <- file.path("results", id)
  if (parameterization %in% c("centered", "ratio")) base else paste0(base, "_", parameterization)
}
comparison_out_dir <- function(id) file.path("results", "comparisons", id)

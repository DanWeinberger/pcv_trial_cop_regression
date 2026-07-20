# =====================================================================
# Model-spec registry for the COP regression project.
#
# This is the ONLY file you edit to add a new outcome/immunogenicity source,
# a new model spec, or a comparison. It declares three layers, from reusable
# building blocks up to what actually gets fit:
#
#   OUTCOME_SOURCES   - named outcome datasets (IPD case counts). Each becomes
#                       one outcome STUDY (its own a[k,s] baseline per
#                       serotype it reports) when included in a MODEL_SPEC.
#   IMMUNO_SOURCES    - named immunogenicity datasets (GMC + 95% CI). Each
#                       one included in a MODEL_SPEC is POOLED onto the SAME
#                       shared latent log-GMC per serotype/arm.
#   MODEL_SPECS       - a model to fit: a set of outcome_sources (by id) +
#                       a set of immuno_sources (by id). The production
#                       default is a single spec ("pooled") that pools
#                       EVERY registered outcome study and EVERY registered
#                       immunogenicity source into one fit.
#   MODEL_COMPARISONS - named sets of model spec ids to rank against each
#                       other by WAIC (see R/compare_waic.R) -- empty until a
#                       second spec exists to compare "pooled" against.
#
# The run driver (R/run_models.R) and the WAIC comparison (R/compare_waic.R)
# read this file; neither hard-codes a study name.
#
# ---- How to add an alternative model spec to compare against "pooled" --
# 1. If your outcome/immunogenicity dataset isn't already registered, add it
#    to OUTCOME_SOURCES / IMMUNO_SOURCES (a data_file + label; immuno sources
#    also need the `study` value selecting the right Study-column subset).
# 2. Add a list() to MODEL_SPECS: an id (-> results/<id>/ output folder), the
#    outcome ids to pool (normally the SAME set "pooled" uses, so WAIC is
#    comparing like with like), and the immuno ids to pool.
# 3. Add a MODEL_COMPARISONS entry listing the spec ids to rank (include
#    "pooled" as the baseline).
# 4. `Rscript R/run_models.R <new spec id>` then
#    `Rscript R/compare_waic.R <comparison id>`.
# =====================================================================

DEFAULT_DATA_FILE   <- file.path("data", "siber_whitney_merged.csv")
ANDREWS_DATA_FILE   <- file.path("data", "siber_andrews_merged.csv")
WISSPAR_DATA_FILE   <- file.path("data", "wisspar_andrews_merged.csv")
VDL_PCV7_DATA_FILE  <- file.path("data", "siber_vanderlinden_pcv7_merged.csv")
VDL_PCV13_DATA_FILE <- file.path("data", "wisspar_vanderlinden_pcv13_merged.csv")

# ---- Outcome sources (IPD case counts) -------------------------------
# PCV7 serotypes (4, 6B, 9V, 14, 18C, 19F, 23F):
#   whitney       - Whitney IPD (NCKP surveillance, US)
#   andrews       - Andrews 2019 PCV7 IPD (England & Wales, >=1 dose)
#   vdl_pcv7      - van der Linden 2016 PCV7 IPD (Germany, at least one dose)
# PCV13-additional serotypes (1, 3, 6A, 7F, 19A):
#   andrews_pcv13 - Andrews 2019 PCV13-additional-serotype IPD (England &
#                   Wales, >=1 dose), broadcast across the WISSPAR trials
#   vdl_pcv13     - van der Linden 2016 PCV13-additional-serotype IPD (Germany)
OUTCOME_SOURCES <- list(
  whitney = list(
    data_file = DEFAULT_DATA_FILE,
    label     = "Whitney IPD case counts (NCKP surveillance)"
  ),
  andrews = list(
    data_file = ANDREWS_DATA_FILE,
    label     = "Andrews 2019 PCV7 IPD case counts (England & Wales, >=1 dose)"
  ),
  vdl_pcv7 = list(
    data_file = VDL_PCV7_DATA_FILE,
    label     = "van der Linden 2016 PCV7 IPD (Germany, at least one dose)"
  ),
  andrews_pcv13 = list(
    data_file = WISSPAR_DATA_FILE,
    label     = "Andrews 2019 PCV13 additional-serotype IPD (England & Wales, >=1 dose)"
  ),
  vdl_pcv13 = list(
    data_file = VDL_PCV13_DATA_FILE,
    label     = "van der Linden 2016 PCV13 additional-serotype IPD (Germany, at least one dose)"
  )
)

# ---- Immunogenicity sources (GMC + 95% CI) ---------------------------
# PCV7 serotypes -- three SIBER predictor sources, selected by the `Study`
# column of DEFAULT_DATA_FILE:
#   nckp, navajo, south_africa
# PCV13-additional serotypes -- WISSPAR head-to-head PCV13-vs-PCV7 GMC, one
# per trial (NCT00205803 excluded: no reported 95% CI, which the EIV model
# requires):
#   wisspar_de (Germany), wisspar_tw (Taiwan), wisspar_kr (South Korea)
IMMUNO_SOURCES <- list(
  nckp = list(
    data_file = DEFAULT_DATA_FILE, study = "NCKP",
    label     = "NCKP immunogenicity (US, Whitney/Kaiser)"
  ),
  navajo = list(
    data_file = DEFAULT_DATA_FILE, study = "Am_Indian",
    label     = "Navajo / American Indian immunogenicity"
  ),
  south_africa = list(
    data_file = DEFAULT_DATA_FILE, study = "South_Africa",
    label     = "South Africa immunogenicity"
  ),
  wisspar_de = list(
    data_file = WISSPAR_DATA_FILE, study = "NCT00366340",
    label     = "WISSPAR NCT00366340 PCV13/PCV7 GMC (Germany)"
  ),
  wisspar_tw = list(
    data_file = WISSPAR_DATA_FILE, study = "NCT00688870",
    label     = "WISSPAR NCT00688870 PCV13/PCV7 GMC (Taiwan)"
  ),
  wisspar_kr = list(
    data_file = WISSPAR_DATA_FILE, study = "NCT00689351",
    label     = "WISSPAR NCT00689351 PCV13/PCV7 GMC (South Korea)"
  )
)

# ---- Model specs: each pools a set of outcome studies + a set of --------
# ---- immunogenicity ("predictor") sources into one multistudy fit. ------
# Production default: ONE model, pooling every registered outcome study
# (PCV7 and PCV13-additional together -- there's no "reference study" or
# connectedness requirement in the multistudy design, see
# JAGS/cop_eiv_model_multistudy.jags, so the two serotype blocks link through
# the shared mu_a/tau_a and mu_b1/tau_b1 hyperpriors even though no single
# outcome study reports both) and every registered immunogenicity source
# onto the shared latent log-GMC.
#
# A spec is still the unit R/compare_waic.R ranks: if you later want to test
# an alternative predictor set (e.g. drop a source, add a new one) against
# this pooled baseline, add a second list() here with the SAME `outcomes` and
# a different `immuno`, then add a MODEL_COMPARISONS entry listing both ids.
MODEL_SPECS <- list(
  list(
    id       = "pooled",
    outcomes = c("whitney", "andrews", "vdl_pcv7", "andrews_pcv13", "vdl_pcv13"),
    immuno   = c("nckp", "navajo", "south_africa", "wisspar_de", "wisspar_tw", "wisspar_kr")
  )
)

# Named comparisons. `specs` lists model spec ids (in display order);
# `reference` (optional) is the id used as the baseline for pairwise slope
# differences (defaults to the first listed). Empty for now -- there is only
# one spec in MODEL_SPECS, so nothing to rank yet. R/run_models.R still
# writes waic.csv for "pooled" so a future alternative spec has something to
# be compared against; add a MODEL_COMPARISONS entry once a second spec
# exists.
MODEL_COMPARISONS <- list()

# ---- Registry helpers (no need to edit below) -----------------------
# OUTCOME_SOURCES / IMMUNO_SOURCES already carry names from their literal
# definitions above; only the list-of-lists registries need id -> name.
names(MODEL_SPECS)       <- vapply(MODEL_SPECS,       function(s) s$id, character(1))
names(MODEL_COMPARISONS) <- vapply(MODEL_COMPARISONS, function(c) c$id, character(1))

get_outcome_source <- function(id) {
  if (!id %in% names(OUTCOME_SOURCES)) {
    stop("Unknown outcome source id '", id, "'. Known: ", paste(names(OUTCOME_SOURCES), collapse = ", "))
  }
  OUTCOME_SOURCES[[id]]
}

get_immuno_source <- function(id) {
  if (!id %in% names(IMMUNO_SOURCES)) {
    stop("Unknown immunogenicity source id '", id, "'. Known: ", paste(names(IMMUNO_SOURCES), collapse = ", "))
  }
  IMMUNO_SOURCES[[id]]
}

get_model_spec <- function(id) {
  if (!id %in% names(MODEL_SPECS)) {
    stop("Unknown model spec id '", id, "'. Known: ", paste(names(MODEL_SPECS), collapse = ", "))
  }
  MODEL_SPECS[[id]]
}

get_model_comparison <- function(id) {
  if (!id %in% names(MODEL_COMPARISONS)) {
    stop("Unknown model comparison id '", id, "'. Known: ", paste(names(MODEL_COMPARISONS), collapse = ", "))
  }
  MODEL_COMPARISONS[[id]]
}

# Resolve a model spec's outcome/immuno ids into the outcome_sources /
# immuno_sources lists that prepare_cop_data_multistudy() expects.
spec_outcome_sources <- function(spec) {
  lapply(spec$outcomes, function(id) {
    src <- get_outcome_source(id)
    list(data_file = src$data_file, study_id = id, label = src$label)
  })
}
spec_immuno_sources <- function(spec) {
  lapply(spec$immuno, function(id) {
    src <- get_immuno_source(id)
    list(data_file = src$data_file, study = src$study, label = src$label)
  })
}

spec_out_dir <- function(id, predictor_error = "se") {
  base <- file.path("results", id)
  if (predictor_error == "sd") paste0(base, "_sd") else base
}
comparison_out_dir <- function(id) file.path("results", "comparisons", id)

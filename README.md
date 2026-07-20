# PCV correlate-of-protection regression

Bayesian error-in-variables Poisson regression estimating a pneumococcal
correlate of protection (COP) across multiple outcome studies and multiple
immunogenicity (GMC) sources. Serotype-specific IPD case counts, from one or
more outcome studies, are regressed on a single shared latent absolute
log-GMC per serotype/arm, pooled across one or more immunogenicity sources
and measured with error. The global slope `mu_b1` is the pooled COP, shared
across every outcome study; baselines vary freely by study x serotype.

## Layout

| Path | Role |
|------|------|
| `JAGS/cop_eiv_model_multistudy.jags` | The Bayesian model (the only one this project carries forward). |
| `R/cop_model.R` | Reusable engine: `prepare_cop_data_multistudy()`, `fit_cop_multistudy()`, `compute_waic()`, `slope_summary()`, `study_summary()`, `plot_cop_absolute_multistudy()`. No model spec is hard-coded here. |
| `R/config.R` | **The only file you edit to add a model spec or comparison.** Registry of `OUTCOME_SOURCES`, `IMMUNO_SOURCES`, `MODEL_SPECS` (which outcomes + which predictor sources to fit), and `MODEL_COMPARISONS` (sets of specs to rank by WAIC). |
| `R/run_models.R` | Driver: fit + summarise one or more model specs. |
| `R/compare_waic.R` | Rank model specs by WAIC and compare their global slope `mu_b1`. |
| `results/<id>/` | Per-spec outputs (git-ignored; regenerated on demand). |
| `results/comparisons/<id>/` | Per-comparison outputs (git-ignored). |

## Running

```sh
Rscript R/run_models.R           # fit every model spec in the registry (currently just "pooled")
Rscript R/run_models.R pooled    # fit a spec by id
Rscript R/compare_waic.R         # run every WAIC comparison (none defined yet -- see below)
```

Each `results/<id>/` gets: `posterior_summary.csv`, `mcmc.rds`,
`diagnostics.pdf`, `waic.csv`, `slope_summary.csv`, `study_summary.csv`,
`cop_scatter_absolute_risk_gmc_multistudy.png`.

Each `results/comparisons/<id>/` gets a **WAIC leaderboard**
(`waic_comparison.csv`, `waic_comparison.pdf`/`.png`) ranking the member
specs, and the **global slope** `mu_b1` side by side (`slopes_by_spec.csv`,
`slope_pairwise_diffs.csv`, `slope_forest.pdf`/`.png`,
`slope_overlay.pdf`/`.png`) — WAIC tells you which predictor combination the
outcome data prefers; the slope forest tells you whether that preference
actually changes the COP estimate.

## Model specs and comparisons

A **model spec** (`MODEL_SPECS` in `R/config.R`) pools a fixed set of outcome
studies with a chosen set of immunogenicity ("predictor") sources into one
`cop_eiv_model_multistudy.jags` fit.

The production default is a single spec, `pooled`, which pools **every**
registered outcome study and **every** registered immunogenicity source into
one fit: PCV7 serotypes (4, 6B, 9V, 14, 18C, 19F, 23F) and PCV13-additional
serotypes (1, 3, 6A, 7F, 19A) together, spanning 12 serotypes total.

- Outcomes: Whitney IPD (`data/siber_whitney_merged.csv`), Andrews 2019 PCV7
  (`data/siber_andrews_merged.csv`, built by `R/build_andrews_merged.R`), van
  der Linden 2016 PCV7 (`data/siber_vanderlinden_pcv7_merged.csv`, built by
  `R/build_vanderlinden_pcv7_merged.R`), Andrews 2019 PCV13-additional
  (`data/wisspar_andrews_merged.csv`, built by
  `R/build_wisspar_andrews_merged.R`), van der Linden 2016 PCV13-additional
  (`data/wisspar_vanderlinden_pcv13_merged.csv`, built by
  `R/build_vanderlinden_pcv13_merged.R`).
- Immunogenicity: NCKP, Navajo/American Indian, South Africa (PCV7
  serotypes), plus the three WISSPAR head-to-head PCV13-vs-PCV7 trials —
  Germany, Taiwan, South Korea (PCV13-additional serotypes).

There's no "reference study" or connectedness requirement in the multistudy
design (see `JAGS/cop_eiv_model_multistudy.jags`), so the PCV7 and
PCV13-additional blocks link through the shared `mu_a`/`tau_a` and
`mu_b1`/`tau_b1` hyperpriors even though no single outcome study reports
both sets of serotypes.

`MODEL_COMPARISONS` is empty for now — there's only one spec, so nothing to
rank. If you register an alternative predictor set later (e.g. drop a
source, add a new one) as a second spec sharing `pooled`'s outcomes, add a
`MODEL_COMPARISONS` entry and `R/compare_waic.R` will rank it against
`pooled` by WAIC, with the global slope `mu_b1` shown side by side — WAIC
tells you which predictor combination the outcome data prefers; the slope
forest tells you whether that preference actually changes the COP estimate.

## WAIC methodology

WAIC only compares models fit to the *same* data. Swapping which
immunogenicity sources are pooled changes the number of immunogenicity rows
fed into the model but never the outcome rows — so `compute_waic()` (in
`R/cop_model.R`) computes WAIC from the **outcome (case-count) likelihood
only**, monitored per outcome row as `log_lik[m]` in
`cop_eiv_model_multistudy.jags`, never from the immunogenicity
measurement-error likelihood. Each row (a study x serotype's paired
unimmunized/immunized case counts) is one WAIC observation unit, since both
arms share the same `b1[s]`/`x_u[s]`/`x_i[s]` and aren't meaningfully
separable. WAIC itself is hand-rolled from those pointwise log-likelihoods
(standard formulas: `lppd`, `p_waic` from the pointwise posterior variance,
`elpd_waic = lppd - p_waic`, `waic = -2 * elpd_waic`) rather than pulled in
via the `loo` package, since nothing else in this project needs it. Requires
JAGS >= 4.3 (for `logdensity.pois()`).

## Adding a model spec or comparison

1. If the outcome/immunogenicity dataset isn't already registered, add it to
   `OUTCOME_SOURCES` / `IMMUNO_SOURCES` in `R/config.R` (a `data_file` +
   `label`; immuno sources also need the `study` value selecting the right
   `Study`-column subset). If it's a new outcome, also add its id to
   `pooled`'s `outcomes` so it's covered by the production model.
2. To test an alternative predictor set against `pooled`: add a
   `list(id=, outcomes=, immuno=)` to `MODEL_SPECS` (normally with the SAME
   `outcomes` as `pooled`, so WAIC is comparing like with like), then a
   `MODEL_COMPARISONS` entry listing both spec ids (include `pooled` as the
   baseline).
3. `Rscript R/run_models.R <new spec id>` then
   `Rscript R/compare_waic.R <comparison id>`.

## WISSPAR immunogenicity import

Two build scripts pull immunogenicity from the [WISSPAR](https://github.com/PneumococcalCapsules/WISSPAR_v2)
database (Worldwide Index of Serotype Specific Pneumococcal Antibody Responses):

1. `R/build_wisspar_head2head.R` — downloads the WISSPAR export
   (`data/wisspar_export.json`, ~4 MB, git-ignored, re-fetched on demand) and
   writes the first-pass slice `data/wisspar_pcv7_pcv13_child_postprimary.csv`:
   head-to-head **PCV7 vs PCV13 (Pfizer)** trials in **children** at the
   **post-primary** timepoint, both **GMC and OPA** assays. "Head-to-head" =
   a single trial that measured both products in the same population/timepoint.
2. `R/build_wisspar_andrews_merged.R` — turns the GMC head-to-head slice into
   the COP predictor sources described above (`data/wisspar_andrews_merged.csv`).

```sh
Rscript R/build_wisspar_head2head.R        # 1. import the slice (GMC + OPA)
Rscript R/build_wisspar_andrews_merged.R   # 2. build the COP merged CSV
Rscript R/run_models.R pooled              # 3. refit the pooled model with it
```

Scope notes for the first pass:

- **OPA is imported but not modelled.** The head-to-head OPA trials cover no
  PCV13-additional serotype (only shared PCV7 serotypes), so the COP predictor
  source is GMC-only. OPA stays in the descriptive slice CSV.
- **NCT00205803 (USA) and NCT00452790 are excluded** from the COP fits — the
  former reports no 95% CI (the error-in-variables model needs it); the latter
  measured only shared PCV7 serotypes, none of the PCV13-additional set.
- Serotype **5** has a head-to-head GMC but no Andrews PCV13 outcome; **6C** has
  the outcome but no GMC. Both drop out, leaving the COP set {1, 3, 6A, 7F, 19A}.

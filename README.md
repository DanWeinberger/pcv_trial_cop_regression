# PCV correlate-of-protection regression

Bayesian error-in-variables Poisson regression estimating a pneumococcal
correlate of protection (COP): serotype-specific IPD case counts regressed on
centered log-GMC immunogenicity (measured with error), with hierarchical
serotype intercepts and slopes. The global slope `mu_b1` is the pooled COP.

## Layout

| Path | Role |
|------|------|
| `JAGS/cop_eiv_model.jags` | The Bayesian model (unchanged across analyses). |
| `R/cop_model.R` | Reusable engine: `prepare_cop_data()`, `fit_cop()`, `plot_cop_scatter()`, `slope_summary()`. No analysis is hard-coded here. |
| `R/config.R` | **The only file you edit to add a comparison.** Registry of `ANALYSES` (outcome + predictor) and `COMPARISONS` (sets of analyses whose slopes are contrasted). |
| `R/run_analysis.R` | Driver: fit + plot + summarise one or more analyses. |
| `R/compare_slopes.R` | Contrast the global slope `mu_b1` across analyses in a comparison. |
| `results/<id>/` | Per-analysis outputs. |
| `results/comparisons/<id>/` | Per-comparison outputs. |

## Running

```sh
Rscript R/run_analysis.R              # fit every analysis in the registry
Rscript R/run_analysis.R nckp navajo  # fit selected analyses by id
Rscript R/compare_slopes.R            # run every comparison
Rscript R/compare_slopes.R predictor_source   # one comparison by id
```

Each `results/<id>/` gets: `posterior_summary.csv`, `mcmc.rds`,
`diagnostics.pdf`, `cop_scatter_gmr_rr.png`, `slope_summary.csv`.
Each `results/comparisons/<id>/` gets, at the **global** slope level:
`slopes_by_analysis.csv`, `slope_pairwise_diffs.csv`, `slope_forest.pdf`/`.png`,
`slope_overlay.pdf`/`.png`; and at the **per-serotype** slope level:
`serotype_slopes_by_analysis.csv`, `serotype_slope_pairwise_diffs.csv`,
`serotype_slope_forest.pdf`/`.png` (forest faceted by serotype).

## Current analyses

Two outcomes, each regressed on the same three immunogenicity predictor
sources (selected by the `Study` column: `NCKP`, `Am_Indian`, `South_Africa`).

**Whitney IPD outcome** — `data/siber_whitney_merged.csv` (case counts identical
across studies):

- `nckp` — NCKP immunogenicity (US, Whitney/Kaiser)
- `navajo` — Navajo / American Indian immunogenicity (`Study == "Am_Indian"`)
- `south_africa` — South Africa immunogenicity

**Andrews 2019 PCV7 outcome** — `data/siber_andrews_merged.csv`, built by
`R/build_andrews_merged.R`. Reuses the SIBER immunogenicity predictors but swaps
the outcome to the Andrews 2019 PCV7 trial (Table 2, `>=1 dose` schedule,
England & Wales; Vaccinated→Immunized, Unvaccinated→Unimmunized) over the seven
PCV7 serotypes (4, 6B, 9V, 14, 18C, 19F, 23F):

- `nckp_andrews`, `navajo_andrews`, `south_africa_andrews` — same three
  predictors, Andrews outcome.

**WISSPAR head-to-head PCV13-vs-PCV7 outcome** — `data/wisspar_andrews_merged.csv`,
built by `R/build_wisspar_andrews_merged.R` from the WISSPAR immunogenicity
database (see below). Children, post-primary, restricted to the PCV13-additional
serotypes (1, 3, 6A, 7F, 19A), where the PCV7 arm carries no antigen and is a
genuine no-antigen comparator: `PCV13 → Immunized`, `PCV7 → Unimmunized`. The
outcome is the Andrews 2019 **PCV13** IPD data (Table 3, `>=1 dose`), broadcast
identically across trials so slopes are directly comparable. One analysis per
head-to-head trial (predictor = that trial's per-serotype GMC + 95% CI):

- `wisspar_de` — NCT00366340 (Germany)
- `wisspar_tw` — NCT00688870 (Taiwan)
- `wisspar_kr` — NCT00689351 (South Korea)

**van der Linden 2016 outcome (Germany)** — PCV7 `data/siber_vanderlinden_pcv7_merged.csv`
(built by `R/build_vanderlinden_pcv7_merged.R`) and PCV13
`data/wisspar_vanderlinden_pcv13_merged.csv` (built by
`R/build_vanderlinden_pcv13_merged.R`). The IPD outcome is van der Linden et al.
2016 (PLOS ONE, doi:10.1371/journal.pone.0161257), a German indirect-cohort /
screening-method case-control study, Table 1 (PCV7) and Table 2 (PCV13), both
"at least one dose", children under two. Tables were extracted to
`data/vanderlinden_tables_tidy.csv` in the Andrews schema by
`R/build_vanderlinden_tidy.R`, using the same arm-specific denominator as
Andrews (`Total_Cases = serotype cases + controls`, per vaccination arm). PCV7
reuses the three SIBER predictors over the 7 PCV7 serotypes; PCV13 reuses the
three WISSPAR head-to-head GMC predictors over the additional serotypes
(1, 3, 6A, 7F, 19A — serotype 5 had zero German cases and drops out):

- `nckp_vdl_pcv7`, `navajo_vdl_pcv7`, `south_africa_vdl_pcv7` — PCV7 outcome.
- `wisspar_de_vdl_pcv13`, `wisspar_tw_vdl_pcv13`, `wisspar_kr_vdl_pcv13` — PCV13.

Comparisons:

- `predictor_source` — three global slopes on the Whitney outcome.
- `predictor_source_andrews` — three global slopes on the Andrews outcome.
- `outcome_nckp`, `outcome_navajo`, `outcome_south_africa` — Whitney vs Andrews
  PCV7 outcome, holding each predictor source fixed (global + per-serotype).
- `wisspar_by_trial` — COP slope contrasted across the three WISSPAR
  head-to-head trials (consistency check for the additional-serotype COP).
- `predictor_source_vdl_pcv7` — three global slopes on the van der Linden PCV7
  outcome; `vdl_pcv13_by_trial` — three WISSPAR trials on the van der Linden
  PCV13 outcome.
- `pcv7_outcome_{nckp,navajo,south_africa}` — PCV7 COP slope, Andrews (England &
  Wales) vs van der Linden (Germany) outcome, holding each SIBER predictor fixed.
- `pcv13_outcome_{de,tw,kr}` — PCV13 COP slope, Andrews vs van der Linden
  outcome, holding each WISSPAR head-to-head GMC predictor fixed.

Every comparison now reports both the pooled global slope `mu_b1` and the
serotype-specific slopes `b1[s]`, so "compare the Andrews slopes with the
Whitney slopes" is answered both overall and serotype by serotype.

## Adding a comparison

1. Add a `list(...)` to `ANALYSES` in `R/config.R` (give it an `id`, a
   `data_file`, the `predictor_study` value, and human-readable labels).
2. Add or extend a `COMPARISONS` entry listing the analysis ids to contrast.
3. `Rscript R/run_analysis.R <new ids>` then `Rscript R/compare_slopes.R <comparison id>`.

When a future dataset carries a genuinely different **outcome**, build the
merged CSV for it and point a new analysis's `data_file` at it — no engine code
changes needed.

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
Rscript R/run_analysis.R wisspar_de wisspar_tw wisspar_kr
Rscript R/compare_slopes.R wisspar_by_trial
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

# =====================================================================
# Core reusable engine for the pneumococcal correlate-of-protection (COP)
# error-in-variables Poisson regression -- MULTI-STUDY MODEL ONLY
# (cop_eiv_model_multistudy.jags).
#
# This file defines FUNCTIONS only (it fits nothing when sourced).
#
#   read_outcome_source()          : one outcome dataset -> per-serotype
#                                     Cases/Total_Cases (deduped across the
#                                     immunogenicity Study subsets that share
#                                     the same merged CSV)
#   prepare_gmc_predictor()         : one immunogenicity source -> aligned
#                                     log-GMC error-in-variables predictor
#   prepare_cop_data_multistudy()  : assemble the ragged, multi-study JAGS
#                                     data list from a list of outcome
#                                     sources and a list of immunogenicity
#                                     sources
#   fit_cop_multistudy()           : compile + sample + save outputs
#                                     (posterior_summary.csv, mcmc.rds,
#                                     diagnostics.pdf, waic.csv)
#   compute_waic()                 : WAIC from the model's monitored
#                                     per-row outcome log-likelihood
#   slope_summary()                : tidy per-model slope table (global +
#                                     per serotype)
#   study_summary()                : tidy per-outcome-study baseline table
#   plot_cop_absolute_multistudy() : absolute risk vs. absolute GMC, faceted
#                                     by outcome study
#
# Model structure (see JAGS/cop_eiv_model_multistudy.jags for the full
# statistical rationale):
#   Outcome    : serotype x study-specific IPD case counts (Cases), Poisson,
#                offset = Total_Cases (log-exposure)
#   Predictor  : absolute log-GMC, measured with error (precision from 95%
#                CI), pooled across every immunogenicity source onto a single
#                shared latent per serotype/arm
#   Hierarchy  : study x serotype random intercepts a[k,s], serotype-specific
#                slopes b1[s] shared across every outcome study
# =====================================================================

suppressPackageStartupMessages({
  library(rjags)
  library(coda)
})

# ---------------------------------------------------------------------
# Measurement-error helpers (shared by fit and plot so they never drift)
# ---------------------------------------------------------------------
# SE of a normally-distributed quantity from a symmetric 95% CI already on
# the scale the quantity is normal on, assuming (upper - lower) spans 2*1.96 SE.
ci_se <- function(lower, upper) (upper - lower) / (2 * qnorm(0.975))

# Relative asymmetry of that CI (0 = symmetric). Large values flag
# misspecification of the normal-on-this-scale assumption.
ci_asym <- function(lower, mid, upper) {
  lo <- mid - lower
  hi <- upper - mid
  2 * (hi - lo) / (hi + lo)
}

# SE(log GMC) from a 95% CI reported on the natural (GMC) scale.
log_se <- function(lower, upper) ci_se(log(lower), log(upper))

# Relative asymmetry of a natural-scale CI, checked on the log scale (the
# normal-on-log EIV prior assumes symmetry there; large values flag
# misspecification, usually from aggressive rounding of small GMCs).
log_asym <- function(lower, gmc, upper) ci_asym(log(lower), log(gmc), log(upper))

SE_FLOOR <- 1e-3   # guard against zero-width (rounded) CIs -> finite precision
ASYM_TOL <- 0.10   # flag |asymmetry| > 10%

# ---------------------------------------------------------------------
# prepare_gmc_predictor(): read one merged CSV, subset to `predictor_study`,
# align the unimmunized/immunized arms by serotype, and build the log-GMC
# error-in-variables predictor (lgmc_u/i, prec_u/i) for that one source.
#
#   data_file       : merged CSV (serotype x arm x study rows)
#   predictor_study : value of the `Study` column supplying the immunogenicity
#                     GMCs used as the predictor (e.g. "NCKP", "Am_Indian").
#   predictor_error : "se" (default) uses the SE of the mean log-GMC, i.e. the
#                     reported 95% CI taken at face value - the EIV precision
#                     reflects how well the GROUP MEAN is known. "sd" instead
#                     scales that SE up to the individual-level population SD
#                     (SD = SE * sqrt(N_immuno), the analyzed immunogenicity
#                     sample size), so the EIV layer absorbs subject-to-subject
#                     titer heterogeneity rather than just sampling error of
#                     the mean - this is the more conservative choice for a
#                     COP claimed at the individual level, and it widens the
#                     predictor error a lot (heavier attenuation of the slope).
#                     Requires an `N_immuno` column in data_file; not all
#                     merged CSVs carry one yet.
#   quiet           : suppress the diagnostic messages (asymmetry / SE floor).
# ---------------------------------------------------------------------
prepare_gmc_predictor <- function(data_file, predictor_study,
                                  predictor_error = c("se", "sd"), quiet = FALSE) {
  predictor_error <- match.arg(predictor_error)
  d <- read.csv(data_file, stringsAsFactors = FALSE)
  d <- subset(d, Study == predictor_study)
  if (nrow(d) == 0) {
    stop(sprintf("No rows with Study == '%s' in %s", predictor_study, data_file))
  }
  if (predictor_error == "sd" && is.null(d$N_immuno)) {
    stop("predictor_error = 'sd' requires an N_immuno column (analyzed ",
         "immunogenicity sample size) in ", data_file)
  }

  u <- subset(d, Vaccine_Arm == "Unimmunized")
  i <- subset(d, Vaccine_Arm == "Immunized")

  serotypes <- intersect(u$Serotype, i$Serotype)
  u <- u[match(serotypes, u$Serotype), ]
  i <- i[match(serotypes, i$Serotype), ]
  S <- length(serotypes)
  stopifnot(nrow(u) == S, nrow(i) == S)

  se_u <- log_se(u$Lower_CL, u$Upper_CL)
  se_i <- log_se(i$Lower_CL, i$Upper_CL)

  # Flag CIs that are not symmetric on the log scale.
  asym_u <- log_asym(u$Lower_CL, u$GMC, u$Upper_CL)
  asym_i <- log_asym(i$Lower_CL, i$GMC, i$Upper_CL)
  flagged <- data.frame(
    Serotype = c(u$Serotype, i$Serotype),
    Arm      = c(rep("Unimmunized", S), rep("Immunized", S)),
    asym_pct = round(100 * c(asym_u, asym_i), 1)
  )
  flagged <- flagged[!is.na(flagged$asym_pct) & abs(flagged$asym_pct) > 100 * ASYM_TOL, ]
  if (!quiet && nrow(flagged)) {
    message("Note: ", nrow(flagged), " CI(s) not symmetric on the log scale ",
            "(|asymmetry| > ", 100 * ASYM_TOL, "%); normal-on-log EIV prior is ",
            "misspecified for these:")
    print(flagged[order(-abs(flagged$asym_pct)), ], row.names = FALSE)
  }

  # Floor zero/near-zero-width CIs so precision stays finite.
  n_floored <- sum(se_u <= SE_FLOOR) + sum(se_i <= SE_FLOOR)
  if (!quiet && n_floored > 0) {
    message("Note: floored ", n_floored, " zero/near-zero-width CI(s) at SE = ", SE_FLOOR)
  }
  se_u <- pmax(se_u, SE_FLOOR)
  se_i <- pmax(se_i, SE_FLOOR)

  # se_u/se_i (reported CI -> SE of the mean) stay as-is for plotting - the
  # measurement-error bars always show what was actually reported. The EIV
  # precision fed to JAGS uses a separate, possibly-scaled sd_u/sd_i so the
  # two never get confused with each other.
  sd_u <- se_u
  sd_i <- se_i
  if (predictor_error == "sd") {
    if (any(is.na(u$N_immuno)) || any(is.na(i$N_immuno))) {
      stop("N_immuno is NA for one or more serotypes in ", data_file)
    }
    sd_u <- se_u * sqrt(u$N_immuno)
    sd_i <- se_i * sqrt(i$N_immuno)
    if (!quiet) {
      message("Note: predictor_error = 'sd' - EIV precision uses SE * sqrt(N_immuno) ",
              "(population spread of individual log-titers, not the SE of the mean). ",
              "Plotted measurement-error bars still show the reported 95% CI.")
    }
  }

  list(
    jags_data = list(
      lgmc_u = log(u$GMC), prec_u = 1 / sd_u^2,
      lgmc_i = log(i$GMC), prec_i = 1 / sd_i^2
    ),
    serotypes = serotypes, S = S,
    u = u, i = i,
    se_u = se_u, se_i = se_i,
    predictor_error = predictor_error,
    x_obs = log(i$GMC) - log(u$GMC),        # observed centered log-GMC (immunized)
    flagged = flagged
  )
}

# ---------------------------------------------------------------------
# read_outcome_source(): read ONE merged CSV as an outcome dataset. Cases/
# Total_Cases are identical across every Study subset within a merged file
# (the Study column there only distinguishes immunogenicity GMC sources), so
# this collapses to one row per (Vaccine_Arm, Serotype) via unique() rather
# than picking a specific Study -- the two are interchangeable by
# construction of the merged CSVs in data/. Errors out loud if that
# assumption doesn't hold for a given file (see anyDuplicated() check below).
# ---------------------------------------------------------------------
read_outcome_source <- function(data_file, study_id, label = study_id) {
  d <- read.csv(data_file, stringsAsFactors = FALSE)
  d$Serotype <- as.character(d$Serotype)
  d <- unique(d[, c("Vaccine_Arm", "Serotype", "Cases", "Total_Cases")])

  u <- subset(d, Vaccine_Arm == "Unimmunized")
  i <- subset(d, Vaccine_Arm == "Immunized")
  if (anyDuplicated(u$Serotype) || anyDuplicated(i$Serotype)) {
    stop("read_outcome_source('", data_file, "'): Cases/Total_Cases are not ",
         "constant across Study subsets for at least one serotype/arm -- the ",
         "dedup assumption behind this function does not hold for this file.")
  }

  serotypes <- intersect(u$Serotype, i$Serotype)
  u <- u[match(serotypes, u$Serotype), ]
  i <- i[match(serotypes, i$Serotype), ]

  list(study_id = study_id, label = label,
       serotypes = serotypes, S = length(serotypes), u = u, i = i)
}

# ---------------------------------------------------------------------
# prepare_cop_data_multistudy(): assemble ragged, multi-study JAGS data for
# cop_eiv_model_multistudy.jags.
#
#   outcome_sources : list of list(data_file=, study_id=, label= [optional]),
#                     one entry per outcome dataset (merged CSV). study_id
#                     must be unique. No entry is a "reference study" -- every
#                     study gets its own a[k,s] baseline per serotype it
#                     reports, shrunk toward the shared mu_a/tau_a hyperprior.
#   immuno_sources  : list of list(data_file=, study=, label= [optional]),
#                     one entry per (file, Study) immunogenicity source, as
#                     already selected by prepare_gmc_predictor(). ALL
#                     entries pool onto the SAME shared x_u[s]/x_i[s] -- this
#                     is the "predictor set" being tested/compared via WAIC.
#   predictor_error : "se" or "sd", applied identically to every
#                     immunogenicity source (see prepare_gmc_predictor()).
#   quiet           : suppress per-source diagnostic messages.
#
# Returns a list with jags_data (ragged vectors/index arrays), the canonical
# serotype universe (serotypes, S), and the per-source bookkeeping
# (outcome_sources, immuno_sources) needed for study_summary() and
# plot_cop_absolute_multistudy().
# ---------------------------------------------------------------------
prepare_cop_data_multistudy <- function(outcome_sources, immuno_sources,
                                        predictor_error = c("se", "sd"),
                                        quiet = FALSE) {
  predictor_error <- match.arg(predictor_error)

  out_list <- lapply(outcome_sources, function(o) {
    read_outcome_source(o$data_file, o$study_id,
                        if (is.null(o$label)) o$study_id else o$label)
  })
  study_ids <- vapply(out_list, `[[`, character(1), "study_id")
  if (anyDuplicated(study_ids)) {
    stop("Duplicate outcome study_id(s): ",
         paste(study_ids[duplicated(study_ids)], collapse = ", "))
  }

  imm_list <- lapply(immuno_sources, function(p) {
    gp <- prepare_gmc_predictor(p$data_file, p$study, predictor_error, quiet)
    list(study = p$study, label = if (is.null(p$label)) p$study else p$label, gp = gp)
  })

  serotypes <- sort(unique(c(
    unlist(lapply(out_list, `[[`, "serotypes")),
    unlist(lapply(imm_list, function(x) x$gp$serotypes))
  )))
  S <- length(serotypes)
  sero_index <- setNames(seq_len(S), serotypes)
  K <- length(out_list)

  study_out <- integer(0); sero_out <- integer(0)
  cases_u <- numeric(0); total_u <- numeric(0)
  cases_i <- numeric(0); total_i <- numeric(0)
  for (k in seq_len(K)) {
    o <- out_list[[k]]
    if (o$S == 0) next
    study_out <- c(study_out, rep(k, o$S))
    sero_out  <- c(sero_out, unname(sero_index[o$serotypes]))
    cases_u   <- c(cases_u, o$u$Cases); total_u <- c(total_u, o$u$Total_Cases)
    cases_i   <- c(cases_i, o$i$Cases); total_i <- c(total_i, o$i$Total_Cases)
  }

  sero_imm_u <- integer(0); lgmc_u <- numeric(0); prec_imm_u <- numeric(0)
  sero_imm_i <- integer(0); lgmc_i <- numeric(0); prec_imm_i <- numeric(0)
  for (im in imm_list) {
    gp <- im$gp
    if (gp$S == 0) next
    sero_imm_u <- c(sero_imm_u, unname(sero_index[gp$serotypes]))
    lgmc_u     <- c(lgmc_u, gp$jags_data$lgmc_u)
    prec_imm_u <- c(prec_imm_u, gp$jags_data$prec_u)
    sero_imm_i <- c(sero_imm_i, unname(sero_index[gp$serotypes]))
    lgmc_i     <- c(lgmc_i, gp$jags_data$lgmc_i)
    prec_imm_i <- c(prec_imm_i, gp$jags_data$prec_i)
  }

  jags_data <- list(
    S = S, K = K,
    N_out = length(sero_out), study_out = study_out, sero_out = sero_out,
    cases_u = cases_u, total_u = total_u, cases_i = cases_i, total_i = total_i,
    N_imm_u = length(sero_imm_u), sero_imm_u = sero_imm_u,
    lgmc_u = lgmc_u, prec_imm_u = prec_imm_u,
    N_imm_i = length(sero_imm_i), sero_imm_i = sero_imm_i,
    lgmc_i = lgmc_i, prec_imm_i = prec_imm_i
  )

  list(jags_data = jags_data, serotypes = serotypes, S = S, K = K,
       study_labels = vapply(out_list, `[[`, character(1), "label"),
       outcome_sources = out_list, immuno_sources = imm_list,
       predictor_error = predictor_error)
}

# ---------------------------------------------------------------------
# fit_cop_multistudy(): compile, burn in, sample, summarise, and save
# outputs for the multi-study model (cop_eiv_model_multistudy.jags).
#
#   prep     : output of prepare_cop_data_multistudy()
#   out_dir  : results folder for this model spec (created if missing)
#   n.*      : MCMC controls
#   seed     : master seed (per-chain RNG seeds are derived from it)
#
# Writes into out_dir: posterior_summary.csv, mcmc.rds, diagnostics.pdf,
# waic.csv (see compute_waic()).
# Returns (invisibly) a list with the samples, convergence stats, and WAIC.
# ---------------------------------------------------------------------
MULTISTUDY_MODEL_FILE <- file.path("JAGS", "cop_eiv_model_multistudy.jags")

MULTISTUDY_PARAMS <- c("mu_b1", "b1", "rr", "rr_global",
                       "a", "mu_a", "sigma_a", "sigma_b1",
                       "mu_x_u", "sigma_x_u", "x_u",
                       "mu_x_i", "sigma_x_i", "x_i",
                       "log_lik")

# Data-driven inits: per-(study, serotype) baseline computed from whichever
# outcome row actually reports that combination (combinations with no row --
# i.e. that study didn't report that serotype -- fall back to the overall
# mean; the hierarchical prior on a[k,s] does the rest). Immunogenicity
# latents still init from per-serotype means across whichever immuno rows
# reference that serotype.
make_inits_multistudy <- function(prep) {
  jd <- prep$jags_data
  S <- prep$S; K <- prep$K

  mean_by_sero <- function(sero_idx, value) {
    out <- vapply(seq_len(S), function(s) {
      rows <- which(sero_idx == s)
      if (length(rows)) mean(value[rows]) else NA_real_
    }, numeric(1))
    out[is.na(out)] <- mean(out, na.rm = TRUE)
    out
  }

  a_start <- matrix(NA_real_, K, S)
  a_start[cbind(jd$study_out, jd$sero_out)] <- log(jd$cases_u / jd$total_u)
  a_start[is.na(a_start)] <- mean(a_start, na.rm = TRUE)

  x_u_start <- mean_by_sero(jd$sero_imm_u, jd$lgmc_u)
  x_i_start <- mean_by_sero(jd$sero_imm_i, jd$lgmc_i)

  list(
    a = a_start, b1 = rep(-0.5, S),
    x_u = x_u_start, x_i = x_i_start,
    mu_a = mean(a_start), mu_b1 = -0.5,
    mu_x_u = mean(x_u_start), mu_x_i = mean(x_i_start)
  )
}

# ---------------------------------------------------------------------
# compute_waic(): WAIC computed from the model's monitored per-outcome-row
# log-likelihood ("log_lik[1..N_out]" -- see JAGS/cop_eiv_model_multistudy.jags
# for why only the outcome likelihood is used). Hand-rolled (no `loo`
# dependency): standard pointwise WAIC formulas, log-sum-exp for numerical
# stability in lppd.
#
#   samp    : mcmc.list from fit_cop_multistudy() (must monitor "log_lik")
#   out_dir : if given, writes a one-row waic.csv
#
# Returns list(waic, elpd_waic, p_waic, se_waic, n_obs).
# ---------------------------------------------------------------------
compute_waic <- function(samp, out_dir = NULL) {
  M <- as.matrix(samp)
  ll_cols <- grep("^log_lik\\[", colnames(M))
  if (length(ll_cols) == 0) {
    stop("compute_waic() requires 'log_lik' to be monitored in the fit.")
  }
  LL <- M[, ll_cols, drop = FALSE]     # draws x N_out
  n_draws <- nrow(LL)

  logmeanexp <- function(x) {
    mx <- max(x)
    mx + log(mean(exp(x - mx)))
  }

  lppd_i   <- apply(LL, 2, logmeanexp)
  p_waic_i <- apply(LL, 2, var)
  elpd_i   <- lppd_i - p_waic_i

  elpd_waic <- sum(elpd_i)
  p_waic    <- sum(p_waic_i)
  waic      <- -2 * elpd_waic
  se_waic   <- sqrt(length(elpd_i) * var(elpd_i)) * sqrt(2)

  out <- list(waic = waic, elpd_waic = elpd_waic, p_waic = p_waic,
              se_waic = se_waic, n_obs = length(elpd_i), n_draws = n_draws)

  if (!is.null(out_dir)) {
    write.csv(as.data.frame(out), file.path(out_dir, "waic.csv"), row.names = FALSE)
  }
  out
}

fit_cop_multistudy <- function(prep, out_dir,
                               n.adapt = 2000, n.burnin = 5000, n.iter = 20000, thin = 5,
                               n.chains = 3, seed = 42, quiet = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  jd <- prep$jags_data

  if (!quiet) {
    cat("Multi-study model (", MULTISTUDY_MODEL_FILE, ")\n", sep = "")
    cat("Serotypes (", prep$S, "):", paste(prep$serotypes, collapse = ", "), "\n")
    cat("Outcome studies (", prep$K, ", each with its own a[k,s] baseline):\n", sep = "")
    print(setNames(prep$study_labels, seq_along(prep$study_labels)))
    cat("Immunogenicity sources pooled onto shared x_u/x_i:", length(prep$immuno_sources), "\n")
  }

  base_inits <- make_inits_multistudy(prep)
  make_inits <- function(chain_seed) {
    c(base_inits, list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = chain_seed))
  }
  inits <- lapply(seq_len(n.chains), function(k) make_inits(11 * k))

  model <- jags.model(MULTISTUDY_MODEL_FILE, data = jd,
                      inits = inits,
                      n.chains = n.chains, n.adapt = n.adapt,
                      quiet = quiet)
  update(model, n.burnin)
  samp <- coda.samples(model, variable.names = MULTISTUDY_PARAMS,
                       n.iter = n.iter, thin = thin)

  summ  <- summary(samp)
  stats <- cbind(summ$statistics[, c("Mean", "SD")], summ$quantiles)

  M <- as.matrix(samp)
  mu_b1_draws <- M[, "mu_b1"]
  rr_draws    <- M[, "rr_global"]

  es <- effectiveSize(samp)

  gd  <- gelman.diag(samp, multivariate = FALSE)
  max_psrf <- max(gd$psrf[, "Point est."], na.rm = TRUE)
  min_ess  <- min(es)

  waic <- compute_waic(samp, out_dir = out_dir)

  if (!quiet) {
    cat("\n--- Global correlate of protection ---\n")
    cat(sprintf("mu_b1 (log RR per unit absolute log-GMC): %.3f  (95%% CrI %.3f, %.3f)\n",
                median(mu_b1_draws), quantile(mu_b1_draws, .025), quantile(mu_b1_draws, .975)))
    cat(sprintf("rr_global (rate ratio per unit absolute log-GMC): %.3f  (95%% CrI %.3f, %.3f)\n",
                median(rr_draws), quantile(rr_draws, .025), quantile(rr_draws, .975)))
    cat(sprintf("P(mu_b1 < 0) = %.3f\n", mean(mu_b1_draws < 0)))
    cat("\n--- Convergence ---\n")
    cat("Max PSRF (all monitored):", round(max_psrf, 3), "\n")
    cat("Min effective sample size:", round(min_ess), "\n")
    cat("\n--- WAIC (outcome likelihood only) ---\n")
    cat(sprintf("WAIC = %.2f (SE %.2f), elpd_waic = %.2f, p_waic = %.2f, n_obs = %d\n",
                waic$waic, waic$se_waic, waic$elpd_waic, waic$p_waic, waic$n_obs))
  }

  # posterior_summary.csv would otherwise carry N_out rows of log_lik[m] --
  # drop those before writing (they belong in waic.csv, not the parameter table).
  keep <- !grepl("^log_lik\\[", rownames(stats))
  write.csv(round(stats[keep, , drop = FALSE], 4), file.path(out_dir, "posterior_summary.csv"))
  saveRDS(samp, file.path(out_dir, "mcmc.rds"))

  diag_params <- c("mu_b1", "sigma_b1", "mu_a", "sigma_a")
  pdf(file.path(out_dir, "diagnostics.pdf"), width = 9, height = 6)
  plot(samp[, diag_params])
  dev.off()

  if (!quiet) cat("\nSaved fit outputs to ", out_dir, "\n", sep = "")

  invisible(list(samp = samp, stats = stats,
                 max_psrf = max_psrf, min_ess = min_ess, waic = waic))
}

# ---------------------------------------------------------------------
# slope_summary(): tidy per-model slope table (global + per serotype),
# derived from an mcmc.rds sample. One row per level; written to
# out_dir/slope_summary.csv and returned as a data frame.
# ---------------------------------------------------------------------
slope_summary <- function(samp, prep, analysis_id, predictor_label, outcome_label,
                          out_dir = NULL) {
  M <- as.matrix(samp)
  qtiles <- c(0.025, 0.5, 0.975)

  row_from <- function(level, b1_draws) {
    q <- quantile(b1_draws, qtiles)
    data.frame(
      analysis_id     = analysis_id,
      predictor_label = predictor_label,
      outcome_label   = outcome_label,
      level           = level,
      b1_mean         = mean(b1_draws),
      b1_median       = q[["50%"]],
      b1_lo           = q[["2.5%"]],
      b1_hi           = q[["97.5%"]],
      rr_median       = exp(q[["50%"]]),
      rr_lo           = exp(q[["2.5%"]]),
      rr_hi           = exp(q[["97.5%"]]),
      p_negative      = mean(b1_draws < 0),
      stringsAsFactors = FALSE
    )
  }

  rows <- list(row_from("global", M[, "mu_b1"]))
  for (s in seq_len(prep$S)) {
    rows[[length(rows) + 1]] <- row_from(prep$serotypes[s], M[, sprintf("b1[%d]", s)])
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  if (!is.null(out_dir)) {
    write.csv(out, file.path(out_dir, "slope_summary.csv"), row.names = FALSE)
  }
  out
}

# ---------------------------------------------------------------------
# study_summary(): tidy per-outcome-study baseline table. There is no single
# per-study fixed effect (a[k,s] varies by serotype within a study too), so
# this summarises each study's OWN average baseline -- mean_s(a[k,s]) over
# just the serotypes that study reports -- and how that compares to the
# global hypermean mu_a (the pooled baseline across every study/serotype).
# rate_ratio_vs_global > 1 means this study's serotypes run hotter (higher
# case ascertainment / incidence) than the pooled average; < 1 means cooler.
# Written to out_dir/study_summary.csv (if out_dir given) and returned as a
# data frame.
# ---------------------------------------------------------------------
study_summary <- function(samp, prep, out_dir = NULL) {
  M <- as.matrix(samp)
  qtiles <- c(0.025, 0.5, 0.975)

  rows <- lapply(seq_len(prep$K), function(k) {
    o <- prep$outcome_sources[[k]]
    sero_idx <- match(o$serotypes, prep$serotypes)
    a_cols <- sprintf("a[%d,%d]", k, sero_idx)
    a_bar_draws <- if (length(a_cols) > 1) rowMeans(M[, a_cols, drop = FALSE]) else M[, a_cols]
    rr_draws <- exp(a_bar_draws - M[, "mu_a"])

    q    <- quantile(a_bar_draws, qtiles)
    rr_q <- quantile(rr_draws, qtiles)
    data.frame(
      study_id       = o$study_id,
      label          = prep$study_labels[k],
      n_serotypes    = o$S,
      a_bar_mean     = mean(a_bar_draws),
      a_bar_median   = q[["50%"]],
      a_bar_lo       = q[["2.5%"]],
      a_bar_hi       = q[["97.5%"]],
      rate_ratio_vs_global_median = rr_q[["50%"]],
      rate_ratio_vs_global_lo     = rr_q[["2.5%"]],
      rate_ratio_vs_global_hi     = rr_q[["97.5%"]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  if (!is.null(out_dir)) {
    write.csv(out, file.path(out_dir, "study_summary.csv"), row.names = FALSE)
  }
  out
}

# ---------------------------------------------------------------------
# plot_cop_absolute_multistudy(): correlate-of-protection relationship in
# ABSOLUTE coordinates, faceted by outcome study:
#
#   x = POOLED absolute GMC per serotype/arm -- the model's posterior
#       estimate of the SHARED latent x_u[s]/x_i[s] (exponentiated median),
#       i.e. the combined estimate across every immunogenicity source. There
#       is no longer a single "observed" GMC once multiple sources are
#       pooled, so the x-axis (and its horizontal error bar, the posterior
#       95% credible interval of that latent) is a model estimate, not a raw
#       observation.
#   y = observed absolute risk (case proportion), per OUTCOME STUDY -- one
#       facet per outcome study k, each facet showing only the serotypes
#       that study reports, plotted against the SAME shared/pooled x-axis
#       (so a point's horizontal position never changes across facets).
#
# Overlay per facet: each serotype's own fitted line
#   log(risk) = a[k,s] + b1[s] * log(GMC)
# across its own pooled absolute-GMC range, and that study's own
# population-average curve
#   log(risk) = mean_s(a[k,s]) + mu_b1 * log(GMC)
# with a 95% credible band -- i.e. the SAME global COP slope (mu_b1) in every
# facet, but each facet's intercept is that study's OWN average baseline
# (mean of a[k,s] over just the serotypes it reports -- see study_summary()),
# since a[k,s] no longer decomposes into a single global a[s] plus a per-study
# shift.
#
# Saves cop_scatter_absolute_risk_gmc_multistudy.png.
# ---------------------------------------------------------------------
plot_cop_absolute_multistudy <- function(prep, samp, out_dir, title_suffix = "") {
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  serotypes <- prep$serotypes
  S <- prep$S
  K <- prep$K
  z <- qnorm(0.975)

  M <- as.matrix(samp)
  a_pairs <- do.call(rbind, lapply(seq_len(K), function(k) {
    o <- prep$outcome_sources[[k]]
    if (o$S == 0) return(NULL)
    data.frame(k = k, s = match(o$serotypes, serotypes))
  }))
  need <- c("mu_a", "mu_b1", sprintf("a[%d,%d]", a_pairs$k, a_pairs$s),
           sprintf("b1[%d]", seq_len(S)),
           sprintf("x_u[%d]", seq_len(S)), sprintf("x_i[%d]", seq_len(S)))
  if (!all(need %in% colnames(M))) {
    stop("plot_cop_absolute_multistudy() requires a fit from fit_cop_multistudy() -- ",
         "missing: ", paste(setdiff(need, colnames(M)), collapse = ", "))
  }

  # ---- Pooled absolute log-GMC per serotype/arm: posterior of the shared --
  # ---- latent x_u[s]/x_i[s], combining every immunogenicity source. -------
  x_u_med <- vapply(seq_len(S), function(s) median(M[, sprintf("x_u[%d]", s)]), numeric(1))
  x_i_med <- vapply(seq_len(S), function(s) median(M[, sprintf("x_i[%d]", s)]), numeric(1))
  x_u_lo  <- vapply(seq_len(S), function(s) quantile(M[, sprintf("x_u[%d]", s)], 0.025), numeric(1))
  x_u_hi  <- vapply(seq_len(S), function(s) quantile(M[, sprintf("x_u[%d]", s)], 0.975), numeric(1))
  x_i_lo  <- vapply(seq_len(S), function(s) quantile(M[, sprintf("x_i[%d]", s)], 0.025), numeric(1))
  x_i_hi  <- vapply(seq_len(S), function(s) quantile(M[, sprintf("x_i[%d]", s)], 0.975), numeric(1))

  b1_hat <- vapply(seq_len(S), function(s) mean(M[, sprintf("b1[%d]", s)]), numeric(1))

  # ---- Per-facet (outcome study) observed risk points, count sampling ----
  # ---- error (Haldane 0.5 correction for zero-count arms) ----------------
  pts_list <- lapply(seq_len(K), function(k) {
    o <- prep$outcome_sources[[k]]
    if (o$S == 0) return(NULL)
    idx <- match(o$serotypes, serotypes)
    zero_u <- o$u$Cases == 0; zero_i <- o$i$Cases == 0
    cu <- ifelse(zero_u, o$u$Cases + 0.5, o$u$Cases)
    ci <- ifelse(zero_i, o$i$Cases + 0.5, o$i$Cases)
    risk_u <- cu / o$u$Total_Cases
    risk_i <- ci / o$i$Total_Cases
    se_risk_u <- sqrt(1 / cu)
    se_risk_i <- sqrt(1 / ci)
    data.frame(
      Study = prep$study_labels[k], study_idx = k,
      Serotype = rep(o$serotypes, 2),
      Arm      = rep(c("Unimmunized", "Immunized"), each = o$S),
      gmc      = exp(c(x_u_med[idx], x_i_med[idx])),
      gmc_lo   = exp(c(x_u_lo[idx],  x_i_lo[idx])),
      gmc_hi   = exp(c(x_u_hi[idx],  x_i_hi[idx])),
      risk     = c(risk_u, risk_i),
      risk_lo  = c(risk_u * exp(-z * se_risk_u), risk_i * exp(-z * se_risk_i)),
      risk_hi  = c(risk_u * exp( z * se_risk_u), risk_i * exp( z * se_risk_i)),
      zero     = c(zero_u, zero_i),
      stringsAsFactors = FALSE
    )
  })
  pts <- do.call(rbind, pts_list)

  # ---- Per-serotype fitted line, per facet: a[k,s] + b1[s]*log(GMC) -------
  sero_lines_list <- lapply(seq_len(K), function(k) {
    o <- prep$outcome_sources[[k]]
    if (o$S == 0) return(NULL)
    idx <- match(o$serotypes, serotypes)
    do.call(rbind, lapply(seq_along(idx), function(j) {
      s <- idx[j]
      a_ks <- mean(M[, sprintf("a[%d,%d]", k, s)])
      xg <- seq(min(x_u_med[s], x_i_med[s]), max(x_u_med[s], x_i_med[s]), length.out = 50)
      data.frame(Study = prep$study_labels[k], study_idx = k,
                Serotype = o$serotypes[j], gmc = exp(xg),
                risk = exp(a_ks + b1_hat[s] * xg),
                stringsAsFactors = FALSE)
    }))
  })
  sero_lines <- do.call(rbind, sero_lines_list)

  # ---- Per-study population-average curve + 95% band, per facet ----------
  # Intercept = mean_s(a[k,s]) over just the serotypes THIS study reports
  # (same quantity as study_summary()'s a_bar), not the global mu_a -- a[k,s]
  # no longer decomposes into a shared a[s] plus a per-study shift.
  xr  <- range(c(x_u_med, x_i_med))
  pad <- 0.1 * diff(xr)
  xg  <- seq(xr[1] - pad, xr[2] + pad, length.out = 200)
  fit_glob_list <- lapply(seq_len(K), function(k) {
    o <- prep$outcome_sources[[k]]
    if (o$S == 0) return(NULL)
    idx <- match(o$serotypes, serotypes)
    a_cols <- sprintf("a[%d,%d]", k, idx)
    a_bar_draws <- if (length(a_cols) > 1) rowMeans(M[, a_cols, drop = FALSE]) else M[, a_cols]
    line_mat <- outer(a_bar_draws, rep(1, length(xg))) + outer(M[, "mu_b1"], xg)
    data.frame(
      Study = prep$study_labels[k], study_idx = k,
      gmc = exp(xg),
      med = exp(apply(line_mat, 2, median)),
      lo  = exp(apply(line_mat, 2, quantile, 0.025)),
      hi  = exp(apply(line_mat, 2, quantile, 0.975))
    )
  })
  fit_glob <- do.call(rbind, fit_glob_list)

  study_order <- prep$study_labels
  pts$Study       <- factor(pts$Study,       levels = study_order)
  sero_lines$Study <- factor(sero_lines$Study, levels = study_order)
  fit_glob$Study  <- factor(fit_glob$Study,  levels = study_order)

  serotypes_used <- sort(unique(pts$Serotype))
  slope_labels <- setNames(sprintf("%s  (post. b1 = %+.2f)", serotypes, b1_hat),
                          serotypes)[serotypes_used]

  main_title <- "Correlate of protection: absolute risk vs pooled absolute GMC"
  if (nzchar(title_suffix)) main_title <- paste0(main_title, "\n", title_suffix)
  subt <- paste0("Points = observed serotype x arm per outcome study (circle = unimmunized, triangle = immunized); ",
                "x = pooled posterior GMC estimate (all immunogenicity sources combined); ",
                "dotted lines = each serotype's own fitted risk-GMC relationship; ",
                "black line = population-average fit per study (that study's own mean a[k,s], mu_b1) with 95% credible band")
  cap <- if (any(pts$zero))
    "Open symbols: an arm has zero cases; risk uses a 0.5 continuity correction." else NULL

  if (has_ggplot) {
    library(ggplot2)
    p <- ggplot() +
      geom_vline(xintercept = 0.35, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
      geom_ribbon(data = fit_glob, aes(x = gmc, ymin = lo, ymax = hi),
                  fill = "grey60", alpha = 0.25) +
      geom_line(data = fit_glob, aes(x = gmc, y = med), colour = "black", linewidth = 1.2) +
      geom_line(data = sero_lines, aes(x = gmc, y = risk, colour = Serotype),
               linetype = "dotted", alpha = 0.6) +
      geom_errorbarh(data = pts, aes(y = risk, xmin = gmc_lo, xmax = gmc_hi, colour = Serotype),
                     height = 0, alpha = 0.7) +
      geom_errorbar(data = pts, aes(x = gmc, ymin = risk_lo, ymax = risk_hi, colour = Serotype),
                    width = 0, alpha = 0.7) +
      geom_point(data = pts, aes(x = gmc, y = risk, colour = Serotype, shape = Arm), size = 2.8) +
      scale_shape_manual(values = c(Unimmunized = 16, Immunized = 17)) +
      scale_x_log10() + scale_y_log10() +
      scale_colour_discrete(labels = slope_labels) +
      facet_wrap(~ Study, ncol = 1) +
      labs(
        x = "Pooled absolute GMC (posterior estimate)", y = "Absolute risk (case proportion)",
        colour = "Serotype (slope b1)", shape = "Arm",
        title = main_title, subtitle = subt, caption = cap
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right", strip.text = element_text(face = "bold"))
    ggsave(file.path(out_dir, "cop_scatter_absolute_risk_gmc_multistudy.png"), p,
          width = 9, height = 3.3 * K + 1.5, dpi = 150, limitsize = FALSE)
  } else {
    cols <- setNames(rainbow(S), serotypes)
    xlim <- range(c(pts$gmc_lo, pts$gmc_hi, fit_glob$gmc))
    ylim <- range(c(pts$risk_lo, pts$risk_hi, fit_glob$lo, fit_glob$hi))
    png(file.path(out_dir, "cop_scatter_absolute_risk_gmc_multistudy.png"),
       width = 9, height = 3.3 * K + 1.5, units = "in", res = 150)
    par(mfrow = c(K, 1), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (k in seq_len(K)) {
      lab <- study_order[k]
      fg <- fit_glob[fit_glob$Study == lab, ]
      sl <- sero_lines[sero_lines$Study == lab, ]
      pk <- pts[pts$Study == lab, ]
      plot(NA, xlim = xlim, ylim = ylim, log = "xy",
           xlab = "Pooled absolute GMC", ylab = "Absolute risk",
           main = lab)
      abline(v = 0.35, lty = 2, col = "grey40")
      polygon(c(fg$gmc, rev(fg$gmc)), c(fg$lo, rev(fg$hi)),
              col = adjustcolor("grey60", 0.25), border = NA)
      lines(fg$gmc, fg$med, col = "black", lwd = 2)
      for (s in unique(sl$Serotype)) {
        ss <- sl[sl$Serotype == s, ]
        lines(ss$gmc, ss$risk, col = adjustcolor(cols[s], alpha.f = 0.6), lty = 3)
      }
      pch_arm <- ifelse(pk$zero, ifelse(pk$Arm == "Unimmunized", 1, 2),
                        ifelse(pk$Arm == "Unimmunized", 16, 17))
      arrows(pk$gmc_lo, pk$risk, pk$gmc_hi, pk$risk, code = 3, angle = 90, length = 0.03,
             col = cols[pk$Serotype])
      arrows(pk$gmc, pk$risk_lo, pk$gmc, pk$risk_hi, code = 3, angle = 90, length = 0.03,
             col = cols[pk$Serotype])
      points(pk$gmc, pk$risk, col = cols[pk$Serotype], pch = pch_arm, cex = 1.3)
      if (k == 1) {
        legend("topleft", legend = slope_labels, col = cols[names(slope_labels)], pch = 16,
               bty = "n", cex = 0.7, title = "Serotype (slope b1)")
      }
    }
    mtext(main_title, outer = TRUE, cex = 1, font = 2)
    dev.off()
  }
  invisible(NULL)
}

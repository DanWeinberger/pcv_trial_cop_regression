# =====================================================================
# Core reusable engine for the pneumococcal correlate-of-protection (COP)
# error-in-variables Poisson regression.
#
# This file defines FUNCTIONS only (it fits nothing when sourced). It is the
# single place where the modelling logic lives, so a new comparison never
# requires copy-pasting the analysis - you just point the functions at a
# different predictor/outcome via R/config.R and R/run_analysis.R.
#
#   prepare_cop_data() : build the JAGS data list for one (outcome, predictor)
#   fit_cop()          : compile + sample the JAGS model, save summaries
#   plot_cop()         : draw the fitted case-proportion curves
#   slope_summary()    : tidy per-analysis slope table (global + per serotype)
#
# Model structure (unchanged from the original single-study analysis):
#   Outcome    : serotype-specific IPD case counts (Cases), Poisson
#   Offset     : Total_Cases (log-exposure)
#   Predictor  : centered log-GMC, measured with error (precision from 95% CI)
#   Hierarchy  : serotype-specific random intercepts a[s] and slopes b1[s]
#   Centering  : in-model latent reference per serotype -> the unimmunized true
#                log-GMC is the reference (x = 0); its measurement error
#                propagates into the immunized predictor.
# =====================================================================

suppressPackageStartupMessages({
  library(rjags)
  library(coda)
})

# ---------------------------------------------------------------------
# Parameterizations
# ---------------------------------------------------------------------
# Four JAGS models are supported: "centered" and "ratio" are algebraically-
# equivalent reparameterizations of the same case-count outcome; "logor" is a
# distinct model where the outcome is a directly-supplied log-OR instead of
# case counts (see cop_eiv_model_logor.jags); "RE" is ALSO a distinct model --
# it regresses each arm's ABSOLUTE log-GMC against its absolute log case-rate
# (random effects on the absolute GMC scale, by arm), rather than the
# within-serotype paired ratio (see cop_eiv_model_RE.jags for why its b1[s] is
# still comparable to the other three). The rest of the pipeline (plots,
# slope tables, comparisons) speaks a single CANONICAL set of variable names;
# each parameterization declares how its own JAGS variables map onto those
# canonical names. fit_cop() monitors the model's OWN names and then renames
# the posterior samples back to canonical, so nothing downstream needs to know
# which parameterization produced a fit (mcmc.rds is always canonical).
#
#   canonical  meaning                                centered   ratio      logor      RE
#   ---------  -------------------------------------  --------   --------   --------   --------
#   a[s]       serotype baseline log case-proportion  a[s]       logp0[s]   (none)     a[s]
#   mu_a       its hypermean                          mu_a       mu_p0      (none)     mu_a
#   sigma_a    its hyper-SD                           sigma_a    sigma_p0   (none)     sigma_a
#   x[s]       latent log-GMC ratio (predictor)       x[s]       lrgmc[s]   lrgmc[s]   (none; see x_u/x_i)
#   mu_x       its hypermean                          mu_x       mu_lr      mu_lr      (none)
#   sigma_x    its hyper-SD                           sigma_x    sigma_lr   sigma_lr   (none)
#   ref[s]     latent unimmunized true log-GMC        ref[s]     g0[s]      g0[s]      (none; see x_u)
#   rr[s]      ratio per unit predictor                rr[s]      rr[s]      or[s]      rr[s]
#   rr_global  pooled ratio                            rr_global  rr_global  or_global  rr_global
#   b1[s], mu_b1, sigma_b1  --- identical in all four ---
#
# The "logor" parameterization has no baseline (a/mu_a/sigma_a): an odds ratio
# already nets out the serotype's baseline risk, so there is nothing to
# estimate there. "RE" has no single centered "x"/"ref": the absolute log-GMC
# of each arm is its own latent (x_u[s], x_i[s], with hypermeans mu_x_u/mu_x_i
# and hyper-SDs sigma_x_u/sigma_x_i) -- these are ADDITIONAL RE-only canonical
# names, not aliases of x/mu_x/sigma_x, because there is no ratio predictor to
# rename onto. Each spec's own `params` and `inits` entries account for this
# (see below); fit_cop() no longer assumes a fixed monitor/inits list.
PARAMETERIZATIONS <- list(
  centered = list(
    model_file = file.path("JAGS", "cop_eiv_model.jags"),
    map = c(a = "a", mu_a = "mu_a", sigma_a = "sigma_a",
            x = "x", mu_x = "mu_x", sigma_x = "sigma_x", ref = "ref")
  ),
  ratio = list(
    model_file = file.path("JAGS", "cop_eiv_model_ratio.jags"),
    map = c(a = "logp0", mu_a = "mu_p0", sigma_a = "sigma_p0",
            x = "lrgmc", mu_x = "mu_lr", sigma_x = "sigma_lr", ref = "g0"),
    # sigma_p0/sigma_b1/sigma_lr are LOGICAL nodes here (derived from a
    # dgamma-on-precision prior, unlike centered's dnorm-on-SD prior) -- JAGS
    # errors if you supply an init for a logical node, so this parameterization
    # gets its own inits builder that omits them (tau_* auto-initializes from
    # its own prior).
    make_inits = function(prep) {
      a_start <- log(prep$u$Cases / prep$u$Total_Cases)
      list(ref = log(prep$u$GMC), x = prep$x_obs,
           a = a_start, b1 = rep(-0.5, prep$S),
           mu_a = mean(a_start), mu_b1 = -0.5,
           mu_x = mean(prep$x_obs))
    }
  ),
  logor = list(
    model_file = file.path("JAGS", "cop_eiv_model_logor.jags"),
    map = c(x = "lrgmc", mu_x = "mu_lr", sigma_x = "sigma_lr", ref = "g0",
            rr = "or", rr_global = "or_global"),
    # No baseline (a/mu_a/sigma_a) to monitor for this parameterization.
    params = c("mu_b1", "b1", "rr", "rr_global",
               "sigma_b1", "mu_x", "sigma_x", "x", "ref"),
    diag_params = c("mu_b1", "sigma_b1", "mu_x", "sigma_x"),
    # Data-driven inits built from prep (no Cases/GMC baseline available).
    # sigma_b1/sigma_x are LOGICAL nodes (derived from a dgamma-on-precision
    # prior) -- omitted here for the same reason as the "ratio" entry above.
    make_inits = function(prep) {
      list(ref = log(prep$u$GMC), x = prep$x_obs,
           b1 = rep(-0.5, prep$S),
           mu_b1 = -0.5, mu_x = mean(prep$x_obs))
    }
  ),
  RE = list(
    model_file = file.path("JAGS", "cop_eiv_model_RE.jags"),
    # A GENUINELY DIFFERENT model, not an algebraic reparameterization: a[s]/
    # b1[s]/mu_a/mu_b1/rr/rr_global are already literal canonical names (no
    # renaming needed), but the predictor is split into two absolute-scale
    # latents (x_u[s], x_i[s], each its own random effect) instead of a single
    # centered/ratio "x"/"ref" pair -- see cop_eiv_model_RE.jags for why this
    # still yields a b1[s] directly comparable to the other parameterizations,
    # while ALSO supporting absolute risk-vs-GMC prediction (plot_cop_absolute()).
    map = c(),
    params = c("mu_b1", "b1", "rr", "rr_global",
               "a", "mu_a", "sigma_a", "sigma_b1",
               "mu_x_u", "sigma_x_u", "x_u", "mu_x_i", "sigma_x_i", "x_i"),
    # sigma_a/sigma_b1/sigma_x_u/sigma_x_i are LOGICAL nodes (derived from a
    # dgamma-on-precision prior) -- omitted from inits for the same reason as
    # the "ratio" entry above.
    make_inits = function(prep) {
      a_start <- log(prep$u$Cases / prep$u$Total_Cases)
      list(x_u = log(prep$u$GMC), x_i = log(prep$i$GMC),
           a = a_start, b1 = rep(-0.5, prep$S),
           mu_a = mean(a_start), mu_b1 = -0.5,
           mu_x_u = mean(log(prep$u$GMC)), mu_x_i = mean(log(prep$i$GMC)))
    }
  )
)

# Default monitor list / diagnostics panel used by "centered" and "ratio"
# (both have the full baseline + predictor + slope hierarchy).
DEFAULT_PARAMS <- c("mu_b1", "b1", "rr", "rr_global",
                    "a", "mu_a", "sigma_a", "sigma_b1",
                    "mu_x", "sigma_x", "x", "ref")
DEFAULT_DIAG_PARAMS <- c("mu_b1", "sigma_b1", "mu_a", "sigma_a")

# Default data-driven inits builder for "centered" and "ratio": uses the
# observed unimmunized case-proportion as the baseline start value.
default_make_inits <- function(prep) {
  list(ref = log(prep$u$GMC), x = prep$x_obs,
       a = log(prep$u$Cases / prep$u$Total_Cases), b1 = rep(-0.5, prep$S),
       mu_a = mean(log(prep$u$Cases / prep$u$Total_Cases)), sigma_a = 0.3,
       mu_b1 = -0.5, sigma_b1 = 0.3,
       mu_x = mean(prep$x_obs), sigma_x = 0.5)
}

# Rename inits keys from canonical -> model-specific (leaving .RNG.* / unmapped
# entries untouched).
canon_to_model_inits <- function(inits, map) {
  nm  <- names(inits)
  hit <- nm %in% names(map)
  nm[hit] <- map[nm[hit]]
  names(inits) <- nm
  inits
}

# Map canonical monitor names -> model-specific names.
canon_to_model_params <- function(params, map) {
  ifelse(params %in% names(map), map[params], params)
}

# Rename mcmc.list columns from model-specific -> canonical, preserving any
# [index] suffix (e.g. logp0[3] -> a[3]).
model_to_canon_samples <- function(samp, map) {
  rev_map <- setNames(names(map), unname(map))    # model-specific -> canonical
  rename_one <- function(v) {
    base <- sub("\\[.*$", "", v)
    idx  <- substring(v, nchar(base) + 1L)        # "" or "[i]" / "[i,j]"
    ifelse(base %in% names(rev_map), paste0(rev_map[base], idx), v)
  }
  as.mcmc.list(lapply(samp, function(ch) {
    colnames(ch) <- rename_one(colnames(ch))
    ch
  }))
}

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
# prepare_cop_data(): assemble the aligned two-arm, by-serotype JAGS data.
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
#
# The outcome (Cases / Total_Cases) travels in the same rows; subsetting by
# `predictor_study` selects one immunogenicity source while keeping whatever
# case counts that study's rows carry. Returns a list with the JAGS data plus
# the aligned per-arm data frames and derived quantities used for plotting.
# ---------------------------------------------------------------------
prepare_cop_data <- function(data_file, predictor_study,
                             predictor_error = c("se", "sd"), quiet = FALSE) {
  gp <- prepare_gmc_predictor(data_file, predictor_study, predictor_error, quiet)
  u <- gp$u; i <- gp$i

  jags_data <- c(
    list(S = gp$S,
         cases_u = u$Cases,    total_u = u$Total_Cases,
         cases_i = i$Cases,    total_i = i$Total_Cases),
    gp$jags_data
  )

  c(list(jags_data = jags_data), gp[setdiff(names(gp), "jags_data")])
}

# ---------------------------------------------------------------------
# prepare_gmc_predictor(): shared helper -- reads the merged CSV, subsets to
# `predictor_study`, aligns the unimmunized/immunized arms by serotype, and
# builds the log-GMC error-in-variables predictor (lgmc_u/i, prec_u/i). Used
# by both prepare_cop_data() (case-count outcome) and prepare_cop_data_logor()
# (directly-supplied log-OR outcome) so the predictor logic never drifts
# between the two. See prepare_cop_data() above for argument meanings.
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
# prepare_cop_data_logor(): assemble JAGS data for the "logor" parameterization
# (cop_eiv_model_logor.jags), where the outcome is a directly-supplied log-OR
# per serotype instead of case counts.
#
#   data_file       : merged CSV (serotype x arm x study rows) supplying the
#                     immunogenicity GMCs -- same format/columns as
#                     prepare_cop_data() (Cases/Total_Cases columns, if
#                     present, are simply ignored).
#   predictor_study : value of the `Study` column selecting the GMC predictor.
#   or_file         : CSV with one row per serotype giving the outcome:
#                       Serotype    matches the Serotype column in data_file
#                       logOR       observed log odds ratio
#                       logOR_Lower, logOR_Upper   95% CI, ALREADY ON THE LOG
#                                   ODDS-RATIO SCALE (not exponentiated)
#   predictor_error : see prepare_cop_data() / prepare_gmc_predictor().
#   quiet           : suppress the diagnostic messages (asymmetry / SE floor).
#
# The log-OR CI is checked for symmetry the same way the GMC CIs are (0 =
# symmetric; the normal-on-log-OR EIV likelihood assumes it), and its SE is
# computed directly (no additional log transform, since logOR_Lower/Upper are
# already on the log scale): SE = (logOR_Upper - logOR_Lower) / (2 * 1.96).
# ---------------------------------------------------------------------
prepare_cop_data_logor <- function(data_file, predictor_study, or_file,
                                   predictor_error = c("se", "sd"), quiet = FALSE) {
  gp <- prepare_gmc_predictor(data_file, predictor_study, predictor_error, quiet)
  serotypes <- gp$serotypes

  or_d <- read.csv(or_file, stringsAsFactors = FALSE)
  or_d <- or_d[match(serotypes, or_d$Serotype), ]
  if (any(is.na(or_d$Serotype))) {
    stop("or_file '", or_file, "' is missing rows for serotype(s): ",
         paste(serotypes[is.na(or_d$Serotype)], collapse = ", "))
  }

  # Flag log-OR CIs that are not symmetric (already on the log scale).
  asym_or <- ci_asym(or_d$logOR_Lower, or_d$logOR, or_d$logOR_Upper)
  flagged_or <- data.frame(
    Serotype = serotypes,
    asym_pct = round(100 * asym_or, 1)
  )
  flagged_or <- flagged_or[!is.na(flagged_or$asym_pct) & abs(flagged_or$asym_pct) > 100 * ASYM_TOL, ]
  if (!quiet && nrow(flagged_or)) {
    message("Note: ", nrow(flagged_or), " log-OR CI(s) not symmetric ",
            "(|asymmetry| > ", 100 * ASYM_TOL, "%); normal-on-log-OR EIV ",
            "likelihood is misspecified for these:")
    print(flagged_or[order(-abs(flagged_or$asym_pct)), ], row.names = FALSE)
  }

  se_or <- ci_se(or_d$logOR_Lower, or_d$logOR_Upper)
  n_floored <- sum(se_or <= SE_FLOOR)
  if (!quiet && n_floored > 0) {
    message("Note: floored ", n_floored, " zero/near-zero-width log-OR CI(s) at SE = ", SE_FLOOR)
  }
  se_or <- pmax(se_or, SE_FLOOR)

  jags_data <- c(
    list(S = gp$S,
         logor = or_d$logOR, prec_or = 1 / se_or^2),
    gp$jags_data
  )

  c(list(jags_data = jags_data,
         logOR = or_d$logOR, se_or = se_or, flagged_or = flagged_or),
    gp[setdiff(names(gp), "jags_data")])
}

# ---------------------------------------------------------------------
# fit_cop(): compile, burn in, sample, summarise, and save outputs.
#
#   prep     : output of prepare_cop_data()
#   out_dir  : results folder for this analysis (created if missing)
#   n.*      : MCMC controls (defaults match the original analysis)
#   seed     : master seed (per-chain RNG seeds are derived from it)
#
# Writes into out_dir: posterior_summary.csv, mcmc.rds, diagnostics.pdf.
# Returns (invisibly) a list with the samples and key convergence stats.
# ---------------------------------------------------------------------
fit_cop <- function(prep, out_dir,
                    n.adapt = 2000, n.burnin = 5000, n.iter = 20000, thin = 5,
                    n.chains = 3, seed = 42, quiet = FALSE,
                    parameterization = "centered") {
  if (!parameterization %in% names(PARAMETERIZATIONS)) {
    stop("Unknown parameterization '", parameterization, "'. Known: ",
         paste(names(PARAMETERIZATIONS), collapse = ", "))
  }
  spec <- PARAMETERIZATIONS[[parameterization]]

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  jd <- prep$jags_data

  if (!quiet) {
    cat("Parameterization:", parameterization, "(", spec$model_file, ")\n")
    cat("Serotypes (order):", paste(prep$serotypes, collapse = ", "), "\n")
    cat("Observed centered log-GMC (immunized - unimmunized):\n")
    print(setNames(round(prep$x_obs, 3), prep$serotypes))
  }

  # Data-driven inits + per-chain RNG for reproducibility. Each spec supplies
  # its own canonical-named inits builder (falls back to the shared default
  # for parameterizations with the full baseline hierarchy).
  build_inits <- if (!is.null(spec$make_inits)) spec$make_inits else default_make_inits
  base_inits  <- build_inits(prep)
  make_inits <- function(chain_seed) {
    # built with CANONICAL names, then renamed to the model's own names.
    canon_to_model_inits(c(base_inits,
      list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = chain_seed)
    ), spec$map)
  }
  inits <- lapply(seq_len(n.chains), function(k) make_inits(11 * k))

  params <- if (!is.null(spec$params)) spec$params else DEFAULT_PARAMS

  model <- jags.model(spec$model_file, data = jd,
                      inits = inits,
                      n.chains = n.chains, n.adapt = n.adapt,
                      quiet = quiet)
  update(model, n.burnin)
  samp <- coda.samples(model, variable.names = canon_to_model_params(params, spec$map),
                       n.iter = n.iter, thin = thin)
  # Back to canonical names so all downstream code is parameterization-agnostic.
  samp <- model_to_canon_samples(samp, spec$map)

  summ  <- summary(samp)
  stats <- cbind(summ$statistics[, c("Mean", "SD")], summ$quantiles)

  M <- as.matrix(samp)
  mu_b1_draws <- M[, "mu_b1"]
  rr_draws    <- M[, "rr_global"]

  gd  <- gelman.diag(samp, multivariate = FALSE)
  max_psrf <- max(gd$psrf[, "Point est."], na.rm = TRUE)
  min_ess  <- min(effectiveSize(samp))

  if (!quiet) {
    cat("\n--- Global correlate of protection ---\n")
    cat(sprintf("mu_b1 (log RR per unit centered log-GMC): %.3f  (95%% CrI %.3f, %.3f)\n",
                median(mu_b1_draws), quantile(mu_b1_draws, .025), quantile(mu_b1_draws, .975)))
    cat(sprintf("rr_global (rate ratio per unit centered log-GMC): %.3f  (95%% CrI %.3f, %.3f)\n",
                median(rr_draws), quantile(rr_draws, .025), quantile(rr_draws, .975)))
    cat(sprintf("P(mu_b1 < 0) = %.3f\n", mean(mu_b1_draws < 0)))
    cat("\n--- Convergence ---\n")
    cat("Max PSRF (all monitored):", round(max_psrf, 3), "\n")
    cat("Min effective sample size:", round(min_ess), "\n")
  }

  write.csv(round(stats, 4), file.path(out_dir, "posterior_summary.csv"))
  saveRDS(samp, file.path(out_dir, "mcmc.rds"))

  diag_params <- if (!is.null(spec$diag_params)) spec$diag_params else DEFAULT_DIAG_PARAMS
  pdf(file.path(out_dir, "diagnostics.pdf"), width = 9, height = 6)
  plot(samp[, diag_params])
  dev.off()

  if (!quiet) cat("\nSaved fit outputs to ", out_dir, "\n", sep = "")

  invisible(list(samp = samp, stats = stats,
                 max_psrf = max_psrf, min_ess = min_ess))
}

# ---------------------------------------------------------------------
# slope_summary(): tidy per-analysis slope table (global + per serotype),
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
# plot_cop_scatter(): the correlate-of-protection relationship in its natural
# "ratio" coordinates - one point per serotype:
#
#   x = observed log GMC-ratio   log(GMC1 / GMC0)   (immunized vs unimmunized)
#   y = observed log rate-ratio  log(p1  / p0 )     (case proportions by arm)
#
# Both axes carry a 95% uncertainty interval:
#   x: GMC measurement error, SE = sqrt(se_i^2 + se_u^2) on the log scale
#   y: sampling error in the two case counts, SE = sqrt(1/cases_i + 1/cases_u)
#      (delta method for a log ratio of Poisson counts; offsets treated as fixed)
# Serotypes with a zero-count arm use a 0.5 continuity correction for the plotted
# y (open symbol) so the point/interval remain finite.
#
# Overlaid: the fitted GLOBAL relationship log RR = mu_b1 * log GMR - a line
# THROUGH THE ORIGIN (the serotype intercept a[s] cancels in the ratio), with a
# 95% credible band from the posterior of mu_b1. Saves cop_scatter_gmr_rr.png.
# Consumes only canonical names, so it is parameterization-agnostic -- EXCEPT
# for the y-axis source: if `prep` carries a `logOR` field (i.e. it came from
# prepare_cop_data_logor()), the observed log-OR and its own CI-derived SE are
# plotted directly instead of being derived from case counts.
# ---------------------------------------------------------------------
plot_cop_scatter <- function(prep, samp, out_dir, title_suffix = "") {
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  serotypes <- prep$serotypes
  S <- prep$S
  u <- prep$u; i <- prep$i
  z <- qnorm(0.975)

  # x-axis: observed log GMC-ratio with combined measurement error.
  logGMR <- prep$x_obs
  se_gmr <- sqrt(prep$se_i^2 + prep$se_u^2)

  if (!is.null(prep$logOR)) {
    # y-axis: directly-supplied observed log-OR, with its own reported SE.
    logRR <- prep$logOR
    se_rr <- prep$se_or
    zero  <- rep(FALSE, S)
  } else {
    # y-axis: observed log rate-ratio with count sampling error (Haldane 0.5
    # correction where an arm has zero cases).
    ci <- i$Cases; cu <- u$Cases
    zero <- ci == 0 | cu == 0
    cic <- ifelse(zero, ci + 0.5, ci)
    cuc <- ifelse(zero, cu + 0.5, cu)
    logRR <- log((cic / i$Total_Cases) / (cuc / u$Total_Cases))
    se_rr <- sqrt(1 / cic + 1 / cuc)
  }

  M <- as.matrix(samp)
  b1_hat <- vapply(seq_len(S), function(s) mean(M[, sprintf("b1[%d]", s)]), numeric(1))
  mu_b1  <- mean(M[, "mu_b1"])

  # Global fitted line through the origin: log RR = mu_b1 * log GMR.
  xr <- range(c(0, logGMR))
  xg <- seq(xr[1] - 0.15, xr[2] + 0.15, length.out = 200)
  line_mat <- outer(M[, "mu_b1"], xg)
  fit_glob <- data.frame(
    x   = xg,
    med = apply(line_mat, 2, median),
    lo  = apply(line_mat, 2, quantile, 0.025),
    hi  = apply(line_mat, 2, quantile, 0.975)
  )

  pts <- data.frame(
    Serotype = serotypes, x = logGMR, y = logRR,
    xlo = logGMR - z * se_gmr, xhi = logGMR + z * se_gmr,
    ylo = logRR  - z * se_rr,  yhi = logRR  + z * se_rr,
    zero = zero, stringsAsFactors = FALSE
  )
  slope_labels <- setNames(sprintf("%s  (post. b1 = %+.2f)", serotypes, b1_hat), serotypes)

  main_title <- "Correlate of protection: log rate-ratio vs log GMC-ratio"
  if (nzchar(title_suffix)) main_title <- paste0(main_title, "\n", title_suffix)
  subt <- sprintf(paste0("Points = observed serotypes; bars = 95%% CI on both axes; ",
                         "dotted lines = each serotype's own observed slope (logRR/logGMR); ",
                         "black line = global fit through origin, b1 = %+.2f (95%% CrI)"), mu_b1)
  cap <- if (any(zero))
    "Open symbols: an arm has zero cases; log RR uses a 0.5 continuity correction." else NULL

  if (has_ggplot) {
    library(ggplot2)
    p <- ggplot() +
      geom_ribbon(data = fit_glob, aes(x = x, ymin = lo, ymax = hi),
                  fill = "grey60", alpha = 0.25) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      geom_line(data = fit_glob, aes(x = x, y = med), colour = "black", linewidth = 1.2) +
      geom_segment(data = pts, aes(x = 0, y = 0, xend = x, yend = y, colour = Serotype),
                   linetype = "dotted", alpha = 0.5, linewidth = 0.5) +
      geom_errorbarh(data = pts, aes(y = y, xmin = xlo, xmax = xhi, colour = Serotype),
                     height = 0, alpha = 0.7) +
      geom_errorbar(data = pts, aes(x = x, ymin = ylo, ymax = yhi, colour = Serotype),
                    width = 0, alpha = 0.7) +
      geom_point(data = pts, aes(x = x, y = y, colour = Serotype, shape = zero), size = 2.8) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
      scale_colour_discrete(labels = slope_labels) +
      labs(
        x = "Observed log GMC-ratio   log(GMC1 / GMC0)",
        y = "Observed log rate-ratio   log(p1 / p0)",
        colour = "Serotype (slope b1)",
        title = main_title, subtitle = subt, caption = cap
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right")
    ggsave(file.path(out_dir, "cop_scatter_gmr_rr.png"), p, width = 9, height = 6, dpi = 150)
  } else {
    cols <- setNames(rainbow(S), serotypes)
    png(file.path(out_dir, "cop_scatter_gmr_rr.png"), width = 9, height = 6, units = "in", res = 150)
    xlim <- range(c(pts$xlo, pts$xhi, xg))
    ylim <- range(c(pts$ylo, pts$yhi, fit_glob$lo, fit_glob$hi))
    plot(NA, xlim = xlim, ylim = ylim,
         xlab = "Observed log GMC-ratio   log(GMC1 / GMC0)",
         ylab = "Observed log rate-ratio   log(p1 / p0)",
         main = main_title)
    polygon(c(fit_glob$x, rev(fit_glob$x)), c(fit_glob$lo, rev(fit_glob$hi)),
            col = adjustcolor("grey60", 0.25), border = NA)
    abline(h = 0, v = 0, lty = 3, col = "grey50")
    lines(fit_glob$x, fit_glob$med, col = "black", lwd = 2)
    segments(0, 0, pts$x, pts$y, col = adjustcolor(cols, alpha.f = 0.5), lty = 3)
    arrows(pts$xlo, pts$y, pts$xhi, pts$y, code = 3, angle = 90, length = 0.03, col = cols)
    arrows(pts$x, pts$ylo, pts$x, pts$yhi, code = 3, angle = 90, length = 0.03, col = cols)
    points(pts$x, pts$y, col = cols, pch = ifelse(pts$zero, 1, 16), cex = 1.3)
    legend("topleft", legend = slope_labels, col = cols, pch = 16, bty = "n",
           cex = 0.8, title = "Serotype (slope b1)")
    if (!is.null(cap)) mtext(cap, side = 1, line = 4, cex = 0.7)
    dev.off()
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------
# plot_cop_absolute(): the correlate-of-protection relationship in ABSOLUTE
# coordinates - one point per (serotype, arm):
#
#   x = observed absolute GMC     (natural units, log scale on the plot)
#   y = observed absolute risk    Cases / Total_Cases  (case proportion)
#
# This is only meaningful for a fit whose a[s]/mu_a/b1[s]/mu_b1 live on the
# ABSOLUTE log-GMC scale, i.e. the "RE" parameterization (see
# cop_eiv_model_RE.jags) - the centered/ratio/logor models have no serotype
# intercept defined at absolute log-GMC = 0, so plugging their a[s]/b1[s] into
# this plot would silently mispredict. `samp` must therefore come from an RE fit.
#
# Both axes carry a 95% uncertainty interval, same conventions as
# plot_cop_scatter(): x from the reported log-GMC CI, y from Poisson count
# sampling error (delta method, 0.5 continuity correction for zero counts).
#
# Overlaid: each serotype's own fitted line log(risk) = a[s] + b1[s]*log(GMC)
# across its own observed absolute-GMC range (dotted, colour = serotype), and
# the GLOBAL population-average curve log(risk) = mu_a + mu_b1*log(GMC) with a
# 95% credible band (black), spanning the full observed absolute-GMC range.
# Saves cop_scatter_absolute_risk_gmc.png.
# ---------------------------------------------------------------------
plot_cop_absolute <- function(prep, samp, out_dir, title_suffix = "") {
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  serotypes <- prep$serotypes
  S <- prep$S
  u <- prep$u; i <- prep$i
  z <- qnorm(0.975)

  M <- as.matrix(samp)
  need <- c("mu_a", "mu_b1", sprintf("a[%d]", seq_len(S)), sprintf("b1[%d]", seq_len(S)))
  if (!all(need %in% colnames(M))) {
    stop("plot_cop_absolute() requires a fit with canonical a[]/mu_a/b1[]/mu_b1 on the ",
         "ABSOLUTE log-GMC scale (fit with parameterization = 'RE') -- missing: ",
         paste(setdiff(need, colnames(M)), collapse = ", "))
  }

  # x-axis: observed absolute log-GMC, by arm, with reported measurement error.
  lgmc_u <- log(u$GMC); lgmc_i <- log(i$GMC)

  # y-axis: observed absolute risk (case proportion), by arm, with count
  # sampling error (Haldane 0.5 correction where an arm has zero cases).
  zero_u <- u$Cases == 0; zero_i <- i$Cases == 0
  cu <- ifelse(zero_u, u$Cases + 0.5, u$Cases)
  ci <- ifelse(zero_i, i$Cases + 0.5, i$Cases)
  risk_u <- cu / u$Total_Cases
  risk_i <- ci / i$Total_Cases
  se_risk_u <- sqrt(1 / cu)
  se_risk_i <- sqrt(1 / ci)

  a_hat  <- vapply(seq_len(S), function(s) mean(M[, sprintf("a[%d]", s)]),  numeric(1))
  b1_hat <- vapply(seq_len(S), function(s) mean(M[, sprintf("b1[%d]", s)]), numeric(1))

  # Per-serotype fitted line across ITS OWN observed absolute log-GMC range.
  sero_lines <- do.call(rbind, lapply(seq_len(S), function(s) {
    xg <- seq(min(lgmc_u[s], lgmc_i[s]), max(lgmc_u[s], lgmc_i[s]), length.out = 50)
    data.frame(Serotype = serotypes[s], gmc = exp(xg),
               risk = exp(a_hat[s] + b1_hat[s] * xg), stringsAsFactors = FALSE)
  }))

  # Global fitted line + 95% credible band: log(risk) = mu_a + mu_b1*log(GMC),
  # spanning the full observed absolute log-GMC range across all serotypes/arms.
  xr  <- range(c(lgmc_u, lgmc_i))
  pad <- 0.1 * diff(xr)
  xg  <- seq(xr[1] - pad, xr[2] + pad, length.out = 200)
  line_mat <- outer(M[, "mu_a"], rep(1, length(xg))) + outer(M[, "mu_b1"], xg)
  fit_glob <- data.frame(
    gmc = exp(xg),
    med = exp(apply(line_mat, 2, median)),
    lo  = exp(apply(line_mat, 2, quantile, 0.025)),
    hi  = exp(apply(line_mat, 2, quantile, 0.975))
  )

  pts <- data.frame(
    Serotype = rep(serotypes, 2),
    Arm      = rep(c("Unimmunized", "Immunized"), each = S),
    gmc      = exp(c(lgmc_u, lgmc_i)),
    gmc_lo   = exp(c(lgmc_u - z * prep$se_u, lgmc_i - z * prep$se_i)),
    gmc_hi   = exp(c(lgmc_u + z * prep$se_u, lgmc_i + z * prep$se_i)),
    risk     = c(risk_u, risk_i),
    risk_lo  = c(risk_u * exp(-z * se_risk_u), risk_i * exp(-z * se_risk_i)),
    risk_hi  = c(risk_u * exp( z * se_risk_u), risk_i * exp( z * se_risk_i)),
    zero     = c(zero_u, zero_i),
    stringsAsFactors = FALSE
  )
  slope_labels <- setNames(sprintf("%s  (post. b1 = %+.2f)", serotypes, b1_hat), serotypes)

  main_title <- "Correlate of protection: absolute risk vs absolute GMC"
  if (nzchar(title_suffix)) main_title <- paste0(main_title, "\n", title_suffix)
  subt <- paste0("Points = observed serotype x arm (circle = unimmunized, triangle = immunized); ",
                "dotted lines = each serotype's own fitted risk-GMC relationship; ",
                "black line = population-average fit (mu_a, mu_b1) with 95% credible band")
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
      labs(
        x = "Absolute GMC", y = "Absolute risk (case proportion)",
        colour = "Serotype (slope b1)", shape = "Arm",
        title = main_title, subtitle = subt, caption = cap
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right")
    ggsave(file.path(out_dir, "cop_scatter_absolute_risk_gmc.png"), p, width = 9, height = 6, dpi = 150)
  } else {
    cols <- setNames(rainbow(S), serotypes)
    png(file.path(out_dir, "cop_scatter_absolute_risk_gmc.png"), width = 9, height = 6, units = "in", res = 150)
    xlim <- range(c(pts$gmc_lo, pts$gmc_hi, fit_glob$gmc))
    ylim <- range(c(pts$risk_lo, pts$risk_hi, fit_glob$lo, fit_glob$hi))
    plot(NA, xlim = xlim, ylim = ylim, log = "xy",
         xlab = "Absolute GMC", ylab = "Absolute risk (case proportion)",
         main = main_title)
    abline(v = 0.35, lty = 2, col = "grey40")
    polygon(c(fit_glob$gmc, rev(fit_glob$gmc)), c(fit_glob$lo, rev(fit_glob$hi)),
            col = adjustcolor("grey60", 0.25), border = NA)
    lines(fit_glob$gmc, fit_glob$med, col = "black", lwd = 2)
    for (s in seq_len(S)) {
      sl <- sero_lines[sero_lines$Serotype == serotypes[s], ]
      lines(sl$gmc, sl$risk, col = adjustcolor(cols[s], alpha.f = 0.6), lty = 3)
    }
    pch_arm <- ifelse(pts$zero, ifelse(pts$Arm == "Unimmunized", 1, 2),
                      ifelse(pts$Arm == "Unimmunized", 16, 17))
    arrows(pts$gmc_lo, pts$risk, pts$gmc_hi, pts$risk, code = 3, angle = 90, length = 0.03,
           col = cols[pts$Serotype])
    arrows(pts$gmc, pts$risk_lo, pts$gmc, pts$risk_hi, code = 3, angle = 90, length = 0.03,
           col = cols[pts$Serotype])
    points(pts$gmc, pts$risk, col = cols[pts$Serotype], pch = pch_arm, cex = 1.3)
    legend("topleft", legend = slope_labels, col = cols, pch = 16, bty = "n",
           cex = 0.8, title = "Serotype (slope b1)")
    legend("bottomright", legend = c("Unimmunized", "Immunized"), pch = c(16, 17),
           bty = "n", cex = 0.8, title = "Arm")
    if (!is.null(cap)) mtext(cap, side = 1, line = 4, cex = 0.7)
    dev.off()
  }
  invisible(NULL)
}

# =====================================================================
# MULTI-STUDY GENERALIZATION (cop_eiv_model_multistudy.jags)
#
# Generalizes the single-outcome/single-predictor design above to MANY
# outcome studies (ragged over serotypes, sharing b1[s] but with a baseline
# a[k,s] that varies FREELY by both study and serotype, shrunk toward one
# shared hypermean/hypervariance mu_a/tau_a) and MANY immunogenicity sources
# (pooled onto a single shared latent true log-GMC per serotype/arm, under
# the assumption that the population distribution of log-titers for a given
# serotype and vaccine is the same across studies). See the JAGS file header
# for the full statistical rationale.
#
#   read_outcome_source()          : one outcome dataset -> per-serotype
#                                     Cases/Total_Cases (deduped across the
#                                     immunogenicity Study subsets that share
#                                     the same merged CSV)
#   prepare_cop_data_multistudy()  : assemble the ragged, multi-study JAGS
#                                     data list from a list of outcome
#                                     sources and a list of immunogenicity
#                                     sources
#   fit_cop_multistudy()           : compile + sample + save outputs
#   study_summary()                : tidy per-outcome-study baseline (mean of
#                                     a[k,s] across that study's own
#                                     serotypes, vs. the global hypermean
#                                     mu_a) table
#
# slope_summary() above is generic over prep$S/prep$serotypes and canonical
# mcmc column names, so it works UNMODIFIED against a multi-study fit's
# global a[k,s]/b1[s]. plot_cop_scatter()/plot_cop_absolute() do NOT apply
# here unmodified -- they assume a single u/i pair per serotype, which no
# longer exists once outcomes are ragged across studies; use
# prep$outcome_sources / prep$immuno_sources for any per-source plotting
# instead.
# =====================================================================

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
#                     entries pool onto the SAME shared x_u[s]/x_i[s].
#   predictor_error : "se" or "sd", applied identically to every
#                     immunogenicity source (see prepare_gmc_predictor()).
#   quiet           : suppress per-source diagnostic messages.
#
# Returns a list with jags_data (ragged vectors/index arrays), the canonical
# serotype universe (serotypes, S), and the per-source bookkeeping
# (outcome_sources, immuno_sources) needed for study_summary() and any
# future per-source plotting.
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
# outputs for the multi-study model (cop_eiv_model_multistudy.jags). Mirrors
# fit_cop() above, but the data/inits/params are specific to this model's
# ragged, multi-study structure rather than the single-outcome/single-
# predictor PARAMETERIZATIONS registry.
#
#   prep     : output of prepare_cop_data_multistudy()
#   out_dir  : results folder for this analysis (created if missing)
#   n.*      : MCMC controls (defaults match fit_cop())
#   seed     : master seed (per-chain RNG seeds are derived from it)
#
# Writes into out_dir: posterior_summary.csv, mcmc.rds, diagnostics.pdf.
# Returns (invisibly) a list with the samples and key convergence stats.
# ---------------------------------------------------------------------
MULTISTUDY_MODEL_FILE <- file.path("JAGS", "cop_eiv_model_multistudy.jags")

MULTISTUDY_PARAMS <- c("mu_b1", "b1", "rr", "rr_global",
                       "a", "mu_a", "sigma_a", "sigma_b1",
                       "mu_x_u", "sigma_x_u", "x_u",
                       "mu_x_i", "sigma_x_i", "x_i")

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
  }

  write.csv(round(stats, 4), file.path(out_dir, "posterior_summary.csv"))
  saveRDS(samp, file.path(out_dir, "mcmc.rds"))

  diag_params <- c("mu_b1", "sigma_b1", "mu_a", "sigma_a")
  pdf(file.path(out_dir, "diagnostics.pdf"), width = 9, height = 6)
  plot(samp[, diag_params])
  dev.off()

  if (!quiet) cat("\nSaved fit outputs to ", out_dir, "\n", sep = "")

  invisible(list(samp = samp, stats = stats,
                 max_psrf = max_psrf, min_ess = min_ess))
}

# ---------------------------------------------------------------------
# study_summary(): tidy per-outcome-study baseline table. There is no longer
# a single per-study fixed effect (a[k,s] now varies by serotype within a
# study too), so this summarises each study's OWN average baseline --
# mean_s(a[k,s]) over just the serotypes that study reports -- and how that
# compares to the global hypermean mu_a (the pooled baseline across every
# study/serotype). rate_ratio_vs_global > 1 means this study's serotypes run
# hotter (higher case ascertainment / incidence) than the pooled average;
# < 1 means cooler. Written to out_dir/study_summary.csv (if out_dir given)
# and returned as a data frame.
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
# plot_cop_absolute_multistudy(): analogous to plot_cop_absolute() above, but
# for a multi-study fit (prepare_cop_data_multistudy() / fit_cop_multistudy()):
#
#   x = POOLED absolute GMC per serotype/arm -- the model's posterior
#       estimate of the SHARED latent x_u[s]/x_i[s] (exponentiated median),
#       i.e. the combined estimate across every immunogenicity source. There
#       is no longer a single "observed" GMC once multiple sources are
#       pooled, so the x-axis (and its horizontal error bar, the posterior
#       95% credible interval of that latent) is a model estimate, not a raw
#       observation -- unlike plot_cop_absolute()'s x-axis.
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

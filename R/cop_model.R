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
# Two algebraically-equivalent JAGS models are supported. The rest of the
# pipeline (plots, slope tables, comparisons) speaks a single CANONICAL set of
# variable names; each parameterization declares how its own JAGS variables map
# onto those canonical names. fit_cop() monitors the model's OWN names and then
# renames the posterior samples back to canonical, so nothing downstream needs
# to know which parameterization produced a fit (mcmc.rds is always canonical).
#
#   canonical  meaning                                centered   ratio
#   ---------  -------------------------------------  --------   --------
#   a[s]       serotype baseline log case-proportion  a[s]       logp0[s]
#   mu_a       its hypermean                          mu_a       mu_p0
#   sigma_a    its hyper-SD                           sigma_a    sigma_p0
#   x[s]       latent log-GMC ratio (predictor)       x[s]       lrgmc[s]
#   mu_x       its hypermean                          mu_x       mu_lr
#   sigma_x    its hyper-SD                           sigma_x    sigma_lr
#   ref[s]     latent unimmunized true log-GMC        ref[s]     g0[s]
#   b1[s], mu_b1, sigma_b1, rr[s], rr_global  --- identical in both ---
PARAMETERIZATIONS <- list(
  centered = list(
    model_file = file.path("JAGS", "cop_eiv_model.jags"),
    map = c(a = "a", mu_a = "mu_a", sigma_a = "sigma_a",
            x = "x", mu_x = "mu_x", sigma_x = "sigma_x", ref = "ref")
  ),
  ratio = list(
    model_file = file.path("JAGS", "cop_eiv_model_ratio.jags"),
    map = c(a = "logp0", mu_a = "mu_p0", sigma_a = "sigma_p0",
            x = "lrgmc", mu_x = "mu_lr", sigma_x = "sigma_lr", ref = "g0")
  )
)

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
# SE(log GMC) from a 95% CI, assuming (log Upper - log Lower) spans 2*1.96 SE.
log_se <- function(lower, upper) (log(upper) - log(lower)) / (2 * qnorm(0.975))

# Relative asymmetry of the CI on the log scale (0 = symmetric). The normal-on-
# log EIV prior assumes symmetry; large values flag misspecification (usually
# from aggressive rounding of small GMCs).
log_asym <- function(lower, gmc, upper) {
  lo <- log(gmc) - log(lower)
  hi <- log(upper) - log(gmc)
  2 * (hi - lo) / (hi + lo)
}

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

  jags_data <- list(
    S       = S,
    cases_u = u$Cases,    total_u = u$Total_Cases,
    cases_i = i$Cases,    total_i = i$Total_Cases,
    lgmc_u  = log(u$GMC), prec_u  = 1 / sd_u^2,
    lgmc_i  = log(i$GMC), prec_i  = 1 / sd_i^2
  )

  list(
    jags_data = jags_data,
    serotypes = serotypes, S = S,
    u = u, i = i,
    se_u = se_u, se_i = se_i,
    predictor_error = predictor_error,
    x_obs = log(i$GMC) - log(u$GMC),        # observed centered log-GMC (immunized)
    flagged = flagged
  )
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
  S  <- prep$S
  u  <- prep$u
  i  <- prep$i

  if (!quiet) {
    cat("Parameterization:", parameterization, "(", spec$model_file, ")\n")
    cat("Serotypes (order):", paste(prep$serotypes, collapse = ", "), "\n")
    cat("Observed centered log-GMC (immunized - unimmunized):\n")
    print(setNames(round(prep$x_obs, 3), prep$serotypes))
  }

  # Data-driven inits + per-chain RNG for reproducibility.
  a_start   <- log(u$Cases / u$Total_Cases)
  x_start   <- prep$x_obs
  ref_start <- log(u$GMC)
  make_inits <- function(chain_seed) {
    # built with CANONICAL names, then renamed to the model's own names.
    canon_to_model_inits(list(
      ref = ref_start, x = x_start,
      a = a_start, b1 = rep(-0.5, S),
      mu_a = mean(a_start), sigma_a = 0.3,
      mu_b1 = -0.5, sigma_b1 = 0.3,
      mu_x = mean(x_start), sigma_x = 0.5,
      .RNG.name = "base::Mersenne-Twister", .RNG.seed = chain_seed
    ), spec$map)
  }
  inits <- lapply(seq_len(n.chains), function(k) make_inits(11 * k))

  params <- c("mu_b1", "b1", "rr", "rr_global",
              "a", "mu_a", "sigma_a", "sigma_b1",
              "mu_x", "sigma_x", "x", "ref")

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

  pdf(file.path(out_dir, "diagnostics.pdf"), width = 9, height = 6)
  plot(samp[, c("mu_b1", "sigma_b1", "mu_a", "sigma_a")])
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
# plot_cop(): fitted case-proportion curves vs centered log-GMC.
#   y = exp(a[s] + b1[s]*x) = expected proportion of cases (log scale)
#   x = centered log-GMC; x = 0 is the unimmunized reference
# Coloured lines = serotype-specific fits, black line = global (pooled) fit
# with 95% credible band; points = observed proportions with 95% GMC
# measurement-error bars. Saves cop_fit_proportion.pdf (+ .png if ggplot2).
# ---------------------------------------------------------------------
plot_cop <- function(prep, samp, out_dir, title_suffix = "") {
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  serotypes <- prep$serotypes
  S <- prep$S
  u <- prep$u; i <- prep$i
  x_obs  <- prep$x_obs
  prop_u <- u$Cases / u$Total_Cases
  prop_i <- i$Cases / i$Total_Cases
  # unimmunized point at x=0 with its own SE; immunized combines both arms.
  se_x <- c(prep$se_u, sqrt(prep$se_i^2 + prep$se_u^2))

  M <- as.matrix(samp)
  post_mean <- function(par) mean(M[, par])
  a_hat  <- vapply(seq_len(S), function(s) post_mean(sprintf("a[%d]",  s)), numeric(1))
  b1_hat <- vapply(seq_len(S), function(s) post_mean(sprintf("b1[%d]", s)), numeric(1))
  mu_b1  <- post_mean("mu_b1")

  slope_labels <- setNames(sprintf("%s  (b1 = %+.2f)", serotypes, b1_hat), serotypes)

  xr <- range(c(0, x_obs))
  xg <- seq(xr[1] - 0.2, xr[2] + 0.2, length.out = 200)

  fit_s <- data.frame(
    Serotype = rep(serotypes, each = length(xg)),
    x        = rep(xg, times = S),
    prop     = as.vector(sapply(seq_len(S), function(s) exp(a_hat[s] + b1_hat[s] * xg)))
  )
  xg_mat    <- outer(M[, "mu_a"], rep(1, length(xg))) + outer(M[, "mu_b1"], xg)
  glob_draw <- exp(xg_mat)
  fit_glob <- data.frame(
    x   = xg,
    med = apply(glob_draw, 2, median),
    lo  = apply(glob_draw, 2, quantile, 0.025),
    hi  = apply(glob_draw, 2, quantile, 0.975)
  )
  obs <- data.frame(
    Serotype = rep(serotypes, 2),
    Arm      = rep(c("Unimmunized", "Immunized"), each = S),
    x        = c(rep(0, S), x_obs),
    prop     = c(prop_u, prop_i),
    xlo      = c(rep(0, S), x_obs) - qnorm(0.975) * se_x,
    xhi      = c(rep(0, S), x_obs) + qnorm(0.975) * se_x
  )

  main_title <- "Correlate of protection: fitted case proportion vs centered log-GMC"
  if (nzchar(title_suffix)) main_title <- paste0(main_title, "\n", title_suffix)

  if (has_ggplot) {
    library(ggplot2)
    p <- ggplot() +
      geom_ribbon(data = fit_glob, aes(x = x, ymin = lo, ymax = hi),
                  fill = "grey60", alpha = 0.25) +
      geom_line(data = fit_s, aes(x = x, y = prop, colour = Serotype), size = 0.7) +
      geom_line(data = fit_glob, aes(x = x, y = med), colour = "black", size = 1.3) +
      geom_errorbarh(data = obs, aes(y = prop, xmin = xlo, xmax = xhi, colour = Serotype),
                     height = 0, alpha = 0.7) +
      geom_point(data = obs, aes(x = x, y = prop, colour = Serotype, shape = Arm), size = 2.6) +
      scale_shape_manual(values = c(Unimmunized = 1, Immunized = 16)) +
      scale_colour_discrete(labels = slope_labels) +
      scale_y_log10() +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      labs(
        x = "Centered log-GMC  (immunized - unimmunized reference)",
        y = "Expected proportion of cases  ( lambda / Total_Cases, log scale )",
        colour = "Serotype (slope b1)",
        title = main_title,
        subtitle = sprintf(paste0("Coloured lines = serotype-specific fits; ",
                                  "black line = global (pooled) fit, b1 = %+.2f, 95%% CrI; ",
                                  "bars = 95%% GMC measurement error"), mu_b1)
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right")
    ggsave(file.path(out_dir, "cop_fit_proportion.pdf"), p, width = 9, height = 6)
    ggsave(file.path(out_dir, "cop_fit_proportion.png"), p, width = 9, height = 6, dpi = 150)
  } else {
    cols <- setNames(rainbow(S), serotypes)
    pdf(file.path(out_dir, "cop_fit_proportion.pdf"), width = 9, height = 6)
    ylim <- range(c(fit_s$prop, obs$prop, fit_glob$hi))
    ylim[1] <- max(ylim[1], min(obs$prop[obs$prop > 0]) / 2)
    plot(NA, xlim = range(xg), ylim = ylim, log = "y",
         xlab = "Centered log-GMC (immunized - unimmunized reference)",
         ylab = "Expected proportion of cases (lambda / Total_Cases, log scale)",
         main = main_title)
    polygon(c(fit_glob$x, rev(fit_glob$x)), c(fit_glob$lo, rev(fit_glob$hi)),
            col = adjustcolor("grey60", 0.25), border = NA)
    for (s in seq_len(S)) lines(xg, exp(a_hat[s] + b1_hat[s] * xg), col = cols[s], lwd = 1.5)
    lines(fit_glob$x, fit_glob$med, col = "black", lwd = 3)
    abline(v = 0, lty = 3, col = "grey50")
    arrows(obs$xlo, obs$prop, obs$xhi, obs$prop, code = 3, angle = 90,
           length = 0.03, col = rep(cols, 2))
    points(0 * seq_len(S), prop_u, col = cols, pch = 1, cex = 1.3)
    points(x_obs, prop_i, col = cols, pch = 16, cex = 1.3)
    legend("topright", legend = slope_labels, col = cols, lwd = 1.5,
           title = "Serotype (slope b1)", bty = "n", cex = 0.8)
    legend("top", legend = c("Global (pooled)", "Unimmunized (x=0)", "Immunized"),
           col = c("black", "grey30", "grey30"), lwd = c(3, NA, NA),
           pch = c(NA, 1, 16), bty = "n", cex = 0.8)
    dev.off()
  }
  invisible(NULL)
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
# 95% credible band from the posterior of mu_b1. Saves cop_scatter_gmr_rr.pdf
# (+ .png if ggplot2). Consumes only canonical names, so it is
# parameterization-agnostic.
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

  # y-axis: observed log rate-ratio with count sampling error (Haldane 0.5
  # correction where an arm has zero cases).
  ci <- i$Cases; cu <- u$Cases
  zero <- ci == 0 | cu == 0
  cic <- ifelse(zero, ci + 0.5, ci)
  cuc <- ifelse(zero, cu + 0.5, cu)
  logRR <- log((cic / i$Total_Cases) / (cuc / u$Total_Cases))
  se_rr <- sqrt(1 / cic + 1 / cuc)

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
    ggsave(file.path(out_dir, "cop_scatter_gmr_rr.pdf"), p, width = 9, height = 6)
    ggsave(file.path(out_dir, "cop_scatter_gmr_rr.png"), p, width = 9, height = 6, dpi = 150)
  } else {
    cols <- setNames(rainbow(S), serotypes)
    pdf(file.path(out_dir, "cop_scatter_gmr_rr.pdf"), width = 9, height = 6)
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

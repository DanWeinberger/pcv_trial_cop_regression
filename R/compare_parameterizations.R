# =====================================================================
# Compare the COP slope across MODEL PARAMETERIZATIONS (not across analyses
# -- see R/compare_slopes.R for that). For one analysis <id> this fits
# "centered", "ratio", and "RE" (see R/cop_model.R PARAMETERIZATIONS) on the
# SAME data and compares their mu_b1 (global) and b1[s] (per-serotype)
# posteriors. "ratio" is algebraically equivalent to "centered" and is
# included mainly as a sanity check that the two agree; "RE" is a genuinely
# different model (cop_eiv_model_RE.jags: absolute-scale random effects on
# each arm's GMC, rather than the within-serotype paired ratio) whose b1[s]
# is nonetheless directly comparable (see that file's header for why). "logor"
# is excluded here: it needs a separately-supplied log-OR outcome file that
# isn't wired into the R/config.R analysis registry.
#
# Usage (from the project root):
#   Rscript R/compare_parameterizations.R                  # analysis "nckp"
#   Rscript R/compare_parameterizations.R navajo south_africa
#
# For each analysis <id> this fits every entry of PARAM_IDS into
# results/<id>_<parameterization>/ (see analysis_out_dir() in R/config.R;
# "centered" stays unsuffixed at results/<id>/) and writes
# results/comparisons/param_<id>/:
#   GLOBAL slope (mu_b1):
#     slopes_by_parameterization.csv
#     slope_pairwise_diffs.csv        (reference = "centered")
#     slope_forest.pdf (+ .png), slope_overlay.pdf (+ .png)
#   SEROTYPE-SPECIFIC slope (b1[s]):
#     serotype_slopes_by_parameterization.csv
#     serotype_slope_pairwise_diffs.csv
#     serotype_slope_forest.pdf (+ .png)
#   ABSOLUTE-SCALE figure (from the RE fit only -- see plot_cop_absolute()):
#     cop_scatter_absolute_risk_gmc.png
# =====================================================================

source(file.path("R", "config.R"))
source(file.path("R", "cop_model.R"))
suppressPackageStartupMessages(library(coda))
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

PARAM_IDS <- c("centered", "ratio", "RE")
REFERENCE <- "centered"

args <- commandArgs(trailingOnly = TRUE)
ids  <- if (length(args) == 0) "nckp" else args

q <- c(0.025, 0.5, 0.975)

# Numeric-aware serotype ordering: 4, 6A, 6B, 9V, 14, 18C, 19A, 19F, 23F.
sero_order <- function(x) {
  num <- suppressWarnings(as.numeric(gsub("[^0-9].*$", "", x)))
  x[order(num, x)]
}

# Summarise a vector of slope draws (log-RR) to median/CrI + RR + P(<0).
summarise_draws <- function(d) {
  qs <- quantile(d, q)
  data.frame(
    b1_mean = mean(d), b1_median = qs[["50%"]],
    b1_lo = qs[["2.5%"]], b1_hi = qs[["97.5%"]],
    rr_median = exp(qs[["50%"]]), rr_lo = exp(qs[["2.5%"]]), rr_hi = exp(qs[["97.5%"]]),
    p_negative = mean(d < 0), stringsAsFactors = FALSE
  )
}

# Difference posterior (a - reference) with index-wise pairing on common n.
diff_draws <- function(d_a, d_ref) {
  n  <- min(length(d_a), length(d_ref))
  dd <- d_a[seq_len(n)] - d_ref[seq_len(n)]
  qs <- quantile(dd, q)
  data.frame(
    diff_mean = mean(dd), diff_median = qs[["50%"]],
    diff_lo = qs[["2.5%"]], diff_hi = qs[["97.5%"]],
    p_gt_0 = mean(dd > 0), p_lt_0 = mean(dd < 0), stringsAsFactors = FALSE
  )
}

run_one <- function(id) {
  cfg     <- get_analysis(id)
  out_dir_cmp <- comparison_out_dir(paste0("param_", id))
  dir.create(out_dir_cmp, recursive = TRUE, showWarnings = FALSE)

  cat("\n=====================================================================\n")
  cat(sprintf("Parameterization comparison for analysis '%s'\n", id))
  cat(sprintf("  predictor: %s  (Study == '%s')\n", cfg$predictor_label, cfg$predictor_study))
  cat(sprintf("  outcome  : %s\n", cfg$outcome_label))
  cat(sprintf("  comparing: %s\n", paste(PARAM_IDS, collapse = ", ")))
  cat("=====================================================================\n")

  prep <- prepare_cop_data(cfg$data_file, cfg$predictor_study, quiet = TRUE)

  fits <- setNames(lapply(PARAM_IDS, function(p) {
    out_dir <- analysis_out_dir(id, p)
    cat(sprintf("\n--- Fitting '%s' -> %s ---\n", p, out_dir))
    fit <- fit_cop(prep, out_dir, parameterization = p, quiet = TRUE)
    slope_summary(fit$samp, prep, id, cfg$predictor_label, cfg$outcome_label, out_dir = out_dir)
    cat(sprintf("  max PSRF %.3f, min ESS %d\n", fit$max_psrf, round(fit$min_ess)))
    fit
  }), PARAM_IDS)

  # ================= GLOBAL slope (mu_b1) =================
  draws <- lapply(fits, function(f) as.matrix(f$samp)[, "mu_b1"])
  others <- setdiff(PARAM_IDS, REFERENCE)

  by_param <- do.call(rbind, lapply(PARAM_IDS, function(p) {
    cbind(data.frame(analysis_id = id, parameterization = p, stringsAsFactors = FALSE),
          setNames(summarise_draws(draws[[p]]),
                   c("mu_b1_mean", "mu_b1_median", "mu_b1_lo", "mu_b1_hi",
                     "rr_median", "rr_lo", "rr_hi", "p_negative")))
  }))
  write.csv(by_param, file.path(out_dir_cmp, "slopes_by_parameterization.csv"), row.names = FALSE)
  cat("\nGlobal slope (mu_b1, log-RR per unit centered log-GMC), by parameterization:\n")
  print(format(by_param[, c("parameterization", "mu_b1_median", "mu_b1_lo", "mu_b1_hi",
                           "rr_median", "p_negative")], digits = 3), row.names = FALSE)

  diffs <- do.call(rbind, lapply(others, function(p) {
    cbind(data.frame(contrast = sprintf("%s - %s", p, REFERENCE), stringsAsFactors = FALSE),
          diff_draws(draws[[p]], draws[[REFERENCE]]))
  }))
  write.csv(diffs, file.path(out_dir_cmp, "slope_pairwise_diffs.csv"), row.names = FALSE)
  cat(sprintf("\nPairwise GLOBAL slope differences (reference = %s):\n", REFERENCE))
  print(format(diffs, digits = 3), row.names = FALSE)

  plot_title <- sprintf("COP slope by parameterization\nPredictor: %s | Outcome: %s",
                        cfg$predictor_label, cfg$outcome_label)

  if (has_ggplot) {
    library(ggplot2)
    pf <- by_param
    pf$parameterization <- factor(pf$parameterization, levels = rev(PARAM_IDS))
    forest <- ggplot(pf, aes(x = mu_b1_median, y = parameterization)) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      geom_errorbarh(aes(xmin = mu_b1_lo, xmax = mu_b1_hi), height = 0.15) +
      geom_point(size = 3) +
      labs(x = "Global COP slope  mu_b1  (log rate-ratio per unit centered log-GMC)",
           y = NULL, title = plot_title,
           subtitle = "Points = posterior median; bars = 95% credible interval") +
      theme_bw(base_size = 12)
    ggsave(file.path(out_dir_cmp, "slope_forest.pdf"), forest, width = 8, height = 4)
    ggsave(file.path(out_dir_cmp, "slope_forest.png"), forest, width = 8, height = 4, dpi = 150)

    dens <- do.call(rbind, lapply(PARAM_IDS, function(p) {
      data.frame(parameterization = p, mu_b1 = draws[[p]])
    }))
    dens$parameterization <- factor(dens$parameterization, levels = PARAM_IDS)
    overlay <- ggplot(dens, aes(x = mu_b1, fill = parameterization, colour = parameterization)) +
      geom_density(alpha = 0.3) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      labs(x = "Global COP slope  mu_b1  (log rate-ratio per unit centered log-GMC)",
           y = "Posterior density", fill = "Parameterization", colour = "Parameterization",
           title = plot_title) +
      theme_bw(base_size = 12)
    ggsave(file.path(out_dir_cmp, "slope_overlay.pdf"), overlay, width = 8, height = 5)
    ggsave(file.path(out_dir_cmp, "slope_overlay.png"), overlay, width = 8, height = 5, dpi = 150)
  } else {
    ord <- rev(seq_along(PARAM_IDS))
    pdf(file.path(out_dir_cmp, "slope_forest.pdf"), width = 8, height = 4)
    plot(by_param$mu_b1_median, ord, xlim = range(by_param[, c("mu_b1_lo", "mu_b1_hi")]),
         ylim = c(0.5, length(ord) + 0.5), yaxt = "n", pch = 19,
         xlab = "Global COP slope mu_b1 (log RR per unit centered log-GMC)", ylab = "",
         main = plot_title)
    abline(v = 0, lty = 3, col = "grey50")
    arrows(by_param$mu_b1_lo, ord, by_param$mu_b1_hi, ord, code = 3, angle = 90, length = 0.03)
    axis(2, at = ord, labels = PARAM_IDS, las = 1)
    dev.off()

    pdf(file.path(out_dir_cmp, "slope_overlay.pdf"), width = 8, height = 5)
    ds <- lapply(draws, density)
    xlim <- range(sapply(ds, function(z) range(z$x)))
    ylim <- range(sapply(ds, function(z) range(z$y)))
    cols <- seq_along(ds)
    plot(NA, xlim = xlim, ylim = ylim, xlab = "mu_b1", ylab = "Posterior density",
         main = plot_title)
    for (k in seq_along(ds)) lines(ds[[k]], col = cols[k], lwd = 2)
    abline(v = 0, lty = 3, col = "grey50")
    legend("topright", legend = PARAM_IDS, col = cols, lwd = 2, bty = "n")
    dev.off()
  }

  # ================= SEROTYPE-SPECIFIC slope (b1[s]) =================
  all_sero <- sero_order(prep$serotypes)
  sero_sum <- do.call(rbind, lapply(PARAM_IDS, function(p) {
    M <- as.matrix(fits[[p]]$samp)
    do.call(rbind, lapply(seq_len(prep$S), function(s) {
      cbind(data.frame(parameterization = p, serotype = prep$serotypes[s], stringsAsFactors = FALSE),
            summarise_draws(M[, sprintf("b1[%d]", s)]))
    }))
  }))
  sero_sum <- sero_sum[order(match(sero_sum$serotype, all_sero),
                             match(sero_sum$parameterization, PARAM_IDS)), ]
  write.csv(sero_sum, file.path(out_dir_cmp, "serotype_slopes_by_parameterization.csv"),
           row.names = FALSE)

  sero_diffs <- do.call(rbind, lapply(others, function(p) {
    Mp   <- as.matrix(fits[[p]]$samp)
    Mref <- as.matrix(fits[[REFERENCE]]$samp)
    do.call(rbind, lapply(seq_len(prep$S), function(s) {
      cbind(data.frame(contrast = sprintf("%s - %s", p, REFERENCE),
                       serotype = prep$serotypes[s], stringsAsFactors = FALSE),
            diff_draws(Mp[, sprintf("b1[%d]", s)], Mref[, sprintf("b1[%d]", s)]))
    }))
  }))
  write.csv(sero_diffs, file.path(out_dir_cmp, "serotype_slope_pairwise_diffs.csv"),
           row.names = FALSE)

  cat("\nSerotype-specific slope (b1[s], RR per unit centered log-GMC), by parameterization:\n")
  wide <- reshape(sero_sum[, c("serotype", "parameterization", "rr_median")],
                  idvar = "serotype", timevar = "parameterization", direction = "wide")
  names(wide) <- sub("^rr_median\\.", "", names(wide))
  wide <- wide[match(all_sero, wide$serotype), ]
  print(format(wide, digits = 3), row.names = FALSE)

  if (has_ggplot) {
    library(ggplot2)
    pf <- sero_sum
    pf$serotype        <- factor(pf$serotype, levels = all_sero)
    pf$parameterization <- factor(pf$parameterization, levels = rev(PARAM_IDS))
    n_facet <- length(all_sero)
    g <- ggplot(pf, aes(x = b1_median, y = parameterization)) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      geom_errorbarh(aes(xmin = b1_lo, xmax = b1_hi), height = 0.15) +
      geom_point(size = 2) +
      facet_wrap(~ serotype) +
      labs(x = "Serotype-specific COP slope  b1[s]  (log RR per unit centered log-GMC)",
           y = NULL, title = paste0(plot_title, " - by serotype"),
           subtitle = "Points = posterior median; bars = 95% credible interval") +
      theme_bw(base_size = 11)
    h <- 2 + 1.6 * ceiling(n_facet / 3)
    ggsave(file.path(out_dir_cmp, "serotype_slope_forest.pdf"), g, width = 9, height = h)
    ggsave(file.path(out_dir_cmp, "serotype_slope_forest.png"), g, width = 9, height = h, dpi = 150)
  } else {
    np <- length(PARAM_IDS)
    yv <- setNames(seq_along(PARAM_IDS), PARAM_IDS)
    pdf(file.path(out_dir_cmp, "serotype_slope_forest.pdf"), width = 9,
        height = 2 + 0.5 * length(all_sero) * np)
    op <- par(mfrow = c(ceiling(length(all_sero) / 3), 3), mar = c(4, 5, 2, 1))
    for (s in all_sero) {
      sub <- sero_sum[sero_sum$serotype == s, ]
      yy  <- yv[sub$parameterization]
      plot(sub$b1_median, yy, xlim = range(sub[, c("b1_lo", "b1_hi")]),
           ylim = c(0.5, np + 0.5), yaxt = "n", pch = 19,
           xlab = "b1[s] (log RR)", ylab = "", main = paste("Serotype", s))
      abline(v = 0, lty = 3, col = "grey50")
      arrows(sub$b1_lo, yy, sub$b1_hi, yy, code = 3, angle = 90, length = 0.03)
      axis(2, at = yv, labels = names(yv), las = 1, cex.axis = 0.7)
    }
    par(op)
    dev.off()
  }

  # ================= Absolute risk vs absolute GMC (RE fit only) =================
  plot_cop_absolute(prep, fits[["RE"]]$samp, out_dir_cmp,
                    title_suffix = sprintf("Predictor: %s", cfg$predictor_label))

  cat(sprintf("\nSaved parameterization-comparison outputs to %s\n", out_dir_cmp))
  invisible(list(by_param = by_param, diffs = diffs,
                sero_sum = sero_sum, sero_diffs = sero_diffs))
}

invisible(lapply(ids, run_one))
cat("\nAll requested parameterization comparisons complete.\n")

# =====================================================================
# Compare model specs (candidate predictor-set combinations) by WAIC, and
# show their global COP slope (mu_b1) side by side.
#
# Currently a no-op: R/config.R's MODEL_SPECS has a single spec ("pooled")
# and MODEL_COMPARISONS is empty, since there's nothing to rank yet. This
# script is scaffolding for when you register an alternative predictor set
# and want to rank it against "pooled" -- see the "How to add an alternative
# model spec" instructions at the top of R/config.R.
#
# Usage (from the project root), once a comparison exists:
#   Rscript R/compare_waic.R                  # run every comparison
#   Rscript R/compare_waic.R <comparison id>  # a specific comparison
#
# For a comparison <id> it reads results/<spec>/waic.csv and mcmc.rds for
# each member spec (see R/config.R MODEL_COMPARISONS) and writes
# results/comparisons/<id>/:
#   waic_comparison.csv          WAIC leaderboard, ranked best-first
#   waic_comparison.pdf (+ .png) forest plot of elpd_waic with its SE
#   slopes_by_spec.csv           per-spec global mu_b1 (log-RR) + rr_global
#   slope_pairwise_diffs.csv     posterior of (spec - reference) slope diffs
#   slope_forest.pdf (+ .png)    forest plot of mu_b1 with 95% CrI
#   slope_overlay.pdf (+ .png)   overlaid posterior densities of mu_b1
#
# All specs in a comparison must share the SAME outcome studies (WAIC only
# compares models fit to the same data -- see R/cop_model.R compute_waic()
# and the header of JAGS/cop_eiv_model_multistudy.jags).
#
# Because each spec is an independent MCMC run, the slope difference
# posterior is formed by pairing draws across chains index-wise (any pairing
# is valid under independence); runs are truncated to a common draw count.
# =====================================================================

source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(coda))
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
ids  <- if (length(args) == 0) names(MODEL_COMPARISONS) else args

if (length(ids) == 0) {
  cat("No MODEL_COMPARISONS defined in R/config.R yet -- there's only one model\n",
      "spec ('pooled') right now, so nothing to rank. Add a second spec + a\n",
      "MODEL_COMPARISONS entry when you have an alternative to compare it against.\n",
      sep = "")
  quit(save = "no", status = 0)
}

q <- c(0.025, 0.5, 0.975)

load_spec_outputs <- function(id) {
  dir <- spec_out_dir(id)
  wf  <- file.path(dir, "waic.csv")
  mf  <- file.path(dir, "mcmc.rds")
  if (!file.exists(wf) || !file.exists(mf)) {
    stop("Missing outputs for spec '", id, "' - run R/run_models.R ", id, " first.")
  }
  list(waic = read.csv(wf, stringsAsFactors = FALSE),
       mu_b1 = as.matrix(readRDS(mf))[, "mu_b1"])
}

summarise_draws <- function(d) {
  qs <- quantile(d, q)
  data.frame(
    b1_mean = mean(d), b1_median = qs[["50%"]],
    b1_lo = qs[["2.5%"]], b1_hi = qs[["97.5%"]],
    rr_median = exp(qs[["50%"]]), rr_lo = exp(qs[["2.5%"]]), rr_hi = exp(qs[["97.5%"]]),
    p_negative = mean(d < 0), stringsAsFactors = FALSE
  )
}

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

run_comparison <- function(id) {
  cmp     <- get_model_comparison(id)
  out_dir <- comparison_out_dir(id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n=====================================================================\n")
  cat(sprintf("Comparison '%s': %s\n", id, cmp$label))
  cat(sprintf("  specs : %s\n", paste(cmp$specs, collapse = ", ")))
  cat("=====================================================================\n")

  loaded <- lapply(cmp$specs, load_spec_outputs)
  names(loaded) <- cmp$specs
  ref    <- if (!is.null(cmp$reference)) cmp$reference else cmp$specs[1]
  others <- setdiff(cmp$specs, ref)

  # ================= WAIC leaderboard =================
  waic_tab <- do.call(rbind, lapply(cmp$specs, function(s) {
    cbind(data.frame(spec_id = s, stringsAsFactors = FALSE), loaded[[s]]$waic)
  }))
  best <- waic_tab$spec_id[which.min(waic_tab$waic)]
  waic_tab$delta_waic <- waic_tab$waic - min(waic_tab$waic)
  rel_lik <- exp(-0.5 * waic_tab$delta_waic)
  waic_tab$waic_weight <- rel_lik / sum(rel_lik)
  waic_tab <- waic_tab[order(waic_tab$waic), ]
  write.csv(waic_tab, file.path(out_dir, "waic_comparison.csv"), row.names = FALSE)

  cat(sprintf("\nWAIC leaderboard (lower WAIC = better fit to the outcome data; best = '%s'):\n", best))
  print(format(waic_tab[, c("spec_id", "waic", "se_waic", "delta_waic", "p_waic", "waic_weight")],
              digits = 4), row.names = FALSE)

  waic_tab$spec_id <- factor(waic_tab$spec_id, levels = rev(cmp$specs))
  if (has_ggplot) {
    library(ggplot2)
    g <- ggplot(waic_tab, aes(x = elpd_waic, y = spec_id)) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      geom_errorbarh(aes(xmin = elpd_waic - se_waic, xmax = elpd_waic + se_waic), height = 0.15) +
      geom_point(size = 3) +
      labs(x = "elpd_waic (higher = better fit to the outcome data)", y = NULL,
           title = cmp$label, subtitle = "Points = elpd_waic; bars = +/- 1 SE") +
      theme_bw(base_size = 12)
    ggsave(file.path(out_dir, "waic_comparison.pdf"), g, width = 8, height = 4)
    ggsave(file.path(out_dir, "waic_comparison.png"), g, width = 8, height = 4, dpi = 150)
  } else {
    ord <- rev(seq_along(cmp$specs))
    pdf(file.path(out_dir, "waic_comparison.pdf"), width = 8, height = 4)
    plot(waic_tab$elpd_waic, ord,
         xlim = range(c(waic_tab$elpd_waic - waic_tab$se_waic, waic_tab$elpd_waic + waic_tab$se_waic)),
         ylim = c(0.5, length(ord) + 0.5), yaxt = "n", pch = 19,
         xlab = "elpd_waic (higher = better)", ylab = "", main = cmp$label)
    abline(v = 0, lty = 3, col = "grey50")
    arrows(waic_tab$elpd_waic - waic_tab$se_waic, ord, waic_tab$elpd_waic + waic_tab$se_waic, ord,
           code = 3, angle = 90, length = 0.03)
    axis(2, at = ord, labels = as.character(waic_tab$spec_id)[order(ord)], las = 1)
    dev.off()
  }

  # ================= Global slope (mu_b1) =================
  draws <- lapply(loaded, `[[`, "mu_b1")

  by_spec <- do.call(rbind, lapply(cmp$specs, function(s) {
    cbind(data.frame(spec_id = s, stringsAsFactors = FALSE),
          setNames(summarise_draws(draws[[s]]),
                   c("mu_b1_mean","mu_b1_median","mu_b1_lo","mu_b1_hi",
                     "rr_median","rr_lo","rr_hi","p_negative")))
  }))
  write.csv(by_spec, file.path(out_dir, "slopes_by_spec.csv"), row.names = FALSE)
  cat("\nGlobal slope (mu_b1, log-RR per unit absolute log-GMC):\n")
  print(format(by_spec[, c("spec_id", "mu_b1_median", "mu_b1_lo", "mu_b1_hi",
                           "rr_median", "p_negative")], digits = 3), row.names = FALSE)

  diffs <- do.call(rbind, lapply(others, function(s) {
    cbind(data.frame(contrast = sprintf("%s - %s", s, ref), stringsAsFactors = FALSE),
          diff_draws(draws[[s]], draws[[ref]]))
  }))
  write.csv(diffs, file.path(out_dir, "slope_pairwise_diffs.csv"), row.names = FALSE)
  cat(sprintf("\nPairwise slope differences (reference = %s):\n", ref))
  print(format(diffs, digits = 3), row.names = FALSE)

  if (has_ggplot) {
    library(ggplot2)
    pf <- by_spec
    pf$spec_id <- factor(pf$spec_id, levels = rev(cmp$specs))
    forest <- ggplot(pf, aes(x = mu_b1_median, y = spec_id)) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      geom_errorbarh(aes(xmin = mu_b1_lo, xmax = mu_b1_hi), height = 0.15) +
      geom_point(size = 3) +
      labs(x = "Global COP slope  mu_b1  (log rate-ratio per unit absolute log-GMC)",
           y = NULL, title = cmp$label,
           subtitle = "Points = posterior median; bars = 95% credible interval") +
      theme_bw(base_size = 12)
    ggsave(file.path(out_dir, "slope_forest.pdf"), forest, width = 8, height = 4)
    ggsave(file.path(out_dir, "slope_forest.png"), forest, width = 8, height = 4, dpi = 150)

    dens <- do.call(rbind, lapply(cmp$specs, function(s) {
      data.frame(spec_id = s, mu_b1 = draws[[s]])
    }))
    dens$spec_id <- factor(dens$spec_id, levels = cmp$specs)
    overlay <- ggplot(dens, aes(x = mu_b1, fill = spec_id, colour = spec_id)) +
      geom_density(alpha = 0.3) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      labs(x = "Global COP slope  mu_b1  (log rate-ratio per unit absolute log-GMC)",
           y = "Posterior density", fill = "Model spec", colour = "Model spec",
           title = cmp$label) +
      theme_bw(base_size = 12)
    ggsave(file.path(out_dir, "slope_overlay.pdf"), overlay, width = 8, height = 5)
    ggsave(file.path(out_dir, "slope_overlay.png"), overlay, width = 8, height = 5, dpi = 150)
  } else {
    ord <- rev(seq_along(cmp$specs))
    pdf(file.path(out_dir, "slope_forest.pdf"), width = 8, height = 4)
    plot(by_spec$mu_b1_median, ord, xlim = range(by_spec[, c("mu_b1_lo", "mu_b1_hi")]),
         ylim = c(0.5, length(ord) + 0.5), yaxt = "n", pch = 19,
         xlab = "Global COP slope mu_b1 (log RR per unit absolute log-GMC)", ylab = "",
         main = cmp$label)
    abline(v = 0, lty = 3, col = "grey50")
    arrows(by_spec$mu_b1_lo, ord, by_spec$mu_b1_hi, ord, code = 3, angle = 90, length = 0.03)
    axis(2, at = ord, labels = cmp$specs, las = 1)
    dev.off()

    pdf(file.path(out_dir, "slope_overlay.pdf"), width = 8, height = 5)
    ds <- lapply(draws, density)
    xlim <- range(sapply(ds, function(z) range(z$x)))
    ylim <- range(sapply(ds, function(z) range(z$y)))
    cols <- seq_along(ds)
    plot(NA, xlim = xlim, ylim = ylim, xlab = "mu_b1", ylab = "Posterior density",
         main = cmp$label)
    for (k in seq_along(ds)) lines(ds[[k]], col = cols[k], lwd = 2)
    abline(v = 0, lty = 3, col = "grey50")
    legend("topright", legend = cmp$specs, col = cols, lwd = 2, bty = "n")
    dev.off()
  }

  cat(sprintf("\nSaved comparison outputs to %s\n", out_dir))
  invisible(list(waic_tab = waic_tab, by_spec = by_spec, diffs = diffs))
}

invisible(lapply(ids, run_comparison))
cat("\nAll requested comparisons complete.\n")

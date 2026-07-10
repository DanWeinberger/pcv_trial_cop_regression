# =====================================================================
# Compare the COP slope across analyses - GLOBAL (mu_b1) and per SEROTYPE.
#
# Usage (from the project root):
#   Rscript R/compare_slopes.R                    # run every comparison
#   Rscript R/compare_slopes.R outcome_nckp       # a specific comparison id
#
# For a comparison <id> it reads results/<analysis>/mcmc.rds (+ slope_summary.csv
# for the serotype labels) for each member analysis and writes
# results/comparisons/<id>/:
#   GLOBAL slope (mu_b1):
#     slopes_by_analysis.csv       per-analysis mu_b1 (log-RR) + rr_global
#     slope_pairwise_diffs.csv     posterior of (analysis - reference) slope diffs
#     slope_forest.pdf (+ .png)    forest plot of mu_b1 with 95% CrI
#     slope_overlay.pdf (+ .png)   overlaid posterior densities of mu_b1
#   SEROTYPE-SPECIFIC slope (b1[s]):
#     serotype_slopes_by_analysis.csv    per-analysis x serotype b1 + RR summary
#     serotype_slope_pairwise_diffs.csv  per-serotype (analysis - reference) diffs
#     serotype_slope_forest.pdf (+ .png) forest plot faceted by serotype
#
# Because each analysis is an independent MCMC run on the SAME outcome, the
# difference posterior is formed by pairing draws across chains index-wise
# (any pairing is valid under independence); runs are truncated to a common
# draw count first. The same pairing is used at the global and serotype level.
# =====================================================================

source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(coda))
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
ids  <- if (length(args) == 0) names(COMPARISONS) else args

q <- c(0.025, 0.5, 0.975)

# Numeric-aware serotype ordering: 4, 6A, 6B, 9V, 14, 18C, 19A, 19F, 23F.
sero_order <- function(x) {
  num <- suppressWarnings(as.numeric(gsub("[^0-9].*$", "", x)))
  x[order(num, x)]
}

# Load mu_b1 draws + per-serotype b1 draws (labelled by serotype) for one
# analysis. The non-"global" rows of slope_summary.csv are in b1[] index order,
# so they map the b1[k] columns to their serotype labels.
load_analysis_draws <- function(analysis_id) {
  f  <- file.path(analysis_out_dir(analysis_id), "mcmc.rds")
  ss <- file.path(analysis_out_dir(analysis_id), "slope_summary.csv")
  if (!file.exists(f))
    stop("Missing ", f,  " - run R/run_analysis.R ", analysis_id, " first.")
  if (!file.exists(ss))
    stop("Missing ", ss, " - run R/run_analysis.R ", analysis_id, " first.")
  M    <- as.matrix(readRDS(f))
  sero <- subset(read.csv(ss, stringsAsFactors = FALSE), level != "global")$level
  by_sero <- setNames(
    lapply(seq_along(sero), function(k) M[, sprintf("b1[%d]", k)]),
    sero
  )
  list(mu_b1 = M[, "mu_b1"], by_sero = by_sero)
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

run_comparison <- function(id) {
  cmp     <- get_comparison(id)
  out_dir <- comparison_out_dir(id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n=====================================================================\n")
  cat(sprintf("Comparison '%s': %s\n", id, cmp$label))
  cat(sprintf("  analyses : %s\n", paste(cmp$analyses, collapse = ", ")))
  cat("=====================================================================\n")

  labels <- vapply(cmp$analyses, function(a) get_analysis(a)$predictor_label, character(1))
  ad     <- lapply(cmp$analyses, load_analysis_draws)
  names(ad) <- cmp$analyses
  draws  <- lapply(ad, `[[`, "mu_b1")     # global mu_b1 draws, by analysis
  ref    <- if (!is.null(cmp$reference)) cmp$reference else cmp$analyses[1]
  others <- setdiff(cmp$analyses, ref)

  # ================= GLOBAL slope (mu_b1) =================
  by_analysis <- do.call(rbind, lapply(seq_along(cmp$analyses), function(k) {
    cbind(data.frame(analysis_id = cmp$analyses[k], predictor_label = labels[k],
                     stringsAsFactors = FALSE),
          setNames(summarise_draws(draws[[k]]),
                   c("mu_b1_mean","mu_b1_median","mu_b1_lo","mu_b1_hi",
                     "rr_median","rr_lo","rr_hi","p_negative")))
  }))
  write.csv(by_analysis, file.path(out_dir, "slopes_by_analysis.csv"), row.names = FALSE)
  cat("\nGlobal slope (mu_b1, log-RR per unit centered log-GMC):\n")
  print(format(by_analysis[, c("analysis_id", "mu_b1_median", "mu_b1_lo", "mu_b1_hi",
                               "rr_median", "p_negative")], digits = 3),
        row.names = FALSE)

  diffs <- do.call(rbind, lapply(others, function(a) {
    cbind(data.frame(contrast = sprintf("%s - %s", a, ref), stringsAsFactors = FALSE),
          diff_draws(draws[[a]], draws[[ref]]))
  }))
  write.csv(diffs, file.path(out_dir, "slope_pairwise_diffs.csv"), row.names = FALSE)
  cat(sprintf("\nPairwise GLOBAL slope differences (reference = %s):\n", ref))
  print(format(diffs, digits = 3), row.names = FALSE)

  # ---- global forest + overlay ----
  if (has_ggplot) {
    library(ggplot2)
    pf <- by_analysis
    pf$analysis_id <- factor(pf$analysis_id, levels = rev(cmp$analyses))
    forest <- ggplot(pf, aes(x = mu_b1_median, y = analysis_id)) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      geom_errorbarh(aes(xmin = mu_b1_lo, xmax = mu_b1_hi), height = 0.15) +
      geom_point(size = 3) +
      labs(x = "Global COP slope  mu_b1  (log rate-ratio per unit centered log-GMC)",
           y = NULL, title = cmp$label,
           subtitle = "Points = posterior median; bars = 95% credible interval") +
      theme_bw(base_size = 12)
    ggsave(file.path(out_dir, "slope_forest.pdf"), forest, width = 8, height = 4)
    ggsave(file.path(out_dir, "slope_forest.png"), forest, width = 8, height = 4, dpi = 150)

    dens <- do.call(rbind, lapply(seq_along(cmp$analyses), function(k) {
      data.frame(analysis_id = cmp$analyses[k], mu_b1 = draws[[k]])
    }))
    dens$analysis_id <- factor(dens$analysis_id, levels = cmp$analyses)
    overlay <- ggplot(dens, aes(x = mu_b1, fill = analysis_id, colour = analysis_id)) +
      geom_density(alpha = 0.3) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      labs(x = "Global COP slope  mu_b1  (log rate-ratio per unit centered log-GMC)",
           y = "Posterior density", fill = "Analysis", colour = "Analysis",
           title = cmp$label) +
      theme_bw(base_size = 12)
    ggsave(file.path(out_dir, "slope_overlay.pdf"), overlay, width = 8, height = 5)
    ggsave(file.path(out_dir, "slope_overlay.png"), overlay, width = 8, height = 5, dpi = 150)
  } else {
    ord <- rev(seq_along(cmp$analyses))
    pdf(file.path(out_dir, "slope_forest.pdf"), width = 8, height = 4)
    plot(by_analysis$mu_b1_median, ord, xlim = range(by_analysis[, c("mu_b1_lo", "mu_b1_hi")]),
         ylim = c(0.5, length(ord) + 0.5), yaxt = "n", pch = 19,
         xlab = "Global COP slope mu_b1 (log RR per unit centered log-GMC)", ylab = "",
         main = cmp$label)
    abline(v = 0, lty = 3, col = "grey50")
    arrows(by_analysis$mu_b1_lo, ord, by_analysis$mu_b1_hi, ord,
           code = 3, angle = 90, length = 0.03)
    axis(2, at = ord, labels = cmp$analyses, las = 1)
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
    legend("topright", legend = cmp$analyses, col = cols, lwd = 2, bty = "n")
    dev.off()
  }

  # ================= SEROTYPE-SPECIFIC slope (b1[s]) =================
  sero_sum <- do.call(rbind, lapply(cmp$analyses, function(a) {
    bs <- ad[[a]]$by_sero
    do.call(rbind, lapply(names(bs), function(s) {
      cbind(data.frame(analysis_id = a, predictor_label = get_analysis(a)$predictor_label,
                       serotype = s, stringsAsFactors = FALSE),
            summarise_draws(bs[[s]]))
    }))
  }))
  all_sero <- sero_order(unique(sero_sum$serotype))
  sero_sum <- sero_sum[order(match(sero_sum$serotype, all_sero),
                             match(sero_sum$analysis_id, cmp$analyses)), ]
  write.csv(sero_sum, file.path(out_dir, "serotype_slopes_by_analysis.csv"), row.names = FALSE)

  ref_bs <- ad[[ref]]$by_sero
  sero_diffs <- do.call(rbind, lapply(others, function(a) {
    bs     <- ad[[a]]$by_sero
    common <- sero_order(intersect(names(bs), names(ref_bs)))
    if (length(common) == 0) return(NULL)
    do.call(rbind, lapply(common, function(s) {
      cbind(data.frame(contrast = sprintf("%s - %s", a, ref), serotype = s,
                       stringsAsFactors = FALSE),
            diff_draws(bs[[s]], ref_bs[[s]]))
    }))
  }))
  if (!is.null(sero_diffs))
    write.csv(sero_diffs, file.path(out_dir, "serotype_slope_pairwise_diffs.csv"),
              row.names = FALSE)

  cat("\nSerotype-specific slope (b1[s], RR per unit centered log-GMC):\n")
  wide <- reshape(sero_sum[, c("serotype", "analysis_id", "rr_median")],
                  idvar = "serotype", timevar = "analysis_id", direction = "wide")
  names(wide) <- sub("^rr_median\\.", "", names(wide))
  wide <- wide[match(all_sero, wide$serotype), ]
  print(format(wide, digits = 3), row.names = FALSE)

  # ---- serotype-faceted forest ----
  if (has_ggplot) {
    library(ggplot2)
    pf <- sero_sum
    pf$serotype    <- factor(pf$serotype, levels = all_sero)
    pf$analysis_id <- factor(pf$analysis_id, levels = rev(cmp$analyses))
    n_facet <- length(all_sero)
    g <- ggplot(pf, aes(x = b1_median, y = analysis_id)) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
      geom_errorbarh(aes(xmin = b1_lo, xmax = b1_hi), height = 0.15) +
      geom_point(size = 2) +
      facet_wrap(~ serotype) +
      labs(x = "Serotype-specific COP slope  b1[s]  (log RR per unit centered log-GMC)",
           y = NULL, title = paste0(cmp$label, " - by serotype"),
           subtitle = "Points = posterior median; bars = 95% credible interval") +
      theme_bw(base_size = 11)
    h <- 2 + 1.6 * ceiling(n_facet / 3)
    ggsave(file.path(out_dir, "serotype_slope_forest.pdf"), g, width = 9, height = h)
    ggsave(file.path(out_dir, "serotype_slope_forest.png"), g, width = 9, height = h, dpi = 150)
  } else {
    na <- length(cmp$analyses)
    yv <- setNames(seq_along(cmp$analyses), cmp$analyses)
    pdf(file.path(out_dir, "serotype_slope_forest.pdf"), width = 9,
        height = 2 + 0.5 * length(all_sero) * na)
    op <- par(mfrow = c(ceiling(length(all_sero) / 3), 3), mar = c(4, 5, 2, 1))
    for (s in all_sero) {
      sub <- sero_sum[sero_sum$serotype == s, ]
      yy  <- yv[sub$analysis_id]
      plot(sub$b1_median, yy, xlim = range(sub[, c("b1_lo", "b1_hi")]),
           ylim = c(0.5, na + 0.5), yaxt = "n", pch = 19,
           xlab = "b1[s] (log RR)", ylab = "", main = paste("Serotype", s))
      abline(v = 0, lty = 3, col = "grey50")
      arrows(sub$b1_lo, yy, sub$b1_hi, yy, code = 3, angle = 90, length = 0.03)
      axis(2, at = yv, labels = names(yv), las = 1, cex.axis = 0.7)
    }
    par(op)
    dev.off()
  }

  cat(sprintf("\nSaved comparison outputs to %s\n", out_dir))
  invisible(list(by_analysis = by_analysis, diffs = diffs,
                 sero_sum = sero_sum, sero_diffs = sero_diffs))
}

invisible(lapply(ids, run_comparison))
cat("\nAll requested comparisons complete.\n")

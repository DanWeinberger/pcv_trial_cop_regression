# =====================================================================
# Plot the fitted correlate-of-protection curves from the EIV Poisson
# JAGS model (fit_eiv_poisson_cop.R).
#
#   y-axis : expected proportion of cases = lambda / Total_Cases
#            = exp(a[s] + b1[s] * x)               (Poisson rate)
#   x-axis : centered log-GMC (x); x = 0 is the unimmunized reference
#   lines  : one serotype-specific curve per serotype (a[s], b1[s])
#            + one global/pooled curve (mu_a, mu_b1)
#   points : observed proportions
#            - unimmunized arm at x = 0  (cases_u / total_u)
#            - immunized  arm at x = observed centered log-GMC
#                                        (cases_i / total_i)
#
# Run:  Rscript R/plot_eiv_poisson_cop.R   (from the project root)
# =====================================================================

suppressPackageStartupMessages({
  library(coda)
})
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

out_dir <- file.path("results", "eiv_poisson_cop")

# ---------------------------------------------------------------------
# 1. Re-derive the data exactly as in the fitting script
# ---------------------------------------------------------------------
d <- read.csv("data/siber_whitney_merged.csv", stringsAsFactors = FALSE)
d <- subset(d, Study == "NCKP")
u <- subset(d, Vaccine_Arm == "Unimmunized")
i <- subset(d, Vaccine_Arm == "Immunized")
serotypes <- intersect(u$Serotype, i$Serotype)
u <- u[match(serotypes, u$Serotype), ]
i <- i[match(serotypes, i$Serotype), ]
S <- length(serotypes)

# Observed centered log-GMC for the immunized arm (unimmunized = 0 reference).
x_obs   <- log(i$GMC) - log(u$GMC)
prop_u  <- u$Cases / u$Total_Cases      # observed proportion, unimmunized
prop_i  <- i$Cases / i$Total_Cases      # observed proportion, immunized

# Measurement SE of log-GMC from the 95% CI (same derivation as the fit script),
# used to draw horizontal (x-axis) uncertainty on the observed points.
log_se   <- function(lower, upper) (log(upper) - log(lower)) / (2 * qnorm(0.975))
se_floor <- 1e-3
se_u <- pmax(log_se(u$Lower_CL, u$Upper_CL), se_floor)
se_i <- pmax(log_se(i$Lower_CL, i$Upper_CL), se_floor)
# Unimmunized point sits at x = 0 (the reference) with its own measurement SE;
# immunized centered log-GMC combines both arms' errors in quadrature.
se_x <- c(se_u, sqrt(se_i^2 + se_u^2))

# ---------------------------------------------------------------------
# 2. Posterior draws
# ---------------------------------------------------------------------
samp <- readRDS(file.path(out_dir, "mcmc.rds"))
M    <- as.matrix(samp)

post_mean <- function(par) mean(M[, par])
a_hat   <- vapply(seq_len(S), function(s) post_mean(sprintf("a[%d]",  s)), numeric(1))
b1_hat  <- vapply(seq_len(S), function(s) post_mean(sprintf("b1[%d]", s)), numeric(1))
mu_a    <- post_mean("mu_a")
mu_b1   <- post_mean("mu_b1")

# Legend labels annotated with each serotype's posterior-mean slope b1[s]
# (the global/pooled slope mu_b1 is shown for reference).
slope_labels <- setNames(sprintf("%s  (b1 = %+.2f)", serotypes, b1_hat),
                         serotypes)

# ---------------------------------------------------------------------
# 3. Prediction grid + fitted curves
# ---------------------------------------------------------------------
xr <- range(c(0, x_obs))
xg <- seq(xr[1] - 0.2, xr[2] + 0.2, length.out = 200)

# serotype-specific fitted proportions
fit_s <- data.frame(
  Serotype = rep(serotypes, each = length(xg)),
  x        = rep(xg, times = S),
  prop     = as.vector(sapply(seq_len(S),
                              function(s) exp(a_hat[s] + b1_hat[s] * xg)))
)
# global (pooled) fitted proportion, with a posterior credible band
xg_mat   <- outer(M[, "mu_a"], rep(1, length(xg))) +
            outer(M[, "mu_b1"], xg)                       # draws x grid on log scale
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
  xlo      = c(rep(0, S), x_obs) - qnorm(0.975) * se_x,   # 95% x error bar
  xhi      = c(rep(0, S), x_obs) + qnorm(0.975) * se_x
)

# ---------------------------------------------------------------------
# 4. Plot
# ---------------------------------------------------------------------
if (has_ggplot) {
  library(ggplot2)
  p <- ggplot() +
    geom_ribbon(data = fit_glob, aes(x = x, ymin = lo, ymax = hi),
                fill = "grey60", alpha = 0.25) +
    geom_line(data = fit_s, aes(x = x, y = prop, colour = Serotype),
              size = 0.7) +
    geom_line(data = fit_glob, aes(x = x, y = med),
              colour = "black", size = 1.3) +
    geom_errorbarh(data = obs, aes(y = prop, xmin = xlo, xmax = xhi,
                                   colour = Serotype), height = 0, alpha = 0.7) +
    geom_point(data = obs, aes(x = x, y = prop, colour = Serotype,
                               shape = Arm), size = 2.6) +
    scale_shape_manual(values = c(Unimmunized = 1, Immunized = 16)) +
    scale_colour_discrete(labels = slope_labels) +
    scale_y_log10() +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
    labs(
      x = "Centered log-GMC  (immunized - unimmunized reference)",
      y = "Expected proportion of cases  ( lambda / Total_Cases, log scale )",
      colour = "Serotype (slope b1)",
      title = "Correlate of protection: fitted case proportion vs centered log-GMC",
      subtitle = sprintf(paste0("Coloured lines = serotype-specific fits; ",
                                "black line = global (pooled) fit, b1 = %+.2f, 95%% CrI; ",
                                "bars = 95%% GMC measurement error"), mu_b1)
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  ggsave(file.path(out_dir, "cop_fit_proportion.pdf"), p, width = 9, height = 6)
  ggsave(file.path(out_dir, "cop_fit_proportion.png"), p, width = 9, height = 6,
         dpi = 150)
} else {
  # base-R fallback
  cols <- setNames(rainbow(S), serotypes)
  pdf(file.path(out_dir, "cop_fit_proportion.pdf"), width = 9, height = 6)
  ylim <- range(c(fit_s$prop, obs$prop, fit_glob$hi))
  ylim[1] <- max(ylim[1], min(obs$prop[obs$prop > 0]) / 2)   # keep off zero for log
  plot(NA, xlim = range(xg), ylim = ylim, log = "y",
       xlab = "Centered log-GMC (immunized - unimmunized reference)",
       ylab = "Expected proportion of cases (lambda / Total_Cases, log scale)",
       main = "Correlate of protection: fitted case proportion vs centered log-GMC")
  polygon(c(fit_glob$x, rev(fit_glob$x)), c(fit_glob$lo, rev(fit_glob$hi)),
          col = adjustcolor("grey60", 0.25), border = NA)
  for (s in seq_len(S)) {
    lines(xg, exp(a_hat[s] + b1_hat[s] * xg), col = cols[s], lwd = 1.5)
  }
  lines(fit_glob$x, fit_glob$med, col = "black", lwd = 3)
  abline(v = 0, lty = 3, col = "grey50")
  # horizontal 95% GMC measurement-error bars, then points
  arrows(obs$xlo, obs$prop, obs$xhi, obs$prop, code = 3, angle = 90,
         length = 0.03, col = rep(cols, 2))
  # points: open = unimmunized, filled = immunized
  points(0 * seq_len(S), prop_u, col = cols, pch = 1, cex = 1.3)
  points(x_obs, prop_i, col = cols, pch = 16, cex = 1.3)
  legend("topright", legend = slope_labels, col = cols, lwd = 1.5,
         title = "Serotype (slope b1)", bty = "n", cex = 0.8)
  legend("top", legend = c("Global (pooled)", "Unimmunized (x=0)", "Immunized"),
         col = c("black", "grey30", "grey30"), lwd = c(3, NA, NA),
         pch = c(NA, 1, 16), bty = "n", cex = 0.8)
  dev.off()
}

cat("Saved cop_fit_proportion.pdf",
    if (has_ggplot) "and cop_fit_proportion.png" else "",
    "to", out_dir, "\n")

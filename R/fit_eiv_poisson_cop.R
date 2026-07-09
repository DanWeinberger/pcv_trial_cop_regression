# =====================================================================
# Bayesian error-in-variables Poisson regression for a pneumococcal
# correlate of protection (COP).
#
#   Outcome     : serotype-specific IPD case counts (Cases), Poisson
#   Offset      : Total_Cases  (used as log-exposure)
#   Predictor   : centered log-GMC, measured with error (precision from 95% CI)
#   Structure   : serotype-specific random intercepts a[s] ~ N(mu_a, sigma_a)
#                 hierarchical serotype slopes        b1[s] ~ N(mu_b1, sigma_b1)
#   Centering   : in-model latent reference per serotype -> unimmunized true
#                 log-GMC is the reference (=0 on predictor scale); its
#                 measurement error propagates into the immunized predictor.
#   Data        : NCKP immunogenicity only.
#
# Run:  Rscript R/fit_eiv_poisson_cop.R   (from the project root)
# =====================================================================

suppressPackageStartupMessages({
  library(rjags)
  library(coda)
})

set.seed(42)

# ---------------------------------------------------------------------
# 1. Load data and restrict to NCKP immunogenicity
# ---------------------------------------------------------------------
# Output directory: separate results folder with an analysis-specific subfolder.
out_dir <- file.path("results", "eiv_poisson_cop")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv("data/siber_whitney_merged.csv", stringsAsFactors = FALSE)
d <- subset(d, Study == "NCKP")

# Split into the two arms and align on serotype so each serotype gives one
# unimmunized (reference) and one immunized observation.
u <- subset(d, Vaccine_Arm == "Unimmunized")
i <- subset(d, Vaccine_Arm == "Immunized")

serotypes <- intersect(u$Serotype, i$Serotype)
u <- u[match(serotypes, u$Serotype), ]
i <- i[match(serotypes, i$Serotype), ]
S <- length(serotypes)
stopifnot(nrow(u) == S, nrow(i) == S)

# ---------------------------------------------------------------------
# 2. Measurement error of log-GMC from the 95% CI
#    SE(log GMC) = (log(Upper) - log(Lower)) / (2 * 1.96)
#    precision   = 1 / SE^2
# ---------------------------------------------------------------------
log_se <- function(lower, upper) (log(upper) - log(lower)) / (2 * qnorm(0.975))

se_u <- log_se(u$Lower_CL, u$Upper_CL)
se_i <- log_se(i$Lower_CL, i$Upper_CL)

# The normal-on-log-scale measurement-error prior below assumes the CI is
# symmetric about log(GMC): (log GMC - log Lower) ~= (log Upper - log GMC).
# Check this and flag serotypes where it fails materially. Aggressive rounding
# of small GMCs (e.g. unimmunized values ~0.03-0.1 reported to 2 dp) can break
# symmetry outright and even collapse one half-width to zero, in which case the
# derived SE reflects rounding granularity, not the study's true uncertainty.
log_asym <- function(lower, gmc, upper) {
  lo <- log(gmc) - log(lower)          # lower half-width on log scale
  hi <- log(upper) - log(gmc)          # upper half-width on log scale
  2 * (hi - lo) / (hi + lo)            # relative to mean half-width; 0 = symmetric
}
asym_tol <- 0.10                        # flag |asymmetry| > 10%
asym_u <- log_asym(u$Lower_CL, u$GMC, u$Upper_CL)
asym_i <- log_asym(i$Lower_CL, i$GMC, i$Upper_CL)
flagged <- data.frame(
  Serotype = c(u$Serotype, i$Serotype),
  Arm      = c(rep("Unimmunized", S), rep("Immunized", S)),
  asym_pct = round(100 * c(asym_u, asym_i), 1)
)
flagged <- flagged[abs(flagged$asym_pct) > 100 * asym_tol, ]
if (nrow(flagged)) {
  message("Note: ", nrow(flagged), " CI(s) not symmetric on the log scale ",
          "(|asymmetry| > ", 100 * asym_tol, "%); normal-on-log EIV prior is ",
          "misspecified for these:")
  print(flagged[order(-abs(flagged$asym_pct)), ], row.names = FALSE)
}

# Guard against zero-width CIs (rounded bounds that collide): floor the SE so
# precision stays finite. Report if any were adjusted.
se_floor <- 1e-3
if (any(se_u <= se_floor) || any(se_i <= se_floor)) {
  message("Note: floored ", sum(se_u <= se_floor) + sum(se_i <= se_floor),
          " zero/near-zero-width CI(s) at SE = ", se_floor)
}
se_u <- pmax(se_u, se_floor)
se_i <- pmax(se_i, se_floor)

jags_data <- list(
  S       = S,
  cases_u = u$Cases,        total_u = u$Total_Cases,
  cases_i = i$Cases,        total_i = i$Total_Cases,
  lgmc_u  = log(u$GMC),     prec_u  = 1 / se_u^2,
  lgmc_i  = log(i$GMC),     prec_i  = 1 / se_i^2
)

cat("Serotypes (order):", paste(serotypes, collapse = ", "), "\n")
cat("Observed centered log-GMC (immunized - unimmunized):\n")
print(setNames(round(log(i$GMC) - log(u$GMC), 3), serotypes))

# ---------------------------------------------------------------------
# 3. Initial values (data-driven) + per-chain RNG for reproducibility
# ---------------------------------------------------------------------
a_start  <- log(u$Cases / u$Total_Cases)
x_start  <- log(i$GMC) - log(u$GMC)
ref_start <- log(u$GMC)

make_inits <- function(seed) {
  list(
    ref = ref_start, x = x_start,
    a = a_start, b1 = rep(-0.5, S),
    mu_a = mean(a_start), sigma_a = 0.3,
    mu_b1 = -0.5, sigma_b1 = 0.3,
    mu_x = mean(x_start), sigma_x = 0.5,
    .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed
  )
}
inits <- list(make_inits(11), make_inits(22), make_inits(33))

# ---------------------------------------------------------------------
# 4. Compile, burn-in, sample
# ---------------------------------------------------------------------
params <- c("mu_b1", "b1", "rr", "rr_global",
            "a", "mu_a", "sigma_a", "sigma_b1",
            "mu_x", "sigma_x", "x", "ref")

model <- jags.model(file.path("JAGS", "cop_eiv_model.jags"), data = jags_data,
                    inits = inits, n.chains = 3, n.adapt = 2000)
update(model, 5000)                                   # burn-in
samp <- coda.samples(model, variable.names = params,
                     n.iter = 20000, thin = 5)

# ---------------------------------------------------------------------
# 5. Summaries & diagnostics
# ---------------------------------------------------------------------
summ <- summary(samp)
stats <- cbind(summ$statistics[, c("Mean", "SD")], summ$quantiles)
print(round(stats, 3))

cat("\n--- Global correlate of protection ---\n")
mu_b1_draws <- as.matrix(samp)[, "mu_b1"]
cat(sprintf("mu_b1 (log RR per unit centered log-GMC): %.3f  (95%% CrI %.3f, %.3f)\n",
            median(mu_b1_draws), quantile(mu_b1_draws, .025), quantile(mu_b1_draws, .975)))
rr_draws <- as.matrix(samp)[, "rr_global"]
cat(sprintf("rr_global (rate ratio per unit centered log-GMC): %.3f  (95%% CrI %.3f, %.3f)\n",
            median(rr_draws), quantile(rr_draws, .025), quantile(rr_draws, .975)))
cat(sprintf("P(mu_b1 < 0) = %.3f\n", mean(mu_b1_draws < 0)))

cat("\n--- Convergence (Gelman-Rubin, upper CI should be ~1) ---\n")
gd <- gelman.diag(samp, multivariate = FALSE)
print(round(gd$psrf[c("mu_b1", "sigma_b1", "mu_a", "sigma_a"), ], 3))
cat("Max PSRF (all monitored):", round(max(gd$psrf[, "Point est."], na.rm = TRUE), 3), "\n")
cat("Min effective sample size:", round(min(effectiveSize(samp))), "\n")

# ---------------------------------------------------------------------
# 6. Save outputs
# ---------------------------------------------------------------------
write.csv(round(stats, 4), file.path(out_dir, "posterior_summary.csv"))
saveRDS(samp, file.path(out_dir, "mcmc.rds"))

pdf(file.path(out_dir, "diagnostics.pdf"), width = 9, height = 6)
plot(samp[, c("mu_b1", "sigma_b1", "mu_a", "sigma_a")])
dev.off()

cat("\nSaved to ", out_dir, ": posterior_summary.csv, mcmc.rds, diagnostics.pdf\n", sep = "")

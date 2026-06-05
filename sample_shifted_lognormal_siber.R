# Sampling from fitted shifted lognormal distributions
# Model: X = theta + exp(mu + sigma * Z),  Z ~ N(0,1)
# Equivalently: (X - theta) ~ Lognormal(mu, sigma)
# GMC = exp(mu) + theta

# ------------------------------------------------------------
# Parameters from CDF fits (natural log scale)
# ------------------------------------------------------------

params <- data.frame(
  trial   = c("nckp", "nckp", "amer_indian", "amer_indian", "s_africa", "s_africa"),
  group   = c("Immunized", "Un-immunized", "Immunized", "Un-immunized", "Immunized", "Un-immunized"),
  mu      = c( 0.7031, -2.7895,  0.8544, -2.2223,  1.1878, -2.5469),
  sigma   = c( 1.4456,  0.9685,  1.3538,  1.4453,  1.2412,  1.1782),
  theta   = c( 0.0000,  0.0049,  0.0000,  0.0043,  0.0000,  0.0059),
  gmc     = c( 2.02,    0.0664,  2.35,    0.1127,   3.28,   0.0843)
)

# ------------------------------------------------------------
# Core functions
# ------------------------------------------------------------

#' Sample from a shifted lognormal distribution
#'
#' @param n       Number of samples
#' @param mu      Log-scale mean (natural log)
#' @param sigma   Log-scale standard deviation
#' @param theta   Shift parameter (minimum value; set 0 for standard lognormal)
#' @return Numeric vector of length n (antibody concentrations in ug/mL)
rshiftlnorm <- function(n, mu, sigma, theta = 0) {
  theta + rlnorm(n, meanlog = mu, sdlog = sigma)
}

#' CDF of the shifted lognormal
pshiftlnorm <- function(q, mu, sigma, theta = 0) {
  plnorm(q - theta, meanlog = mu, sdlog = sigma)
}

#' Quantile function of the shifted lognormal
qshiftlnorm <- function(p, mu, sigma, theta = 0) {
  theta + qlnorm(p, meanlog = mu, sdlog = sigma)
}

#' Density of the shifted lognormal
dshiftlnorm <- function(x, mu, sigma, theta = 0) {
  dlnorm(x - theta, meanlog = mu, sdlog = sigma)
}

# ------------------------------------------------------------
# Convenience wrapper: sample by trial and group name
# ------------------------------------------------------------

#' Sample antibody concentrations for a given trial and group
#'
#' @param n      Number of samples
#' @param trial  One of "nckp", "amer_indian", "s_africa"
#' @param group  One of "Immunized", "Un-immunized"
#' @return Numeric vector of antibody concentrations (ug/mL)
sample_antibody <- function(n, trial, group) {
  row <- params[params$trial == trial & params$group == group, ]
  if (nrow(row) == 0) stop("Unknown trial/group combination")
  rshiftlnorm(n, mu = row$mu, sigma = row$sigma, theta = row$theta)
}

# ------------------------------------------------------------
# Examples
# ------------------------------------------------------------

set.seed(42)

# Sample 1000 subjects from each group
samples <- lapply(seq_len(nrow(params)), function(i) {
  rshiftlnorm(1000, mu = params$mu[i], sigma = params$sigma[i], theta = params$theta[i])
})
names(samples) <- paste(params$trial, params$group, sep = " - ")

# Check geometric means recover published values
cat("Geometric mean check (sample vs published):\n")
for (nm in names(samples)) {
  gmc_sample   <- exp(mean(log(samples[[nm]])))
  gmc_published <- params$gmc[paste(params$trial, params$group, sep = " - ") == nm]
  cat(sprintf("  %-40s  sample GMC = %.3f  |  published GMC = %.3f\n",
              nm, gmc_sample, gmc_published))
}

# Proportion above a protective concentration threshold
cprot <- c(nckp = 0.20, `amer_indian` = 0.99, `s_africa` = 0.68)

cat("\nProportion of immunized subjects above [C]prot:\n")
for (trial in unique(params$trial)) {
  s    <- sample_antibody(100000, trial, "Immunized")
  prot <- mean(s > cprot[trial])
  row  <- params[params$trial == trial & params$group == "Immunized", ]
  # Exact value from CDF
  prot_exact <- 1 - pshiftlnorm(cprot[trial], row$mu, row$sigma, row$theta)
  cat(sprintf("  %-20s  [C]prot = %.2f  simulated = %.1f%%  exact = %.1f%%\n",
              trial, cprot[trial], prot * 100, prot_exact * 100))
}

# Quick plot (base R)
if (requireNamespace("graphics", quietly = TRUE)) {
  op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  trials <- unique(params$trial)
  cols   <- c(Immunized = "#1a3a6b", `Un-immunized` = "#1a3a6b")
  lty    <- c(Immunized = 1, `Un-immunized` = 2)

  for (trial in trials) {
    x_seq <- exp(seq(log(0.01), log(100), length.out = 500))
    plot(NA, xlim = c(0.01, 100), ylim = c(0, 1),
         log = "x", xlab = "Antibody concentration (ug/mL)",
         ylab = "Proportion above concentration", main = trial)
    abline(v = cprot[trial], col = "grey60", lty = 3)
    for (grp in c("Immunized", "Un-immunized")) {
      row <- params[params$trial == trial & params$group == grp, ]
      y   <- 1 - pshiftlnorm(x_seq, row$mu, row$sigma, row$theta)
      lines(x_seq, y, col = cols[grp], lty = lty[grp], lwd = 2)
    }
    legend("topleft", legend = c("Immunized", "Un-immunized"),
           col = "#1a3a6b", lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  }
  par(op)
}


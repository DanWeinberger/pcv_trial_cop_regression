#nct00366340 serotype 6A

#From Siber:
#NCKP: VE 97.4% (82.7, 99.9)  2.02 (1.90, 2.15) (N=10,940); 0.05 (0.05, 0.06) (N=10,995)
#AI: VE 76.8% (−9.4, 95.1)  2.35 (2.24, 2.46) (N=2974); 0.08 (0.07, 0.08) (N=2818)
#S. Africa VE:90% (29.7, 99.8) 3.28 (3.10, 3.46) (N=18,557);0.12 (0.11, 0.13) (N=18,550)

#From Claude, estimate parameters for a shifted log-normal distribution
"trial,group,mu_nat_log,sigma,theta_ug_per_mL,geometric_mean_ug_per_mL
NCKP,Immunized,-0.994,1.108,0.00962,0.370
NCKP,Un-immunized,-2.795,1.104,0.00373,0.0614
American Indian,Immunized,-0.376,1.304,0.01046,0.686
American Indian,Un-immunized,-2.303,1.382,0.00492,0.100
South African,Immunized,-0.259,1.175,0.000663,0.772
South African,Un-immunized,-2.891,1.232,0.00466,0.0553"


extract_gmc_sd <- function(gmc, gmc_lcl, gmc_ucl, N){
  log_gmc = log
  
  
}

# =====================================================================
# 1. GENERATE DATA
# =====================================================================
set.seed(42)
n_per_group <- 10000

# Assuming log IgG levels are higher in the vaccinated group
conc <- data.frame(
  log_igg = c(rnorm(n_per_group, mean = -1.47, sd = 1.12),   # Unvax
              rnorm(n_per_group, mean = 0.285, sd = 0.939)),  # Vax
  grp = c(rep(1, n_per_group), rep(2, n_per_group))   # 1=Unvax, 2=Vax
)

# =====================================================================
# 2. DEFINE HELPER FUNCTIONS
# =====================================================================
ilogit <- function(x) 1 / (1 + exp(-x))

# Calculates predicted VE given a set of parameters
calc_pred_ve <- function(b0, b1, data) {
  prob_ipd <- ilogit(b0 + b1 * data$log_igg)
  
  risk_unvax <- mean(prob_ipd[data$grp == 1], na.rm = TRUE)
  risk_vax <- mean(prob_ipd[data$grp == 2], na.rm = TRUE)
  
  # Protect against division by zero in extreme parameter proposals
  if (risk_unvax == 0) return(-Inf) 
  
  # Standard VE formula: 1 - (Risk_Vax / Risk_Unvax)
  return(1 - (risk_vax / risk_unvax))
}

# Calculates the Log-Posterior (Likelihood + Prior)
log_posterior <- function(b0, b1, data, target_log_rr, target_rr_sd) {
  
  # 1. Log-Likelihood: How well does pred_ve match the trial target_ve?
  pred_ve <- calc_pred_ve(b0, b1, data)
  if (is.infinite(pred_ve)) return(-Inf)
  
  log_pred_rr = log(1-pred_ve)
  ll <- dnorm(target_log_rr, mean = log_pred_rr, sd = target_rr_sd, log = TRUE)
  
  # 2. Log-Priors
  # Prior for baseline risk: tightly clustered around log(0.01)
  prior_b0 <- dnorm(b0, mean = log(0.01), sd = 0.5, log = TRUE)
  # Prior for slope: weakly informative
  prior_b1 <- dnorm(b1, mean = 0, sd = 5, log = TRUE) 
  
  return(ll + prior_b0 + prior_b1)
}

# =====================================================================
# 3. RUN METROPOLIS-HASTINGS MCMC
# =====================================================================
# Settings
iterations <- 5000
burn_in <- 1000

# Vectors to store results
chain_b0 <- numeric(iterations)
chain_b1 <- numeric(iterations)

# Trial constraints (e.g., VE = 80%, standard error = 0.05)
#andrewa: VE: 98%, 64-99.8 (logRR= -3.91, SD = 1.3)
target_log_rr <- -3.91
target_rr_sd <- 1.3

# Initial values
chain_b0[1] <- log(0.01)
chain_b1[1] <- -0.5 

current_log_post <- log_posterior(chain_b0[1], chain_b1[1], conc, target_log_rr, target_rr_sd)

# MCMC Loop
for (i in 2:iterations) {
  # Propose new parameters (random walk with small step size)
  prop_b0 <- rnorm(1, mean = chain_b0[i-1], sd = 0.1)
  prop_b1 <- rnorm(1, mean = chain_b1[i-1], sd = 0.1)
  
  prop_log_post <- log_posterior(prop_b0, prop_b1, conc, target_log_rr, target_rr_sd)
  
  # Calculate acceptance probability
  # (Using log scale to prevent numerical underflow)
  acceptance_ratio <- exp(prop_log_post - current_log_post)
  
  # Accept or reject
  if (runif(1) < acceptance_ratio) {
    chain_b0[i] <- prop_b0
    chain_b1[i] <- prop_b1
    current_log_post <- prop_log_post
  } else {
    chain_b0[i] <- chain_b0[i-1]
    chain_b1[i] <- chain_b1[i-1]
  }
}

# =====================================================================
# 4. SUMMARIZE RESULTS
# =====================================================================
# Remove burn-in phase
post_b0 <- chain_b0[(burn_in+1):iterations]
post_b1 <- chain_b1[(burn_in+1):iterations]

cat("Median Baseline Log-Odds (b0):", median(post_b0), 
    "\n95% CrI:", quantile(post_b0, c(0.025, 0.975)), "\n\n")

cat("Median Slope (b1):", median(post_b1), 
    "\n95% CrI:", quantile(post_b1, c(0.025, 0.975)), "\n")

#################################################################
#PLOT OUTPUT

library(tidyverse)
# 1. Sample 300 random iterations from your posterior matrix
set.seed(42)
sample_indices <- sample(1:nrow(pred_matrix), 1000)

# 2. Reshape those 300 draws into a "long" format for ggplot
pred_long <- data.frame(
  draw = rep(1:1000, each = length(x_seq)),
  log_igg = rep(x_seq, times = 1000),
  # t() transposes the matrix so as.vector reads it curve-by-curve
  risk = as.vector(t(pred_matrix[sample_indices, ])) 
)

# 3. Plot the individual draws
ggplot() +
  # Plot 300 highly transparent lines
  # The overlapping alpha creates the visual density gradient
  geom_line(data = pred_long, 
            aes(x = log_igg, y = risk, group = draw), 
            color = "#2c7fb8", alpha = 0.03) + 
  # Overlay the median line in solid dark blue
  geom_line(data = pred_summary, 
            aes(x = log_igg, y = median_risk), 
            color = "#081d58", linewidth = 1) +
  geom_rug(data = conc, 
           aes(x = log_igg, color = factor(grp)), 
           alpha = 0.1, sides = "b") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_color_manual(values = c("1" = "darkred", "2" = "darkgreen"), 
                     labels = c("1" = "Unvaccinated", "2" = "Vaccinated"),
                     name = "Group") +
  labs(
    title = "Posterior Dose-Response Curve for IPD Risk",
    subtitle = "Darker bands indicate higher posterior probability density",
    x = "Log IgG Concentration",
    y = "Predicted Probability of IPD"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


plot(post_b0, post_b1, pch = ".", 
     xlab = "b0", ylab = "b1",
     main = "Joint posterior — is it a ridge?")

contraction_b0 <- 1 - var(post_b0) / var(prior_draws_b0)
contraction_b1 <- 1 - var(post_b1) / var(prior_draws_b1)



##DIAGNOSTICS

# =====================================================================
# DIAGNOSTICS: Is the model learning from data, or reflecting priors?
# =====================================================================
# Run this AFTER your main MCMC. Assumes the following objects exist:
#   - conc           (data frame with log_igg, grp)
#   - post_b0, post_b1   (post-burn-in posterior draws)
#   - target_log_rr, target_rr_sd
#   - calc_pred_ve(), ilogit(), log_posterior()  (from your main script)

library(tidyverse)

# Reusable MCMC function so we can swap in different log-posteriors
run_mh <- function(log_post_fn, n_iter = 5000, burn = 1000,
                   init_b0 = log(0.01), init_b1 = -0.5,
                   step_b0 = 0.1, step_b1 = 0.1) {
  b0 <- numeric(n_iter); b1 <- numeric(n_iter)
  b0[1] <- init_b0; b1[1] <- init_b1
  curr <- log_post_fn(b0[1], b1[1])
  for (i in 2:n_iter) {
    p0 <- rnorm(1, b0[i-1], step_b0)
    p1 <- rnorm(1, b1[i-1], step_b1)
    prop <- log_post_fn(p0, p1)
    if (runif(1) < exp(prop - curr)) {
      b0[i] <- p0; b1[i] <- p1; curr <- prop
    } else {
      b0[i] <- b0[i-1]; b1[i] <- b1[i-1]
    }
  }
  list(b0 = b0[(burn+1):n_iter], b1 = b1[(burn+1):n_iter])
}


# =====================================================================
# DIAGNOSTIC 1: PRIOR-ONLY SAMPLING
# =====================================================================
# Sample from the prior alone, then compare prior vs posterior on
# (a) the parameters and (b) the implied VE and dose-response curve.

log_prior_only <- function(b0, b1) {
  dnorm(b0, mean = log(0.01), sd = 0.5, log = TRUE) +
    dnorm(b1, mean = 0,         sd = 5,   log = TRUE)
}

set.seed(123)
prior_chain <- run_mh(log_prior_only, n_iter = 5000, burn = 1000,
                      step_b0 = 0.4, step_b1 = 2.0)  # larger steps for diffuse prior

prior_b0 <- prior_chain$b0
prior_b1 <- prior_chain$b1

# Compare parameter marginals
par(mfrow = c(1, 2))
plot(density(prior_b0), col = "grey50", lwd = 2, main = "b0: prior vs posterior",
     xlab = "b0", ylim = range(density(prior_b0)$y, density(post_b0)$y))
lines(density(post_b0), col = "firebrick", lwd = 2)
legend("topright", c("prior", "posterior"), col = c("grey50", "firebrick"), lwd = 2, bty = "n")

plot(density(prior_b1), col = "grey50", lwd = 2, main = "b1: prior vs posterior",
     xlab = "b1", ylim = range(density(prior_b1)$y, density(post_b1)$y))
lines(density(post_b1), col = "firebrick", lwd = 2)
legend("topright", c("prior", "posterior"), col = c("grey50", "firebrick"), lwd = 2, bty = "n")
par(mfrow = c(1, 1))

# Contraction: 1 - var(posterior)/var(prior). Near 0 = no learning, near 1 = strong.
contraction_b0 <- 1 - var(post_b0) / var(prior_b0)
contraction_b1 <- 1 - var(post_b1) / var(prior_b1)
cat("\nPrior-to-posterior contraction:\n")
cat("  b0:", round(contraction_b0, 3), "\n")
cat("  b1:", round(contraction_b1, 3), "\n")

# Compare implied VE under prior vs posterior
prior_ve <- mapply(calc_pred_ve, prior_b0, prior_b1,
                   MoreArgs = list(data = conc))
post_ve  <- mapply(calc_pred_ve, post_b0,  post_b1,
                   MoreArgs = list(data = conc))

# Drop pathological draws (VE can be < 0 or NaN under diffuse prior)
prior_ve_clean <- prior_ve[is.finite(prior_ve) & prior_ve > -2 & prior_ve < 1]
post_ve_clean  <- post_ve[is.finite(post_ve)]

plot(density(prior_ve_clean), col = "grey50", lwd = 2,
     main = "Implied VE: prior vs posterior",
     xlab = "VE", xlim = c(-1, 1))
lines(density(post_ve_clean), col = "firebrick", lwd = 2)
abline(v = 1 - exp(target_log_rr), col = "blue", lty = 2)
legend("topleft", c("prior", "posterior", "trial target"),
       col = c("grey50", "firebrick", "blue"), lwd = 2, lty = c(1,1,2), bty = "n")


# =====================================================================
# DIAGNOSTIC 2: JOINT POSTERIOR — IDENTIFIABILITY RIDGE
# =====================================================================
plot(post_b0, post_b1, pch = 16, cex = 0.3, col = rgb(0.7, 0.1, 0.1, 0.2),
     xlab = "b0", ylab = "b1",
     main = "Joint posterior — ridge indicates only the combination is identified")
points(prior_b0, prior_b1, pch = 16, cex = 0.3, col = rgb(0.3, 0.3, 0.3, 0.05))
legend("topright", c("posterior", "prior"),
       col = c("firebrick", "grey40"), pch = 16, bty = "n")

cor_b0_b1 <- cor(post_b0, post_b1)
cat("\nPosterior correlation between b0 and b1:", round(cor_b0_b1, 3), "\n")
cat("  (|cor| close to 1 means a ridge — only the linear combination is identified)\n")


# =====================================================================
# DIAGNOSTIC 3: PRIOR SENSITIVITY ANALYSIS
# =====================================================================
# Re-fit under several alternative priors and compare implied curves.

prior_specs <- list(
  baseline   = list(b0_mu = log(0.01), b0_sd = 0.5, b1_mu = 0, b1_sd = 5),
  tight_b1   = list(b0_mu = log(0.01), b0_sd = 0.5, b1_mu = 0, b1_sd = 1),
  loose_b1   = list(b0_mu = log(0.01), b0_sd = 0.5, b1_mu = 0, b1_sd = 20),
  high_base  = list(b0_mu = log(0.02), b0_sd = 0.5, b1_mu = 0, b1_sd = 5),
  low_base   = list(b0_mu = log(0.005),b0_sd = 0.5, b1_mu = 0, b1_sd = 5),
  neg_slope  = list(b0_mu = log(0.01), b0_sd = 0.5, b1_mu = -1, b1_sd = 5)
)

make_logpost <- function(spec) {
  function(b0, b1) {
    pred_ve <- calc_pred_ve(b0, b1, conc)
    if (is.infinite(pred_ve) || pred_ve >= 1) return(-Inf)
    ll <- dnorm(target_log_rr, mean = log(1 - pred_ve), sd = target_rr_sd, log = TRUE)
    ll +
      dnorm(b0, spec$b0_mu, spec$b0_sd, log = TRUE) +
      dnorm(b1, spec$b1_mu, spec$b1_sd, log = TRUE)
  }
}

set.seed(456)
sensitivity_results <- lapply(prior_specs, function(spec) {
  run_mh(make_logpost(spec), n_iter = 5000, burn = 1000)
})

# Summary table
sens_summary <- map_dfr(names(sensitivity_results), function(nm) {
  ch <- sensitivity_results[[nm]]
  ve <- mapply(calc_pred_ve, ch$b0, ch$b1, MoreArgs = list(data = conc))
  ve <- ve[is.finite(ve)]
  tibble(
    prior       = nm,
    b0_median   = median(ch$b0),
    b0_lo       = quantile(ch$b0, 0.025),
    b0_hi       = quantile(ch$b0, 0.975),
    b1_median   = median(ch$b1),
    b1_lo       = quantile(ch$b1, 0.025),
    b1_hi       = quantile(ch$b1, 0.975),
    ve_median   = median(ve),
    ve_lo       = quantile(ve, 0.025),
    ve_hi       = quantile(ve, 0.975)
  )
})
print(sens_summary)

# Plot dose-response curves under each prior
x_seq <- seq(min(conc$log_igg), max(conc$log_igg), length.out = 100)

curves <- map_dfr(names(sensitivity_results), function(nm) {
  ch <- sensitivity_results[[nm]]
  med_b0 <- median(ch$b0); med_b1 <- median(ch$b1)
  # Also compute 95% pointwise band
  risk_mat <- sapply(seq_along(ch$b0), function(i) ilogit(ch$b0[i] + ch$b1[i] * x_seq))
  tibble(
    prior   = nm,
    log_igg = x_seq,
    median  = apply(risk_mat, 1, median),
    lo      = apply(risk_mat, 1, quantile, 0.025),
    hi      = apply(risk_mat, 1, quantile, 0.975)
  )
})

ggplot(curves, aes(x = log_igg, y = median, color = prior, fill = prior)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Dose-response curve under different priors",
       subtitle = "Curves that move with the prior indicate prior-driven inference",
       x = "Log IgG", y = "Predicted IPD risk") +
  theme_minimal()


# =====================================================================
# DIAGNOSTIC 4: POSTERIOR PREDICTIVE CHECK ON THE TARGET
# =====================================================================
# If the likelihood is binding, the posterior predictive log RR should
# concentrate near target_log_rr = -3.91.

pred_log_rr <- log(1 - post_ve_clean)
hist(pred_log_rr, breaks = 40, col = "lightsteelblue", border = "white",
     main = "Posterior predictive log RR",
     xlab = "log RR")
abline(v = target_log_rr, col = "red", lwd = 2)
abline(v = target_log_rr + c(-1.96, 1.96) * target_rr_sd, col = "red", lty = 2)
legend("topright", c("trial target", "trial 95% CI"),
       col = "red", lty = c(1, 2), lwd = 2, bty = "n")

cat("\nPosterior predictive log RR:\n")
cat("  median:", round(median(pred_log_rr), 2),
    "  (trial target:", target_log_rr, ")\n")
cat("  95% CrI:", round(quantile(pred_log_rr, c(0.025, 0.975)), 2), "\n")


# =====================================================================
# UPDATED MODEL: Add baseline rate constraint
# =====================================================================
# Adds a second likelihood term: observed IPD rate in unvaccinated arm
# = 165 per 100,000 person-years.
#
# This makes b0 data-identified rather than prior-anchored.
# =====================================================================

library(tidyverse)

# --- Inputs you can tune ---
baseline_rate <- 165 / 100000   # 165 per 100,000 person-years
# How much person-time backs this estimate? This controls how tight the
# constraint is. Set based on the actual study/surveillance source:
#   small trial (~30k PY total)  -> few cases, weak constraint
#   meta-analysis (~500k PY)     -> moderate
#   national surveillance (5M+)  -> very tight
person_years_unvax <- 100000   # <-- adjust to reflect data source

observed_cases <- round(baseline_rate * person_years_unvax)
cat("Implied observed cases:", observed_cases,
    "over", person_years_unvax, "person-years\n")
cat("Implied SE on log rate:", round(1/sqrt(observed_cases), 3), "\n\n")

# --- Helper functions (re-defined for clarity) ---
ilogit <- function(x) 1 / (1 + exp(-x))

calc_group_risks <- function(b0, b1, data) {
  prob_ipd <- ilogit(b0 + b1 * data$log_igg)
  list(
    unvax = mean(prob_ipd[data$grp == 1], na.rm = TRUE),
    vax   = mean(prob_ipd[data$grp == 2], na.rm = TRUE)
  )
}

# --- New log-posterior with TWO likelihood terms ---
log_posterior_v2 <- function(b0, b1, data,
                             target_log_rr, target_rr_sd,
                             obs_cases, person_years) {
  
  risks <- calc_group_risks(b0, b1, data)
  if (risks$unvax <= 0 || risks$unvax >= 1) return(-Inf)
  if (risks$vax  <  0 || risks$vax  >= 1) return(-Inf)
  
  # Likelihood 1: log RR / VE (as before)
  pred_log_rr <- log(risks$vax / risks$unvax)
  ll_ve <- dnorm(target_log_rr, mean = pred_log_rr,
                 sd = target_rr_sd, log = TRUE)
  
  # Likelihood 2: observed cases in unvax arm
  # cases ~ Poisson(risk_unvax_pred * person_years)
  expected_cases <- risks$unvax * person_years
  ll_baseline <- dpois(obs_cases, lambda = expected_cases, log = TRUE)
  
  # Priors (b0 can be looser now since data identifies it)
  prior_b0 <- dnorm(b0, mean = log(0.01), sd = 2, log = TRUE)
  prior_b1 <- dnorm(b1, mean = 0,         sd = 5, log = TRUE)
  
  ll_ve + ll_baseline + prior_b0 + prior_b1
}

# --- Trial constraints (unchanged) ---
target_log_rr <- -3.91
target_rr_sd  <- 1.3

# --- Reusable MCMC ---
run_mh <- function(log_post_fn, n_iter = 10000, burn = 2000,
                   init_b0 = log(0.01), init_b1 = -0.5,
                   step_b0 = 0.1, step_b1 = 0.1) {
  b0 <- numeric(n_iter); b1 <- numeric(n_iter)
  b0[1] <- init_b0; b1[1] <- init_b1
  curr <- log_post_fn(b0[1], b1[1])
  accept <- 0
  for (i in 2:n_iter) {
    p0 <- rnorm(1, b0[i-1], step_b0)
    p1 <- rnorm(1, b1[i-1], step_b1)
    prop <- log_post_fn(p0, p1)
    if (runif(1) < exp(prop - curr)) {
      b0[i] <- p0; b1[i] <- p1; curr <- prop
      accept <- accept + 1
    } else {
      b0[i] <- b0[i-1]; b1[i] <- b1[i-1]
    }
  }
  cat("Acceptance rate:", round(accept / n_iter, 2), "\n")
  list(b0 = b0[(burn+1):n_iter], b1 = b1[(burn+1):n_iter])
}

# --- Run the new model ---
set.seed(42)
fit_v2 <- run_mh(
  function(b0, b1) log_posterior_v2(b0, b1, conc,
                                    target_log_rr, target_rr_sd,
                                    observed_cases, person_years_unvax),
  n_iter = 20000, burn = 5000,
  step_b0 = 0.05, step_b1 = 0.15   # smaller b0 step now that it's constrained
)

post_b0_v2 <- fit_v2$b0
post_b1_v2 <- fit_v2$b1

# =====================================================================
# COMPARE: Old model (post_b0, post_b1) vs New model (post_b0_v2, post_b1_v2)
# =====================================================================

cat("\n=== b0 (baseline log-odds) ===\n")
cat("Old:  median", round(median(post_b0), 2),
    "  95% CrI", round(quantile(post_b0, c(0.025, 0.975)), 2), "\n")
cat("New:  median", round(median(post_b0_v2), 2),
    "  95% CrI", round(quantile(post_b0_v2, c(0.025, 0.975)), 2), "\n")

cat("\n=== b1 (slope on log_igg) ===\n")
cat("Old:  median", round(median(post_b1), 2),
    "  95% CrI", round(quantile(post_b1, c(0.025, 0.975)), 2), "\n")
cat("New:  median", round(median(post_b1_v2), 2),
    "  95% CrI", round(quantile(post_b1_v2, c(0.025, 0.975)), 2), "\n")

cat("\n=== Correlation b0-b1 ===\n")
cat("Old:", round(cor(post_b0, post_b1), 3), "\n")
cat("New:", round(cor(post_b0_v2, post_b1_v2), 3),
    "  <- should now show a ridge\n")

cat("\n=== Implied baseline risk in unvax ===\n")
old_baseline <- mapply(function(b0, b1) calc_group_risks(b0, b1, conc)$unvax,
                       post_b0, post_b1)
new_baseline <- mapply(function(b0, b1) calc_group_risks(b0, b1, conc)$unvax,
                       post_b0_v2, post_b1_v2)
cat("Old:  median", signif(median(old_baseline), 3),
    "  95% CrI", signif(quantile(old_baseline, c(0.025, 0.975)), 3), "\n")
cat("New:  median", signif(median(new_baseline), 3),
    "  95% CrI", signif(quantile(new_baseline, c(0.025, 0.975)), 3), "\n")
cat("Target: ", baseline_rate, "\n")

# --- Joint posterior plot ---
par(mfrow = c(1, 2))
plot(post_b0, post_b1, pch = ".", col = rgb(0.5, 0.5, 0.5, 0.3),
     xlab = "b0", ylab = "b1",
     main = "Old: no ridge (b0 prior-driven)",
     xlim = range(c(post_b0, post_b0_v2)),
     ylim = range(c(post_b1, post_b1_v2)))
plot(post_b0_v2, post_b1_v2, pch = ".", col = rgb(0.7, 0.1, 0.1, 0.3),
     xlab = "b0", ylab = "b1",
     main = "New: ridge (b0 + b1 jointly identified)",
     xlim = range(c(post_b0, post_b0_v2)),
     ylim = range(c(post_b1, post_b1_v2)))
par(mfrow = c(1, 1))

# --- Dose-response curve comparison ---
x_seq <- seq(min(conc$log_igg), max(conc$log_igg), length.out = 100)

curve_summary <- function(b0_draws, b1_draws, label) {
  risk_mat <- sapply(seq_along(b0_draws),
                     function(i) ilogit(b0_draws[i] + b1_draws[i] * x_seq))
  tibble(
    model   = label,
    log_igg = x_seq,
    median  = apply(risk_mat, 1, median),
    lo      = apply(risk_mat, 1, quantile, 0.025),
    hi      = apply(risk_mat, 1, quantile, 0.975)
  )
}

curves_compare <- bind_rows(
  curve_summary(post_b0,    post_b1,    "Original (VE only)"),
  curve_summary(post_b0_v2, post_b1_v2, "Updated (VE + baseline rate)")
)

ggplot(curves_compare, aes(x = log_igg, y = median,
                           color = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_rug(data = conc, aes(x = log_igg, color = NULL, fill = NULL,
                            linetype = factor(grp)),
           alpha = 0.05, sides = "b", inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Dose-response: with vs. without baseline rate constraint",
       subtitle = paste0("Baseline rate target: ", baseline_rate,
                         " per person-year (",
                         format(person_years_unvax, big.mark = ","),
                         " person-years of data)"),
       x = "Log IgG", y = "Predicted IPD risk") +
  theme_minimal() +
  theme(legend.position = "bottom")


##RELATIVE RISK PLOT

# --- Choose reference point ---
# Mean log IgG in unvaccinated group = baseline for "what does vaccination buy you?"
ref_log_igg <- mean(conc$log_igg[conc$grp == 1])
cat("Reference log IgG (mean unvax):", round(ref_log_igg, 3), "\n")

# Grid of log IgG values to evaluate
x_seq <- seq(min(conc$log_igg), max(conc$log_igg), length.out = 100)

# --- Compute relative risk for each posterior draw ---
# For each draw i: RR(x) = ilogit(b0 + b1*x) / ilogit(b0 + b1*ref)
# Note: b0 cancels out in the rare-disease limit, so this is almost
# entirely a function of b1. That's a feature - the curve shape is
# what b1 identifies, and that's what we're plotting.

rr_matrix <- sapply(seq_along(post_b0_v2), function(i) {
  risk_x   <- ilogit(post_b0_v2[i] + post_b1_v2[i] * x_seq)
  risk_ref <- ilogit(post_b0_v2[i] + post_b1_v2[i] * ref_log_igg)
  risk_x / risk_ref
})

rr_summary <- tibble(
  log_igg = x_seq,
  median  = apply(rr_matrix, 1, median),
  lo50    = apply(rr_matrix, 1, quantile, 0.25),
  hi50    = apply(rr_matrix, 1, quantile, 0.75),
  lo95    = apply(rr_matrix, 1, quantile, 0.025),
  hi95    = apply(rr_matrix, 1, quantile, 0.975)
)

# --- Plot ---
# Sample draws for spaghetti
set.seed(1)
draw_idx <- sample(seq_along(post_b0_v2), 300)
rr_draws <- tibble(
  draw    = rep(draw_idx, each = length(x_seq)),
  log_igg = rep(x_seq, times = length(draw_idx)),
  rr      = as.vector(rr_matrix[, draw_idx])
)

p_rr <- ggplot() +
  geom_line(data = rr_draws,
            aes(x = log_igg, y = rr, group = draw),
            color = "#2c7fb8", alpha = 0.03) +
  geom_line(data = rr_summary,
            aes(x = log_igg, y = median),
            color = "#081d58", linewidth = 1) +
  geom_ribbon(data = rr_summary,
              aes(x = log_igg, ymin = lo95, ymax = hi95),
              fill = NA, color = "#081d58", linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = 1, color = "grey50", linetype = "dotted") +
  geom_vline(xintercept = ref_log_igg, color = "darkred",
             linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = ref_log_igg, y = max(rr_summary$hi95) * 0.95,
           label = "reference\n(mean unvax)", hjust = -0.1,
           size = 3, color = "darkred") +
  geom_rug(data = conc,
           aes(x = log_igg, color = factor(grp)),
           alpha = 0.05, sides = "b", inherit.aes = FALSE) +
  scale_color_manual(values = c("1" = "darkred", "2" = "darkgreen"),
                     labels = c("1" = "Unvaccinated", "2" = "Vaccinated"),
                     name = "Group") +
  scale_y_log10(labels = scales::percent_format(accuracy = 1),
                breaks = c(0.001, 0.01, 0.1, 0.5, 1, 2)) +
  labs(
    title = "Relative IPD risk by log IgG",
    subtitle = paste0("Risk relative to mean unvaccinated IgG (log_igg = ",
                      round(ref_log_igg, 2), "); dashed lines = 95% CrI"),
    x = "Log IgG concentration",
    y = "Relative risk (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_rr)

# --- Summary table at key IgG points ---
key_points <- c(
  "very low (-4)"     = -4,
  "low (-2)"          = -2,
  "reference (unvax mean)" = ref_log_igg,
  "vax mean (0.29)"   = mean(conc$log_igg[conc$grp == 2]),
  "high (2)"          = 2
)

rr_table <- map_dfr(names(key_points), function(label) {
  x <- key_points[[label]]
  rr_at_x <- mapply(function(b0, b1) {
    ilogit(b0 + b1 * x) / ilogit(b0 + b1 * ref_log_igg)
  }, post_b0_v2, post_b1_v2)
  tibble(
    point        = label,
    log_igg      = round(x, 2),
    rr_median    = signif(median(rr_at_x), 3),
    rr_lo        = signif(quantile(rr_at_x, 0.025), 3),
    rr_hi        = signif(quantile(rr_at_x, 0.975), 3),
    pct_reduction = paste0(round(100 * (1 - median(rr_at_x)), 1), "%")
  )
})

cat("\n=== Relative risk at key IgG levels ===\n")
print(rr_table)

# --- Bonus: implied VE if everyone shifted from unvax IgG dist to vax IgG dist ---
# This is essentially recapitulating the VE input - a sanity check
implied_ve_draws <- mapply(function(b0, b1) {
  risks <- calc_group_risks(b0, b1, conc)
  1 - risks$vax / risks$unvax
}, post_b0_v2, post_b1_v2)

cat("\n=== Implied VE under updated model (sanity check, target was 0.98) ===\n")
cat("Median:", round(median(implied_ve_draws), 3),
    "  95% CrI:", round(quantile(implied_ve_draws, c(0.025, 0.975)), 3), "\n")


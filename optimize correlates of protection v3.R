library(tidyverse)

#nct00366340 serotype 6A

#From Siber:
#NCKP: VE 97.4% (82.7, 99.9)  2.02 (1.90, 2.15) (N=10,940); 0.05 (0.05, 0.06) (N=10,995)
#AI: VE 76.8% (−9.4, 95.1)  2.35 (2.24, 2.46) (N=2974); 0.08 (0.07, 0.08) (N=2818)
#S. Africa VE:90% (29.7, 99.8) 3.28 (3.10, 3.46) (N=18,557);0.12 (0.11, 0.13) (N=18,550)


#Use claude to try to estimate the SD for the log_GMC; the values presente din Siber can't be used to back this out
obs_values <- list(
'nckp' = c(
  log_gmc_vax = log(2.02),
  log_gmc_vax_SD = 1.446,
  log_gmc_unvax = log(0.05),
  log_gmc_unvax_SD = 0.969,
  obs_cases_unvax = 39,
  obs_cases_vax = 1,
  person_years_vax = 10940,
  person_years_unvax = 10995
),
'amer_indian' = c(
  log_gmc_vax = log(2.35),
  log_gmc_vax_SD = 1.354,
  log_gmc_unvax = log(0.08),
  log_gmc_unvax_SD = 1.445,
  obs_cases_unvax = 8,
  obs_cases_vax = 2,
  person_years_vax = 2974,
  person_years_unvax = 2818
),
's_africa' = c(
  log_gmc_vax = log(3.28),
  log_gmc_vax_SD = 1.241,
  log_gmc_unvax = log(0.12),
  log_gmc_unvax_SD = 1.178,
    obs_cases_unvax =10 ,
  obs_cases_vax = 1,
  person_years_vax = 18557 ,
  person_years_unvax = 18550
)
)

# =====================================================================
# 1. GENERATE DATA
# =====================================================================
set.seed(42)

source('sample_shifted_lognormal_siber.R')

# =====================================================================
# 2. DEFINE HELPER FUNCTIONS
# =====================================================================
ilogit <- function(x) 1 / (1 + exp(-x))

calc_group_risks <- function(b0, b1, log_gmc_unvax, conc) {
  prob_ipd <- ilogit(b0  + b1 * (conc$log_igg-log_gmc_unvax))
  list(
    unvax = mean(prob_ipd[conc$grp == 1], na.rm = TRUE),
    vax   = mean(prob_ipd[conc$grp == 2], na.rm = TRUE)
  )
}

#####################################################

log_posterior_v3 <- function(b0, b1, log_gmc_unvax,
                             obs_cases_unvax, person_years_unvax,
                             obs_cases_vax,   person_years_vax,
                             conc) {

  risks <- calc_group_risks(b0, b1, log_gmc_unvax, conc)
  if (risks$unvax <= 0 || risks$unvax >= 1) return(-Inf)
  if (risks$vax   <= 0 || risks$vax   >= 1) return(-Inf)

  ll_unvax <- dpois(obs_cases_unvax, lambda = risks$unvax * person_years_unvax, log = TRUE)
  ll_vax   <- dpois(obs_cases_vax,   lambda = risks$vax   * person_years_vax,   log = TRUE)

  baseline_rate <- obs_cases_unvax/person_years_unvax

  prior_b0 <- dnorm(b0, mean = log(baseline_rate), sd = 0.1, log = TRUE)
  prior_b1 <- dnorm(b1, mean = 0, sd = 5, log = TRUE)

  ll_unvax + ll_vax + prior_b0 + prior_b1
}


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

# =====================================================================
# 3. ANALYSIS FUNCTION
# =====================================================================

run_population_analysis <- function(pop, sample_prefix = pop) {
  vals <- obs_values[[pop]]

  log_gmc_unvax      <- vals['log_gmc_unvax']
  log_gmc_unvax_sd   <- vals['log_gmc_unvax_SD']
  log_gmc_vax        <- vals['log_gmc_vax']
  log_gmc_vax_sd     <- vals['log_gmc_unvax_SD']
  obs_cases_unvax    <- vals['obs_cases_unvax']
  obs_cases_vax      <- vals['obs_cases_vax']
  person_years_vax   <- vals['person_years_vax']
  person_years_unvax <- vals['person_years_unvax']

  n_per_group <- length(samples[[paste0(sample_prefix, ' - Un-immunized')]])

  conc <- data.frame(
    log_igg = c(log(samples[[paste0(sample_prefix, ' - Un-immunized')]]),
                log(samples[[paste0(sample_prefix, ' - Immunized')]])),
    grp = c(rep(1, n_per_group), rep(2, n_per_group))
  ) %>%
    mutate(log_igg = if_else(log_igg < log(0.01), log(0.01), log_igg))

  print(ggplot(conc) +
    geom_histogram(aes(x=log_igg, group=grp, color=grp)) +
    labs(title=pop))

  print(conc %>%
    mutate(grp = factor(grp, labels = c("Unvaccinated", "Vaccinated"))) %>%
    arrange(log_igg) %>%
    group_by(grp) %>%
    mutate(prop_above = 1 - (rank(log_igg) - 1) / n()) %>%
    ggplot(aes(x = log_igg, y = prop_above, color = grp)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Unvaccinated" = "steelblue", "Vaccinated" = "tomato")) +
    labs(
      x = "Log IgG",
      y = "Proportion with concentration above x",
      color = "Group",
      title = paste("Reverse Cumulative Distribution of Log IgG by Vaccination Status -", pop)
    ) +
    theme_bw() +
    geom_vline(xintercept=log(0.68)))

  set.seed(42)
  fit <- run_mh(
    function(b0, b1) log_posterior_v3(b0, b1, log_gmc_unvax,
                                      obs_cases_unvax, person_years_unvax,
                                      obs_cases_vax,   person_years_vax,
                                      conc),
    n_iter = 20000, burn = 5000,
    step_b0 = 0.05, step_b1 = 0.15
  )

  plot(fit$b1, fit$b0, main=pop)
  plot(fit$b1, type='l', main=pop)
  hist(fit$b1, main=pop)

  print(quantile(fit$b1, probs=c(0.5, 0.025, 0.975)))

  betas <- cbind(fit$b0, fit$b1)

  log_igg_quants_unvax <- qnorm(seq(0.01, 0.99, by=0.01), mean=log_gmc_unvax, sd=log_gmc_unvax_sd)
  log_igg_quants_vax   <- qnorm(seq(0.01, 0.99, by=0.01), mean=log_gmc_vax,   sd=log_gmc_vax_sd)
  log_igg_quants <- sort(c(log_igg_quants_unvax, log_igg_quants_vax))
  log_igg_quants <- log_igg_quants[log_igg_quants >= log(0.01)]

  X <- cbind(rep(1, length(log_igg_quants)), (log_igg_quants - log_gmc_unvax))
  preds <- X %*% t(betas)

  matplot(X[,2] + log_gmc_unvax, ilogit(preds[,1:500]), type='l', col='gray', main=pop)
  abline(v=log(0.68), lwd=2)
  abline(v=log_gmc_unvax, col='red', lty=2)
  abline(v=log_gmc_vax, col='blue', lty=2)

  matplot(X[,2] + log_gmc_unvax, (preds[,1:500]), type='l', col='gray', main=pop)
  abline(v=log(0.68))
  abline(v=log_gmc_unvax, col='red')
  abline(v=log_gmc_vax, col='blue')

  invisible(list(fit=fit, conc=conc, betas=betas))
}

# =====================================================================
# 4. RUN FOR EACH POPULATION
# =====================================================================
#For SAfrica b1 = -0.54;  -0.91,-0.22 )
#For Navajo: b1 = -0.57, -1.14,-0.14
#For NCKP: b1 = -1.22, -1.72, -0.74

results_nckp        <- run_population_analysis('nckp')
results_amer_indian <- run_population_analysis('amer_indian')
results_s_africa    <- run_population_analysis('s_africa')


#In S Africa, 6B 7/25; 14 7/25; 23F 4/25, 19F 3/25
#In NCKP: 19F 13/49; 14 11/49; 18C 9/49; 23F 6/49; 6B 7/49
#Navaj0: 6B 4/8 ; 1 each of 9V, 14, 19F, 23F
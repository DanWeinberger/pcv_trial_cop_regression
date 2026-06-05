library(tidyverse)
library(boot)
################
simulate_gmc <- function(gmc,
                         ci_lower,
                         ci_upper,
                         n,
                         nsim = n,
                         seed = NULL) {
  
  # Optional reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Mean on log scale
  meanlog <- log(gmc)
  
  # Standard error from 95% CI on log scale
  se_log <- (log(ci_upper) - log(ci_lower)) / (2 * 1.96)
  
  # SD on log scale
  sdlog <- se_log * sqrt(n)
  
  # Generate samples
  x <- rlnorm(nsim,
              meanlog = meanlog,
              sdlog = sdlog)
  
  # Return useful outputs
  list(
    samples = x,
    parameters = list(
      meanlog = meanlog,
      sdlog = sdlog,
      gmc = gmc,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      n = n
    )
  )
}
#########################################
# Reverse cumulative distribution curve
# (a.k.a. survival curve / reverse ECDF)

plot_reverse_cdf <- function(x,
                             xlab = "Value",
                             ylab = "Proportion ≥ x",
                             main = "Reverse Cumulative Distribution") {
  
  # Sort values
  x_sorted <- sort(x)
  
  # Reverse cumulative probabilities
  rev_cdf <- rev(seq_along(x_sorted)) / length(x_sorted)
  
  # Plot
  plot(x_sorted,
       rev_cdf,
       type = "s",
       lwd = 2,
       xlab = xlab,
       ylab = ylab,
       main = main)
  
  invisible(data.frame(
    x = x_sorted,
    proportion_geq = rev_cdf
  ))
}
reverse_cdf_df <- function(x) {
  
  x_sorted <- sort(x)
  
  data.frame(
    x = x_sorted,
    proportion_geq =
      rev(seq_along(x_sorted)) / length(x_sorted)
  )
}

#####################################

#nct00366340 serotype 6A
sim1 <- simulate_gmc(gmc=1.33 , ci_lower=1.18 , ci_upper=1.49, n=249, nsim=1000000, seed=123)
sim2 <- simulate_gmc(gmc=0.23 , ci_lower=0.20 , ci_upper=0.26, n=279, nsim=1000000, seed=123)

# Plot reverse cumulative distribution
df1 <- reverse_cdf_df(sim1$samples) %>%
  mutate(sample=1)

df2 <- reverse_cdf_df(sim2$samples)%>%
  mutate(sample=2)

compare1 <- bind_rows(df1,df2) %>%
  mutate(log_prop)

ggplot(compare1) +
  geom_line(aes(x=x, y=proportion_geq, group=sample, color=sample))+
  theme_classic()

 df1 <- cbind.data.frame('igg' =sim1$samples ,  'sample'=1)
 df2 <- cbind.data.frame('igg' =sim2$samples ,  'sample'=2)
 
 log_mean_unvax = mean(log(df2$igg))
 
 #antibody concentration
 conc <- bind_rows(df1,df2) %>%
   mutate(log_igg = log(igg),
          sample = as.factor(sample),
          ipd =NA_real_,
          log_mean_unvax = log_mean_unvax ### #set this to mean of unvax groi[]
          )
 
 ggplot(conc) +
   geom_histogram(aes(x=log_igg, group=sample, color=sample))+
   theme_classic()
 
 
 ilogit  <- function(x){exp(x)/(1+exp(x))}
 
 #first, assume simple linear association
 
 lin_optim <- function(x, obs_ve, log_risk_unvax) {
   # Calculate probability of outcome for each individual
   conc$prob_ipd <- ilogit(log_risk_unvax + x[1] * conc$log_igg)
   
   # Calculate predicted VE based on the ratio of mean probabilities
   # (Note: See the structural note below regarding grp assignments)
   pred_ve <- 1 - mean(conc$prob_ipd[conc$grp == 1], na.rm = TRUE) / 
     mean(conc$prob_ipd[conc$grp == 2], na.rm = TRUE)
   
   # Return Squared Error so optim has a clear minimum at 0
   sq_err <- (obs_ve - pred_ve)^2
   return(sq_err)
 }

 # Run optim using Brent method for a single parameter, specifying a search interval
 res <- optim(par = c(-0.5), 
              fn = lin_optim, 
              obs_ve = 0.8, 
              log_risk_unvax = log(0.01),
              method = "Brent", 
              lower = -5, 
              upper = 0)
 
 print(res$par)  

  
  

 


# =============================================
# HPAI Ontario Wild Birds – Seasonal Migration
# Full model with estimated initial populations
# =============================================

# Clean workspace
rm(list = ls())

# Load libraries
library(deSolve)
library(mcmc)
library(coda)
library(ggmcmc)
library(dplyr)
library(lubridate)

# ---------------------------
# Read Ontario daily cases CSV
# ---------------------------
data <- read.csv("Ontario-wild-birds-waterfowl-daily_cases.csv",
                 header = TRUE,
                 stringsAsFactors = FALSE)

# Convert CollectionDate to Date
#data$CollectionDate <- mdy(data$CollectionDate)  # MM/DD/YYYY
data$CollectionDate <- as.Date(data$CollectionDate)

# Filter for Nov 2024 – Jul 2025
data_filtered <- data %>%
  filter(CollectionDate >= as.Date("2024-11-01"),
         CollectionDate <= as.Date("2025-03-25")) %>%
  arrange(CollectionDate) %>%
  mutate(t = as.integer(CollectionDate - min(CollectionDate)) + 1) # t = 1 for Nov 1
##################


data2 <- read.csv("Ontario-wild-birds-raptor-daily_cases.csv",
                  header = TRUE,
                  stringsAsFactors = FALSE)
data2$CollectionDate <- as.Date(data2$CollectionDate)

data2_filtered <- data2 %>%
  filter(CollectionDate >= as.Date("2024-11-01"),
         CollectionDate <= as.Date("2025-03-25")) %>%
  arrange(CollectionDate) %>%
  mutate(t = as.integer(CollectionDate - min(CollectionDate)) + 1) # t = 1 for Nov 1
###########################
# ---------------------------
# Fixed parameters
# ---------------------------
fixed <- c(
  mu = 0.001, # natural death rate
  d = 0.3   # disease death rate
)


# ---------------------------
# HPAI ODE model
# ---------------------------
hpai_model <- function(params, times, fixed) {
  
  hpai_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      
      Lambda_t <- function(t, Lambda_amp){
        if(t <= 60) return(Lambda_amp)          # Nov: peak arrivals
        else return(0)            # Feb-Jul: low
      }
      
      delta_t <- function(t, delta_amp){
        if(t <= 120) return(0)          
        else return(delta_amp)       
      }
      
      
      
      # Time-varying migration
      Lambda_now <- Lambda_t(time, Lambda_amp)
      delta_now <- delta_t(time, delta_amp)
      
      # ODEs
      dS_w <- Lambda_now - beta*S_w*I_w - (mu + delta_now)*S_w
      dI_w <- beta*S_w*I_w  - (mu + d + delta_now)*I_w     #I_observation
      dT_w <- d * I_w
      
      dS_r <- -beta_wr *S_r * I_w 
      dI_r <-  beta_wr *S_r * I_w  - d *I_r     #I_observation
      dT_r <- d * I_r
      
      return(list(c(dS_w, dI_w, dT_w, dS_r, dI_r, dT_r)))
    })
  }
  
  # Initial conditions
  S_w0 <- params[[1]] * 100
  I_w0 <- params[[2]] 
  beta_w <- params[[3]]/100
  Lambda_amp_w <- params[[4]]*10
  delta_amp_w  <-  params[[5]] 
  #########
  S_r0 <- params[[6]] * 100
  #  I_r0 <- params[[7]] 
  beta_wrr <- params[[7]]/100
  
  
  initial_values <- c(S_w = S_w0, I_w = I_w0, T_d = 0, S_r = S_r0, I_r = 0, T_r=0)
  
  # Parameters for ODE
  parameters_values <- c(
    beta = beta_w,
    Lambda_amp =  Lambda_amp_w,
    delta_amp  =   delta_amp_w,
    beta_wr <- beta_wrr,
    fixed
  )
  
  ode(initial_values, times, hpai_equations, parameters_values)
}

# ---------------------------
# Prior function
# ---------------------------
prior <- function(parameters){
  if(any(parameters < 0)) return(-Inf)
  
  S_w0_prior         <- dnorm(parameters[[1]], mean = 3, sd = 1.5, log = TRUE)
  I_w0_prior         <- dnorm(parameters[[2]], mean = 1.5, sd = 0.5, log = TRUE)
  beta_w_prior       <- dnorm(parameters[[3]], 0.5, 0.2, log = TRUE)
  Lambda_amp_w_prior <- dnorm(parameters[[4]], mean = 2, sd = 1, log = TRUE)
  delta_amp_w_prior  <- dnorm(parameters[[5]], 0.01, sd = 0.005, log = TRUE)
  S_r0_prior         <- dnorm(parameters[[6]], 1, sd = 0.5, log = TRUE)
  beta_wrr_prior    <- dnorm(parameters[[7]], 0.5, sd = 0.2, log = TRUE)
  
  S_w0_prior + I_w0_prior + beta_w_prior + 
    Lambda_amp_w_prior + delta_amp_w_prior + 
    S_r0_prior +   beta_wrr_prior 
}

# ---------------------------
# Likelihood function
# ---------------------------
likelihood <- function(parameters, fixed, data, data2, times) {
  
  if(any(parameters < 0)) return(-Inf)
  
  model_output <- hpai_model(params = parameters, times = times, fixed = fixed)
  
  observations1 <- data$DailyCases[-1]
  predictions1  <- diff(model_output[,"T_d"])
  
  observations2 <- data2$DailyCases[-1]
  predictions2  <- diff(model_output[,"T_r"])
  
  # Safeguard
  predictions1[!is.finite(predictions1) | predictions1 <= 0] <- 1e-6
  predictions2[!is.finite(predictions2) | predictions2 <= 0] <- 1e-6
  
  
  d1 <- dpois(x = observations1, lambda = predictions1, log = TRUE)
  d2 <- dpois(x = observations2, lambda = predictions2, log = TRUE)
  
  log_likelihood <- sum(d1)+sum(d2)
  
  return(log_likelihood)
  
}

# ---------------------------
# Joint log-probability
# ---------------------------
joint_log_prob <- function(parameters, fixed, data, data2,times) {
  likelihood(parameters, fixed, data, data2, times) + prior(parameters)
}

# ---------------------------
# Initial guesses for MCMC
# ---------------------------
parameters_initial <- c(
  S_w0       = 1.26,   # initial susceptible birds
  I_w0       = 0.3,
  beta_w       = 0.3,
  Lambda_amp_w = 0.07,  # peak daily arrival
  delta_amp_w  = 0.006,   # peak departure rate
  S_r0 = 2,
  beta_wrr = 0.17
)

times <- data_filtered$t

# ---------------------------
# Run MCMC
# ---------------------------
nbatch_runs <- 100000
scale_tune  <- 0.1

mcmc_samples <- metrop(
  obj     = joint_log_prob,
  initial = parameters_initial,
  nbatch  = nbatch_runs,
  blen    = 1,
  scale   = scale_tune,
  data    = data_filtered,
  data2   = data2_filtered,
  fixed   = fixed,
  times   = times
)


mcmc_samples <- metrop(
  obj=mcmc_samples,
  nbatch = 100000,
  blen = 100,
  scale = 0.1,
  data    = data_filtered,
  data2   = data2_filtered,
  fixed = fixed,
  times = times
)

cat("MCMC Acceptance Rate:\n")
print(mcmc_samples$accept) # Should ideally be between 0.2 and 0.5

burn_in <- 20000

posterior_samples <- mcmc_samples$batch[1:nbatch_runs, ]



colnames(posterior_samples) <- c(
  "S_w0",
  "I_w0",
  "beta_w",
  "Lambda_amp_w",
  "delta_amp_w",
  "S_r0",
  "beta_wr"
)

posterior_chain <- as.mcmc(posterior_samples)

# 
ggs_object <- ggs(posterior_chain)
ggmcmc(ggs_object)


cat("\nPosterior Mean:\n")
# mle <- colMeans(posterior_samples)
# print(mle)
post_mean <- colMeans(posterior_samples)
print(post_mean)

# 
# cat("\n95% Credible Intervals:\n")
ci <- apply(posterior_samples, 2, quantile, probs = c(0.025, 0.975))
print(ci)

fit_output <- hpai_model(
  params = post_mean,
  times  = times,
  fixed  = fixed
)

fit_df <- as.data.frame(fit_output)
fit_df$t <- times
fit_df <- fit_df %>%
  mutate(
    daily_waterfowl = c(NA, diff(T_d)),
    daily_raptor    = c(NA, diff(T_r))
  )

plot_data <- bind_rows(
  # Waterfowl
  data.frame(
    t = data_filtered$t,
    DailyCases = data_filtered$DailyCases,
    fitted = fit_df$daily_waterfowl,
    species = "Waterfowl"
  ),
  
  # Raptors
  data.frame(
    t = data2_filtered$t,
    DailyCases = data2_filtered$DailyCases,
    fitted = fit_df$daily_raptor,
    species = "Raptors"
  )
)


ggplot(plot_data, aes(x = t)) +
  geom_line(
    aes(y = DailyCases),
    color = "black",
    linewidth = 0.6
  ) +
  geom_line(
    aes(y = fitted),
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  labs(
    #title = "HPAI: Observed vs Fitted Daily Cases",
    x = "Time (days)",
    y = "Daily cases"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold")
  )

# 
# sample1 = mcmc_samples
# save(sample1, file="sample1_ontario-both-waterfowl-raptor.RData")
#  load(file="sample1_ontario-both-waterfowl-raptor.RData")
# #posterior_chain <- as.mcmc(sample1)
# posterior_samples <- sample1$batch[1:100000, ]
# 
# post_mean <- colMeans(posterior_samples)
# print(post_mean)


############


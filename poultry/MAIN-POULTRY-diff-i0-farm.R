rm(list = ls())

# =========================================================
# LIBRARIES
# =========================================================
library(deSolve)
library(mcmc)
library(coda)
library(stringr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(ggmcmc)
library(scales)

# =========================================================
# HARDCODE DATA
# =========================================================
ontario_farm <- data.frame(
  Outbreak = c("ON-IP50", "ON-IP51", "ON-IP52", "ON-IP53", "ON-IP54",
               "ON-IP55", "ON-IP56", "ON-IP57", "ON-IP59", "ON-IP60", "ON-IP61"),
  Susceptible = c(15695, 4300, 11274, 15901, 16371,
                  24032, 61222, 61963, 54736, 47497, 32062),
  Dead = c(4778, 1831, 3637, 3937, 3263,
           3351, 1050, 6143, 5350, 13127, 17256),
  clinical_date = as.Date(c("2024-12-13", "2024-12-18", "2024-12-19", "2024-12-20", "2024-12-20",
                            "2024-12-25", "2025-01-02", "2025-01-07", "2025-02-26", "2025-03-12", "2025-03-21")),
  depop_date = as.Date(c("2024-12-18", "2024-12-22", "2024-12-21", "2024-12-22", "2024-12-23",
                         "2024-12-28", "2025-01-05", "2025-01-15", "2025-03-02", "2025-03-16", "2025-03-27")),
  species_group = c("Turkey", "Turkey", "Turkey", "Turkey", "Turkey",
                    "Turkey", "Turkey", "Other", "Turkey", "Other", "Turkey")
) %>%
  mutate(
    fraction_dead = Dead / Susceptible,
    delay_days = as.numeric(depop_date - clinical_date),
    flock_group = case_when(
      Susceptible <= 5000 ~ "<=5,000",
      Susceptible <= 20000 ~ "5,001–20,000",
      TRUE ~ ">20,000"
    ),
    flock_group = factor(flock_group, levels = c("<=5,000", "5,001–20,000", ">20,000")),
    I_group = case_when(
      Outbreak %in% c("ON-IP56") ~ "I_p1",
      Outbreak %in% c("ON-IP59") ~ "I_p2",
      Outbreak %in% c("ON-IP57") ~ "I_p3",
      Outbreak %in% c("ON-IP55") ~ "I_p4",
      Outbreak %in% c("ON-IP54") ~ "I_p5",
      Outbreak %in% c("ON-IP53") ~ "I_p6",
      Outbreak %in% c("ON-IP60") ~ "I_p7",
      Outbreak %in% c("ON-IP50") ~ "I_p8",
      Outbreak %in% c("ON-IP52") ~ "I_p9",
      Outbreak %in% c("ON-IP51") ~ "I_p10",
      TRUE ~ "I_p11"
    )
  ) %>%
  filter(delay_days >= 0)

N_p_farms <- ontario_farm$Susceptible
dead_poultry_obs <- ontario_farm$Dead
delay_days_farms <- ontario_farm$delay_days
species_group <- ontario_farm$species_group
F <- length(N_p_farms)

# =========================================================
# SINGLE-FARM ODE MODEL
# Count-based states:
#   S = susceptible birds
#   I = infected birds
#   T = cumulative dead birds
# =========================================================
hpai_model_single_farm <- function(I_p0, beta, d_p, times, N_p) {
  
  hpai_equations <- function(time, state, parameters) {
    S <- state[1]
    I <- state[2]
    T <- state[3]
    
    beta <- parameters[1]
    d_p  <- parameters[2]
    
    dS <- -beta * S * I / N_p
    dI <-  beta * S * I / N_p - d_p * I
    dT <-  d_p * I
    
    list(c(dS, dI, dT))
  }
  
  initial_values <- c(
    S = N_p - I_p0,
    I = I_p0,
    T = 0
  )
  
  parameters_values <- c(beta, d_p)
  
  out <- ode(
    y = initial_values,
    times = times,
    func = hpai_equations,
    parms = parameters_values,
    method = "lsoda",
    rtol = 1e-6,
    atol = 1e-8
  )
  
  colnames(out) <- c("time", "S", "I", "T")
  out
}

# =========================================================
# ASSIGN GROUP-SPECIFIC INITIAL INFECTED COUNTS
# params = c(I_p1, I_p2, I_p3, I_p4, d_p)
# =========================================================
get_initial_I <- function(outbreak_name, params) {
  if (outbreak_name %in% c("ON-IP56")) {
    params[1]   # I_p1
  } else if (outbreak_name %in% c("ON-IP59")) {
    params[2]   # I_p2
  } else if (outbreak_name %in% c("ON-IP57")) {
    params[3]   # I_p3
  } else if (outbreak_name %in% c("ON-IP55")) {
    params[4]   # I_p4
  } else if (outbreak_name %in% c("ON-IP54")) {
    params[5]   # I_p5
  } else if (outbreak_name %in% c("ON-IP53")) {
    params[6]   # I_p6
  } else if (outbreak_name %in% c("ON-IP60")) {
    params[7]   # I_p7
  } else if (outbreak_name %in% c("ON-IP50")) {
    params[8]   # I_p8
  } else if (outbreak_name %in% c("ON-IP52")) {
    params[9]   # I_p9
  } else if (outbreak_name %in% c("ON-IP51")) {
    params[10]   # I_p10
  } else {
    params[11]   # I_p11
  }
}

# =========================================================
# GET PREDICTED DEAD FOR ALL FARMS
# params = c(I_p1, I_p2, I_p3, I_p4, d_p)
# =========================================================
get_predicted_dead <- function(params, N_p_farms, delay_days_farms, species_group, F, outbreak_names) {
  
  d_p <- params[12]
  predicted_T <- numeric(F)
  
  beta <- d_p * 1.93817
  
 # beta <- d_p * 1.290563
  
  for (f in 1:F) {
    N_p <- N_p_farms[f]
    delay <- delay_days_farms[f]
    outbreak_name <- outbreak_names[f]
    
    I_p0 <- get_initial_I(outbreak_name, params)
    
    if (!is.finite(I_p0) || I_p0 <= 0 || I_p0 >= N_p) {
      return(rep(NA_real_, F))
    }
    
    times <- seq(0, delay, by = 0.1)
    
    out <- hpai_model_single_farm(
      I_p0 = I_p0,
      beta = beta,
      d_p = d_p,
      times = times,
      N_p = N_p
    )
    
    predicted_T[f] <- tail(out[, "T"], 1)
  }
  
  predicted_T
}

# =========================================================
# LOG-LIKELIHOOD
# =========================================================
loglik <- function(params, N_p_farms, dead_poultry_obs, delay_days_farms,
                   species_group, F, outbreak_names) {
  
  predicted_T <- get_predicted_dead(
    params, N_p_farms, delay_days_farms, species_group, F, outbreak_names
  )
  
  if (any(!is.finite(predicted_T)) || any(predicted_T <= 0)) return(-1e10)
  
  sum(dpois(dead_poultry_obs, lambda = predicted_T, log = TRUE))
}

# =========================================================
# PRIOR
# params = c(I_p1, I_p2, I_p3, I_p4, d_p)
# =========================================================
prior <- function(parameters) {
  if (any(parameters <= 0)) return(-1e10)
  
  I_p1_prior <- dnorm(parameters[1], mean = 30, sd = 20, log = TRUE) #56
  I_p2_prior <- dnorm(parameters[2], mean = 30,  sd = 20,   log = TRUE) #59
  I_p3_prior <- dnorm(parameters[3], mean = 150, sd = 50,  log = TRUE) #57
  I_p4_prior <- dnorm(parameters[4], mean = 30, sd = 20,  log = TRUE) #55
  I_p5_prior <- dnorm(parameters[5], mean = 30, sd = 20,  log = TRUE) #54
  I_p6_prior <- dnorm(parameters[6], mean = 200, sd = 100,  log = TRUE) #53
  I_p7_prior <- dnorm(parameters[7], mean = 100, sd = 50,  log = TRUE) #60
  I_p8_prior <- dnorm(parameters[8], mean = 200, sd = 100,  log = TRUE) #50
  I_p9_prior <- dnorm(parameters[9], mean = 300, sd = 100,  log = TRUE) #52
  I_p10_prior <- dnorm(parameters[10], mean = 200, sd = 100,  log = TRUE) #51
  I_p11_prior <- dnorm(parameters[11], mean = 2000, sd = 1000,  log = TRUE) #61
  d_p_prior  <- dnorm(parameters[12], mean = 1.5, sd = 1,   log = TRUE)
  
  I_p1_prior + I_p2_prior + I_p3_prior + I_p4_prior + I_p5_prior + I_p6_prior+ I_p7_prior+ I_p8_prior + I_p9_prior
  + I_p10_prior+ I_p11_prior  + d_p_prior
}

# =========================================================
# POSTERIOR
# =========================================================
logpost <- function(params) {
  prior_val <- prior(params)
  if (!is.finite(prior_val)) return(-1e10)
  
  lik_val <- loglik(
    params,
    N_p_farms,
    dead_poultry_obs,
    delay_days_farms,
    species_group,
    F,
    ontario_farm$Outbreak
  )
  if (!is.finite(lik_val)) return(-1e10)
  
  prior_val + lik_val
}

# =========================================================
# RUN MCMC
# =========================================================
initial_params <- c(10,
                    10,
                    100,
                    20,
                    20,
                    120,
                    300,
                    150,
                    200,
                    200,
                    2000,
                    1.5)

#cat("Initial log-posterior:", logpost(initial_params), "\n")
scale_vec <- c(
  0.5,   # I_p1 = 10
  0.5,   # I_p2 = 10
  8,     # I_p3 = 100
  1,     # I_p4 = 20
  1,     # I_p5 = 20
  6,     # I_p6 = 120
  15,    # I_p7 = 300
  10,    # I_p8 = 200
  10,    # I_p9 = 200
  10,    # I_p10 = 200
  80,    # I_p11 = 2000
  0.05   # d_p = 1.5
)

mcmc_res <- metrop(
  logpost,
  initial = initial_params,
  nbatch = 100000,
  blen = 10,
  scale = scale_vec
)

mcmc_res <- metrop(
  obj = mcmc_res,
  nbatch = 100000,
  blen = 10,
  scale = scale_vec
)


sample1 <- mcmc_res
save(sample1, file = "sample1_ontario-new-11-I-y-kx-poultry-farm-2024-2025.RData")

# sample1 <- mcmc_res
# save(sample1, file = "sample1_ontario-new-y-kx-poultry-farm-2024-2025.RData")
# 
# print(mcmc_res$accept)
load("sample1_ontario-new-11-I-y-kx-poultry-farm-2024-2025.RData")

# =========================================================
# POSTERIOR SUMMARY
# =========================================================
burnin <- 0
posterior_samples <- sample1$batch[(burnin + 1):nrow(sample1$batch), ]
colnames(posterior_samples) <- c("I_p1", "I_p2", "I_p3", "I_p4","I_p5", "I_p6","I_p7","I_p8","I_p9","I_p10","I_p11","d_p")

param_estimates <- colMeans(posterior_samples)
print(param_estimates)



# =========================================================
# PREDICTED VS OBSERVED
# =========================================================
predicted_dead <- get_predicted_dead(
  param_estimates,
  N_p_farms,
  delay_days_farms,
  species_group,
  F,
  ontario_farm$Outbreak
)

comparison_df <- data.frame(
  farm = ontario_farm$Outbreak,
  I_group = ontario_farm$I_group,
  species = species_group,
  susceptible = ontario_farm$Susceptible,
  observed = dead_poultry_obs,
  predicted = predicted_dead
) %>%
  mutate(
    abs_error = abs(observed - predicted),
    pct_error = 100 * abs(observed - predicted) / observed,
    observed_fraction = observed / susceptible,
    predicted_fraction = predicted / susceptible
  )

print(comparison_df)

# =========================================================
# PLOT: OBSERVED VS FITTED CUMULATIVE DEAD
# =========================================================
plot_df <- comparison_df %>%
  arrange(observed) %>%
  mutate(farm = factor(farm, levels = farm)) %>%
  select(farm, observed, predicted) %>%
  pivot_longer(
    cols = c(observed, predicted),
    names_to = "type",
    values_to = "dead_birds"
  ) %>%
  mutate(
    type = recode(type,
                  observed = "Observed",
                  predicted = "Fitted")
  )

ggplot(plot_df, aes(x = farm, y = dead_birds, fill = type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.8) +
  scale_fill_manual(values = c("Observed" = "#1f77b4", "Fitted" = "#ff7f0e")) +
  labs(
    title = "Observed vs Model-Fitted Cumulative Dead Birds",
    subtitle = paste(
      "Ontario commercial poultry farms - Clade 2.3.4.4b (H5N1) -",
      min(ontario_farm$clinical_date, na.rm = TRUE), "to",
      max(ontario_farm$clinical_date, na.rm = TRUE)
    ),
    x = "Farm (ordered by observed deaths)",
    y = "Cumulative number of dead birds",
    fill = "Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9.5),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# =========================================================
# PLOT: OBSERVED VS FITTED FRACTION DEAD
# =========================================================
comparison_df <- comparison_df %>%
  arrange(observed_fraction) %>%
  mutate(farm = factor(farm, levels = farm))

fraction_plot_df <- comparison_df %>%
  select(farm, observed_fraction, predicted_fraction) %>%
  pivot_longer(
    cols = c(observed_fraction, predicted_fraction),
    names_to = "type",
    values_to = "fraction_dead"
  ) %>%
  mutate(
    type = recode(
      type,
      observed_fraction = "Observed",
      predicted_fraction = "Fitted"
    )
  )

ggplot(fraction_plot_df, aes(x = farm, y = fraction_dead, fill = type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.8) +
  scale_fill_manual(values = c("Observed" = "#1f77b4", "Fitted" = "#ff7f0e")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Observed vs Fitted Fraction of Dead Birds",
    x = "Farm",
    y = "Fraction dead",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

#############
# =========================================================
# PLOT: OBSERVED VS FITTED CUMULATIVE DEAD
# =========================================================
plot_df <- comparison_df %>%
  arrange(observed) %>%
  mutate(farm = factor(farm, levels = farm)) %>%
  select(farm, observed, predicted) %>%
  pivot_longer(
    cols = c(observed, predicted),
    names_to = "type",
    values_to = "dead_birds"
  ) %>%
  mutate(
    type = recode(
      type,
      observed = "Observed",
      predicted = "Fitted"
    )
  )

# save the farm order from the first plot
farm_order <- levels(plot_df$farm)

ggplot(plot_df, aes(x = farm, y = dead_birds, fill = type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.8) +
  scale_fill_manual(values = c("Observed" = "#1f77b4", "Fitted" = "#ff7f0e")) +
  labs(
    title = "Observed vs Model-Fitted Cumulative Dead Birds",
    subtitle = paste(
      "Ontario commercial poultry farms - Clade 2.3.4.4b (H5N1) -",
      min(ontario_farm$clinical_date, na.rm = TRUE), "to",
      max(ontario_farm$clinical_date, na.rm = TRUE)
    ),
    x = "Farm (ordered by observed deaths)",
    y = "Cumulative number of dead birds",
    fill = "Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9.5),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# =========================================================
# PLOT: OBSERVED VS FITTED FRACTION DEAD
# use the SAME farm order as above
# =========================================================
fraction_plot_df <- comparison_df %>%
  mutate(farm = factor(farm, levels = farm_order)) %>%
  select(farm, observed_fraction, predicted_fraction) %>%
  pivot_longer(
    cols = c(observed_fraction, predicted_fraction),
    names_to = "type",
    values_to = "fraction_dead"
  ) %>%
  mutate(
    type = recode(
      type,
      observed_fraction = "Observed",
      predicted_fraction = "Fitted"
    )
  )

ggplot(fraction_plot_df, aes(x = farm, y = fraction_dead, fill = type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.8) +
  scale_fill_manual(values = c("Observed" = "#1f77b4", "Fitted" = "#ff7f0e")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Observed vs Fitted Fraction of Dead Birds",
    x = "Farm",
    y = "Fraction dead",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# =========================================================
# DIAGNOSTICS
# =========================================================
# posterior_chain <- as.mcmc(mcmc_res$batch)
# print(summary(mcmc_res$batch))
# print(effectiveSize(posterior_chain))
# plot(posterior_chain)
# 
# ggs_object <- ggs(posterior_chain)
# ggmcmc(ggs_object)



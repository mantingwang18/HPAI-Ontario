# =============================================
# 3. SPILLOVER LINKAGE: Wild → Poultry Farms
# Inhomogeneous Poisson process driven by I_w(t) + α I_r(t)
# =============================================
rm(list = ls())

library(deSolve)
library(mcmc)
library(coda)
library(dplyr)
library(ggplot2)
library(lubridate)



###################################3
load("sample1_ontario-new-11-I-y-kx-poultry-farm-2024-2025.RData")  # farm: mcmc_res or sample1

farm_burn <- 50000
farm_post <- as.matrix(sample1$batch[-(1:farm_burn), ])   # or sample1$batch
colnames(farm_post) <- c(paste0("I_p",1:11),"d_p")
farm_mean <- colMeans(farm_post)


########################3
# 2. Re-run wild model at posterior mean to get I_w(t), I_r(t) on fine grid
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

  
  initial_values <- c(S_w = 171.4 , 
                      I_w = 0.605 , 
                      T_d = 0, 
                      S_r =107.2, 
                      I_r = 0, 
                      T_r=0)
  
  # Parameters for ODE
  parameters_values <- c(
    fixed
  )
  
  ode(initial_values, times, hpai_equations, parameters_values)
}

fixed <- c(
  mu = 0.001, # natural death rate
  d = 0.3,   # disease death rate
  beta <- 0.174/100,
  beta_wr <- 0.213 /100,
  Lambda_amp <-1.09,
  delta_amp <-0.00536
)

times_full <- 1:150   # Nov 1 2024 → ~April 2025 (covers all farms)
wild_fit <- hpai_model(times = times_full, fixed = fixed)
Iw_vec <- wild_fit[,"I_w"]   # infectious waterfowl
Ir_vec <- wild_fit[,"I_r"]   # infectious raptors

# 3. Clinical dates → t (days since 2024-11-01 = t=1)
clinical_t <- c(43,48,49,50,50,55,63,68,118,132,141)   # from date calc
farm_names <- c("ON-IP50","ON-IP51","ON-IP52","ON-IP53","ON-IP54",
                "ON-IP55","ON-IP56","ON-IP57","ON-IP59","ON-IP60","ON-IP61")

# 4. Back-calculate t_intro,f for each posterior draw (or mean first)
calc_t_intro <- function(Ip_vec, dp) {   # Ip_vec = 11 values
  r <- 1.93817 * dp - dp                  # early growth rate
  tau_pre <- rep(NA, 11)
  for (f in 1:11) {
    ip <- Ip_vec[f]
    if (ip > 1) {
      tau_pre[f] <- log(ip) / r
    } else {
      tau_pre[f] <- 0   # ON-IP57 special case
    }
  }
  t_intro <- clinical_t - tau_pre
  pmax(t_intro, 1)   # can't be before Nov 1
}

t_intro_mean <- calc_t_intro(farm_mean[1:11], farm_mean["d_p"])

# 5. SPILLOVER MODEL
spillover_model <- function(beta_spill, alpha, Iw, Ir, times, t_events) {
  lambda_t <- beta_spill * (Iw + alpha * Ir)          # force of introduction
  # log-likelihood (inhomogeneous Poisson)
  log_lambda_events <- sum(log(lambda_t[round(t_events)] + 1e-10))
  integral <- sum(lambda_t) * 1.0                     # dt=1 day
  log_lambda_events - integral
}

# 6. Prior + joint log-prob for spillover parameters only
prior_spill <- function(par) {
  beta_spill <- par[1]; alpha <- par[2]
  if (beta_spill <= 0 || alpha < 0) return(-Inf)
  dnorm(beta_spill, mean = 0.001, sd = 0.005, log = TRUE) +
    dnorm(alpha, mean = 0.5, sd = 1, log = TRUE)
}

logpost_spill <- function(par) {
  beta_spill <- par[1]; alpha <- par[2]
  lik <- spillover_model(beta_spill, alpha, Iw_vec, Ir_vec, times_full, t_intro_mean)
  prior_spill(par) + lik
}

# 7. Run MCMC for spillover (very fast)
init_spill <- c(beta_spill = 0.001, alpha = 0.3)
spill_mcmc <- metrop(logpost_spill, initial = init_spill,
                     nbatch = 50000, blen = 10, scale = c(0.0005, 0.2))
spill_mcmc <- metrop(obj = spill_mcmc, nbatch = 50000, blen = 10,
                     scale = c(0.0005, 0.2))

# 8. Posterior summary
burn_sp <- 10000
spill_post <- spill_mcmc$batch[-(1:burn_sp), ]
colnames(spill_post) <- c("beta_spill", "alpha")
spill_mean <- colMeans(spill_post)
spill_ci   <- apply(spill_post, 2, quantile, probs = c(0.025, 0.975))
cat("Spillover posterior:\n")
print(spill_mean)
print(spill_ci)
cat("Acceptance rate:", spill_mcmc$accept, "\n")

# 9. Publication figure (main text Figure X)
lambda_mean <- spill_mean[1] * (Iw_vec + spill_mean[2] * Ir_vec)

df_plot <- data.frame(
  t = times_full,
  Iw = Iw_vec,
  Ir = Ir_vec,
  lambda = lambda_mean,
  wild_dead = c(NA, diff(wild_fit[,"T_d"] + wild_fit[,"T_r"]))   # total wild deaths
) %>%
  mutate(date = as.Date("2024-11-01") + t - 1)

events_df <- data.frame(t_intro = t_intro_mean, farm = farm_names)

ggplot() +
  geom_line(data = df_plot, aes(x = date, y = lambda, color = "Spillover force λ(t)"), linewidth = 1.1) +
  geom_point(data = events_df, aes(x = as.Date("2024-11-01") + t_intro - 1, y = 0),
             color = "red", size = 3, shape = 4) +
  geom_vline(data = events_df, aes(xintercept = as.Date("2024-11-01") + t_intro - 1),
             linetype = "dashed", alpha = 0.6) +
  #scale_y_continuous(sec.axis = sec_axis(~./100, name = "Infectious wild birds")) +
  labs(x = "Date (2024–2025)", y = "Spillover rate λ(t) ") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# 10. Predictive check (simulate 500 introduction sets)
simulate_intros <- function(n_sim = 500) {
  sims <- list()
  for (i in 1:n_sim) {
    par <- spill_post[sample(nrow(spill_post),1), ]
    lambda_sim <- par[1] * (Iw_vec + par[2]*Ir_vec)
    # cumulative intensity
    cum_int <- cumsum(lambda_sim)
    # inhomogeneous sampling
    u <- runif(11) * max(cum_int)
    sim_t <- approx(cum_int, times_full, xout = u)$y
    sims[[i]] <- sim_t
  }
  do.call(rbind, sims)
}
sim_matrix <- simulate_intros()
# overlay on figure or compute coverage % of observed events inside 95% sim bands

# 11. Bayes factor vs null (constant hazard) – quick approximation
loglik_null <- 11 * log(11 / max(times_full)) - 11   # rough
bf <- exp(mean(apply(spill_post,1,function(p) spillover_model(p[1],p[2],Iw_vec,Ir_vec,times_full,t_intro_mean))) - loglik_null)
cat("Approximate Bayes factor (full vs null):", bf, "\n")   # >10 = strong support

# After you have farm_post and spill_post


# Now run spillover MCMC *once* with the mean t_intro (as before) 
# OR (better) run it for every draw of t_intro_post and average the log-posterior.
# Then plot 95% credible bands for λ(t).

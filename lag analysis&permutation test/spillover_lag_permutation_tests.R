# ============================================================
# TEMPORAL TESTS FOR WILD-BIRD -> POULTRY SPILLOVER ALIGNMENT
#
# Test 1: Lag-profile analysis, shifting wild-bird pressure
#         from -14 to +14 days.
#
# Test 2: Randomization analysis with:
#         (a) introduction dates randomized over the study period;
#         (b) circular shifts preserving spacing among introductions.
#
# Positive lag means that poultry-farm introduction pressure follows
# wild-bird infection pressure by the specified number of days.
# ============================================================

rm(list = ls())

library(deSolve)
library(dplyr)
library(ggplot2)
library(lubridate)

# ============================================================
# 1. USER SETTINGS
# ============================================================

farm_posterior_file <-
  "sample1_ontario-new-11-I-y-kx-poultry-farm-2024-2025.RData"

farm_burn <- 50000

study_start <- as.Date("2024-11-01")
times_full <- 1:150

lag_values <- -14:14
number_randomizations <- 1000

# Upper bound used when profiling the relative contribution of raptors.
alpha_upper <- 20

# Small positive value used only for numerical protection.
epsilon <- 1e-12

# ============================================================
# 2. LOAD THE WITHIN-FARM POSTERIOR
# ============================================================

load(farm_posterior_file)

if (!exists("sample1")) {
  stop("The RData file must contain an object named 'sample1'.")
}

if (is.null(sample1$batch)) {
  stop("The object 'sample1' does not contain a $batch matrix.")
}

farm_batch <- as.matrix(sample1$batch)

if (ncol(farm_batch) != 12) {
  stop(
    "Expected 12 posterior columns: I_p1,...,I_p11,d_p; found ",
    ncol(farm_batch), "."
  )
}

if (farm_burn >= nrow(farm_batch)) {
  stop(
    "farm_burn must be smaller than the number of stored posterior rows."
  )
}

farm_post <- farm_batch[(farm_burn + 1):nrow(farm_batch), , drop = FALSE]

colnames(farm_post) <- c(
  paste0("I_p", 1:11),
  "d_p"
)

farm_mean <- colMeans(farm_post)

# ============================================================
# 3. RECONSTRUCT WILD-BIRD INFECTION PRESSURE
# ============================================================

hpai_model <- function(times, fixed) {

  hpai_equations <- function(time, variables, parameters) {

    with(as.list(c(variables, parameters)), {

      Lambda_t <- function(t, Lambda_amp) {
        if (t <= 60) {
          Lambda_amp
        } else {
          0
        }
      }

      delta_t <- function(t, delta_amp) {
        if (t <= 120) {
          0
        } else {
          delta_amp
        }
      }

      Lambda_now <- Lambda_t(time, Lambda_amp)
      delta_now <- delta_t(time, delta_amp)

      dS_w <-
        Lambda_now -
        beta * S_w * I_w -
        (mu + delta_now) * S_w

      dI_w <-
        beta * S_w * I_w -
        (mu + d + delta_now) * I_w

      dT_w <- d * I_w

      dS_r <- -beta_wr * S_r * I_w

      dI_r <-
        beta_wr * S_r * I_w -
        d * I_r

      dT_r <- d * I_r

      list(c(
        dS_w,
        dI_w,
        dT_w,
        dS_r,
        dI_r,
        dT_r
      ))
    })
  }

  initial_values <- c(
    S_w = 171.4,
    I_w = 0.605,
    T_w = 0,
    S_r = 107.2,
    I_r = 0,
    T_r = 0
  )

  ode(
    y = initial_values,
    times = times,
    func = hpai_equations,
    parms = fixed,
    method = "lsoda",
    rtol = 1e-8,
    atol = 1e-10
  )
}

wild_parameters <- c(
  mu = 0.001,
  d = 0.3,
  beta = 0.174 / 100,
  beta_wr = 0.213 / 100,
  Lambda_amp = 1.09,
  delta_amp = 0.00536
)

wild_fit <- hpai_model(
  times = times_full,
  fixed = wild_parameters
)

Iw_vec <- wild_fit[, "I_w"]
Ir_vec <- wild_fit[, "I_r"]

if (any(!is.finite(Iw_vec)) || any(!is.finite(Ir_vec))) {
  stop("The reconstructed wild-bird trajectories contain invalid values.")
}

# ============================================================
# 4. RECONSTRUCT FARM INTRODUCTION TIMES
# ============================================================

farm_names <- c(
  "ON-IP50",
  "ON-IP51",
  "ON-IP52",
  "ON-IP53",
  "ON-IP54",
  "ON-IP55",
  "ON-IP56",
  "ON-IP57",
  "ON-IP59",
  "ON-IP60",
  "ON-IP61"
)

clinical_dates <- as.Date(c(
  "2024-12-13",
  "2024-12-18",
  "2024-12-19",
  "2024-12-20",
  "2024-12-20",
  "2024-12-25",
  "2025-01-02",
  "2025-01-07",
  "2025-02-26",
  "2025-03-12",
  "2025-03-21"
))

clinical_t <- as.numeric(clinical_dates - study_start) + 1

# IMPORTANT:
# These mappings come from get_initial_I() in the within-farm model.
# The I_p parameter order is not the same as the farm/date order.
Ip_by_farm <- c(
  "ON-IP50" = unname(farm_mean["I_p8"]),
  "ON-IP51" = unname(farm_mean["I_p10"]),
  "ON-IP52" = unname(farm_mean["I_p9"]),
  "ON-IP53" = unname(farm_mean["I_p6"]),
  "ON-IP54" = unname(farm_mean["I_p5"]),
  "ON-IP55" = unname(farm_mean["I_p4"]),
  "ON-IP56" = unname(farm_mean["I_p1"]),
  "ON-IP57" = unname(farm_mean["I_p3"]),
  "ON-IP59" = unname(farm_mean["I_p2"]),
  "ON-IP60" = unname(farm_mean["I_p7"]),
  "ON-IP61" = unname(farm_mean["I_p11"])
)

if (!setequal(names(Ip_by_farm), farm_names)) {
  stop(
    "Farm mapping failed. Missing farms: ",
    paste(setdiff(farm_names, names(Ip_by_farm)), collapse = ", "),
    ". Unexpected names: ",
    paste(setdiff(names(Ip_by_farm), farm_names), collapse = ", "),
    "."
  )
}

calc_t_intro <- function(
    Ip_by_farm,
    d_p,
    clinical_t,
    farm_names,
    minimum_time,
    maximum_time
) {

  early_growth_rate <- d_p * (1.93817 - 1)

  if (
    !is.finite(early_growth_rate) ||
    early_growth_rate <= 0
  ) {
    stop("The inferred early growth rate must be positive.")
  }

  initial_infected <- Ip_by_farm[farm_names]

  if (any(!is.finite(initial_infected))) {
    stop("Some farms could not be matched to an I_p parameter.")
  }

  time_before_detection <- ifelse(
    initial_infected > 1,
    log(initial_infected) / early_growth_rate,
    0
  )

  introduction_time <-
    clinical_t -
    time_before_detection

  introduction_time <- pmax(
    pmin(introduction_time, maximum_time),
    minimum_time
  )

  setNames(introduction_time, farm_names)
}

t_intro_mean <- calc_t_intro(
  Ip_by_farm = Ip_by_farm,
  d_p = farm_mean["d_p"],
  clinical_t = clinical_t,
  farm_names = farm_names,
  minimum_time = min(times_full),
  maximum_time = max(times_full)
)

introduction_data <- data.frame(
  farm = farm_names,
  clinical_date = clinical_dates,
  clinical_t = clinical_t,
  initial_infected = as.numeric(Ip_by_farm[farm_names]),
  introduction_t = as.numeric(t_intro_mean),
  introduction_date =
    study_start + as.numeric(t_intro_mean) - 1
)

print(introduction_data)

# ============================================================
# 5. CORE POINT-PROCESS FUNCTIONS
# ============================================================

shift_series <- function(x, lag_days, times) {

  # At calendar time t, use the wild-bird value at t - lag.
  # Therefore, a positive lag makes farm pressure follow wild pressure.
  shifted <- approx(
    x = times,
    y = x,
    xout = times - lag_days,
    rule = 1
  )$y

  shifted[!is.finite(shifted)] <- 0
  shifted
}

constant_hazard_loglik <- function(t_events, times) {

  number_events <- length(t_events)
  study_duration <- length(times)

  constant_rate <- number_events / study_duration

  number_events * log(constant_rate) -
    constant_rate * study_duration
}

profile_spillover_lag <- function(
    lag_days,
    Iw,
    Ir,
    times,
    t_events,
    alpha_upper = 20,
    epsilon = 1e-12
) {

  Iw_shifted <- shift_series(
    Iw,
    lag_days,
    times
  )

  Ir_shifted <- shift_series(
    Ir,
    lag_days,
    times
  )

  profile_loglik_alpha <- function(alpha) {

    relative_pressure <-
      Iw_shifted +
      alpha * Ir_shifted

    integrated_pressure <- sum(relative_pressure)

    if (
      !is.finite(integrated_pressure) ||
      integrated_pressure <= 0
    ) {
      return(-1e100)
    }

    event_pressure <- approx(
      x = times,
      y = relative_pressure,
      xout = t_events,
      rule = 1
    )$y

    if (
      any(!is.finite(event_pressure)) ||
      any(event_pressure <= 0)
    ) {
      return(-1e100)
    }

    number_events <- length(t_events)

    # Profile maximum-likelihood estimate for beta_spill.
    beta_hat <-
      number_events /
      integrated_pressure

    sum(
      log(beta_hat * event_pressure + epsilon)
    ) -
      beta_hat * integrated_pressure
  }

  alpha_fit <- optimize(
    f = profile_loglik_alpha,
    interval = c(0, alpha_upper),
    maximum = TRUE
  )

  alpha_hat <- alpha_fit$maximum

  relative_pressure <-
    Iw_shifted +
    alpha_hat * Ir_shifted

  beta_hat <-
    length(t_events) /
    sum(relative_pressure)

  data.frame(
    lag = lag_days,
    logLik = alpha_fit$objective,
    beta_spill = beta_hat,
    alpha = alpha_hat
  )
}

run_lag_profile <- function(
    t_events,
    Iw,
    Ir,
    times,
    lag_values,
    alpha_upper = 20
) {

  bind_rows(
    lapply(
      lag_values,
      profile_spillover_lag,
      Iw = Iw,
      Ir = Ir,
      times = times,
      t_events = t_events,
      alpha_upper = alpha_upper
    )
  ) %>%
    mutate(
      delta_logLik = logLik - max(logLik),
      relative_likelihood = exp(delta_logLik)
    )
}

max_lag_statistic <- function(
    t_events,
    Iw,
    Ir,
    times,
    lag_values,
    alpha_upper = 20
) {

  lag_profile <- run_lag_profile(
    t_events = t_events,
    Iw = Iw,
    Ir = Ir,
    times = times,
    lag_values = lag_values,
    alpha_upper = alpha_upper
  )

  best_spillover_loglik <-
    max(lag_profile$logLik)

  null_loglik <- constant_hazard_loglik(
    t_events = t_events,
    times = times
  )

  best_spillover_loglik - null_loglik
}

# ============================================================
# 6. TEST 1: LAG-PROFILE ANALYSIS
# ============================================================

observed_events <- as.numeric(t_intro_mean)

lag_results <- run_lag_profile(
  t_events = observed_events,
  Iw = Iw_vec,
  Ir = Ir_vec,
  times = times_full,
  lag_values = lag_values,
  alpha_upper = alpha_upper
)

best_lag_result <- lag_results %>%
  slice_max(
    order_by = logLik,
    n = 1,
    with_ties = FALSE
  )

zero_lag_result <- lag_results %>%
  filter(lag == 0)

observed_statistic <- max_lag_statistic(
  t_events = observed_events,
  Iw = Iw_vec,
  Ir = Ir_vec,
  times = times_full,
  lag_values = lag_values,
  alpha_upper = alpha_upper
)

cat("\n==============================================\n")
cat("LAG-PROFILE RESULTS\n")
cat("==============================================\n")
cat("Best lag:", best_lag_result$lag, "days\n")
cat("Best-lag log-likelihood:", best_lag_result$logLik, "\n")
cat("Zero-lag log-likelihood:", zero_lag_result$logLik, "\n")
cat(
  "Best versus zero-lag improvement:",
  best_lag_result$logLik - zero_lag_result$logLik,
  "\n"
)
cat("Best-lag beta_spill:", best_lag_result$beta_spill, "\n")
cat("Best-lag alpha:", best_lag_result$alpha, "\n")
cat(
  "Best-lag improvement over constant hazard:",
  observed_statistic,
  "\n"
)

lag_plot <- ggplot(
  lag_results,
  aes(x = lag, y = delta_logLik)
) +
  geom_hline(
    yintercept = -2,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "grey45"
  ) +
  geom_line(
    linewidth = 1,
    color = "#2166AC"
  ) +
  geom_point(
    size = 2,
    color = "#2166AC"
  ) +
  geom_point(
    data = best_lag_result,
    aes(
      x = lag,
      y = logLik - max(lag_results$logLik)
    ),
    color = "#B2182B",
    size = 3
  ) +
  annotate(
    "text",
    x = best_lag_result$lag,
    y = 0,
    label = paste0(
      "Best lag = ",
      best_lag_result$lag,
      " days"
    ),
    hjust = ifelse(best_lag_result$lag >= 0, 1.05, -0.05),
    vjust = -0.8,
    color = "#B2182B"
  ) +
  scale_x_continuous(
    breaks = seq(-14, 14, by = 2)
  ) +
  labs(
    x = "Lag between wild-bird pressure and farm introduction (days)",
    y = expression(Delta * " log-likelihood")
    #caption = paste(
    # "Positive lags indicate that poultry-farm introduction pressure",
    #  "follows wild-bird infection pressure."
   # )
  ) +
  theme_minimal(base_size = 13)

print(lag_plot)

# ============================================================
# 7. TEST 2A: UNIFORM-DATE RANDOMIZATION TEST
# ============================================================

uniform_statistics <- numeric(number_randomizations)

for (permutation_index in seq_len(number_randomizations)) {

  randomized_dates <- sort(
    sample(
      times_full,
      size = length(observed_events),
      replace = FALSE
    )
  )

  uniform_statistics[permutation_index] <-
    max_lag_statistic(
      t_events = randomized_dates,
      Iw = Iw_vec,
      Ir = Ir_vec,
      times = times_full,
      lag_values = lag_values,
      alpha_upper = alpha_upper
    )

  if (permutation_index %% 100 == 0) {
    cat(
      "Completed",
      permutation_index,
      "of",
      number_randomizations,
      "uniform randomizations\n"
    )
  }
}

uniform_p_value <- (
  1 +
  sum(uniform_statistics >= observed_statistic)
) / (
  number_randomizations + 1
)

cat("\n==============================================\n")
cat("UNIFORM-DATE RANDOMIZATION RESULTS\n")
cat("==============================================\n")
cat("Observed statistic:", observed_statistic, "\n")
cat("Randomization p-value:", uniform_p_value, "\n")

uniform_results <- data.frame(
  statistic = uniform_statistics,
  randomization = "Uniform dates"
)

uniform_plot <- ggplot(
  uniform_results,
  aes(x = statistic)
) +
  geom_histogram(
    bins = 35,
    fill = "#9ECAE1",
    color = "white"
  ) +
  geom_vline(
    xintercept = observed_statistic,
    color = "#B2182B",
    linewidth = 1.1
  ) +
  annotate(
    "text",
    x = observed_statistic,
    y = Inf,
    label = paste0(
      "Observed\np = ",
      format.pval(uniform_p_value, digits = 3)
    ),
    hjust = ifelse(
      observed_statistic >
        median(uniform_statistics),
      1.05,
      -0.05
    ),
    vjust = 1.2,
    color = "#B2182B"
  ) +
  labs(
    x = "Maximum log-likelihood improvement over constant hazard",
    y = "Number of randomized datasets",
    title = "Uniform-date randomization test"
  ) +
  theme_minimal(base_size = 13)

print(uniform_plot)

# ============================================================
# 8. TEST 2B: CIRCULAR-SHIFT RANDOMIZATION TEST
# ============================================================

circular_shift_dates <- function(
    t_events,
    shift_days,
    times
) {

  start_time <- min(times)
  study_duration <- length(times)

  (
    (t_events - start_time + shift_days) %%
      study_duration
  ) +
    start_time
}

circular_statistics <- numeric(number_randomizations)
circular_shifts <- numeric(number_randomizations)

for (permutation_index in seq_len(number_randomizations)) {

  random_shift <- sample(
    1:(length(times_full) - 1),
    size = 1
  )

  circular_shifts[permutation_index] <- random_shift

  shifted_dates <- circular_shift_dates(
    t_events = observed_events,
    shift_days = random_shift,
    times = times_full
  )

  circular_statistics[permutation_index] <-
    max_lag_statistic(
      t_events = shifted_dates,
      Iw = Iw_vec,
      Ir = Ir_vec,
      times = times_full,
      lag_values = lag_values,
      alpha_upper = alpha_upper
    )

  if (permutation_index %% 100 == 0) {
    cat(
      "Completed",
      permutation_index,
      "of",
      number_randomizations,
      "circular shifts\n"
    )
  }
}

circular_p_value <- (
  1 +
  sum(circular_statistics >= observed_statistic)
) / (
  number_randomizations + 1
)

cat("\n==============================================\n")
cat("CIRCULAR-SHIFT RANDOMIZATION RESULTS\n")
cat("==============================================\n")
cat("Observed statistic:", observed_statistic, "\n")
cat("Circular-shift p-value:", circular_p_value, "\n")

circular_results <- data.frame(
  statistic = circular_statistics,
  shift_days = circular_shifts,
  randomization = "Circular shift"
)

circular_plot <- ggplot(
  circular_results,
  aes(x = statistic)
) +
  geom_histogram(
    bins = 35,
    fill = "#A1D99B",
    color = "white"
  ) +
  geom_vline(
    xintercept = observed_statistic,
    color = "#B2182B",
    linewidth = 1.1
  ) +
  annotate(
    "text",
    x = observed_statistic,
    y = Inf,
    label = paste0(
      "Observed\np = ",
      format.pval(circular_p_value, digits = 3)
    ),
    hjust = ifelse(
      observed_statistic >
        median(circular_statistics),
      1.05,
      -0.05
    ),
    vjust = 1.2,
    color = "#B2182B"
  ) +
  labs(
    x = "Maximum log-likelihood improvement over constant hazard",
    y = "Number of circularly shifted datasets",
    title = "Circular-shift randomization test"
  ) +
  theme_minimal(base_size = 13)

print(circular_plot)

# ============================================================
# 9. COMBINED RANDOMIZATION FIGURE
# ============================================================

combined_randomization_results <- bind_rows(
  uniform_results,
  circular_results %>%
    select(statistic, randomization)
)

combined_plot <- ggplot(
  combined_randomization_results,
  aes(x = statistic, fill = randomization)
) +
  geom_histogram(
    bins = 35,
    color = "white",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = observed_statistic,
    color = "#B2182B",
    linewidth = 1.1
  ) +
  facet_wrap(
    ~randomization,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Uniform dates" = "#9ECAE1",
      "Circular shift" = "#A1D99B"
    ),
    guide = "none"
  ) +
  labs(
    x = "Maximum log-likelihood improvement over constant hazard",
    y = "Number of randomized datasets",
    caption = paste0(
      "Observed statistic = ",
      round(observed_statistic, 3),
      "; uniform p = ",
      format.pval(uniform_p_value, digits = 3),
      "; circular-shift p = ",
      format.pval(circular_p_value, digits = 3),
      "."
    )
  ) +
  theme_minimal(base_size = 13)

print(combined_plot)

# ============================================================
# 10. SAVE RESULTS
# ============================================================

summary_results <- data.frame(
  best_lag_days = best_lag_result$lag,
  best_lag_logLik = best_lag_result$logLik,
  zero_lag_logLik = zero_lag_result$logLik,
  best_vs_zero_logLik =
    best_lag_result$logLik -
    zero_lag_result$logLik,
  best_beta_spill = best_lag_result$beta_spill,
  best_alpha = best_lag_result$alpha,
  observed_improvement_over_constant = observed_statistic,
  uniform_randomization_p = uniform_p_value,
  circular_shift_p = circular_p_value,
  number_randomizations = number_randomizations
)

write.csv(
  introduction_data,
  "reconstructed_farm_introduction_dates.csv",
  row.names = FALSE
)

write.csv(
  lag_results,
  "spillover_lag_profile_results.csv",
  row.names = FALSE
)

write.csv(
  uniform_results,
  "spillover_uniform_randomization_results.csv",
  row.names = FALSE
)

write.csv(
  circular_results,
  "spillover_circular_shift_results.csv",
  row.names = FALSE
)

write.csv(
  summary_results,
  "spillover_temporal_test_summary.csv",
  row.names = FALSE
)

ggsave(
  filename = "spillover_lag_profile.png",
  plot = lag_plot,
  width = 8,
  height = 5.2,
  units = "in",
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = "spillover_uniform_randomization.png",
  plot = uniform_plot,
  width = 7.5,
  height = 5.2,
  units = "in",
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = "spillover_circular_shift.png",
  plot = circular_plot,
  width = 7.5,
  height = 5.2,
  units = "in",
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = "spillover_randomization_combined.png",
  plot = combined_plot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 320,
  bg = "white"
)

cat("\n==============================================\n")
cat("FINAL SUMMARY\n")
cat("==============================================\n")
print(summary_results)

cat("\nSaved files:\n")
cat("  reconstructed_farm_introduction_dates.csv\n")
cat("  spillover_lag_profile_results.csv\n")
cat("  spillover_uniform_randomization_results.csv\n")
cat("  spillover_circular_shift_results.csv\n")
cat("  spillover_temporal_test_summary.csv\n")
cat("  spillover_lag_profile.png\n")
cat("  spillover_uniform_randomization.png\n")
cat("  spillover_circular_shift.png\n")
cat("  spillover_randomization_combined.png\n")

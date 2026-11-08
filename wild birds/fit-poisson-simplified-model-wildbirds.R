# =============================================================================
# HPAI in Ontario wild birds, 2024-11-01 to 2025-03-25
# Seasonal migration model, waterfowl and raptors
# Poisson observation model, two-chain random-walk Metropolis
#
# Natural mortality (mu) is not in the model. The only removal by death is
# the fixed disease-induced rate d.
#
# This file is the analysis used for the manuscript, reconstructed as one
# top-to-bottom script from the original interactive session. Seeds, proposal
# scales, burn-in, and which draws are reported match that session.
#
# Run from this directory (the folder that contains data/):
#   Rscript fit_poisson_simplified_waterfowl_raptor.R
#
# The production lengths are long: every proposal solves an ODE, and the
# recorded chains are 20,000 + 100,000 x 50 + pilots + 50,000 + 150,000
# + 150,000 iterations, twice. Expect days, not minutes.
#
# Smoke test (checks that data, the ODE, and the MCMC plumbing run; the
# draws are NOT manuscript results):
#   HPAI_SMOKE_TEST=1 Rscript fit_poisson_simplified_waterfowl_raptor.R
#
# Reported posterior
#   The manuscript posterior is the LAST production block only
#   (chain1_c, chain2_c), after discarding the first 20% of that block.
#   That is the original sample5 file. Earlier blocks are kept so the
#   tuning path can be audited; they are not pooled into the reported draws.
#
# Parameter scale
#   The sampler runs on scaled parameters theta. The ODE uses
#     S_w(0)       = theta[1] * 100
#     I_w(0)       = theta[2]
#     beta         = theta[3] / 100
#     Lambda_amp   = theta[4] * 10
#     delta_amp    = theta[5]
#     S_r(0)       = theta[6] * 100
#     beta_wr      = theta[7] / 100
#   posterior_samples and print(post_mean) are on the theta scale, which is
#   what the original session printed. The ODE-scale table is written beside
#   it so the two are not mixed up in the paper.
#
# Original local paths, now replaced by data/ and results/ next to this file:
#   ~/Desktop/HPAI-MODEL/Ontario-wild-birds-waterfowl-daily_cases.csv
#   ~/Desktop/HPAI-MODEL/Ontario-wild-birds-raptor-daily_cases.csv
#   ~/Desktop/HPAI-MODEL-REVISION/wild-birds-more chain-results/
# =============================================================================

# ----------------------- user settings ---------------------------------------
# Production settings below are the manuscript run. HPAI_SMOKE_TEST=1
# overrides SMOKE_TEST without editing this file.
SMOKE_TEST <- FALSE
if (identical(Sys.getenv("HPAI_SMOKE_TEST"), "1")) SMOKE_TEST <- TRUE

# ggmcmc on the full final chains is what the original session wrote.
# Set FALSE if that PDF step runs out of memory. Rhat and ESS do not use it.
SAVE_GGMCMC <- TRUE

START_DATE <- as.Date("2024-11-01")
END_DATE   <- as.Date("2025-03-25")

# Immigration is on through day 60. Departure is on only after day 120.
# t = 1 is the first waterfowl date in the window (intended: 2024-11-01).
T_ARRIVAL_END  <- 60
T_DEPART_START <- 120

# Disease-induced mortality only. Fixed, not estimated.
# Mean time to disease death is 1/d if the bird does not depart.
fixed <- c(d = 0.3)

BURN_FRAC <- 0.2
FLOOR_LAMBDA <- 1e-6          # floor on non-positive Poisson means
SCALAR_SCALE <- 0.05          # phase-1 proposal sd (independent components)
SCALE_MULTIPLIER <- 0.35      # production multiplier on the covariance factor
PILOT_MULTIPLIERS <- c(0.7, 0.35)
N_COV_DRAWS <- 4000           # thinned phase-1 draws used to build Sigma

if (SMOKE_TEST) {
  message("SMOKE_TEST is on. Chains are shortened. Do not use these draws.")
  N_ADAPT <- 40
  N_SAMPLE <- 80
  BLEN_SAMPLE <- 2
  N_PILOT <- 20
  BLEN_PILOT <- 2
  N_BLOCK1 <- 80
  N_BLOCK2 <- 80
  N_BLOCK3 <- 80
} else {
  N_ADAPT <- 20000
  N_SAMPLE <- 100000
  BLEN_SAMPLE <- 50
  N_PILOT <- 2000
  BLEN_PILOT <- 50
  N_BLOCK1 <- 50000   # original sample3
  N_BLOCK2 <- 150000  # original sample4
  N_BLOCK3 <- 150000  # original sample5, this is the reported posterior
}
# -----------------------------------------------------------------------------

needed <- c(
  "deSolve", "mcmc", "coda", "ggmcmc",
  "dplyr", "lubridate", "ggplot2", "posterior"
)
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Install required packages first:\n  install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "),
    "))",
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(deSolve)
  library(mcmc)
  library(coda)
  library(ggmcmc)
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(posterior)
})

script_dir <- (function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", hit))))
  }
  getwd()
})()

data_dir <- file.path(script_dir, "data")
out_dir  <- file.path(script_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(out_dir, "run_log.txt")
cat("", file = log_file)  # truncate
log_msg <- function(...) {
  txt <- paste0(..., collapse = "")
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_file, append = TRUE)
}
log_obj <- function(label, x, digits = NULL) {
  log_msg(label)
  printed <- if (is.null(digits)) {
    capture.output(print(x))
  } else {
    capture.output(print(x, digits = digits))
  }
  log_msg(paste(printed, collapse = "\n"))
}

if (!interactive()) {
  options(error = function() {
    msg <- geterrmessage()
    cat(msg, file = log_file, append = TRUE)
    cat(msg, file = stderr())
    quit(status = 1)
  })
}
options(warn = 1)

log_msg("R: ", R.version.string)
log_msg("RNG: ", paste(RNGkind(), collapse = " / "))
for (p in needed) {
  log_msg(sprintf("Package %s %s", p, as.character(utils::packageVersion(p))))
}
log_msg("SMOKE_TEST = ", SMOKE_TEST)
log_msg("Output directory: ", out_dir)

# ----------------------------- data ------------------------------------------
# Both files need columns CollectionDate and DailyCases.
# The likelihood matches the two series row by row, so the dates must agree.
# Days that are absent are NOT filled with zero here: that would change the
# likelihood relative to the fitted model. If the warning below fires, stop
# and decide whether missing dates are real gaps or unreported zeros.

waterfowl_file <- file.path(
  data_dir, "Ontario-wild-birds-waterfowl-daily_cases.csv"
)
raptor_file <- file.path(
  data_dir, "Ontario-wild-birds-raptor-daily_cases.csv"
)
if (!file.exists(waterfowl_file) || !file.exists(raptor_file)) {
  stop(
    "Case files not found. Expected:\n  ",
    waterfowl_file, "\n  ", raptor_file, "\n",
    "Put both CSVs in data/. Columns: CollectionDate, DailyCases.",
    call. = FALSE
  )
}

parse_dates <- function(x) {
  x <- trimws(as.character(x))
  d <- suppressWarnings(as.Date(x))
  if (anyNA(d)) {
    d2 <- suppressWarnings(lubridate::ymd(x, quiet = TRUE))
    d[is.na(d)] <- d2[is.na(d)]
  }
  if (anyNA(d)) {
    d2 <- suppressWarnings(lubridate::mdy(x, quiet = TRUE))
    d[is.na(d)] <- d2[is.na(d)]
  }
  if (anyNA(d)) {
    stop("Could not parse CollectionDate. Use YYYY-MM-DD.", call. = FALSE)
  }
  d
}

read_cases <- function(path) {
  raw <- utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE)
  need <- c("CollectionDate", "DailyCases")
  if (!all(need %in% names(raw))) {
    stop(path, " must contain columns: ", paste(need, collapse = ", "),
         call. = FALSE)
  }
  out <- raw[, need, drop = FALSE]
  out$CollectionDate <- parse_dates(out$CollectionDate)
  out$DailyCases <- as.numeric(out$DailyCases)
  if (anyNA(out$DailyCases)) {
    stop(path, ": DailyCases has non-numeric values.", call. = FALSE)
  }
  if (any(out$DailyCases < 0 | out$DailyCases != floor(out$DailyCases))) {
    stop(path, ": DailyCases must be non-negative integers.", call. = FALSE)
  }
  if (anyDuplicated(out$CollectionDate)) {
    warning(path, ": duplicate dates; summing DailyCases within each date.")
    out <- out %>%
      group_by(CollectionDate) %>%
      summarise(DailyCases = sum(DailyCases), .groups = "drop")
  }
  out %>%
    filter(CollectionDate >= START_DATE, CollectionDate <= END_DATE) %>%
    arrange(CollectionDate)
}

data_filtered  <- read_cases(waterfowl_file)
data2_filtered <- read_cases(raptor_file)

if (nrow(data_filtered) < 3 || nrow(data2_filtered) < 3) {
  stop("Need at least 3 observation days in each series.", call. = FALSE)
}
if (!identical(data_filtered$CollectionDate, data2_filtered$CollectionDate)) {
  stop(
    "Waterfowl and raptor files do not have the same CollectionDate values ",
    "inside ", START_DATE, " to ", END_DATE, ". The likelihood pairs them ",
    "row by row. Align the dates (put 0 on days with no detections) and rerun.",
    call. = FALSE
  )
}

origin <- min(data_filtered$CollectionDate)
data_filtered$t  <- as.integer(data_filtered$CollectionDate - origin) + 1L
data2_filtered$t <- data_filtered$t
times <- data_filtered$t

if (origin != START_DATE) {
  warning(
    "First included date is ", origin, ", not ", START_DATE,
    ". Arrival/departure switches are t = ", T_ARRIVAL_END,
    " and t = ", T_DEPART_START, " on this clock."
  )
}
if (any(diff(data_filtered$CollectionDate) != 1)) {
  warning(
    "Dates are not a complete daily sequence. The original likelihood still ",
    "matches diff(T) across each gap to that row's single DailyCases value."
  )
}

log_msg(
  "Observations: ", nrow(data_filtered), " days, ",
  origin, " to ", max(data_filtered$CollectionDate),
  " (t = ", min(times), " to ", max(times), ")."
)
log_msg(
  "Case totals: waterfowl ", sum(data_filtered$DailyCases),
  ", raptor ", sum(data2_filtered$DailyCases), "."
)

# ----------------------------- ODE -------------------------------------------
# State order: S_w, I_w, T_w, S_r, I_r, T_r.
# T_w and T_r are cumulative disease-induced deaths. They are what the
# Poisson mean increments. An earlier draft named T_w as T_d; same state.
#
#   dS_w/dt = Lambda(t) - beta * S_w * I_w - delta(t) * S_w
#   dI_w/dt = beta * S_w * I_w - (d + delta(t)) * I_w
#   dT_w/dt = d * I_w
#
#   dS_r/dt = -beta_wr * S_r * I_w
#   dI_r/dt =  beta_wr * S_r * I_w - d * I_r
#   dT_r/dt = d * I_r
#
#   Lambda(t) = Lambda_amp  if t <= 60, else 0
#   delta(t)  = 0           if t <= 120, else delta_amp
#
# Force of infection is mass action (beta * S * I), not frequency dependent.
# Infected waterfowl leave by disease death or by departure. Only disease
# death is counted in T. Raptors do not migrate: no Lambda and no delta.
#
# Lambda and delta jump at t = 60 and t = 120. The fitted model integrates
# across those jumps in one deSolve::ode call (lsoda defaults). Splitting
# the integral at the jumps would be stabler, but it is not what was run.

hpai_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    Lambda_now <- if (time <= T_ARRIVAL_END) Lambda_amp else 0
    delta_now  <- if (time <= T_DEPART_START) 0 else delta_amp
    
    dS_w <- Lambda_now - beta * S_w * I_w - delta_now * S_w
    dI_w <- beta * S_w * I_w - (d + delta_now) * I_w
    dT_w <- d * I_w
    
    dS_r <- -beta_wr * S_r * I_w
    dI_r <-  beta_wr * S_r * I_w - d * I_r
    dT_r <- d * I_r
    
    list(c(dS_w, dI_w, dT_w, dS_r, dI_r, dT_r))
  })
}

hpai_model <- function(params, times, fixed) {
  if (length(params) != 7) stop("hpai_model expects 7 parameters.")
  
  y0 <- c(
    S_w = params[[1]] * 100,
    I_w = params[[2]],
    T_w = 0,
    S_r = params[[6]] * 100,
    I_r = 0,
    T_r = 0
  )
  parms <- c(
    beta       = params[[3]] / 100,
    Lambda_amp = params[[4]] * 10,
    delta_amp  = params[[5]],
    beta_wr    = params[[7]] / 100,
    fixed
  )
  ode(y = y0, times = times, func = hpai_equations, parms = parms)
}

# ----------------------------- prior / likelihood ----------------------------
# Priors are independent normals on the scaled parameters theta:
#   S_w0          N(3,    1.5^2)     ODE initial susceptible waterfowl / 100
#   I_w0          N(1.5,  0.5^2)     ODE initial infected waterfowl
#   beta_w        N(0.5,  0.2^2)     ODE beta * 100
#   Lambda_amp_w  N(2,    1^2)       ODE immigration rate / 10
#   delta_amp_w   N(0.01, 0.005^2)   ODE departure rate
#   S_r0          N(1,    0.5^2)     ODE initial susceptible raptors / 100
#   beta_wr       N(0.5,  0.2^2)     ODE beta_wr * 100
#
# theta < 0 has posterior density zero. The normal log-density is not
# renormalized on [0, Inf). That is the prior that was fitted. Do not
# replace it with a truncated-normal density without refitting.

par_names <- c(
  "S_w0", "I_w0", "beta_w", "Lambda_amp_w", "delta_amp_w",
  "S_r0", "beta_wr"
)

prior <- function(parameters) {
  if (length(parameters) != 7 || any(!is.finite(parameters)) ||
      any(parameters < 0)) {
    return(-Inf)
  }
  sum(
    dnorm(parameters[[1]], mean = 3,    sd = 1.5,   log = TRUE),
    dnorm(parameters[[2]], mean = 1.5,  sd = 0.5,   log = TRUE),
    dnorm(parameters[[3]], mean = 0.5,  sd = 0.2,   log = TRUE),
    dnorm(parameters[[4]], mean = 2,    sd = 1,     log = TRUE),
    dnorm(parameters[[5]], mean = 0.01, sd = 0.005, log = TRUE),
    dnorm(parameters[[6]], mean = 1,    sd = 0.5,   log = TRUE),
    dnorm(parameters[[7]], mean = 0.5,  sd = 0.2,   log = TRUE)
  )
}

# Daily counts, dropping the first day, against increments of cumulative
# disease deaths. Non-positive or non-finite means are floored at 1e-6
# so a numerical undershoot does not return NaN. A structural zero mean
# would otherwise make any positive count impossible.

likelihood <- function(parameters, fixed, data, data2, times) {
  if (length(parameters) != 7 || any(!is.finite(parameters)) ||
      any(parameters < 0)) {
    return(-Inf)
  }
  
  # Solver failures are rejected. Without this, one stiff proposal aborts a
  # multi-day run. It does not change the log posterior at parameters the
  # solver can integrate.
  model_output <- tryCatch(
    hpai_model(params = parameters, times = times, fixed = fixed),
    error = function(e) NULL
  )
  if (is.null(model_output) || any(!is.finite(model_output))) return(-Inf)
  
  observations1 <- data$DailyCases[-1]
  predictions1  <- diff(model_output[, "T_w"])
  observations2 <- data2$DailyCases[-1]
  predictions2  <- diff(model_output[, "T_r"])
  
  if (length(predictions1) != length(observations1) ||
      length(predictions2) != length(observations2)) {
    return(-Inf)
  }
  
  predictions1[!is.finite(predictions1) | predictions1 <= 0] <- FLOOR_LAMBDA
  predictions2[!is.finite(predictions2) | predictions2 <= 0] <- FLOOR_LAMBDA
  
  d1 <- dpois(observations1, lambda = predictions1, log = TRUE)
  d2 <- dpois(observations2, lambda = predictions2, log = TRUE)
  if (any(!is.finite(d1)) || any(!is.finite(d2))) return(-Inf)
  sum(d1) + sum(d2)
}

joint_log_prob <- function(parameters, fixed, data, data2, times) {
  lp <- prior(parameters)
  if (!is.finite(lp)) return(-Inf)
  ll <- likelihood(parameters, fixed, data, data2, times)
  if (!is.finite(ll)) return(-Inf)
  lp + ll
}

# ----------------------------- MCMC helpers ----------------------------------
# mcmc::metrop is fixed-scale random-walk Metropolis, not adaptive Metropolis.
# When obj is a previous metrop fit, the new run starts at that fit's end
# state and does not keep the old batches. Continuation also copies
# obj$final.seed back into .Random.seed, so a set.seed() on the line above
# is overwritten. Only set.seed(1) and set.seed(2), which start a chain
# from the log-posterior function rather than from a previous fit, actually
# initialize the random streams. Later set.seed calls are kept so this file
# matches the original session; they are not extra knobs.
# Proposal: current + scale %*% rnorm(7). A 7 x 7 scale has proposal
# covariance scale %*% t(scale).

metrop_call <- function(obj, nbatch, blen, scale, initial = NULL) {
  args <- list(
    obj    = obj,
    nbatch = nbatch,
    blen   = blen,
    scale  = scale,
    data   = data_filtered,
    data2  = data2_filtered,
    fixed  = fixed,
    times  = times
  )
  if (!is.null(initial)) args$initial <- initial
  do.call(mcmc::metrop, args)
}

drop_burnin <- function(batch, burn_frac = BURN_FRAC) {
  n <- nrow(batch)
  batch[(floor(burn_frac * n) + 1):n, , drop = FALSE]
}

name_draws <- function(batch) {
  colnames(batch) <- par_names
  batch
}

# Phase-1 batch means -> proposal factor.
# Sigma_prop = (2.38^2 / d) * (0.95 * cov + 0.05 * diagonal),
# the Roberts-Rosenthal scaling, with a diagonal shrink so chol() is stable.
# metrop wants the matrix square root, not Sigma itself:
#   scale_mat <- t(chol(Sigma_prop))
# Production then uses 0.35 * scale_mat. That 0.35 was chosen from the
# pilots below and is fixed during every production draw. Because it is
# not updated inside the production chains, those chains are ordinary
# random-walk Metropolis.

proposal_factor <- function(draws) {
  draws <- as.matrix(draws)
  n_keep <- min(N_COV_DRAWS, nrow(draws))
  post <- draws[seq(1, nrow(draws), length.out = n_keep), , drop = FALSE]
  Sigma <- stats::cov(post)
  Sigma <- 0.95 * Sigma + 0.05 * diag(diag(Sigma))
  d <- ncol(Sigma)
  Sigma_prop <- (2.38^2 / d) * Sigma
  factor <- tryCatch(t(chol(Sigma_prop)), error = function(e) NULL)
  if (is.null(factor)) {
    if (!SMOKE_TEST) {
      stop(
        "chol() failed on the phase-1 covariance. The scalar-scale chains ",
        "may not have moved (check their acceptance rates) ",
        "or a parameter margin has zero variance.",
        call. = FALSE
      )
    }
    warning("chol() failed in the smoke test; using a diagonal proposal.")
    factor <- diag(SCALAR_SCALE, d)
  }
  factor
}

diagnose_pair <- function(s1, s2, label) {
  stopifnot(identical(dim(s1), dim(s2)), identical(colnames(s1), par_names))
  mcmc_list <- mcmc.list(as.mcmc(s1), as.mcmc(s2))
  
  # Smoke-test chains are deliberately too short for some coda diagnostics.
  # A failure there is reported and skipped. A failure on the real run stops.
  run_diag <- function(fun, what) {
    tryCatch(fun(), error = function(e) {
      if (!SMOKE_TEST) stop(e)
      warning(what, " failed in the smoke test: ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
  }
  
  log_msg("\n===== ", label, " =====")
  gelman <- run_diag(function() gelman.diag(mcmc_list, autoburnin = FALSE),
                     "gelman.diag")
  ess <- run_diag(function() effectiveSize(mcmc_list), "effectiveSize")
  gew1 <- run_diag(function() geweke.diag(as.mcmc(s1)), "geweke chain 1")
  gew2 <- run_diag(function() geweke.diag(as.mcmc(s2)), "geweke chain 2")
  if (!is.null(gelman)) log_obj("Gelman-Rubin (Rhat; want this near 1):", gelman)
  if (!is.null(ess)) log_obj("Effective sample size (coda):", ess)
  if (!is.null(gew1)) log_obj("Geweke z (chain 1):", gew1)
  if (!is.null(gew2)) log_obj("Geweke z (chain 2):", gew2)
  
  posterior_samples <- rbind(s1, s2)
  post_mean <- colMeans(posterior_samples)
  ci <- apply(posterior_samples, 2, quantile, probs = c(0.025, 0.975))
  log_obj("Posterior mean (scaled parameters):", post_mean)
  log_obj("95% credible intervals (scaled parameters):", ci)
  
  a <- array(
    NA_real_,
    dim = c(nrow(s1), 2, ncol(s1)),
    dimnames = list(
      iteration = NULL,
      chain = c("chain1", "chain2"),
      variable = colnames(s1)
    )
  )
  a[, 1, ] <- s1
  a[, 2, ] <- s2
  diagnostics <- run_diag(function() {
    posterior::summarise_draws(
      posterior::as_draws_array(a),
      "mean", "sd", "rhat",
      "ess_bulk", "ess_tail", "mcse_mean"
    )
  }, "summarise_draws")
  if (!is.null(diagnostics)) {
    log_obj("posterior::summarise_draws:", as.data.frame(diagnostics), digits = 6)
  }
  
  list(
    mcmc_list = mcmc_list,
    posterior_samples = posterior_samples,
    post_mean = post_mean,
    ci = ci,
    diagnostics = diagnostics
  )
}

save_block <- function(path, object_names) {
  save(list = object_names, file = path, envir = parent.frame())
  log_msg("Saved ", path)
}

# Initial values are the ones used in the session, not the prior means.
# Names are the scaled parameters. The 7th initial was called beta_wrr
# in the original vector; it is beta_wr. Same number, same position.
initial_chain1 <- c(1.26, 0.3, 0.3, 0.07, 0.006, 2, 0.17)
initial_chain2 <- c(2.5, 0.8, 0.6, 1.5, 0.012, 1.2, 0.35)
names(initial_chain1) <- par_names
names(initial_chain2) <- par_names

lp1 <- joint_log_prob(
  initial_chain1, fixed, data_filtered, data2_filtered, times
)
lp2 <- joint_log_prob(
  initial_chain2, fixed, data_filtered, data2_filtered, times
)
log_obj("Log posterior at the two initial points:", c(chain1 = lp1, chain2 = lp2))
if (!is.finite(lp1) || !is.finite(lp2)) {
  stop(
    "Initial log posterior is not finite. Check the case files and dates.",
    call. = FALSE
  )
}

# ===================== Phase 1: scalar proposal ==============================
# Warm-up is fixed-scale (scale = 0.05), not adaptive. One seed covers the
# warm-up and the following sample, matching the original run_chain().
# The sample stores batch means of length BLEN_SAMPLE (50 in production),
# not every iteration. Those batch means are used only to build a proposal
# and to see that the scalar scale mixes poorly. They are not the reported
# posterior.

log_msg("\nPhase 1: scalar random walk, scale = ", SCALAR_SCALE)

set.seed(1)
warmup1 <- metrop_call(
  obj = joint_log_prob, nbatch = N_ADAPT, blen = 1,
  scale = SCALAR_SCALE, initial = initial_chain1
)
chain1 <- metrop_call(
  obj = warmup1, nbatch = N_SAMPLE, blen = BLEN_SAMPLE, scale = SCALAR_SCALE
)

set.seed(2)
warmup2 <- metrop_call(
  obj = joint_log_prob, nbatch = N_ADAPT, blen = 1,
  scale = SCALAR_SCALE, initial = initial_chain2
)
chain2 <- metrop_call(
  obj = warmup2, nbatch = N_SAMPLE, blen = BLEN_SAMPLE, scale = SCALAR_SCALE
)

log_obj(
  "Phase-1 acceptance (warm-up, then batch-mean sample; target was about 0.2-0.4):",
  c(
    warmup1 = warmup1$accept, chain1 = chain1$accept,
    warmup2 = warmup2$accept, chain2 = chain2$accept
  )
)

s1 <- name_draws(drop_burnin(chain1$batch))
s2 <- name_draws(drop_burnin(chain2$batch))
phase1 <- diagnose_pair(s1, s2, "Phase 1, after 20% burn-in of the batch means")
posterior_samples <- phase1$posterior_samples
post_mean <- phase1$post_mean

save_block(
  file.path(out_dir, "sample1_poisson-simplified-waterfowl-raptor.RData"),
  c("chain1", "chain2", "posterior_samples", "post_mean")
)

# ===================== Phase 2: covariance proposal ==========================
scale_mat <- proposal_factor(rbind(s1, s2))
log_msg("Built ", nrow(scale_mat), " x ", ncol(scale_mat), " proposal factor.")

# Pilots start from the end of chain 1 and are then discarded.
# They exist to confirm the acceptance rate. Production does not start
# from the pilot state; it starts again from chain1 and chain2.
# Target acceptance on these pilots was about 0.15-0.30.

set.seed(31)
pilot <- metrop_call(
  obj = chain1, nbatch = N_PILOT, blen = BLEN_PILOT,
  scale = PILOT_MULTIPLIERS[1] * scale_mat
)
set.seed(32)
pilot2 <- metrop_call(
  obj = pilot, nbatch = N_PILOT, blen = BLEN_PILOT,
  scale = PILOT_MULTIPLIERS[2] * scale_mat
)
log_obj(
  "Pilot acceptance:",
  c(scale_0.70 = pilot$accept, scale_0.35 = pilot2$accept)
)

scale_use <- SCALE_MULTIPLIER * scale_mat
log_msg(
  "Production proposal is ", SCALE_MULTIPLIER,
  " * scale_mat (chosen from the pilots and then frozen)."
)

# ===================== Production block 1 (sample3) ==========================
# blen = 1 from here on: every iteration is stored, not a batch mean.

set.seed(41)
chain1_new <- metrop_call(
  obj = chain1, nbatch = N_BLOCK1, blen = 1, scale = scale_use
)
set.seed(42)
chain2_new <- metrop_call(
  obj = chain2, nbatch = N_BLOCK1, blen = 1, scale = scale_use
)
log_obj(
  "Block 1 acceptance:",
  c(chain1 = chain1_new$accept, chain2 = chain2_new$accept)
)

s1 <- name_draws(drop_burnin(chain1_new$batch))
s2 <- name_draws(drop_burnin(chain2_new$batch))
block1 <- diagnose_pair(s1, s2, "Production block 1 (sample3), 20% burn-in")
posterior_samples <- block1$posterior_samples
post_mean <- block1$post_mean
ci <- block1$ci

save_block(
  file.path(out_dir, "sample3_poisson-simplified-waterfowl-raptor.RData"),
  c("chain1_new", "chain2_new", "s1", "s2",
    "posterior_samples", "post_mean", "ci", "scale_use")
)

# ===================== Production block 2 (sample4) ==========================
set.seed(61)
chain1_b <- metrop_call(
  obj = chain1_new, nbatch = N_BLOCK2, blen = 1, scale = scale_use
)
set.seed(62)
chain2_b <- metrop_call(
  obj = chain2_new, nbatch = N_BLOCK2, blen = 1, scale = scale_use
)
log_obj(
  "Block 2 acceptance:",
  c(chain1 = chain1_b$accept, chain2 = chain2_b$accept)
)

s1 <- name_draws(drop_burnin(chain1_b$batch))
s2 <- name_draws(drop_burnin(chain2_b$batch))
block2 <- diagnose_pair(s1, s2, "Production block 2 (sample4), 20% burn-in")
posterior_samples <- block2$posterior_samples
post_mean <- block2$post_mean
ci <- block2$ci

save_block(
  file.path(out_dir, "sample4_poisson-simplified-waterfowl-raptor.RData"),
  c("chain1_b", "chain2_b", "s1", "s2",
    "posterior_samples", "post_mean", "ci", "scale_use")
)

# ===================== Production block 3 (sample5) ==========================
# Seeds 61 and 62 are set again because the original session did that.
# They do not replay block 2 and they do not reset the proposal noise:
# metrop restores .Random.seed from chain*_b$final.seed before drawing.
# See the note above metrop_call.

set.seed(61)
chain1_c <- metrop_call(
  obj = chain1_b, nbatch = N_BLOCK3, blen = 1, scale = scale_use
)
set.seed(62)
chain2_c <- metrop_call(
  obj = chain2_b, nbatch = N_BLOCK3, blen = 1, scale = scale_use
)
log_obj(
  "Block 3 acceptance (reported run):",
  c(chain1 = chain1_c$accept, chain2 = chain2_c$accept)
)

s1 <- name_draws(drop_burnin(chain1_c$batch))
s2 <- name_draws(drop_burnin(chain2_c$batch))
block3 <- diagnose_pair(s1, s2, "Production block 3 (sample5), REPORTED posterior")
posterior_samples <- block3$posterior_samples
post_mean <- block3$post_mean
ci <- block3$ci
mcmc_list <- block3$mcmc_list
diagnostics <- block3$diagnostics

accept_rates <- c(
  phase1_warmup1 = warmup1$accept,
  phase1_warmup2 = warmup2$accept,
  phase1_chain1  = chain1$accept,
  phase1_chain2  = chain2$accept,
  pilot_0.70     = pilot$accept,
  pilot_0.35     = pilot2$accept,
  block1_chain1  = chain1_new$accept,
  block1_chain2  = chain2_new$accept,
  block2_chain1  = chain1_b$accept,
  block2_chain2  = chain2_b$accept,
  block3_chain1  = chain1_c$accept,
  block3_chain2  = chain2_c$accept
)
log_obj("Acceptance rates, all stages:", accept_rates)

save_block(
  file.path(out_dir, "sample5_poisson-simplified-waterfowl-raptor.RData"),
  c("chain1_c", "chain2_c", "s1", "s2",
    "posterior_samples", "post_mean", "ci", "scale_use",
    "accept_rates", "diagnostics", "scale_mat")
)

# Scaled theta and ODE-scale summaries. Quantiles transform because every
# map from theta to the ODE parameter is linear and increasing.
ode_multiplier <- c(
  S_w0 = 100, I_w0 = 1, beta_w = 1 / 100, Lambda_amp_w = 10,
  delta_amp_w = 1, S_r0 = 100, beta_wr = 1 / 100
)
ode_parameter <- c(
  S_w0 = "S_w0", I_w0 = "I_w0", beta_w = "beta",
  Lambda_amp_w = "Lambda_amp", delta_amp_w = "delta_amp",
  S_r0 = "S_r0", beta_wr = "beta_wr"
)

paper_table <- data.frame(
  variable = par_names,
  ode_parameter = unname(ode_parameter[par_names]),
  ode_multiplier = unname(ode_multiplier[par_names]),
  mean = unname(post_mean[par_names]),
  q025_scaled = unname(ci["2.5%", par_names]),
  q975_scaled = unname(ci["97.5%", par_names]),
  stringsAsFactors = FALSE
)
paper_table$mean_ode <- paper_table$mean * paper_table$ode_multiplier
paper_table$q025_ode <- paper_table$q025_scaled * paper_table$ode_multiplier
paper_table$q975_ode <- paper_table$q975_scaled * paper_table$ode_multiplier

if (!is.null(diagnostics)) {
  diag_df <- as.data.frame(diagnostics, stringsAsFactors = FALSE)
  # posterior::summarise_draws names the parameter column "variable".
  # ess_* columns can be vctrs/pillar types; write.csv needs plain numbers.
  if (!"variable" %in% names(diag_df) && ncol(diag_df) >= 1) {
    names(diag_df)[1] <- "variable"
  }
  diag_df$variable <- as.character(diag_df$variable)
  for (nm in setdiff(names(diag_df), "variable")) {
    diag_df[[nm]] <- suppressWarnings(as.numeric(diag_df[[nm]]))
  }
  keep_cols <- intersect(
    c("variable", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"),
    names(diag_df)
  )
  paper_table <- merge(paper_table, diag_df[, keep_cols, drop = FALSE],
                       by = "variable", all.x = TRUE, sort = FALSE)
  paper_table <- paper_table[match(par_names, paper_table$variable), ]
  if (!SMOKE_TEST && "rhat" %in% names(paper_table) &&
      any(paper_table$rhat > 1.01, na.rm = TRUE)) {
    warning(
      "At least one Rhat is above 1.01. Read results/run_log.txt before ",
      "putting these numbers in the manuscript."
    )
  }
}

d_row <- paper_table[1, ]
for (nm in names(d_row)) d_row[[nm]] <- NA
d_row$variable <- "d"
d_row$ode_parameter <- "d"
d_row$mean <- unname(fixed[["d"]])
d_row$q025_scaled <- unname(fixed[["d"]])
d_row$q975_scaled <- unname(fixed[["d"]])
d_row$mean_ode <- unname(fixed[["d"]])
d_row$q025_ode <- unname(fixed[["d"]])
d_row$q975_ode <- unname(fixed[["d"]])
if ("sd" %in% names(d_row)) d_row$sd <- 0
paper_table <- rbind(paper_table, d_row)

summary_path <- file.path(out_dir, "posterior_summary.csv")
utils::write.csv(paper_table, summary_path, row.names = FALSE)
log_msg("Saved ", summary_path)
scaled_cols <- intersect(
  c("variable", "mean", "q025_scaled", "q975_scaled", "rhat",
    "ess_bulk", "ess_tail"),
  names(paper_table)
)
log_obj(
  "Paper table (mean / CrI on the scaled theta used by the sampler):",
  paper_table[, scaled_cols, drop = FALSE]
)
log_obj(
  "Same posterior on the ODE scale (S*100, beta/100, Lambda*10):",
  paper_table[, c("ode_parameter", "mean_ode", "q025_ode", "q975_ode"),
              drop = FALSE]
)

if (SAVE_GGMCMC) {
  pdf_path <- file.path(out_dir, "mcmc_diagnostics_poisson_two_chains.pdf")
  log_msg("Writing ggmcmc report to ", pdf_path)
  ggs_object <- ggs(mcmc_list)
  ggmcmc(ggs_object, file = pdf_path)
  log_msg("Saved ", pdf_path)
}

log_msg("\nSession info:")
log_msg(paste(capture.output(sessionInfo()), collapse = "\n"))
log_msg(
  "\nDone. Reported draws: results/sample5_poisson-simplified-waterfowl-raptor.RData"
)
log_msg("Scaled posterior matrix: posterior_samples (rows = both chains).")
log_msg("ODE-scale numbers for the paper: results/posterior_summary.csv")

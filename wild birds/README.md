# Ontario wild-bird HPAI: migration model, Poisson likelihood

Companion code for the manuscript fit. One script,
[`fit_poisson_simplified_waterfowl_raptor.R`](fit_poisson_simplified_waterfowl_raptor.R),
runs the whole path that was originally done interactively: scalar-scale
pilot chains, a frozen covariance proposal, then the production chains.

Natural mortality (`mu`) is not in the model. Disease-induced mortality `d = 0.3`
is fixed. Waterfowl immigrate through day 60 and depart after day 120.
Raptors do not migrate. Daily counts are Poisson increments of cumulative
disease deaths.

## Layout

```
data/Ontario-wild-birds-waterfowl-daily_cases.csv
data/Ontario-wild-birds-raptor-daily_cases.csv
fit_poisson_simplified_waterfowl_raptor.R
results/                  created by the script
```

Each CSV needs columns `CollectionDate` and `DailyCases`. The two files must
contain the same dates from 2024-11-01 through 2025-03-25. The script does
not fill missing days with zeros, because that would change the likelihood
that was fitted. If a day had no detections, put that day in both files with
`DailyCases = 0`.

## Packages

```r
install.packages(c(
  "deSolve", "mcmc", "coda", "ggmcmc",
  "dplyr", "lubridate", "ggplot2", "posterior"
))
```

## Run

From this directory:

```bash
Rscript fit_poisson_simplified_waterfowl_raptor.R
```

That is the manuscript run. It is long (days on a laptop): every proposal
solves the ODE, and there are two chains through a scalar-scale warm-up
(20,000), a scalar-scale sample (100,000 batches of length 50), two short
pilots, then production blocks of 50,000, 150,000, and 150,000 iterations.

Check the plumbing without producing usable draws:

```bash
HPAI_SMOKE_TEST=1 Rscript fit_poisson_simplified_waterfowl_raptor.R
```

## What to report

The posterior in the paper is the **last block only**, after a 20% burn-in:

* `results/sample5_poisson-simplified-waterfowl-raptor.RData`
* `posterior_samples` in that file (both chains, scaled parameters)
* `results/posterior_summary.csv` (scaled and ODE-scale mean and 95% CrI, plus Rhat and ESS)
* `results/mcmc_diagnostics_poisson_two_chains.pdf`
* `results/run_log.txt` (acceptance, Gelman–Rubin, ESS, Geweke, session info)

Earlier files `sample1`, `sample3`, and `sample4` are the tuning path
(scalar proposal, then the first two production blocks). They are not
pooled into the reported posterior. There was no `sample2`.

The sampler's parameters are scaled. `print(post_mean)` in the original
session, and `posterior_samples` here, are on that scaled scale:

| Scaled (MCMC) | In the ODE |
| --- | --- |
| `S_w0` | initial susceptible waterfowl = `S_w0 * 100` |
| `I_w0` | initial infected waterfowl |
| `beta_w` | `beta = beta_w / 100` |
| `Lambda_amp_w` | immigration rate = `Lambda_amp_w * 10` |
| `delta_amp_w` | departure rate |
| `S_r0` | initial susceptible raptors = `S_r0 * 100` |
| `beta_wr` | cross-species transmission = `beta_wr / 100` |

Use `mean_ode` / `q025_ode` / `q975_ode` in `posterior_summary.csv` for
ODE-scale numbers. `d` is fixed at 0.3.

## MCMC, briefly

`mcmc::metrop` is random-walk Metropolis with a **fixed** proposal. The
20,000-iteration "warm-up" does not adapt the scale.

Phase 1 uses independent normal proposals with sd 0.05. Its kept batch means
are thinned to 4,000 rows and turned into a multivariate proposal

```
Sigma_prop = (2.38^2 / 7) * (0.95 * cov + 0.05 * diagonal)
scale_mat  = t(chol(Sigma_prop))
```

Pilots with `0.7 * scale_mat` and `0.35 * scale_mat` were used only to
check the acceptance rate. Production uses `0.35 * scale_mat`, frozen, so
the production chains are ordinary random-walk Metropolis.

`set.seed(1)` and `set.seed(2)` start the two chains. Later `set.seed`
calls are in the file because they were in the original session, but
`mcmc::metrop` restores the previous fit's random-number state whenever a
chain is continued, so those later seeds do not change the draws.

Priors are independent normals on the scaled parameters, with density zero
for negative values. The normal density is **not** renormalized on the
positive orthant. Poisson means that are not strictly positive are floored
at `1e-6`. Both choices match the fitted model; do not "correct" them
without refitting.

## Re-running the numbers

A clean re-run will follow the same seeds and the same proposal
construction, so it should reproduce this script. It will not bit-match an
old `.RData` if the original session had a different R version, a different
`RNGkind`, or a different `deSolve` version (the likelihood is an ODE).
`results/run_log.txt` records those versions.

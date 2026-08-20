# HPAI-Ontario
From Wild-Bird Mortality to Poultry Farm Introductions: HPAI Spillover in Ontario, Canada
# Code and processed data supporting the manuscript

## Overview

This archive contains R code and processed data used for analyses of highly
pathogenic avian influenza (HPAI) in wild birds and poultry farms in Ontario,
Canada. 
The two CSV files are processed daily case-count series derived from a
third-party public dashboard. The poultry-farm observations used in the model
are transcribed in the poultry analysis script from a separate third-party
public dashboard. This archive does not contain the complete underlying
national datasets.

## Directory and file descriptions

### `wild birds/`

- `Ontario-wild-birds-waterfowl-daily_cases.csv`: processed daily counts of
  confirmed HPAI cases in Ontario waterfowl, covering 2022-03-15 through
  2025-03-25 (1,107 observations, excluding the header).
- `Ontario-wild-birds-raptor-daily_cases.csv`: processed daily counts of
  confirmed HPAI cases in Ontario raptors, covering 2022-03-22 through
  2025-04-29 (1,135 observations, excluding the header).
- `ontario-DEAD-CASES-waterfowl-raptor-2024-2025.R`: prepares the two wild-bird
  daily series and fits the wild-bird transmission model.

### `poultry/`

- `MAIN-POULTRY-diff-i0-farm.R`: contains the Ontario poultry-farm observations,
  fits the within-farm model, and writes the posterior object
  `sample1_ontario-new-11-I-y-kx-poultry-farm-2024-2025.RData` to the current
  working directory.

### `wildbirds-farm/`

- `sample1_ontario-new-11-I-y-kx-poultry-farm-2024-2025.RData`: saved R object
  named `sample1`, produced by the poultry analysis. Its `batch` component is a
  100,000 by 12 MCMC matrix. The first 11 columns represent farm-specific
  initial infected-bird parameters (`I_p1` through `I_p11`), and the final
  column represents the poultry mortality-rate parameter (`d_p`).
- `link-wildbirds-farm.R`: uses the saved within-farm posterior and the modeled
  wild-bird infection pressure to estimate wild-bird-to-poultry spillover.

### `lag analysis&permutation test/`

- `spillover_lag_permutation_tests.R`: performs the lag-profile analysis and
  two randomization analyses (uniform randomization of introduction dates and
  circular shifts). It writes five CSV result tables and four PNG figures to
  the current working directory.

## Data dictionary

Both wild-bird CSV files have the same columns:

| Column | Type | Unit/format | Description |
|---|---|---|---|
| `CollectionDate` | date | `YYYY-MM-DD` (calendar day) | Date on which the animal was sampled/collected. |
| `DailyCases` | integer | confirmed animals per day | Number of confirmed HPAI-positive animals in the specified bird group and Ontario on that date. Zeroes denote dates with no cases in the processed series. |

The poultry data frame embedded in `MAIN-POULTRY-diff-i0-farm.R` contains:

| Column | Type | Unit/format | Description |
|---|---|---|---|
| `Outbreak` | character | CFIA premises identifier | Identifier for the infected premises. |
| `Susceptible` | integer | birds | Number of susceptible birds on the premises. |
| `Dead` | integer | birds | Number of dead birds used in the analysis. |
| `clinical_date` | date | `YYYY-MM-DD` | Date of clinical signs/reporting used by the analysis. |
| `depop_date` | date | `YYYY-MM-DD` | Premises depopulation date. |
| `species_group` | categorical | `Turkey` or `Other` | Poultry species grouping used by the analysis. |

The script derives `fraction_dead` (proportion), `delay_days` (days),
`flock_group` (categorical flock-size band), and `I_group` (model parameter
mapping) from these columns.

## Third-party data provenance

### Wild birds

The full source data are provided through the Canadian Wildlife Health
Cooperative (CWHC) Highly Pathogenic Avian Influenza Wildlife Dashboard:

https://www.cwhc-rcsf.ca/avian_influenza.php



### Poultry farms

The full source data are provided through the Canadian Food Inspection Agency
(CFIA) Highly Pathogenic Avian Influenza Dashboards, under “Domestic birds”:

https://inspection.canada.ca/en/animal-health/terrestrial-animals/diseases/reportable/avian-influenza/latest-bird-flu-situation/hpai-dashboards

The dashboard reports Canadian infected premises and related outbreak
information. The Ontario premises-level observations used here are transcribed
in `poultry/MAIN-POULTRY-diff-i0-farm.R`. 

The original dashboard providers retain rights in the source material. Verify
their applicable terms and confirm that redistribution of these processed
subsets is permitted before applying an open licence to the complete deposit.

## Software requirements

The analyses use R. The package was inspected with R 4.6.1; replace this with
the R version actually used for the manuscript analyses if different.

Required R packages across the scripts are: `deSolve`, `mcmc`, `coda`,
`stringr`, `dplyr`, `ggplot2`, `lubridate`, `tidyr`, `ggmcmc`, and `scales`.
Package versions were not recorded in the supplied files. For reproducibility,
add the output of `sessionInfo()` from the analysis environment or create an
`renv.lock` file before deposit.

## Suggested run order

The scripts use paths relative to the current working directory.

1. Set the working directory to `wild birds/` and run
   `ontario-DEAD-CASES-waterfowl-raptor-2024-2025.R`.
2. Set the working directory to `poultry/` and run
   `MAIN-POULTRY-diff-i0-farm.R`. This MCMC analysis may take substantial time.
3. Copy the generated poultry posterior RData file into `wildbirds-farm/`, or
   use the included saved posterior, then set the working directory to
   `wildbirds-farm/` and run `link-wildbirds-farm.R`.
4. Set the working directory to `lag analysis&permutation test/`. Copy the
   included poultry posterior RData file from `wildbirds-farm/` into this
   directory, then run `spillover_lag_permutation_tests.R`. The script performs
   1,000 randomizations by default and may take substantial time.

## Outputs of the lag and permutation analysis

- `reconstructed_farm_introduction_dates.csv`
- `spillover_lag_profile_results.csv`
- `spillover_uniform_randomization_results.csv`
- `spillover_circular_shift_results.csv`
- `spillover_temporal_test_summary.csv`
- `spillover_lag_profile.png`
- `spillover_uniform_randomization.png`
- `spillover_circular_shift.png`
- `spillover_randomization_combined.png`



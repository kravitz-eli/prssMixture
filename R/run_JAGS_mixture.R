# Load packages ---------------------------
library(here) # relative file paths
library(rjags) # R interface for JAGS
library(purrr) # better apply() functions
library(tidybayes) # tools to work with MCMC output

# Load functions I've written  -----------------------
custom_functions <- c(
  "normal_mixture_jags.R", # JAGS code to run mixture prior model
  "make_cov_matrix.R" # JAGS code to run mixture prior model
)

walk(custom_functions, ~source(here("R", .x)))

# Set parameters in bivariate normal likelihood -------------------------------
likelihood_params <- list(
  mu = log(c(0.76, 0.84)),
  sigma_1 = 1 / 48 + 1 / 37,
  sigma_2 = 1 / 31 + 1 / 32,
  rho = 0.55
)

# Make covariance matrix from variance and covariance terms
likelihood_params[["sigma"]] <- with(
  likelihood_params,
  make_cov_matrix(sigma_1, sigma_2, rho)
)

# Set prior parameters ---------------------------------
# Mixing proportion
p_mix <- 0.50

# First component
comp1_prior_params <- list(
  mu = log(c(0.68, 0.74)),
  sigma_1 = 0.2,
  sigma_2 = 0.2,
  rho = 0.30,
  alpha = 1,
  beta = 1
)

# Make covariance matrix from variance and covariance terms
comp1_prior_params[["sigma"]] <- with(
  comp1_prior_params,
  make_cov_matrix(sigma_1, sigma_2, rho)
)

# Second component
comp2_prior_params <- list(
  mu = log(c(1, 1)),
  sigma_1 = 1,
  sigma_2 = 1,
  rho = 0.5,
  alpha = 1,
  beta = 1
)

# Make covariance matrix from variance and covariance terms
comp2_prior_params[["sigma"]] <- with(
  comp2_prior_params,
  make_cov_matrix(sigma_1, sigma_2, rho)
)

# Run JAGS model ------------------------------------------------
jags_results <- normal_mixture_jags(
  likelihood_params,
  p_mix = p_mix,
  comp1_prior_params,
  comp2_prior_params,
  n.chains = 1,
  burn_in = 100,
  n_iter = 1000
)


jags_results %>%
  tidy_draws %>%
  median_qi()

jags_results %>%
  tidy_draws %>%
  group_by(.chain) %>%
  median_qi()


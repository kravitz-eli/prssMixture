# Load packages ---------------------------
# General purpose packages
library(here) # relative file paths
library(magrittr) #pipes
library(purrr) # better apply() functions
library(dplyr) # data.frame manipulation
library(ggplot2)

# Bayes packages
library(rjags) # R interface for JAGS
library(tidybayes) # tools to work with MCMC output
library(bayesplot)

# Load functions I've written  -----------------------
source(here("R", "run_jags.R"))

# Set model parameters -----------------------------------------------------
likelihood_params <- list(
  y_obs = log(c(0.76, 0.84)),
  sigma_1 = 1 / 48 + 1 / 37,
  sigma_2 = 1 / 31 + 1 / 32,
  rho = 0.55
)

# Set prior parameters ---------------------------------
# Mixing proportion
p_mix <- 0.50

# First component
comp1_prior_params <- list(
  mu = log(c(0.50, 0.65)),
  sigma_1 = 0.2,
  sigma_2 = 0.2,
  rho = 0.30,
  alpha = 1,
  beta = 1
)

# Second component
comp2_prior_params <- list(
  mu = log(c(1, 1)),
  sigma_1 = 10^2,
  sigma_2 = 10^2,
  rho = 0,
  alpha = 1,
  beta = 1
)

# Combine into a list to pass to jags
# parameters are prefixed with a
data = c("L" = likelihood_params,
         "p_mix" = p_mix,
         "P1" = comp1_prior_params,
         "P2" = comp2_prior_params)


# Run JAGS model ------------------------------------------------
jags_results <- run_jags(
  model_file = "normal_mixture",
  data = data,
  track_variable_names = c("theta", "ind"),
  iter = 1e5,
  burn = 1e3,
  chains = 1,
  progress.bar = "text")

# Run with prior on P
data[["p_mix"]] = NULL
jags_results_prior_p_mix <- run_jags(
  model_file = "normal_mixture_prior_p",
  data = data,
  track_variable_names = c("theta", "ind", "weight"),
  iter = 1e5,
  burn = 1e3,
  chains = 1,
  progress.bar = "text"
)

# Run with prior on P
data[["p_mix"]] = NULL
jags_results_prior_p_mix <- run_jags(
  model_file = "normal_mixture_prior_p",
  data = data,
  track_variable_names = c("theta", "ind", "weight"),
  iter = 1e5,
  burn = 1e3,
  chains = 1,

  progress.bar = "text"
)

# See posterior of mixing proportion
jags_results_prior_p_mix %>%
  spread_draws(ind) %>%
  summarise(mean(ind), median(ind), min(ind), max(ind))

  post_weight = jags_results_prior_p_mix %>%
    spread_draws(weight, ind) %>%
    pull(weight)

ggplot(data.frame(post_weight), aes(x = post_weight)) +
  geom_histogram(aes(y =..density..), colour="black", fill="white")+
  geom_density(alpha = 0.40, fill="red") +
  theme_minimal()


prior_var = c(1e-2, 1e-1, 1, 1e2, 1e3, 1e4, 1e5, 1e6)
jags_results = vector("list", length(prior_var))

for (i in seq_along(prior_var)) {

  data[c("P2.sigma_1", "P2.sigma_2")] = prior_var[[i]]

  jags_results <- run_jags(
    model_file = "normal_mixture_prior_p",
    data = data,
    track_variable_names = c("theta", "ind", "weight"),
    iter = 1e5,
    burn = 1e3,
    chains = 1,
    progress.bar = "text"
  )




}

# Get individual posterior distributions ------------------------------------
weight_post = jags_results %>%
  map(spread_draws, weight, ind) %>%
  map(pull, weight)

ind_post = jags_results %>%
  map(spread_draws, weight, ind) %>%
  map(pull, ind)

theta_post = jags_results %>%
  map(as.matrix) %>%
  map(~.x[, c("theta[1]", "theta[2]")])

theta_1_post = theta_post %>%
  map(~.x[, "theta[1]"])

theta_2_post = theta_post %>%
  map(~.x[, "theta[2]"])

# Get posterior means and medians --------------------------------
weights_post_means = map_dbl(weight_post, mean)
ind_post_mean = map_dbl(ind_post, mean)
theta1_post_mean = map_dbl(theta_1_post, mean)
theta2_post_mean = map_dbl(theta_2_post, mean)


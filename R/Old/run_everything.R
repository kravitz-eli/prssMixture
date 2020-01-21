library("here")
source(here("R", "posterior_mixture.R"))
source(here("R", "get_mixing_prop.R"))
source(here("R", "make_cov_matrix.R"))
source(here("R", "normal_normal_conjugate.R"))

# Set parameters in bivariate normal likelihood -------------------------------
likelihood_params = list(
  y_obs = log(c(0.76, 0.84)),
  sigma_L_11 = 1 / 48 + 1/ 37,
  sigma_L_22 = 1/31 + 1/32,
  rho_L = 0.55
)
#
# # Make covariance matrix from variance and covariance terms
# likelihood_params[["sigma_L"]] = with(likelihood_params,
#                                       make_cov_matrix(sigma_L_11, sigma_L_22, rho_L))

# Set prior parameters ---------------------------------
# Mixing proportion
p_mix = 0.50

# First component
comp1_prior_params = list(
  mu_P = log(c(0.68, 0.74)),
  sigma_P_11 = 0.2,
  sigma_P_22 = 0.2,
  rho_P = 0.30
)

# # Make covariance matrix from variance and covariance terms
# comp1_prior_params[["sigma_P"]] = with(comp1_prior_params,
#                                        make_cov_matrix(sigma_P_11, sigma_P_22, rho_P))

# Second component
comp2_prior_params = list(
  mu_P = log(c(1, 1)),
  sigma_P_11 = 1,
  sigma_P_22 = 1,
  rho_P = 0.5
)

# # Make covariance matrix from variance and covariance terms
# comp2_prior_params[["sigma_P"]] = with(comp2_prior_params,
#                                        make_cov_matrix(sigma_P_11, sigma_P_22, rho_P))

# Get posterior parameter
post_closed = posterior_mixture(likelihood_params,
                  p_mix,
                  comp1_prior_params,
                  comp2_prior_params)

# Run JAGS model to test against
library("rjags")
source(here("R", "normal_mixture_jags_no_cor_prior.R"))
jags_results = normal_mixture_jags(
  likelihood_params,
  p_mix = p_mix,
  comp1_prior_params,
  comp2_prior_params
)

# See if posterior parameters match JAGS --------------
# Get posterior parameters
p_mix_post = post_closed$p_mix_post # Mixing param
mu_1_post = post_closed$mu_post_1 # Mean of first component
mu_2_post = post_closed$mu_post_2 # Mean of second component
mu_comb_post = mu_1_post * p_mix_post + (1 - p_mix_post) * mu_2_post # Combined mean
sigma_mix_post = p_mix_post^2 * post_closed$sigma_post_1 + #Combined variance
  (1 - p_mix_post)^2 * post_closed$sigma_post_2

# Compare closed form to JAGS
error_p_mix = p_mix_post - jags_results$p_mix_post
error_mu = mu_comb_post - jags_results$mu_post
error_sigma = sigma_mix_post - jags_results$sigma_post





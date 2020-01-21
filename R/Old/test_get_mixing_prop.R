
# Likelihood -------------------
y_obs = log(c(0.76, 0.84))
sigma_L_11 = 1 / 48 + 1/ 37
sigma_L_22 = 1/31 + 1/32
rho_L = 0.55
sigma_L = matrix(c(sigma_L_11, rho_L * sqrt(sigma_L_11 * sigma_L_22),
                   rho_L * sqrt(sigma_L_11 * sigma_L_22), sigma_L_22),
                 nrow = 2)

# Prior -------------------
p_mix = 0.50
# First component
mu_P_1 = log(c(0.68, 0.74))
sigma_P_11_1 = 1
sigma_P_22_1 = 1
rho_P_1 = 0.30

sigma_P_1 = matrix(c(sigma_P_11_1, rho_P_1 * sqrt(sigma_P_11_1 * sigma_P_22_1),
                     rho_P_1 * sqrt(sigma_P_11_1 * sigma_P_22_1), sigma_P_22_1),
                   nrow = 2)

# Second component
mu_P_2 = log(c(1, 1))
sigma_P_11_2 = 1
sigma_P_22_2 = 1
rho_P_2 = 0.5

sigma_P_2 = matrix(c(sigma_P_11_2, rho_P_2 * sqrt(sigma_P_11_2 * sigma_P_22_2),
                     rho_P_2 * sqrt(sigma_P_11_2 * sigma_P_22_2), sigma_P_22_2),
                   nrow = 2)


# Get hyperparameters of both components in mixture ----------------------------
# Component 1
comp1_closed_form = normal_normal_conjugate(
  y_obs = y_obs,
  sigma_L_11 = sigma_L_11,
  sigma_L_22 = sigma_L_22,
  rho_L = rho_L,
  mu_P = mu_P_1,
  sigma_P_11 = sigma_P_11_1,
  sigma_P_22 = sigma_P_22_1,
  rho_P = rho_P_1
)

mu_post_1 = comp1_closed_form$mu_post
sigma_post_1 = comp1_closed_form$sigma_post

# Component 2
comp2_closed_form  = normal_normal_conjugate(
  y_obs = y_obs,
  sigma_L_11 = sigma_L_11 ,
  sigma_L_22 = sigma_L_22 ,
  rho_L = rho_L,
  mu_P = mu_P_2,
  sigma_P_11 = sigma_P_11_2,
  sigma_P_22 = sigma_P_22_2,
  rho_P = rho_P_2
)

mu_post_2 = comp2_closed_form$mu_post
sigma_post_2 = comp2_closed_form$sigma_post


# Calculate posterir mixing weight
f1 = mvtnorm::dmvnorm(y_obs, mean = mu_P_1, sigma = sigma_L + sigma_P_1)
f2 = mvtnorm::dmvnorm(y_obs, mean = mu_P_2, sigma = sigma_L + sigma_P_2)

p_mix_post = p_mix * f1 / (p_mix * f1 + (1-p_mix) * f2)
attributes(p_mix_post) = NULL


library("rjags")

jags_results = normal_mixture_jags(
  y_obs = y_obs,
  tau_L = tau_L,
  mu_P_1 = mu_P_1,
  mu_P_2 = mu_P_2,
  sigma_P_11_1 = sigma_P_11_1,
  sigma_P_22_1 = sigma_P_22_1,
  sigma_P_11_2 = sigma_P_11_2,
  sigma_P_22_1 = sigma_P_22_1,
  rho_P_1 = rho_P_1,
  rho_P_2 = rho_P_2,
  p_mix = p_mix,
)

p_mix_post_jags = jags_results$p_mix_post

c(p_mix_post, p_mix_post_jags, abs(p_mix_post - p_mix_post_jags))

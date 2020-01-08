
# Likelihood -------------------
y_obs = log(c(0.65, 0.88))
sigma_L_11 = 4 / 31
sigma_L_22 = 4 / 12
rho_L = 0.30
sigma_L = matrix(c(sigma_L_11, rho_L * sqrt(sigma_L_11 + sigma_L_22),
                   rho_L * sqrt(sigma_L_11 + sigma_L_22), sigma_L_22),
                 nrow = 2)

# Prior -------------------
p_mix = 0.50
# First component
mu_P_1 = log(c(0.63, 0.90))
sigma_P_11_1 = 0.10
sigma_P_22_1 = 0.20
rho_P_1 = 0.30

# Second component
mu_P_2 = log(c(1.2, 1.5))
sigma_P_11_2 = 100
sigma_P_22_2 = 100
rho_P_2 = 0.50




# Get hyperparameters of both components in mixture ----------------------------
# Component 1
comp1_closed_form  = normal_normal_conjugate(
  y_obs = y_obs,
  sigma_L_11 = sigma_L_11 ,
  sigma_L_22 = sigma_L_22 ,
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

f1 = mvtnorm::pmvnorm(y_obs, mean = mu_post_1, sigma = sigma_L + sigma_post_1)
f2 = mvtnorm::pmvnorm(y_obs, mean = mu_post_2, sigma = sigma_L + sigma_post_2)


p_mix_post = p_mix * f1 / (p_mix * f1 + (1-p_mix) * f2)



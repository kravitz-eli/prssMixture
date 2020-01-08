closed_form  = normal_normal_conjugate(
  y_obs = log(c(0.65, 0.88)),
  sigma_L_11 = 4 / 31,
  sigma_L_22 = 4 / 6,
  rho_L = 0.65,
  mu_P = c(0, 0),
  sigma_P_11 = 0.10,
  sigma_P_22 = 0.20,
  rho_P = 0.5
)

jags_results = normal_normal_jags(
  y_obs = log(c(0.65, 0.88)),
  sigma_L_11 = 4 / 31,
  sigma_L_22 = 4 / 6,
  rho_L = 0.65,
  mu_P = c(0, 0),
  sigma_P_11 = 0.10,
  sigma_P_22 = 0.20,
  rho_P = 0.5
)


# See how different the results are
diffs = c(
  "mu" = closed_form$mu_post - jags_results$mu_post,
  "sigma" = solve(closed_form$tau_post) - jags_results$sigma_post
)

max(diffs)

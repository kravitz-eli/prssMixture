
normal_normal_conjugate = function(
  y_obs,
  sigma_L_11,
  sigma_L_22,
  rho_L,
  mu_P,
  sigma_P_11,
  sigma_P_22,
  rho_P
){

  # Likelihood covaraince and precision
  sigma_L = matrix(c(sigma_L_11, rho_L * sqrt(sigma_L_11 * sigma_L_22),
                     rho_L * sqrt(sigma_L_11 * sigma_L_22), sigma_L_22),
                   nrow = 2)
  tau_L = solve(sigma_L)

  # Prior covariance and precision
  # Likelihood covaraince and precision
  sigma_P = matrix(c(sigma_P_11, rho_P * sqrt(sigma_P_11 * sigma_P_22),
                     rho_P * sqrt(sigma_P_11 * sigma_P_22), sigma_P_22),
                   nrow = 2)

  tau_P = solve(sigma_P)

  # Posterior parameters ----------------------------------------
  # Precision
  tau_post = tau_P + tau_L
  # Mean
  mu_post = solve(tau_post) %*% (tau_P %*% mu_P + tau_L %*% y_obs)
  mu_post = as.vector(mu_post)


  return(list(
    "mu_post" = mu_post,
    "tau_post" = tau_post,
    "sigma_post" = solve(tau_post)
  ))

}

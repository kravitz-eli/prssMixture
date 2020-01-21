
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

  # Likelihood covariance and precision ----------------------------------------
  sigma_L = make_cov_matrix(var_1 = sigma_L_11,
                            var_2 = sigma_L_22,
                            rho = rho_L)
  tau_L = solve(sigma_L)

  # Prior covariance and precision ---------------------------------------------
  sigma_P = make_cov_matrix(var_1 = sigma_P_11,
                            var_2 = sigma_P_22,
                            rho = rho_P)
  tau_P = solve(sigma_P)

  # Posterior parameters ------------------------------------------------------
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

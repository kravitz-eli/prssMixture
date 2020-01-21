normal_normal_jags = function(
  y_obs,
  sigma_L_11,
  sigma_L_22,
  rho_L,
  mu_P,
  sigma_P_11,
  sigma_P_22,
  rho_P,
  n.iter = 1e5,
  burn_in = 1e3){


  # Delare model
  bivariate_norm_model = "model{
  # Likelihood  -----------------------------------
  y_obs ~ dmnorm(theta[1:2], tau_L)

  # Priors -----------------------------------------
  # Multivariate normal
  theta ~ dmnorm(mu_P[1:2], tau_P)
  tau_P <-inverse(sigma_P)
  sigma_P[1,1] <- sigma_P_11
  sigma_P[1,2] <- rho_P * sqrt(sigma_P_11 * sigma_P_22)
  sigma_P[2,1] <- rho_P * sqrt(sigma_P_11 * sigma_P_22)
  sigma_P[2,2] <- sigma_P_22
}"

  # Put in parameters
  ## Likelihood covaraince and precision
  sigma_L = matrix(c(sigma_L_11, rho_L * sqrt(sigma_L_11 * sigma_L_22),
                     rho_L * sqrt(sigma_L_11 * sigma_L_22), sigma_L_22),
                   nrow = 2)
  tau_L = solve(sigma_L)

  data = list(
    y_obs = y_obs,
    tau_L = tau_L,
    mu_P = mu_P,
    sigma_P_11 = sigma_P_11,
    sigma_P_22 = sigma_P_22,
    rho_P = rho_P
  )

  model = jags.model(
    file = textConnection(bivariate_norm_model),
    data = data
  )

  # Do a burn-in
  update(model, n.iter = burn_in)

  # Run model to (hopefully) convergence
  samples <- coda.samples(model,
                          variable.names = c("theta"),
                          n.iter = n.iter)

  # Get posterior samples for parameters
  mu_post_samples = as.matrix(samples[ ,c("theta[1]", "theta[2]")])
  colnames(mu_post_samples) = NULL


  mu_post = colMeans(mu_post_samples)
  sigma_post = cov(mu_post_samples)
  tau_post = solve(sigma_post)

  return(list(
    "mu_post" = mu_post,
    "tau_post" = tau_post,
    "sigma_post" = sigma_post
  ))

}

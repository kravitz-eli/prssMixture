normal_mixture_jags = function(
  likelihood_params,
  comp1_prior_params,
  comp2_prior_params,
  p_mix,
  n.iter = 1e5,
  burn_in = 1e3
){


  # Delare model
  normal_mixture = "model{
    # Likelihood (phase II OR and PFS) -----------------------------------
    y_obs ~ dmnorm(theta[1:2], tau_L)
    tau_L <-inverse(sigma_L)

    theta <- ind * component_1 + (1-ind) * component_2
    ind ~ dbern(p_mix)

    # Priors ------------------------------------------------------------------------------------------------------------------------
    # Prior for first component in normal mixture.
    component_1 ~ dmnorm(mu_P_1[1:2], tau_P_1)
    tau_P_1 <-inverse(sigma_P_1)
    sigma_P_1[1,1] <- sigma_P_11_1
    sigma_P_1[1,2] <- rho_P_1 * sqrt(sigma_P_11_1 * sigma_P_22_1)
    sigma_P_1[2,1] <- rho_P_1 * sqrt(sigma_P_11_1 * sigma_P_22_1)
    sigma_P_1[2,2] <- sigma_P_22_1

    # Prior for first component in normal mixture.
    component_2 ~ dmnorm(mu_P_2[1:2], tau_P_2)
    tau_P_2 <- inverse(sigma_P_2)
    sigma_P_2[1,1] <- sigma_P_11_2
    sigma_P_2[1,2] <- rho_P_2 * sqrt(sigma_P_11_2 * sigma_P_22_2)
    sigma_P_2[2,1] <- rho_P_2 * sqrt(sigma_P_11_2 * sigma_P_22_2)
    sigma_P_2[2,2] <- sigma_P_22_2

  }"


  # Put in parameters
  ## Likelihood covaraince and precision
  sigma_L = with(likelihood_params,
                 matrix(c(sigma_L_11, rho_L * sqrt(sigma_L_11 * sigma_L_22),
                          rho_L * sqrt(sigma_L_11 * sigma_L_22), sigma_L_22),
                        nrow = 2))

  data = list(
    y_obs = likelihood_params$y_obs,
    sigma_L = sigma_L,
    p_mix = p_mix,
    #First component
    mu_P_1 = comp1_prior_params$mu_P,
    sigma_P_11_1 = comp1_prior_params$sigma_P_11,
    sigma_P_22_1 = comp1_prior_params$sigma_P_22,
    rho_P_1 = comp1_prior_params$rho_P,
    # Second component
    mu_P_2 = comp2_prior_params$mu_P,
    sigma_P_11_2 = comp2_prior_params$sigma_P_11,
    sigma_P_22_2 = comp2_prior_params$sigma_P_22,
    rho_P_2 = comp2_prior_params$rho_P
  )

  model = jags.model(
    file = textConnection(normal_mixture),
    data = data
  )

  # Do a burn-in
  update(model, n.iter = burn_in)

  # Run model to (hopefully) convergence
  samples <- coda.samples(model,
                          variable.names = c("theta", "ind"),
                          n.iter = n.iter)

  # Get posterior samples for parameters
  mu_post_samples = as.matrix(samples[ ,c("theta[1]", "theta[2]")])
  colnames(mu_post_samples) = NULL

  p_mix_post_samples = as.matrix(samples[, "ind"])


  mu_post = colMeans(mu_post_samples)
  sigma_post = cov(mu_post_samples)
  tau_post = solve(sigma_post)
  p_mix_post = mean(p_mix_post_samples)

  return(list(
    "mu_post" = mu_post,
    "tau_post" = tau_post,
    "sigma_post" = sigma_post,
    "p_mix_post" = p_mix_post
  ))

}

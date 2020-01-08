normal_mixture_jags = function(...,
                               n.iter = 1e5,
                               burn_in = 1e3){


  # Delare model

  normal_mixture = "model{
    # Likelihood (phase II OR and PFS) -----------------------------------
    y_obs ~ dmnorm(theta[1:2], tau_L)
    tau_L <-inverse(sigma_L)

    mu <- ind * component_1 + (1-ind) * component_2
    ind ~ dbern(mix_prop)

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
    tau_P_2 <-inverse(sigma_P_2)
    sigma_P_2[1,1] <- sigma_P_11_2
    sigma_P_2[1,2] <- rho_P_2 * sqrt(sigma_P_11_2 * sigma_P_22_2)
    sigma_P_2[2,1] <- rho_P_2 * sqrt(sigma_P_11_2 * sigma_P_22_2)

    #Prior for correlation ----------------------------------
    rho~dbeta(1,1)


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



# Likelihood -------------------
likelihood_params = list(
  y_obs = log(c(0.65, 0.88)),
  sigma_L_11 = 4 / 31,
  sigma_L_22 = 4 / 12,
  rho_L = 0.30,
  sigma_L = matrix(c(sigma_L_11, rho_L * sqrt(sigma_L_11 + sigma_L_22),
                     rho_L * sqrt(sigma_L_11 + sigma_L_22), sigma_L_22),
                   nrow = 2)
)

# Prior -------------------
p_mix = 0.50
prior_1_params = list(
  mu_P_1 = log(c(0.63, 0.90)),
  sigma_P_11_1 = 0.10,
  sigma_P_22_1 = 0.20,
  rho_P_1 = 0.30
)

# Second component
prior_2_params = list(
  mu_P_2 = log(c(1.2, 1.5)),
  sigma_P_11_2 = 100,
  sigma_P_22_2 = 100,
  rho_P_2 = 0.50
)

all_params = c(
  likelihood_params,
  p_mix,
  prior_1_params,
  prior_2_params
)

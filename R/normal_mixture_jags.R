normal_mixture_jags = function(
  likelihood_params,
  comp1_prior_params,
  comp2_prior_params,
  p_mix,
  n_iter = 1e5,
  burn_in = 1e3,
  n.chains = 1
){


  normal_mixture = "model{

    # Likelihood (phase II OR and PFS) -----------------------------------
    y_obs ~ dmnorm(theta[], tau_L)
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
    rho_P_1 ~ dbeta(a_1, b_1)

    # Prior for first component in normal mixture.
    component_2 ~ dmnorm(mu_P_2[1:2], tau_P_2)
    tau_P_2 <- inverse(sigma_P_2)
    sigma_P_2[1,1] <- sigma_P_11_2
    sigma_P_2[1,2] <- rho_P_2 * sqrt(sigma_P_11_2 * sigma_P_22_2)
    sigma_P_2[2,1] <- rho_P_2 * sqrt(sigma_P_11_2 * sigma_P_22_2)
    sigma_P_2[2,2] <- sigma_P_22_2
    rho_P_2 ~ dbeta(a_2, b_2)

  }"



  data = list(
    y_obs = likelihood_params$mu,
    sigma_L = likelihood_params$sigma,
    p_mix = p_mix,
    # First component
    mu_P_1 = comp1_prior_params$mu,
    sigma_P_11_1 = comp1_prior_params$sigma_1,
    sigma_P_22_1 = comp1_prior_params$sigma_2,
    a_1 = comp1_prior_params$alpha,
    b_1 = comp1_prior_params$beta,
    # Second component
    mu_P_2 = comp2_prior_params$mu,
    sigma_P_11_2 = comp2_prior_params$sigma_1,
    sigma_P_22_2 = comp2_prior_params$sigma_2,
    rho_P_2 = comp2_prior_params$rho,
    a_2 = comp2_prior_params$alpha,
    b_2 = comp2_prior_params$beta
  )

  model <- jags.model(
    file = textConnection(normal_mixture),
    data = data,
    n.chains = n.chains
  )

  # Do a burn-in
  update(model, n.iter = burn_in)

  # Run model to (hopefully) convergence
  coda.samples(model,
               variable.names = c("theta", "ind"),
               n.iter = n_iter)


}
#
#
#
#
# # Set parameters in bivariate normal likelihood -------------------------------
# likelihood_params = list(
#   y_obs = log(c(0.76, 0.84)),
#   sigma_L_11 = 1 / 48 + 1/ 37,
#   sigma_L_22 = 1/31 + 1/32,
#   rho_L = 0.55
# )
#
# # Set prior parameters ---------------------------------
# # Mixing proportion
# p_mix = 0.50
#
# # First component
# comp1_prior_params = list(
#   mu_P = log(c(0.68, 0.74)),
#   sigma_P_11 = 0.2,
#   sigma_P_22 = 0.2,
#   alpha = 1,
#   beta = 1
# )
#
# # Second component
# comp2_prior_params = list(
#   mu_P = log(c(1, 1)),
#   sigma_P_11 = 1,
#   sigma_P_22 = 1,
#   alpha = 1,
#   beta = 1
# )
#
# library("rjags")
# library("tidybayes")
#
# jags_results = normal_mixture_jags(
#   likelihood_params,
#   p_mix = p_mix,
#   comp1_prior_params,
#   comp2_prior_params,
#   n.iter = 1000
# )
#
#
#
# jags_results %>%
#   tidy_draws %>%
#   median_qi()
#
# jags_results %>%
#   tidy_draws %>%
#   group_by(.chain) %>%
#   median_qi()
#

posterior_mixture = function(
  likelihood_params,
  p_mix,
  comp1_prior_params,
  comp2_prior_params
){

  # Get posterior mean and variance for each component
  # Component 1
  comp1_closed_form  = do.call(normal_normal_conjugate,
                               c(comp1_prior_params, likelihood_params))
  # Component 2
  comp2_closed_form  = do.call(normal_normal_conjugate,
                               c(comp2_prior_params, likelihood_params))

  # Get mixing proportion
  p_mix_post = get_mixing_prop(
    p_mix = p_mix,
    y_obs = likelihood_params$y_obs,
    mu_P_1 = comp1_prior_params$mu_P,
    mu_P_2 = comp2_prior_params$mu_P,
    sigma_L = with(likelihood_params, make_cov_matrix(sigma_L_11, sigma_L_22, rho_L)),
    sigma_P_1 = with(comp1_prior_params, make_cov_matrix(sigma_P_11, sigma_P_22, rho_P)),
    sigma_P_2 = with(comp2_prior_params, make_cov_matrix(sigma_P_11, sigma_P_22, rho_P))
  )

  return(list(
    mu_post_1 = comp1_closed_form$mu_post,
    sigma_post_1 = comp1_closed_form$sigma_post,
    mu_post_2 = comp2_closed_form$mu_post,
    sigma_post_2 = comp2_closed_form$sigma_post,
    p_mix_post = p_mix_post
  ))

}

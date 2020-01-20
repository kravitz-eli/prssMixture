get_mixing_prop = function(
  likelihood_params,
  p_mix,
  comp1_prior_params,
  comp2_prior_params
){

  browser()

  # Seperately calculate posterior parameters of each component
  # Component 1
  comp1_closed_form  = do.call(normal_normal_conjugate,
                               c(likelihood_params, comp1_prior_params))

  mu_post_1 = comp1_closed_form$mu_post
  sigma_post_1 = comp1_closed_form$sigma_post

  # Component 2
  comp2_closed_form  = do.call(normal_normal_conjugate,
                               c(likelihood_params, comp2_prior_params))

  mu_post_2 = comp2_closed_form$mu_post
  sigma_post_2 = comp2_closed_form$sigma_post

  # Calculate posterir mixing weight
  f1 = mvtnorm::dmvnorm(y_obs, mean = mu_P_1, sigma = sigma_L + sigma_P_1)
  f2 = mvtnorm::dmvnorm(y_obs, mean = mu_P_2, sigma = sigma_L + sigma_P_2)

  p_mix_post = p_mix * f1 / (p_mix * f1 + (1-p_mix) * f2)
  attributes(p_mix_post) = NULL

  return(list(
    mu_post_1 = mu_post_1,
    sigma_post_1 = sigma_post_1,
    mu_post_2 = mu_post_2,
    sigma_post_2 = sigma_post_2,
    p_mix_post = p_mix_post
  ))

}

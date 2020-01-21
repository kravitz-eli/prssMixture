get_mixing_prop = function(
  p_mix,
  y_obs,
  mu_P_1,
  mu_P_2,
  sigma_L,
  sigma_P_1,
  sigma_P_2
){

  # Get density of prior predictive at observed data for each component
  f1 = mvtnorm::dmvnorm(y_obs,
                        mean = mu_P_1,
                        sigma = sigma_L + sigma_P_1)
  f2 = mvtnorm::dmvnorm(y_obs,
                        mean = mu_P_2,
                        sigma = sigma_L + sigma_P_2)


  p_mix_post = p_mix * f1 / (p_mix * f1 + (1-p_mix) * f2)
  attributes(p_mix_post) = NULL

  return(p_mix_post)


}

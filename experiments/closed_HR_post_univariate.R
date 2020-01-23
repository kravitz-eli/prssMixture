# # Scratch work to confirm I'm calculate conjugates correctly
# prior_mean = 0
# n_P = 32
# prior_var = 4 / n_P
#
# # Data values
# L_mean = log(2.25)
# n_L = 45
# L_var = 4/45
#
# # Get posterior mean and var with conjugacy
# post_var = 4 / (n_L + n_P)
# post_sd = sqrt(post_var)
# post_mean = (n_L * L_mean + n_P * prior_mean) / (n_L + n_P)

closed_HR_post = function(prior_N, prior_mean, data_N, data_mean){

  post_var = 4 / (prior_N + data_N)
  post_mean = (prior_N * prior_mean + data_N * data_mean) / (prior_N + data_N)

  return(list("post_mean" = post_mean,
              "post_var" = post_var))
}

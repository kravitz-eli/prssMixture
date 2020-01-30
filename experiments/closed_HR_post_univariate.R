
closed_HR_post = function(prior_N, prior_mean, data_N, data_mean){

  post_var = 4 / (prior_N + data_N)
  post_mean = (prior_N * prior_mean + data_N * data_mean) / (prior_N + data_N)

  return(list("post_mean" = post_mean,
              "post_var" = post_var))
}

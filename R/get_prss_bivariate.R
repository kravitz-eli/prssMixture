get_prss_bivariate = function(
  log_HR_PFS = log(0.75),
  log_HR_OS = log(0.83),
  mu1 = c(0,0),
  var11 = 0.20,
  var22 = 0.20,
  mu2 = c(0, 0),
  var2 = 0.5,
  rho_L = 0.5, # correlation between PFS and OS in the likelihood
  rho_P = 0, # correlation in prior
  mixing_prop = 0.5,
  nE2_PFS = 25,
  nE2_OS = 10,
  nE3_PFS = 90,
  nE3_OS = 70,
  alpha = 0.05,
  Npost = 1e3,
  Ntrial = 1e4
) {

  # Likelihood parameters
  # log_HR_PFS = log(HR2)
  # log_HR_OS = log(PFS2)
  var_PFS = 4/nE2_PFS
  var_OS = 4/nE2_OS
  sigma_likelihood = matrix(
    c(var_PFS, rho_L * sqrt(var_PFS * var_OS),
      rho_L * sqrt(var_PFS * var_OS), var_OS),
    nrow = 2, byrow = TRUE)
  tau_likelihood = solve(sigma_likelihood)

  # Prior Parameters -------------------------------------
  mixing_prop = 1/2
  # First component
  sigma1 =  matrix(
    c(var_PFS, rho_P * sqrt(var_PFS * var_OS),
      rho_L * sqrt(var_PFS * var_OS), var_OS),
    nrow = 2, byrow = TRUE)
  tau1 = solve(sigma1)

  # Second component
  sigma2 =


  # Posterior Parameters ------------------------------
  tau2_post = tau2 + precision_likelihood #posterior precision for second component
  mu2_post = (tau2 * mu2 + precision_likelihood * log_HR_P2)/(tau2 + precision_likelihood) #posterior mean for second component

  # Update mean and precision as a normal conjugate
  tau1_post = tau1 + precision_likelihood #posterior precision for first component
  mu1_post = (tau1 * mu1 + precision_likelihood * log_HR_P2)/(tau1 + precision_likelihood) #posterior mean for first component
  # Calculate the part of the new weight that is updated
  c1 = dnorm(log_HR_P2, mean = mu1, sd = sqrt(1/precision_likelihood + 1/tau1)) #used to calculate posterior mixing proprtion for 1st component

  # used to calculate posterior mixing proprtion for 2nd component
  c2 = dnorm(log_HR_P2, mean = mu2, sd = sqrt(1/precision_likelihood + 1/tau2))
  # posterior mixing proportion
  p_post = p * c1 / (p * c1 + (1 - p) * c2)


  # Sample from posterior --------------------------------
  # Draw random variables for assignments to each mixture
  components = sample(1:2, prob = c(p_post, 1 - p_post), size = Npost, replace = TRUE)
  mus = c(mu1_post, mu2_post)
  sds = c(tau1_post, tau2_post)^(-1/2)
  posterior_HR_samples = rnorm(n = Npost, mean = mus[components], sd = sds[components])


  # Simulate trials  -------------------
  simulated_HR3 = rnorm(n = Ntrial, mean = posterior_HR_samples, sd = sqrt(4/nE3))

  # See which trials were successful ------------------
  CV <- qnorm(alpha, 0, sqrt(4 / nE3))
  PrSS = mean(simulated_HR3 < CV)

  return(PrSS)


}

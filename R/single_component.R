# Phase 2 information -----
log_HR_PFS = log(0.65)
log_HR_OS = log(0.88)

n = 99
nE2_PFS = 31
nE2_OS = 06

# Likelihood parameters -----
y_obs = c(log_HR_PFS, log_HR_OS)
var_PFS = 4/nE2_PFS
var_OS = 4/nE2_OS
rho_L = 0.65
sigma_likelihood = matrix(
  c(var_PFS, rho_L * sqrt(var_PFS * var_OS),
    rho_L * sqrt(var_PFS * var_OS), var_OS),
  nrow = 2, byrow = TRUE)
tau_likelihood = solve(sigma_likelihood)

# Prior Parameters -------------------------------------
mu1 = c(0 , 0)
rho_P = 0.5
sigma_11 = 0.10
sigma_22 = 0.20

sigma1 =  matrix(
  c(sigma_11, rho_P * sqrt(sigma_11 * sigma_22),
    rho_P * sqrt(sigma_11 * sigma_22), sigma_22),
  nrow = 2, byrow = TRUE)

tau1 = solve(sigma1)

# Posterior parameters ----------------------------------------
# Precision
tau1_post = tau1 + tau_likelihood
# Mean
mu1_post = solve(tau1_post) %*% (tau1 %*% mu1 + tau_likelihood %*% y_obs)
mu1_post = as.vector(mu1_post)

# See if this answer matches JAGs -------------------------------------------
# Setup R jags
library("rjags")

n.iter = 1e6
burn_in = 1e4

bivariate_norm_model = "model{
  # Likelihood  -----------------------------------
  y ~ dmnorm(theta[1:2], tau_L)

  # Priors -----------------------------------------
  # Multivariate normal
  theta ~ dmnorm(mu[1:2], tau_P)
  tau_P <-inverse(sigma_P)
  sigma_P[1,1] <- sigma_11
  sigma_P[1,2] <- rho_P * sqrt(sigma_11 * sigma_22)
  sigma_P[2,1] <- rho_P * sqrt(sigma_11 * sigma_22)
  sigma_P[2,2] <- sigma_22
}"


data = list(
  y = y_obs,
  tau_L = tau_likelihood,
  mu = mu1,
  sigma_11 = sigma_11,
  sigma_22 = sigma_22,
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
                        variable.names = c("theta", "tau_P"),
                        n.iter = n.iter)

# Get posterior samples for parameters
theta = as.matrix(samples[ ,c("theta[1]", "theta[2]")])
tau_11 = as.matrix(samples[ , "tau_P[1,1]"])
tau_22 = as.matrix(samples[ , "tau_P[2,2]"])
tau_12 = as.matrix(samples[ , "tau_P[1,2]"])


post_mean_jags = colMeans(theta)
post_prec_PFS_jags = mean(tau_11)
post_prec_OS_jags = mean(tau_22)
post_prec_off_jags = mean(tau_12)


# See if parameters match
diffs = c(
  "post_mean_diff" = mu1_post - post_mean_jags,
  "post_prec_PFS_diff" = tau1[1, 1] - post_prec_PFS_jags,
  "post_prec_OS_diff" = tau1[2, 2] - post_prec_OS_jags,
  "post_prec_off_diff" = tau1[1, 2] - post_prec_off_jags
)

print(diffs, digits = 3)
max(abs(diffs))


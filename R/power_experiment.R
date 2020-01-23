# Set constants for phase II
n_P2 = 40
theta_P2 = log(0.5)


# Set constants for phase III -------------------------
n_P3 = 88
theta_0 = log(1)
# theta_1 = log(0.5)
theta_1 = theta_P2

alpha = 0.025
crit_val = qnorm(alpha,
                 mean = theta_0,
                 sd = sqrt(4 / n_P3))

# Calculate regular power ----------------------
# probability a N(theta_1, sigma^2) < critical value
power_analytical = pnorm(crit_val,
      mean = theta_1,
      sd = sqrt(4/88))

print(power_analytical)

# Calculate power by simulation -----------------------
set.seed(3)
n_samp = 1e4
power_sample = rnorm(n_samp, theta_1, sqrt(4/88))
power_simulation = mean(power_sample < crit_val)

std_error = sqrt(power_simulation * (1 - power_simulation) / n_samp)

print(c("estimated power" = power_simulation, "standard error" = std_error))

# Calculate PrSS by simulation -----------------------------------------------
source(here::here("experiments", "closed_HR_post_univariate.R"))
prior_mean = 0
n_prior = 10
prior_var = 4 / n_prior

closed = closed_HR_post(prior_N = n_prior,
               prior_mean = prior_mean,
               data_N = n_P2 ,
               data_mean = theta_P2)

post_mean = closed$post_mean
post_var = closed$post_var

# Prss Calculaton
n_samp = 1e5
n_p3 = 88
theta_0 = 0
crit_val = qnorm(alpha,
                 mean = theta_0,
                 sd = sqrt(4/n_p3))

post_samples = rnorm(n_samp, post_mean, sqrt(post_var))
fake_data = rnorm(n_samp, post_samples, sqrt(4/n_p3))
mean(fake_data < crit_val)

# Try analytical version of posterior
n_post = 4 / post_var
prss_mean = sqrt(n_post / n_p3) * post_mean
# Also can be 4/n_p3
prss_var = n_post / n_p3 * post_var

# Doesn't quite work.
pnorm(crit_val,
      mean = prss_mean,
      sd = prss_var)

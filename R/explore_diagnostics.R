# See posterior of mixing proportion
jags_results_prior_p_mix %>%
  spread_draws(ind) %>%
  summarise(mean(ind), median(ind), min(ind), max(ind))

post_weight = jags_results_prior_p_mix %>%
  spread_draws(weight, ind) %>%
  pull(weight)

ggplot(data.frame(post_weight), aes(x = post_weight)) +
  geom_histogram(aes(y =..density..), colour="black", fill="white")+
  geom_density(alpha = 0.40, fill="red") +
  theme_minimal()


prior_var = c(1e-1, 1, 1e2, 1e3, 1e4, 1e5, 1e6)
jags_results = vector("list", length(prior_var))
data[["p_mix"]] = NULL
for (i in seq_along(prior_var)) {

  data[c("P2.sigma_1", "P2.sigma_2")] = prior_var[[i]]

  jags_results[[i]] <- run_jags(
    model_file = "normal_mixture_prior_p",
    data = data,
    track_variable_names = c("theta", "ind", "weight"),
    iter = 1e4,
    burn = 1e3,
    chains = 2,
    progress.bar = "text"
  )

}

# Get individual posterior distributions ------------------------------------
weight_post = jags_results %>%
  map(spread_draws, weight, ind) %>%
  map(pull, weight)

ind_post = jags_results %>%
  map(spread_draws, weight, ind) %>%
  map(pull, ind)

theta_post = jags_results %>%
  map(as.matrix) %>%
  map(~.x[, c("theta[1]", "theta[2]")])

theta_1_post = theta_post %>%
  map(~.x[, "theta[1]"])

theta_2_post = theta_post %>%
  map(~.x[, "theta[2]"])

# Get posterior means and medians --------------------------------
weights_post_means = map_dbl(weight_post, mean)

ind_post_mean = map_dbl(ind_post, mean)
theta1_post_mean = map_dbl(theta_1_post, mean)
theta2_post_mean = map_dbl(theta_2_post, mean)


par(mfrow = c(1,2))
plot(
  x = log(prior_var, 10),
  y = weights_post_means,
  type = "b",
  ylab = "Posterior Mean",
  xlab = "log10 Prior Variance",
  main = "Posterior Mixing Success Prop",
  lwd = 3
)

plot(
  x = log(prior_var, 10),
  y = ind_post_mean,
  type = "b",
  ylab = "Posterior Mean",
  xlab = "log10 Prior Variance",
  main = "Posterior Mixing Success Prop",
  lwd = 3
)




library(coda)
traceplot(jags_results[[5]])

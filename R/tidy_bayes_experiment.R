data[["P1.rho"]] = 0.5
data[["P2.rho"]] = 0.5

jags_results <- run_jags(
  model_file = "normal_mixture_fixed_rho",
  data = data,
  track_variable_names = c("theta", "ind",
                           "P1.mu", "P2.mu"),
  iter = 1e5,
  burn = 1e3,
  chains = 1,
  progress.bar = "text"
)


summary(jags_result)

draws = jags_results %>%
  spread_draws(ind, P1.mu[1])

colnames(draws)[5] = "mu_p1"

draws %<>%
  mutate(
    b = map_chr(mu_p1, paste, collapse = ",")
  ) %>%
  tidyr::separate(
    b,
    c("p1_mu1", "p1_mu2"),
    ","
  ) %>%
  mutate(
    p1_mu1 = as.numeric(p1_mu1),
    p1_mu2 = as.numeric(p1_mu2)
  ) %>%
  select(-mu_p1)


colMeans(draws)

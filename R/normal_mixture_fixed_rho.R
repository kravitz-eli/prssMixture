model{
  # Likelihood (phase II OR and PFS) -----------------------------------
  L.y_obs ~ dmnorm(theta[1:2], L.tau)
  L.tau <-inverse(L.sigma)
  L.sigma[1,1] <- L.sigma_1
  L.sigma[1,2] <- L.rho * sqrt(L.sigma_1 * L.sigma_2)
  L.sigma[2,1] <- L.rho * sqrt(L.sigma_1 * L.sigma_2)
  L.sigma[2,2] <- L.sigma_2

  # Mixture Structure  ------------------------------------------------------
  theta <- ind * component_1 + (1-ind) * component_2
  ind ~ dbern(p_mix)

  # Priors -------------------------------------------------------------------
  # First component in normal mixture --------------------------
  component_1 ~ dmnorm(P1.mu[1:2], P1.tau)
  P1.tau <-inverse(P1.sigma)
  P1.sigma[1,1] <- P1.sigma_1
  P1.sigma[1,2] <- P1.rho * sqrt(P1.sigma_1 * P1.sigma_2)
  P1.sigma[2,1] <- P1.rho * sqrt(P1.sigma_1 * P1.sigma_2)
  P1.sigma[2,2] <- P1.sigma_2
  P1.rho ~ dbeta(P1.alpha, P1.beta)

  # Second component in normal mixture -------------------
  component_2 ~ dmnorm(P2.mu[1:2], P2.tau)
  P2.tau <-inverse(P2.sigma)
  P2.sigma[1,1] <- P2.sigma_1
  P2.sigma[1,2] <- P2.rho * sqrt(P2.sigma_1 * P2.sigma_2)
  P2.sigma[2,1] <- P2.rho * sqrt(P2.sigma_1 * P2.sigma_2)
  P2.sigma[2,2] <- P2.sigma_2
  P2.rho ~ dbeta(P2.alpha, P2.beta)
}

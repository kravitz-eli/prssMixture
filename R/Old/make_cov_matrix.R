# Take the variance terms and correlation and make a covariance matrix
make_cov_matrix <- function(var_1, var_2, rho) {

  matrix(c(
    var_1, rho * sqrt(var_1 * var_2),
    rho * sqrt(var_1 * var_2), var_2
  ),
  nrow = 2, ncol = 2
  )
}

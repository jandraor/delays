mean_generation_time <- function(j, sigma_inv, gamma_inv) {
  sigma_inv + ( gamma_inv * (j + 1) / (2 * j))
}

calculate_infectious_period <- function(tau, sigma_inv, j) {
  j_factor <- (j + 1) / (2 *j)
  
  (tau - sigma_inv) / j_factor
}

calculate_beta <- function(R0, tau, sigma_inv, j) {
  
  j_factor <- (j + 1) / (2 *j)
  
  R0 * (j_factor/( tau - sigma_inv))
}
#' @keywords internal
c_plateau_gaussian_density <- function(delta, sigma) {
  sqrt_2pi <- sqrt(2 * pi)
  (sqrt_2pi * sigma / 2) + delta
}

#' Computes the normalised plateau–Gaussian density used for y = 0
#' @param x vector of x-values
#' @param mu mean location of the Gaussian tail
#' @param sigma transition spread
#' @param L lower truncation
#' @export
plateau_gaussian_density_left <- function(x, mu, sigma, L) {
  sigma_sq <- sigma^2
  C <- c_plateau_gaussian_density(abs(mu - L), sigma)

  density <- numeric(length(x))

  below_L <- x < L
  below_mu <- (x >= L) & (x < mu)
  above_mu <- x >= mu

  density[below_L] <- 0
  density[below_mu] <- 1
  density[above_mu] <- exp(-0.5 * ((x[above_mu] - mu)^2 / sigma_sq))

  density / C
}

#' Computes the normalised plateau–Gaussian density used for y = 1
#' @param x vector of x-values
#' @param mu mean of the Gaussian tail
#' @param sigma transition spread
#' @param U upper truncation
#' @export
plateau_gaussian_density_right <- function(x, mu, sigma, U) {
  sigma_sq <- sigma^2
  C <- c_plateau_gaussian_density(abs(mu - U), sigma)

  density <- numeric(length(x))

  above_U <- x > U
  above_mu <- (x <= U) & (x > mu)
  below_mu <- x <= mu

  density[above_U] <- 0
  density[above_mu] <- 1
  density[below_mu] <- exp(-0.5 * ((x[below_mu] - mu)^2 / sigma_sq))

  density / C
}


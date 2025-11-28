# internal helpers mirroring Stan functions ----

.uniform_gaussian_normalising_constant <- function(delta, sigma) {
  sqrt_2pi <- sqrt(2 * pi)
  (sqrt_2pi * sigma / 2) + delta
}

.uniform_left_gaussian_log_density_unnormalised <- function(x, mu, sigma, L) {
  ifelse(
    x <= mu,
    0,
    -0.5 * ((x - mu) / sigma)^2
  )
}

.uniform_right_gaussian_log_density_unnormalised <- function(x, mu, sigma, U) {
  ifelse(
    x >= mu,
    0,
    -0.5 * ((x - mu) / sigma)^2
  )
}

#' Predict class probabilities from a STAGE model
#'
#' @param object A `stage_fit` object as returned by [fit_stage()].
#' @param newdata Numeric vector of new x values.
#' @param group Optional group for each new x (length 1 or `length(newdata)`).
#'   If `NULL` and there are multiple populations in the model, population 1
#'   is used for all predictions.
#' @param type `"prob"` for P(y = 1 | x), `"class"` for 0/1 classification.
#' @param ... Currently ignored.
#'
#' @return Numeric vector of probabilities or classes.
#' @export
predict.stage_fit <- function(object, newdata, group = NULL,
                              type = c("prob", "class"), ...) {
  type <- match.arg(type)
  x_new <- as.numeric(newdata)
  N_new <- length(x_new)

  # cmdstanr::summary() output: columns include "variable", "mean", ...
  summ <- object$fit$summary(
    variables = c("mu_m50", "m50_pop", "d", "sigma_x")
  )

  # Convenience helper to grab means by variable name
  get_mean <- function(name) {
    summ$mean[summ$variable == name]
  }

  J <- object$data$J
  L <- object$data$L
  U <- object$data$U

  d       <- get_mean("d")
  sigma_x <- get_mean("sigma_x")

  # Storage
  lp0 <- numeric(N_new)
  lp1 <- numeric(N_new)

  if (J == 1L) {
    # Single-population model: use global mu_m50
    m50 <- get_mean("mu_m50")
    mu0 <- m50 - d / 2
    mu1 <- m50 + d / 2

    C0 <- .uniform_gaussian_normalising_constant(mu0 - L, sigma_x)
    C1 <- .uniform_gaussian_normalising_constant(U - mu1, sigma_x)

    lp0 <- .uniform_left_gaussian_log_density_unnormalised(x_new, mu0, sigma_x, L) - log(C0)
    lp1 <- .uniform_right_gaussian_log_density_unnormalised(x_new, mu1, sigma_x, U) - log(C1)

  } else {
    # Multi-population model: need a population index for each new x
    if (is.null(group)) {
      pop_new <- rep(1L, N_new)
    } else {
      pop_new <- as.integer(as.factor(group))
      if (length(pop_new) == 1L && N_new > 1L) {
        pop_new <- rep(pop_new, N_new)
      }
      if (length(pop_new) != N_new) {
        stop("Length of `group` must be 1 or length(newdata).")
      }
    }

    for (j in 1:J) {
      idx <- which(pop_new == j)
      if (length(idx) == 0L) next

      m50_j <- get_mean(sprintf("m50_pop[%d]", j))
      mu0_j <- m50_j - d / 2
      mu1_j <- m50_j + d / 2

      C0_j <- .uniform_gaussian_normalising_constant(mu0_j - L, sigma_x)
      C1_j <- .uniform_gaussian_normalising_constant(U - mu1_j, sigma_x)

      lp0[idx] <- .uniform_left_gaussian_log_density_unnormalised(
        x_new[idx], mu0_j, sigma_x, L
      ) - log(C0_j)

      lp1[idx] <- .uniform_right_gaussian_log_density_unnormalised(
        x_new[idx], mu1_j, sigma_x, U
      ) - log(C1_j)
    }
  }

  # Convert log-densities to probabilities via 2-class softmax
  max_lp <- pmax(lp0, lp1)
  p1 <- exp(lp1 - max_lp) / (exp(lp0 - max_lp) + exp(lp1 - max_lp))

  if (type == "prob") return(p1)
  as.integer(p1 > 0.5)
}

#' Construct prior hyperparameters for the STAGE model
#'
#' @param x,y Optional data vectors to set sensible defaults.
#' @param prior_mu_m50_mu Mean of the prior on the global m50.
#' @param prior_mu_m50_tau SD of the prior on the global m50.
#' @param prior_d_mu Mean of the prior on d.
#' @param prior_d_tau SD of the prior on d.
#' @param prior_sigma_x_tau SD scale for sigma_x.
#' @param prior_sigma_alpha_tau SD scale for sigma_alpha (between-pop).
#'
#' @return A list with elements matching the Stan data block.
#' @export
stage_priors <- function(
    x = NULL,
    y = NULL,
    prior_mu_m50_mu  = if (!is.null(x)) stats::median(x) else 0,
    prior_mu_m50_tau = if (!is.null(x)) diff(range(x)) / 6 else 100,
    prior_d_mu       = if (!is.null(x)) diff(range(x)) / 4 else 0,
    prior_d_tau      = if (!is.null(x)) diff(range(x)) / 4 else 100,
    prior_sigma_x_tau      = if (!is.null(x)) stats::sd(x) else 100,
    prior_sigma_alpha_tau  = if (!is.null(x)) diff(range(x)) / 6 else 100
) {
  list(
    prior_mu_m50_mu  = prior_mu_m50_mu,
    prior_mu_m50_tau = prior_mu_m50_tau,
    prior_d_mu       = prior_d_mu,
    prior_d_tau      = prior_d_tau,
    prior_sigma_x_tau     = prior_sigma_x_tau,
    prior_sigma_alpha_tau = prior_sigma_alpha_tau
  )
}

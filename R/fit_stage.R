#' Fit the STAGE (uniform + Gaussian) generative model
#'
#' @param x Numeric vector of predictor values (e.g. length).
#' @param y Binary vector: 0 = immature, 1 = mature.
#' @param group Optional grouping variable (population). If NULL, all data
#'   are treated as one population (J = 1).
#' @param L,U Truncation bounds for x. If NULL, use min(x) and max(x).
#' @param priors A list of prior hyperparameters created by [stage_priors()].
#' @param chains,iter,... Passed to [cmdstanr::sample()].
#'
#' @return An object of class "stage_fit".
#' @export
fit_stage <- function(
    x,
    y,
    group   = NULL,
    L = NULL,
    U = NULL,
    priors  = stage_priors(x, y),
    chains  = 4,
    iter    = 2000,
    ...
) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg cmdstanr} is required for fit_stage().")
  }

  x <- as.numeric(x)
  y <- as.integer(y)

  if (!all(y %in% c(0L, 1L))) {
    cli::cli_abort("{.arg y} must be binary (0/1).")
  }

  N <- length(x)
  if (length(y) != N) {
    cli::cli_abort("x and y must have the same length.")
  }


  # Truncation bounds
  if (is.null(L)) L <- min(x)
  if (is.null(U)) U <- max(x)
  if (L >= U) cli::cli_abort("`L` must be < `U`.")

  # Group / populations
  if (is.null(group)) {
    pop <- rep(1L, N)
    J   <- 1L
  } else {
    pop_f <- as.factor(group)
    pop   <- as.integer(pop_f)
    J     <- length(levels(pop_f))
  }

  # Build Stan data list: all arguments explicit
  data_list <- c(
    list(
      N   = N,
      J   = J,
      L   = L,
      U   = U,
      x   = x,
      y   = y,
      pop = pop
    ),
    priors   # adds prior_mu_m50_mu, prior_mu_m50_tau, etc.
  )

  mod <- get_cmdstan_model("stage_Jpop.stan")

  fit <- mod$sample(
    data          = data_list,
    chains        = chains,
    iter_sampling = iter,
    ...
  )

  out <- list(
    fit   = fit,
    data  = data_list,
    call  = match.call()
  )
  class(out) <- "stage_fit"
  out
}

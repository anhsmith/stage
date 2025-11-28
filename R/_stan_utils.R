#' Get compiled cmdstanr model
#' @keywords internal
get_cmdstan_model <- function(stan_file) {
  stopifnot(requireNamespace("cmdstanr", quietly = TRUE))

  # stan_file is relative to inst/stan
  stan_path <- system.file("stan", stan_file, package = "stage", mustWork = TRUE)

  cmdstanr::cmdstan_model(stan_path)
}

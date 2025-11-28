#' Get a compiled cmdstanr model from inst/stan
#'
#' @keywords internal
get_cmdstan_model <- function(stan_file) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr must be installed to compile Stan models.")
  }

  # location inside inst/stan
  stan_path <- system.file("stan", stan_file, package = "stage", mustWork = TRUE)

  cmdstanr::cmdstan_model(stan_path)
}

#' Generic transition point extractor
#'
#' @param object A model fit object.
#' @param ... Further arguments passed to methods.
#'
#' @return Typically a numeric vector or list summarising the transition point.
#' @export
transition_point <- function(object, ...) {
  UseMethod("transition_point")
}

#' Posterior summaries of STAGE transition points (m50)
#'
#' Returns posterior summaries of \eqn{m_{50}}, the value of \eqn{x} where
#' \eqn{P(y=1 | x) = 0.5}. Under equal class priors, \eqn{m_{50}} is the
#' direct sampling parameter of the STAGE model.
#'
#' @param object A `stage_fit` object as returned by [fit_stage()].
#' @param population Optional population index (integer) or label (character).
#'   If `NULL`, return global and all populations.
#' @param ... Currently ignored. Included for method compatibility.
#'
#' @return A named list. Always contains `global` (posterior summary of the
#'   global `m50`). For J > 1 models, also contains `pop1`, `pop2`, etc. for
#'   each population-specific transition point. Each summary is a named numeric
#'   vector with `mean`, `median`, `q2.5`, `q97.5`.
#' @export
transition_point.stage_fit <- function(object, population = NULL, ...) {
  J <- object$data$J

  draws <- object$fit$draws(
    variables = c("m50", "m50_pop"),
    format = "draws_df"
  )

  posterior_summary <- function(z) {
    c(
      mean   = mean(z),
      median = stats::median(z),
      q2.5   = stats::quantile(z, 0.025),
      q97.5  = stats::quantile(z, 0.975)
    )
  }

  # specific population requested
  if (!is.null(population)) {
    if (is.character(population)) {
      cli::cli_abort("Character population labels not yet implemented; use integer index.")
    }
    j <- as.integer(population)
    if (J == 1L) {
      return(list(global = posterior_summary(draws$m50)))
    }
    return(list(
      global = posterior_summary(draws$m50),
      pop    = posterior_summary(draws[[paste0("m50_pop[", j, "]")]])
    ))
  }

  # No population specified: global + all populations (J > 1)
  out <- list(global = posterior_summary(draws$m50))
  if (J > 1L) {
    for (j in seq_len(J)) {
      out[[paste0("pop", j)]] <- posterior_summary(draws[[paste0("m50_pop[", j, "]")]])
    }
  }
  out
}

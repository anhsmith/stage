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
#' @param object A `stage_fit` object as returned by [fit_stage()].
#' @param population Optional population index (integer) or label (character).
#'   If `NULL`, return global and all populations.
#' @param ... Currently ignored. Included for method compatibility.
#'
#' @return A named list or numeric vector with posterior summary stats.
#' @export
transition_point.stage_fit <- function(object, population = NULL, ...) {
  draws <- object$fit$draws(
    variables = c("mu_m50", "m50_pop"),
    format    = "draws_df"
  )

  J <- object$data$J

  # helper
  posterior_summary <- function(z) {
    c(
      mean   = mean(z),
      median = stats::median(z),
      q2.5   = stats::quantile(z, 0.025),
      q97.5  = stats::quantile(z, 0.975)
    )
  }

  # global only (or no populations in model)
  if (J == 1L && is.null(population)) {
    return(posterior_summary(draws$mu_m50))
  }

  # specific population requested?
  if (!is.null(population)) {
    if (is.character(population)) {
      # if you later store the group levels, you can map label -> index here
      stop("Character population labels not yet implemented; use integer index.")
    }
    j <- as.integer(population)
    colname <- paste0("m50_pop[", j, "]")
    return(posterior_summary(draws[[colname]]))
  }

  # otherwise: global + all populations
  out <- list(
    global = posterior_summary(draws$mu_m50)
  )
  for (j in 1:J) {
    colname <- paste0("m50_pop[", j, "]")
    out[[paste0("pop", j)]] <- posterior_summary(draws[[colname]])
  }
  out
}

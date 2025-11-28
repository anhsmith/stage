#' @export
print.stage_fit <- function(x, ...) {
  cat("STAGE model fit\n\n")
  print(x$fit$summary(
    variables = c("mu_m50", "d", "sigma_x", "sigma_alpha"),
    ...
  ))
  invisible(x)
}

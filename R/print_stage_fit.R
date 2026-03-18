#' @export
print.stage_fit <- function(x, ...) {
  cat("STAGE model fit\n\n")
  J <- x$data$J
  vars <- if (J == 1L) {
    c("m50", "d", "sigma_x")
  } else {
    c("m50", "m50_pop", "d", "sigma_x", "sigma_alpha")
  }
  print(x$fit$summary(variables = vars, ...))
  invisible(x)
}

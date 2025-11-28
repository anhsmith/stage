test_that("fit_stage runs and returns sensible object", {
  set.seed(123)
  N <- 100
  L <- 1000
  U <- 1500
  x <- runif(N, L, U)
  y <- rbinom(N, 1, plogis((x - 1250)/50))

  fit <- fit_stage(x, y, L = L, U = U, chains = 1, iter = 500)

  expect_s3_class(fit, "stage_fit")
  tp <- transition_point(fit)
  expect_true(is.numeric(tp["median"]))
  p <- predict(fit, x, type = "prob")
  expect_length(p, N)
  expect_true(all(p >= 0 & p <= 1))
})

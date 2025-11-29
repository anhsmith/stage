# Fit the STAGE (uniform + Gaussian) generative model

Fit the STAGE (uniform + Gaussian) generative model

## Usage

``` r
fit_stage(
  x,
  y,
  group = NULL,
  L = NULL,
  U = NULL,
  priors = stage_priors(x, y),
  chains = 4,
  iter = 2000,
  ...
)
```

## Arguments

- x:

  Numeric vector of predictor values (e.g. length).

- y:

  Binary vector: 0 = immature, 1 = mature.

- group:

  Optional grouping variable (population). If NULL, all data are treated
  as one population (J = 1).

- L, U:

  Truncation bounds for x. If NULL, use min(x) and max(x).

- priors:

  A list of prior hyperparameters created by
  [`stage_priors()`](https://anhsmith.github.io/stage/reference/stage_priors.md).

- chains, iter, ...:

  Passed to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).

## Value

An object of class "stage_fit".

# Construct prior hyperparameters for the STAGE model

Construct prior hyperparameters for the STAGE model

## Usage

``` r
stage_priors(
  x = NULL,
  y = NULL,
  prior_mu_m50_mu = if (!is.null(x)) stats::median(x) else 0,
  prior_mu_m50_tau = if (!is.null(x)) diff(range(x))/6 else 100,
  prior_d_mu = if (!is.null(x)) diff(range(x))/4 else 0,
  prior_d_tau = if (!is.null(x)) diff(range(x))/4 else 100,
  prior_sigma_x_tau = if (!is.null(x)) stats::sd(x) else 100,
  prior_sigma_alpha_tau = if (!is.null(x)) diff(range(x))/6 else 100
)
```

## Arguments

- x, y:

  Optional data vectors to set sensible defaults.

- prior_mu_m50_mu:

  Mean of the prior on the global m50.

- prior_mu_m50_tau:

  SD of the prior on the global m50.

- prior_d_mu:

  Mean of the prior on d.

- prior_d_tau:

  SD of the prior on d.

- prior_sigma_x_tau:

  SD scale for sigma_x.

- prior_sigma_alpha_tau:

  SD scale for sigma_alpha (between-pop).

## Value

A list with elements matching the Stan data block.

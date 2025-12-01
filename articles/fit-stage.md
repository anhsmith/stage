# Fitting the STAGE Model

This vignette shows how to:

- fit the STAGE model to a simple simulated dataset,
- extract posterior transition points, and
- extend to a hierarchical (multi-population) model.

## Simulating data

To demonstrate the workflow, we simulate a simple dataset of lengths in
$\lbrack L,U\rbrack$ with a smooth transition around a true $m_{50}$.

``` r
N <- 200
L <- 1000
U <- 1500
true_m50 <- 1250
true_d   <- 100

# Simulate covariate
x <- runif(N, L, U)

# Smooth transition in P(mature)
p <- plogis((x - true_m50) / (true_d / 4))

# Binary maturity indicator
y <- rbinom(N, 1, p)

head(tibble(x = x, y = y))
#> # A tibble: 6 × 2
#>       x     y
#>   <dbl> <int>
#> 1 1133.     0
#> 2 1186.     0
#> 3 1286.     1
#> 4 1454.     1
#> 5 1101.     0
#> 6 1449.     1
```

## Fitting the STAGE model

We now fit the STAGE model using
[`fit_stage()`](https://anhsmith.github.io/stage/reference/fit_stage.md).  
Because STAGE is a Bayesian model fit via Stan, we keep `chains` and
`iter` small in this vignette to keep runtime acceptable; for real
applications you should increase these.

``` r
fit <- fit_stage(x, y, L = L, U = U, chains = 2, iter = 1000)
fit
```

You can inspect the posterior with standard tools such as
`summary(fit$fit)`, `bayesplot`, or `posterior` utilities if you want
more diagnostics.

### Extracting the transition point

The helper
[`transition_point()`](https://anhsmith.github.io/stage/reference/transition_point.md)
returns posterior summaries of the global transition point $m_{50}$.

``` r
transition_point(fit)
```

Typical output is a small data frame with posterior mean, median, and
credible intervals for `m50`.

### Predicting maturity probabilities

We can also predict $\Pr\left( \text{mature} \mid x \right)$ over a grid
of lengths, and plot the posterior mean transition curve.

``` r
x_grid <- seq(L, U, length.out = 200)

p_hat <- predict(fit, x_grid, type = "prob")

plot(
  x_grid, p_hat, type = "l",
  xlab = "Length (x)",
  ylab = "P(mature)",
  las = 1
)
```

## Hierarchical (multi-population) model

In many applications we have multiple populations (e.g. regions, stocks,
or sexes) and we expect:

- a **shared average transition point** $\mu_{m50}$, and
- **population-specific deviations**.

STAGE implements this using a hierarchical structure
$$m50_{\text{pop}{\lbrack j\rbrack}} = \mu_{m50} + z_{j}\,\sigma_{\alpha},$$
with a non-centred parameterisation for efficient sampling.

A simple two-group example:

``` r
group <- rep(c("A", "B"), each = N/2)

fit2 <- fit_stage(
  x, y,
  group = group,
  L = L, U = U,
  chains = 2, iter = 1000
)

transition_point(fit2)
```

`transition_point(fit2)` returns both:

- a global `m50` (overall mean transition)  
- group-specific `m50_pop[j]` values.

## Summary

In this vignette we:

- introduced the plateau–Gaussian likelihood underlying STAGE,
- visualised the class densities and resulting transition probabilities,
- simulated a simple dataset,
- fit the STAGE model to estimate the transition point, and
- extended the model to a hierarchical multi-population setting.

In practice, you will want to:

- choose $\lbrack L,U\rbrack$ to capture the transition zone plus some
  plateau on each side,
- check posterior diagnostics (divergences, R-hat, effective sample
  sizes),
- and interpret $m_{50}$ and its uncertainty in the biological context.

Future vignettes will expand on:

- comparing STAGE to logistic regression and Gaussian LDA, and
- sensitivity to the choice of $\lbrack L,U\rbrack$ and prior settings.

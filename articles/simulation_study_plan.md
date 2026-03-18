# Simulation Study Design

> **Status: planned.** This vignette describes the design of a
> simulation study that is not yet complete. It serves as a public
> record of the design and will be updated with results as they become
> available. The simulation study will form the basis of a methods paper
> describing and validating the STAGE model.

------------------------------------------------------------------------

## Overview

The goal of the simulation study is to evaluate the performance of the
STAGE model for estimating biological transition points (e.g. length at
50% maturity, $m_{50}$) against competing methods, under a range of
conditions reflecting real-world data collection challenges. The study
should demonstrate where STAGE performs well, where it degrades
gracefully, and where its assumptions are violated.

------------------------------------------------------------------------

## Methods compared

All methods described below are implemented in the `stage` package,
providing a consistent interface and comparable Bayesian inference
across methods. $m_{50}$ is extracted from all methods as the 50%
crossing point under equal class priors ($\pi_{0} = \pi_{1} = 0.5$).

``` r
library(stage)

# STAGE — plateau-Gaussian generative classifier (this package)
fit_stage(x, y, L = L, U = U)

# Bayesian LDA — equal-sigma Gaussian class-conditionals, pi fixed at 0.5
fit_lda(x, y)

# Bayesian logistic regression — standard discriminative model
fit_logistic(x, y)

# Bayesian logistic regression with full inverse-probability weights
fit_logistic(x, y, weights = "full")

# Bayesian logistic regression with square-root inverse-probability weights
fit_logistic(x, y, weights = "sqrt")
```

Frequentist logistic regression (MLE) will also be included as a
baseline, fitted via [`glm()`](https://rdrr.io/r/stats/glm.html).

------------------------------------------------------------------------

## Data-generating process

### Biological truth

Maturity status is assigned probabilistically via a logistic transition
function of length $x$:

$$\Pr(y = 1 \mid x) = \text{logistic}\left( \text{slope} \times \left( x - m_{50} \right) \right)$$

This is a **discriminative truth** — $\Pr(y = 1 \mid x)$ is defined
directly rather than via a generative model. STAGE is therefore tested
under mild model misspecification, which is the relevant real-world
condition.

A secondary set of scenarios will simulate from the STAGE generative
model directly (drawing $x$ from the plateau-Gaussian
class-conditionals) to test performance under the correctly specified
model. These scenarios are clearly labelled.

### Baseline true parameters

``` r
m50_true   <- 1250   # mm — true transition point
slope_true <- 0.08   # logistic slope; roughly equivalent to d=150, sigma=50
L_pop      <- 200    # true lower bound of population length distribution
U_pop      <- 2000   # true upper bound
```

### Sampling model

Real samples are not iid draws from the population length distribution.
Sampling is patchy: fish are taken in discrete events (hauls, trawl
sets, survey stations) from locations that may favour a particular size
class; gear selectivity imposes length-dependent catchability; and
cohort structure means fish from a single event may be highly similar in
length.

A two-level sampling model captures this non-independence:

``` r
simulate_clustered <- function(
    n_events,       # number of fishing events
    n_per_event,    # fish per event (scalar or vector)
    m50, slope,     # true biological parameters
    event_mu_mean,  # mean of event-level length modes
    event_mu_sd,    # SD across events (spatial/cohort variation)
    length_sd,      # within-event SD (individual variation)
    L, U            # analyst-specified truncation bounds
) {
  event_means <- rnorm(n_events, mean = event_mu_mean, sd = event_mu_sd)
  x <- unlist(mapply(
    function(mu, n) rnorm(n, mean = mu, sd = length_sd),
    mu = event_means,
    n  = rep(n_per_event, length.out = n_events)
  ))
  p <- plogis(slope * (x - m50))
  y <- rbinom(length(x), size = 1, prob = p)
  keep <- x >= L & x <= U
  tibble::tibble(x = x[keep], y = y[keep])
}
```

The degree of clustering is controlled by `event_mu_sd` relative to
`length_sd`. When `event_mu_sd` is small, the sample is approximately
iid.

A simpler iid baseline is also used:

``` r
simulate_simple <- function(n0, n1, m50, slope, L, U) {
  x <- c(runif(n0, L, m50), runif(n1, m50, U))
  p <- plogis(slope * (x - m50))
  y <- rbinom(n0 + n1, size = 1, prob = p)
  tibble::tibble(x = x, y = y)
}
```

------------------------------------------------------------------------

## Simulation factors

### Factor 1: Sample size

Total $N$, approximately equal classes at baseline:

``` r
sample_sizes <- c(50, 100, 200, 500, 1000)
```

All methods degrade at low $N$. Bayesian methods should degrade more
gracefully than MLE. Prior sensitivity is most acute at low $N$ and will
be documented explicitly.

------------------------------------------------------------------------

### Factor 2: Class imbalance

Proportion of class 0 (immature), total $N$ fixed at 200. Imbalance here
refers to **sampling imbalance** only — the true $m_{50}$ is fixed
regardless of how many of each class are sampled.

``` r
imbalance_scenarios <- tibble::tribble(
  ~label,      ~prop_class0,
  "50|50",     0.50,
  "75|25",     0.75,
  "90|10",     0.90,
  "10|90",     0.10,
  "95|5",      0.95
)
```

This is the **central claim** of the paper. Logistic regression (MLE and
Bayesian, unweighted) should show systematic bias as imbalance
increases. STAGE and LDA should be largely unaffected because they use
equal class priors and do not condition on the observed class
frequencies. The weighted logistic regression variants should partially
correct for imbalance.

------------------------------------------------------------------------

### Factor 3: Truncation choice

$L$ and $U$ are analyst choices. The true population range is \[200,
2000\]; the true transition zone is approximately
$\left\lbrack m_{50} - 2\sigma,m_{50} + 2\sigma \right\rbrack \approx \lbrack 1150,1350\rbrack$.

``` r
truncation_scenarios <- tibble::tribble(
  ~label,        ~L,    ~U,    ~notes,
  "Full range",  200,   2000,  "Baseline; wide plateaus",
  "Moderate",    500,   1700,  "Well outside transition zone",
  "Tight",       900,   1600,  "Narrower plateaus; still outside",
  "Marginal",    1100,  1400,  "Very narrow plateaus; near mu0/mu1",
  "Stress test", 1150,  1350,  "Inside transition zone; expected failure",
  # Asymmetric — tests L/U invariance claim
  "Tight lower", 1100,  2000,  "L near mu0, U wide",
  "Tight upper", 200,   1400,  "L wide, U near mu1"
)
```

STAGE should be robust to truncation as long as $L < \mu_{0}$ and
$U > \mu_{1}$ (some plateau data present). When bounds encroach on the
Gaussian tails, posterior uncertainty should increase — this is honest
degradation, not bias. Verify using **coverage**, not just interval
width.

Logistic regression has no dependence on $L$ and $U$ and will be
insensitive to truncation choice, but this is not a virtue: it means the
model ignores information about the data-generating bounds entirely.

The asymmetric scenarios are particularly important: they test the claim
that $m_{50}$ is invariant to the placement of $L$ and $U$ relative to
the transition zone.

------------------------------------------------------------------------

### Factor 4: Spatial/cohort clustering

Using the two-level sampling model, varying the degree of within-sample
clustering. Total $N = 200$, 20 events of 10 fish each.

``` r
clustering_scenarios <- tibble::tribble(
  ~label,       ~event_mu_sd, ~length_sd, ~approx_icc,
  "iid",        10,           150,        "~0",
  "Low",        50,           150,        "~0.1",
  "Moderate",   100,          150,        "~0.3",
  "High",       200,          150,        "~0.6"
)
```

None of the methods explicitly model clustering — this is a test of
robustness to a violated iid assumption, which is realistic for field
data. Clustering inflates effective variance and reduces effective
sample size for all methods. The key question is whether clustering
introduces **bias** in any method, not just increased uncertainty.

------------------------------------------------------------------------

### Factor 5: Transition sharpness

Varying the slope of the logistic truth:

``` r
sharpness_scenarios <- tibble::tribble(
  ~label,         ~slope,
  "Very sharp",   0.15,
  "Moderate",     0.08,
  "Diffuse",      0.03,
  "Very diffuse", 0.01
)
```

All methods should perform well on sharp transitions. On diffuse
transitions, separation between classes is low and all methods
accumulate uncertainty. The question is which methods remain unbiased as
transitions become diffuse.

------------------------------------------------------------------------

### Factor 6: Correctly specified model

Simulate from the STAGE generative model (plateau-Gaussian
class-conditionals) rather than a logistic truth. Fit STAGE (correctly
specified) and logistic regression (misspecified):

``` r
# Simulate x from plateau-Gaussian class-conditionals
simulate_stage_truth <- function(n0, n1, m50, d, sigma, L, U) {
  mu0 <- m50 - d / 2
  mu1 <- m50 + d / 2
  # Draw from truncated plateau-Gaussian distributions
  # ... (implementation using rejection sampling or inverse CDF)
}
```

This complements factor 5 by showing performance under correct
specification. It should not be over-emphasised — the logistic truth
scenarios are the primary test.

------------------------------------------------------------------------

### Factor 7: Multi-population hierarchical model

$J$ populations with varying $m_{50}$, using
`fit_stage(..., group = ...)`:

``` r
hierarchical_scenarios <- tibble::tribble(
  ~label,             ~J,  ~sigma_alpha, ~n_per_pop,
  "Few, similar",     5,   20,           100,
  "Few, divergent",   5,   100,          100,
  "Many, similar",    20,  20,           100,
  "Many, unbalanced", 20,  50,           "Poisson(50)"
)
```

Compare hierarchical STAGE against fitting STAGE independently to each
population. The hierarchical model should demonstrate partial pooling:
small populations borrow strength from larger ones.

------------------------------------------------------------------------

## Performance metrics

Across $R = 500$ simulated datasets per scenario:

``` r
compute_metrics <- function(m50_hat, m50_true, lower, upper) {
  tibble::tibble(
    bias     = mean(m50_hat) - m50_true,
    variance = var(m50_hat),
    rmse     = sqrt(mean((m50_hat - m50_true)^2)),
    coverage = mean(lower <= m50_true & m50_true <= upper),
    width    = mean(upper - lower)
  )
}
```

Bias and variance are reported separately, not just RMSE — a method can
have good RMSE by trading bias against variance, and these are different
kinds of failure. For Bayesian methods, `m50_hat` is the posterior
median.

Coverage at the nominal 95% level is the key diagnostic for uncertainty
quantification. With $R = 500$ replicates, the SE on a coverage estimate
of 0.95 is approximately 0.010, giving reasonable precision.

------------------------------------------------------------------------

## Priority

For the methods paper, the primary factors are:

1.  **Class imbalance** (factor 2) — central claim, must be thorough
2.  **Sample size** (factor 1) — necessary context
3.  **Truncation** (factor 3) — robustness claim with honest
    characterisation of boundary conditions
4.  **Clustering** (factor 4) — biological realism, novel test

Factors 5–7 are supplementary.

------------------------------------------------------------------------

## Implementation notes

``` r
# Fix global seed; derive per-replicate seeds from it
set.seed(42)
sim_seeds <- sample.int(1e6, size = 500)

# Precompile Stan models before simulation loop
# (compilation is slow; sampling is fast once compiled)
model_stage    <- stage:::get_cmdstan_model("stage_Jpop.stan")
model_lda      <- stage:::get_cmdstan_model("lda.stan")
model_logistic <- stage:::get_cmdstan_model("logistic.stan")

# Parallelise across replicates
future::plan(future::multisession)
results <- furrr::future_map(sim_seeds, run_one_replicate, .options =
  furrr::furrr_options(seed = TRUE))
```

- Store thinned posterior draws, not just summaries, to allow
  retrospective analysis
- Prior choices are fixed across all scenarios and documented explicitly
  — do not tune priors per scenario
- For truncation scenarios, $L$ and $U$ passed to STAGE are set to the
  truncation bounds used for sampling, not the true population bounds
- Stress-test truncation scenarios are clearly flagged as boundary
  conditions

------------------------------------------------------------------------

## Open questions

1.  Should any scenario use an asymmetric population (true $m_{50}$
    closer to $L_{pop}$ than to $U_{pop}$)? This would provide a more
    stringent test of $L$/$U$ invariance under realistic conditions.

2.  Should the simulation code be included as a separate vignette once
    complete, or integrated into this one?

3.  Should pre-computed results be stored in `inst/extdata/` for display
    on the pkgdown site without requiring re-execution during build?

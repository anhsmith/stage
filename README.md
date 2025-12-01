
# **📦 stage** <img src="man/figures/logo.png" align="right" width="120"/>

The **stage** package implements the **STAGE** model  
(**S**tate **T**ransition **A**nalysis via **G**enerative
**E**stimator),  
a Bayesian generative classifier for estimating **transition points**  
(e.g., *length-at-maturity*) from two groups where the underlying  
distributions consist of:

- a **uniform plateau**,  
- a **Gaussian tail**, and  
- a **shared transition region** centred at **m50** with width **d**.

The goal is **robust transition-point estimation** even with  
highly unbalanced sample sizes and data far from the transition zone.

------------------------------------------------------------------------

# **Introduction**

## Problems with existing methods

Estimating a transition point—such as the length at which individuals
shift from immature to mature—is a simple and common task, but it is
surprisingly difficult. Existing methods each fail in important ways.

### Ad hoc midpoint methods

A common practice is to avoid fitting a model altogether and instead
compute a transition point directly from the data.  
A widespread example is the midpoint between the smallest mature and
largest immature individuals. The advantage of this approach is that it
focuses on the transition zone only. However, this approach (1) has no
probabilistic foundation and so cannot quantify uncertainty, and (2)
uses only two data points, so has very high variance.

### Discriminative models (e.g., logistic regression)

Discriminative models such as logistic regression estimate the
conditional probability $P(y = 1 \mid x)$ *directly*. In other words,
the model focuses on predicting class labels based on the prevalence of
one class vs another at each value of $x$.

While simple and familiar, they suffer from two major issues in this
setting:

1.  **Sensitivity to class imbalance.**  
    Logistic regression (and other discriminative models) estimates the
    probability of maturity in the sample, not the population. Unless
    sampling is independent of both class and size, the estimated curve
    can be seriously biased. If one class has a much greater sample size
    than the other (i.e., they are unbalanced), logistic regression
    effectively places more weight on the more abundant class. This
    pulls the estimated transition curve toward the majority class, even
    though class counts should not determine where the transition
    occurs.

2.  **Influence from biologically irrelevant regions.**  
    Observations far from the transition zone change the fitted slope
    and intercept—even though they contain **no information** about
    where the transition takes place. For example, if the true
    transition is around 1000 mm, individuals at 200–400 mm should not
    influence the estimate of the maturity transition, but in a
    discriminative model they often do.

An advantage of the discriminative-model approach is that inferences
about the transition point are based mainly on data near the transition
zone.

### Gaussian generative models (e.g., Linear Discriminant Analysis)

Generative classifiers model the class-conditional densities

$f(x \mid y = 0), \qquad f(x \mid y = 1)$,

and then obtain $P(y=1 \mid x)$ *via* Bayes’ rule.

The generative classifier approach has the advantage of being less
sensitive to class imbalance, since the transition point estimate
depends on the shapes of the class distributions, not their relative
sample sizes. This is achieved by fitting a density to each class
separately, and then calculating the transition point based on the
relative densities of the two classes. Linear Discriminant Analysis
(LDA) is a widely used generative classifier that fits a Gaussian
(Normal) density to each class.

However, while generative models are less sensitive to class imbalance,
they suffer from another important limitation: *estimation of the
transition point is influenced by all values of $x$, not just those near
the transition zone*. In LDA, the Gaussian tail behaviour near the
transition is determined by the estimated mean and variance of each
class—and those estimates are strongly affected by observations far from
the transition. In the context of length-at-maturity, this means that
very small immature individuals (e.g., neonates) and very large mature
individuals have as much influence on the transition point as
individuals whose lengths fall near the maturity boundary. This is
unappealing, because these distant observations contain almost no
information about the location of the transition, yet they shape the
Gaussian fits and can shift the inferred transition point substantially.
This creates a further structural problem: the two tails that actually
matter for the transition—the upper tail of the immature class and the
lower tail of the mature class—are not modelled directly against each
other. Instead, each tail is constrained to be the mirror image of the
opposite tail of its own class (because of the global Gaussian
assumption).

------------------------------------------------------------------------

## Summary

What is needed is a model for estimating transition points that:

- is probabilistic and can quantify uncertainty,
- robust to class imbalance, and
- focuses inference on data in the overlap region rather than distant
  data.

------------------------------------------------------------------------

# **The STAGE solution**

The STAGE model is a Bayesian generative classifier where each class is
fit with an **asymmetric Uniform–Gaussian mixture**.

- **Immature (y = 0):**  
  uniform plateau → Gaussian decay approaching the transition.
- **Mature (y = 1):**  
  Gaussian rise out of the transition → uniform plateau.

This design has several key consequences:

### Generative modelling

Transition probability arises from comparing **relative densities**  
(Bayes rule), not from discriminative regression.

### Irrelevant regions contribute minimally

Observations far from the transition zone enter via the **uniform**
component,  
so they cannot distort estimates of the transition point.

### Robust to imbalance

Because inference depends on the **shape** of the densities—not the
counts—  
unbalanced datasets do not bias m50.

### Bayesian inference

Posterior uncertainty is natural, and multi-population structure is  
handled with a hierarchical model.

The STAGE likelihood is a structured mixture inspired by the
plateau–Gaussian  
formulation of Lau & Krumscheid (2022), adapted for transition-point
estimation.

------------------------------------------------------------------------

# **Package features**

- `fit_stage()` — fit STAGE model (single or multi-population)  
- `transition_point()` — summarise posteriors for ***m50***  
- `predict.stage_fit()` — predict class probabilities or labels  
- `stage_priors()` — prior helper  
- `get_cmdstan_model()` — internal CmdStan compiler helper

Planned additions:

- Reparameterised logistic/HOF model  
- Bayesian LDA  
- Posterior predictive classification  
- Diagnostic + visualisation tools  
- Full vignette

------------------------------------------------------------------------

# **Installation**

You need **cmdstanr** + CmdStan installed:

``` r
install.packages("cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
)
cmdstanr::install_cmdstan()   # run once
```

Then install from GitHub:

``` r
devtools::install_github("anhsmith/stage")
```

------------------------------------------------------------------------

# **Example**

``` r
library(stage)

set.seed(123)

N <- 200
L <- 1000
U <- 1500
true_m50 <- 1250
true_d   <- 100

# simulate data
x <- runif(N, L, U)
p <- plogis((x - true_m50) / (true_d / 4))
y <- rbinom(N, 1, p)

# fit STAGE model
fit <- fit_stage(x, y, L = L, U = U, chains = 2, iter = 1000)
fit
```

Estimated transition point:

``` r
transition_point(fit)
```

Prediction:

``` r
x_grid <- seq(L, U, length.out = 100)
p_hat <- predict(fit, x_grid, type = "prob")

plot(x_grid, p_hat, type = "l",
     xlab = "x", ylab = "P(mature)", las = 1)
```

------------------------------------------------------------------------

# **Hierarchical model**

Multiple populations:

``` r
group <- rep(c("A","B"), each = N/2)

fit2 <- fit_stage(x, y, group = group, L = L, U = U)

transition_point(fit2)
```

This returns:

- global m50  
- group-specific m50_pop\[j\]

------------------------------------------------------------------------

# **Model summary**

Each class is a **uniform plateau** with a **Gaussian edge**.  
Key parameters:

- `m50` — transition centre  
- `d` — width of transition region  
- `sigma_x` — Gaussian spread

Hierarchical structure:

``` md
m50_pop[j] = mu_m50 + z[j] * sigma_alpha
```

(non-centred parameterisation for efficient sampling)

------------------------------------------------------------------------

# **Development**

``` r
devtools::test()
devtools::document()
devtools::build_readme()
```

------------------------------------------------------------------------

# **License**

MIT © Adam N. H. Smith

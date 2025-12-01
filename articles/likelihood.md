# The STAGE Likelihood

## Overview

This vignette describes the **likelihood and parameterisation** used in
the STAGE model.

STAGE is a **hierarchical Bayesian generative model** designed to
estimate a **transition point** $m_{50}$ on variable $x$ between two
classes—for example, the length at which individuals transition from
immature to mature. The model produces a posterior distribution for key
parameter $m_{50}$, representing the value of $x$ at which the
probabilities $\Pr(y = 0)$ and $\Pr(y = 1)$ are equal. The model can
also be used to derive other useful quantities (such as class
probabilities $\Pr(y = 1 \mid x)$).

Conceptually, STAGE is a **Bayesian generative classifier** for
estimating a transition point (e.g. length at maturity). Each class is
fit with an asymmetric **plateau–Gaussian** density:

- the lower class has a **uniform plateau** followed by a **Gaussian
  tail down into the transition**,  
- the higher class has a **Gaussian tail up out of the transition**
  followed by a **uniform plateau**.

![](likelihood_files/figure-html/visualize-1.png)

This construction focuses inference on the **overlap region** and
prevents distant observations (e.g. very small immatures, very large
matures) from unduly influencing the transition point.

The central idea is to model each class with a **piecewise density**
that combines a **uniform plateau** in regions far from the transition,
and a **Gaussian tail** in the transition region. The left group
($y = 0$) is modelled with a uniform on the left that begins at the
lower truncation value $L$. The distribution remains flat until it
approaches the transition range, at which point the distribution
declines according to a Gaussian. The right group ($y = 1$) is a mirror
image of the left—it increases according to a Gaussian distribution and
then levels off to a uniform, terminating at the upper truncation value
$U$.

The model focuses inference on the **transition zone**, where the values
of $x$ of the two classes overlap, and down-weights observations that
are biologically uninformative because they are far away from the
transition point.

The uniform–Gaussian pdf used here is adapted from a plateau proposal
distribution combining a uniform and two Gaussians described by Lau &
Krumscheid (2022).

## Model summary

We observe pairs $\left( x_{i},y_{i} \right)$ for $i = 1,\ldots,N$,
where

- $x_{i}$ is a continuous covariate (e.g., length), truncated to an
  interval $\lbrack L,U\rbrack$, and  
- $y_{i} \in \{ 0,1\}$ is a binary indicator of state (e.g., immature /
  mature).

The model assumes that the distribution of $x$ conditional on state is
**piecewise**:

- For $y = 0$ (immature individuals), the distribution of $x$ is
  approximately **uniform** up to a point $\mu_{0}$, after which it
  follows a **Gaussian decay**.
- For $y = 1$ (mature individuals), the distribution is **Gaussian**
  below a point $\mu_{1}$, and **uniform** above $\mu_{1}$ up to the
  upper truncation $U$.

The **transition point** $m_{50}$ is defined as the value of $x$ where
the probability of observing either class is equal:

$$\Pr\left( y = 0 \mid x = m_{50} \right) = \Pr\left( y = 1 \mid x = m_{50} \right).$$

Additional parameters govern the **width** and **sharpness** of the
transition:

- $d$ controls the distance between $\mu_{0}$ and $\mu_{1}$ (the width
  of the transition region),
- $\sigma$ (denoted `sigma_x` in the Stan code) controls the spread
  (standard deviation) of the Gaussian tails.

By allowing different behaviours on either side of $m_{50}$, the model
can capture **asymmetric transition patterns**: the immature class can
decline towards the transition at a different rate than the mature class
rises out of it.

This is particularly useful in settings where:

- the boundary between classes is **not sharp**, and  
- the data exhibit **gradual overlap** between the two classes.

## Likelihood

The likelihood is built from two **unnormalised log-density functions**
for the two classes. We work on the log scale for numerical stability
and then subtract log-normalising constants to obtain valid densities.

Let $\ell_{0}(x)$ and $\ell_{1}(x)$ denote the unnormalised
log-densities for $y = 0$ and $y = 1$, respectively. The corresponding
(unnormalised) densities are
${\widetilde{f}}_{k}(x) = \exp\{\ell_{k}(x)\}$.

We restrict attention to the truncated interval $\lbrack L,U\rbrack$;
the densities are zero outside this range.

### Lower class ($y = 0$)

For $y_{i} = 0$, the unnormalised piecewise log-density is

$$\ell_{0}\left( x_{i} \mid \mu_{0},\sigma,L \right) = \begin{cases}
0 & {{\text{if}\mspace{6mu}}L \leq x_{i} \leq \mu_{0},} \\
{- \frac{\left( x_{i} - \mu_{0} \right)^{2}}{2\sigma^{2}}} & {{\text{if}\mspace{6mu}}x_{i} > \mu_{0}.}
\end{cases}$$

Here:

- the plateau region $\left\lbrack L,\mu_{0} \right\rbrack$ has constant
  log-density 0 (density 1),
- the tail region $x_{i} > \mu_{0}$ decays according to a Gaussian
  kernel.

Any constant could be added to $\ell_{0}( \cdot )$ without changing the
posterior; we fix the plateau at 0 so that the plateau density is 1
after exponentiation.

### Higher class ($y = 1$)

For $y_{i} = 1$, the unnormalised log-density is

$$\ell_{1}\left( x_{i} \mid \mu_{1},\sigma,U \right) = \begin{cases}
{- \frac{\left( x_{i} - \mu_{1} \right)^{2}}{2\sigma^{2}}} & {{\text{if}\mspace{6mu}}x_{i} < \mu_{1},} \\
0 & {{\text{if}\mspace{6mu}}\mu_{1} \leq x_{i} \leq U.}
\end{cases}$$

Now:

- the Gaussian tail is **below** $\mu_{1}$, and  
- the plateau covers the upper region
  $\left\lbrack \mu_{1},U \right\rbrack$.

In both cases we have a **plateau–Gaussian** shape: a flat density where
the class “dominates” and a Gaussian tail in the transition region where
the other class is present.

### Normalisation constants

To convert ${\widetilde{f}}_{k}(x) = \exp\{\ell_{k}(x)\}$ into proper
probability density functions on $\lbrack L,U\rbrack$, we divide by
constants $C_{0}$ and $C_{1}$ chosen so that the densities integrate to
1.

For the immature class ($y = 0$):

$$C_{0} = \int_{L}^{\mu_{0}}1\, dx + \int_{\mu_{0}}^{\infty}\frac{1}{\sqrt{2\pi}\sigma}\exp\!\left( - \frac{\left( x - \mu_{0} \right)^{2}}{2\sigma^{2}} \right)dx = \left( \mu_{0} - L \right) + \frac{\sqrt{2\pi}\,\sigma}{2}.$$

For the mature class ($y = 1$):

$$C_{1} = \int_{- \infty}^{\mu_{1}}\frac{1}{\sqrt{2\pi}\sigma}\exp\!\left( - \frac{\left( x - \mu_{1} \right)^{2}}{2\sigma^{2}} \right)dx + \int_{\mu_{1}}^{U}1\, dx = \frac{\sqrt{2\pi}\,\sigma}{2} + \left( U - \mu_{1} \right).$$

Intuitively:

- $\left( \mu_{0} - L \right)$ and $\left( U - \mu_{1} \right)$ are the
  **widths of the uniform plateaus**,  
- $\sqrt{2\pi}\,\sigma/2$ is the **contribution from the half-Gaussian
  tail**.

The corresponding **log-likelihood contribution** for observation $i$ is

$$\log\mathcal{L}_{i} = \begin{cases}
{\ell_{0}\left( x_{i} \right) - \log C_{0}} & {{\text{if}\mspace{6mu}}y_{i} = 0,} \\
{\ell_{1}\left( x_{i} \right) - \log C_{1}} & {{\text{if}\mspace{6mu}}y_{i} = 1.}
\end{cases}$$

Summing over all observations,

$$\log p\left( x,y \mid m_{50},d,\sigma \right) = \sum\limits_{i = 1}^{N}\lbrack\ell_{y_{i}}\left( x_{i} \right) - \log C_{y_{i}}\rbrack,$$

which is exactly the likelihood implemented in the Stan code used by
[`fit_stage()`](https://anhsmith.github.io/stage/reference/fit_stage.md).

## Prior specification

The priors reflect prior beliefs about plausible transition points and
transition widths in typical length-at-maturity problems. In practice,
these are **defaults** that should be adapted to the scale of the data
and prior biological knowledge.

By default, STAGE uses weakly informative normal priors:

- **Transition point** $m_{50}$:

$$m_{50} \sim \mathcal{N}(1300,50),$$

which centres the transition around 1300 (mm) with moderate uncertainty.

- **Transition width** $d$:

$$d \sim \mathcal{N}(100,100),$$

which allows the width of the overlap region to vary around 100 units
with a wide spread. Note that $d$ is the distance between the
class-specific Gaussian centres.

- **Gaussian spread** $\sigma$ (Stan parameter `sigma_x`):

$$\sigma \sim \mathcal{N}(100,50),$$

representing prior beliefs about the scale of variability in the
transition tails.

In practice, you should adapt these priors to your scale (e.g., lengths,
ages) and prior biological knowledge. The STAGE implementation allows
you to pass in hyperparameters via helper functions such as
[`stage_priors()`](https://anhsmith.github.io/stage/reference/stage_priors.md).

## Transformed parameters

The Gaussian “centres” for the two classes are defined in terms of the
transition point $m_{50}$ and the width parameter $d$:

$$\mu_{0} = m_{50} - \frac{d}{2},\qquad\mu_{1} = m_{50} + \frac{d}{2}.$$

Thus:

- $d$ controls the **distance** between the immature and mature Gaussian
  means,  
- $m_{50}$ remains the **midpoint** between them.

This is convenient both conceptually and computationally: it allows the
model to focus directly on $m_{50}$ and $d$, which are often the primary
scientific targets, while $\mu_{0}$ and $\mu_{1}$ are derived
quantities.

## From class densities to class probabilities (Bayes rule)

The likelihood above defines the **class-conditional densities**  
$f_{0}(x) = p(x \mid y = 0)$ and $f_{1}(x) = p(x \mid y = 1)$ on
$\lbrack L,U\rbrack$, obtained by normalising the unnormalised densities
${\widetilde{f}}_{k}(x) = \exp\{\ell_{k}(x)\}$ with the constants
$C_{k}$.

To obtain the probability that an individual with a particular $x$ value
is in the higher class ($y = 1$), we use Bayes rule:

$$\Pr(y = 1 \mid x) = \frac{\pi_{1}f_{1}(x)}{\pi_{0}f_{0}(x) + \pi_{1}f_{1}(x)},$$

and similarly

$$\Pr(y = 0 \mid x) = \frac{\pi_{0}f_{0}(x)}{\pi_{0}f_{0}(x) + \pi_{1}f_{1}(x)},$$

where $\pi_{0}$ and $\pi_{1}$ are the prior class probabilities  
(e.g., prior probabilities of being immature or mature *before* seeing
$x$).

In many applications, such as estimating length-at-maturity, we are
interested in a transition point that is **not driven by sample sizes or
prevalence**, so it is natural to use **equal class priors**,
$\pi_{0} = \pi_{1} = 1/2$. In that case, Bayes rule simplifies to

$$\Pr(y = 1 \mid x) = \frac{f_{1}(x)}{f_{0}(x) + f_{1}(x)},\qquad\Pr(y = 0 \mid x) = \frac{f_{0}(x)}{f_{0}(x) + f_{1}(x)}.$$

These are the class probabilities that the STAGE model uses for
prediction, and they are what you see when you plot $\Pr(y = 1 \mid x)$
as a function of $x$.

### Connection to $m_{50}$

By definition, the transition point $m_{50}$ is the value of $x$ where
the two classes are equally probable:

$$\Pr\left( y = 0 \mid x = m_{50} \right) = \Pr\left( y = 1 \mid x = m_{50} \right) = \frac{1}{2}.$$

Using Bayes rule with equal class priors, this equality is equivalent to

$$\left. \frac{f_{1}\left( m_{50} \right)}{f_{0}\left( m_{50} \right) + f_{1}\left( m_{50} \right)} = \frac{f_{0}\left( m_{50} \right)}{f_{0}\left( m_{50} \right) + f_{1}\left( m_{50} \right)}\quad\Leftrightarrow\quad f_{0}\left( m_{50} \right) = f_{1}\left( m_{50} \right). \right.$$

Thus, $m_{50}$ can be viewed in two equivalent ways:

1.  As a **model parameter** that sets the midpoint between the Gaussian
    centres $\mu_{0}$ and $\mu_{1}$ via  
    $$\mu_{0} = m_{50} - \frac{d}{2},\qquad\mu_{1} = m_{50} + \frac{d}{2},$$
2.  As the **Bayes-optimal classifier cut-point** where the posterior
    class probabilities for immature and mature are both 0.5.

In practice, the STAGE model samples $m_{50}$, $d$, and $\sigma$ from
their posterior distribution. For any posterior draw, we can compute
$\Pr(y = 1 \mid x)$ via Bayes rule and verify that the point where this
curve crosses 0.5 aligns with the sampled value of $m_{50}$. This makes
the interpretation of $m_{50}$ as a “50% maturity” transition point
directly linked to the underlying generative model and Bayes-optimal
classifier.

### Visualising a single STAGE transition

We can visualise a simple STAGE configuration, showing:

- the immature and mature densities, and
- the resulting transition probability $\Pr(y = 1 \mid x)$.

![](likelihood_files/figure-html/visualize-likelihood-1.png)

## Relationship to the `stage` implementation

In the `stage` package:

- The single-population model uses parameters `m50`, `d`, and `sigma_x`
  that correspond to $m_{50}$, $d$, and $\sigma$ in this vignette.
- The multi-population (hierarchical) model introduces
  population-specific transition points $m_{50,j}$ via a random-effects
  structure:
  $$m50_{\text{pop}{\lbrack j\rbrack}} = \mu_{m50} + z_{j}\,\sigma_{\alpha},$$
  where $\mu_{m50}$ is the overall mean transition point, $z_{j}$ are
  standard-normal effects, and $\sigma_{\alpha}$ is the population-level
  standard deviation.

The Stan code used in
[`fit_stage()`](https://anhsmith.github.io/stage/reference/fit_stage.md)
is a direct implementation of the likelihood and parameterisation
described above.

Future vignettes will provide a more detailed walk-through of the Stan
code and extensions such as alternative priors and model comparison.

## References

Lau, F. Din-Houn, and Sebastian Krumscheid. 2022. “Plateau Proposal
Distributions for Adaptive Component-Wise Multiple-Try Metropolis.”
*METRON* 80 (3): 343–70. <https://doi.org/10.1007/s40300-022-00235-y>.

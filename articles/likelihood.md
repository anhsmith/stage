# The STAGE Likelihood

## Overview

This vignette describes the **likelihood and parameterisation** used in
the STAGE model.

STAGE is a **hierarchical Bayesian generative model** designed to
estimate a **transition point** $m_{50}$ on variable $x$ between two
classes—for example, the length at which individuals transition from
immature to mature. The model produces a posterior distribution for key
parameter $m_{50}$, representing the value of $x$ at which the
probabilities $Pr(y = 0)$ and $\Pr(y = 1)$ are equal. The model can also
be used to derive other useful quantities (such as class probabilities
$\Pr(y = 1 \mid x)$).

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
of $x$ of the two classes overlap, thereby limiting the influence of
observations that are biologically uninformative because they are far
away from the transition point.

The PDF of the uniform-Gaussian used here was modified from the PDF
uniform and two Guassians described by Lau & Krumscheid (2022).

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
- $\sigma$ controls the spread (standard deviation) of the Gaussian
  tails.

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
and then subtract log-normalising constants to ensure valid densities.

Let $f_{0}(x)$ and $f_{1}(x)$ denote the unnormalised densities for the
two classes.

### Lower class ($y = 0$)

For $y_{i} = 0$, the unnormalised piecewise log-density is:

$$f\left( x_{i} \mid y_{i} = 0,\mu_{0},\sigma,L \right) = \begin{cases}
0 & {{\text{if}\mspace{6mu}}x_{i} < L} \\
1 & {{\text{if}\mspace{6mu}}L \leq x_{i} \leq \mu_{0}} \\
{- \frac{\left( x_{i} - \mu_{0} \right)^{2}}{2\sigma^{2}}} & {{\text{if}\mspace{6mu}}x_{i} > \mu_{0}}
\end{cases}$$

Here:

- the plateau region $\left\lbrack L,\mu_{0} \right\rbrack$ has constant
  (log-)density 1,  
- the tail region $x_{i} > \mu_{0}$ decays according to a Gaussian
  kernel.

In implementation, this is treated as an **unnormalised log-density**;
the overall density is recovered after subtracting a normalising
constant.

### Higher class ($y = 1$)

For $y_{i} = 1$, the unnormalised piecewise log-density is:

$$f\left( x_{i} \mid y_{i} = 1,\mu_{1},\sigma,U \right) = \begin{cases}
{- \frac{\left( x_{i} - \mu_{1} \right)^{2}}{2\sigma^{2}}} & {{\text{if}\mspace{6mu}}x_{i} < \mu_{1}} \\
1 & {{\text{if}\mspace{6mu}}\mu_{1} \leq x_{i} \leq U} \\
0 & {{\text{if}\mspace{6mu}}x_{i} > U}
\end{cases}$$

Now:

- the Gaussian tail is **below** $\mu_{1}$, and  
- the plateau covers the upper region
  $\left\lbrack \mu_{1},U \right\rbrack$.

In both cases, we are effectively using a **plateau–Gaussian** mixture:
a flat density where the class “dominates” and a Gaussian shape in the
transition region where the other class is present.

### Normalisation constants

To convert the unnormalised densities into proper probability density
functions (pdf), we divide by constants $C_{0}$ and $C_{1}$ chosen such
that the densities integrate to 1 over $\lbrack L,U\rbrack$.

For the immature class, the normalising constant is:

$$C_{0} = \frac{\sqrt{2\pi}\,\sigma}{2} + \left( \mu_{0} - L \right).$$

For the mature class, the normalising constant is:

$$C_{1} = \frac{\sqrt{2\pi}\,\sigma}{2} + \left( U - \mu_{1} \right).$$

Intuitively:

- $\left( \mu_{0} - L \right)$ and $\left( U - \mu_{1} \right)$ are the
  **widths of the uniform plateaus**,  
- $\sqrt{2\pi}\,\sigma/2$ is the **contribution from the half-Gaussian
  tail**.

The **log-likelihood contributions** for each observation are then:

$$\log\mathcal{L}_{i}\; = \;\begin{cases}
{\log f\left( x_{i} \right) - \log C_{0}} & {{\text{if}\mspace{6mu}}y_{i} = 0} \\
{\log f\left( x_{i} \right) - \log C_{1}} & {{\text{if}\mspace{6mu}}y_{i} = 1.}
\end{cases}$$

### Full log-likelihood

Summing over all observations, the joint log-likelihood is:

$$\log p\left( x,y \mid m_{50},d,\sigma \right)\; = \;\sum\limits_{i = 1}^{N}\left\lbrack \log f\left( x_{i} \mid y_{i},\mu_{y_{i}},\sigma \right) - \log C_{y_{i}} \right\rbrack,$$

where

- $f( \cdot )$ is the appropriate unnormalised density for class
  $y_{i}$, and  
- $C_{y_{i}} \in \{ C_{0},C_{1}\}$ is the normalising constant for that
  class.

## Prior specification

The priors reflect prior beliefs about plausible transition points and
transition widths in typical length-at-maturity problems (but can be
changed by the user).

- **Transition point** $m_{50}$:

$$m_{50} \sim \mathcal{N}(1300,50),$$

which centres the transition around 1300 (mm) with moderate uncertainty.

- **Transition width** $d$:

$$d \sim \mathcal{N}(100,100),$$

which allows the width of the overlap region to vary around 100 units
with a wide spread.

- **Gaussian spread** $\sigma$:

$$\sigma \sim \mathcal{N}(100,50),$$

representing prior beliefs about the scale of variability in the
transition tails.

In practice, you should adapt these priors to your scale (e.g., lengths,
ages) and prior biological knowledge. The STAGE implementation allows
you to pass in hyperparameters via helper functions
(e.g. [`stage_priors()`](https://anhsmith.github.io/stage/reference/stage_priors.md)).

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
$\lbrack L,U\rbrack$.

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

#### Connection to $m_{50}$

By definition, the transition point $m_{50}$ is the value of $x$ where
the two classes are equally probable:

$$\Pr\left( y = 0 \mid x = m_{50} \right) = \Pr\left( y = 1 \mid x = m_{50} \right) = \frac{1}{2}.$$

Using Bayes rule with equal class priors, this equality is equivalent to

$$\left. \frac{f_{1}\left( m_{50} \right)}{f_{0}\left( m_{50} \right) + f_{1}\left( m_{50} \right)} = \frac{f_{0}\left( m_{50} \right)}{f_{0}\left( m_{50} \right) + f_{1}\left( m_{50} \right)}\quad\Leftrightarrow\quad f_{0}\left( m_{50} \right) = f_{1}\left( m_{50} \right). \right.$$

Thus, $m_{50}$ can be viewed in two equivalent ways:

1.  As a **model parameter** that sets the midpoint between the Gaussian
    centres $\mu_{0}$ and $\mu_{1}$ via  
    $$\mu_{0} = m_{50} - \frac{d}{2},\qquad\mu_{1} = m_{50} + \frac{d}{2},$$
2.  As the **Bayes classifier cut-point** where the posterior class
    probabilities for immature and mature are both 0.5.

In practice, the STAGE model samples $m_{50}$, $d$, and $\sigma$ from
their posterior distribution. For any posterior draw, we can compute
$\Pr(y = 1 \mid x)$ via Bayes rule and verify that the point where this
curve crosses 0.5 aligns with the sampled value of $m_{50}$. This makes
the interpretation of $m_{50}$ as a “50% maturity” transition point
directly linked to the underlying generative model and Bayes-optimal
classifier.

## Advantages of the STAGE likelihood

1.  **Piecewise densities**

    The plateau–Gaussian form naturally captures situations where
    individuals in one class (e.g., immatures) are approximately uniform
    over a wide region, but then decline gradually into the transition.
    Similarly for the mature class. This is more realistic than assuming
    a global Gaussian for each class.

2.  **Focus on the transition region**

    Because only the Gaussian tails enter the overlap region, the
    estimation of $m_{50}$ is driven by observations near the
    transition. Distant observations (e.g., very small immatures, very
    large matures) mainly contribute through the uniform plateau and
    therefore do not distort the transition point.

3.  **Flexibility via $d$**

    The parameter $d$ controls the width of the transition zone,
    allowing the model to capture anything from very sharp to very
    gradual transitions.

4.  **Bayesian inference**

    A fully Bayesian treatment:

    - incorporates prior knowledge (e.g., plausible ranges for
      $m_{50}$),  
    - provides full posterior uncertainty for all parameters,  
    - extends naturally to **hierarchical models** with multiple
      populations.

5.  **Proper normalisation**

    The explicit normalising constants $C_{0}$ and $C_{1}$ ensure that
    each class-conditional density is a valid pdf. This is crucial for
    stable inference in Stan and for interpreting posterior predictions
    as probabilities.

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

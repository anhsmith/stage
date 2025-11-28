functions {
  // Normalising constant for uniform+Gaussian on one side:
  // delta is the uniform segment length (e.g. mu - L or U - mu).
  real uniform_gaussian_normalising_constant(real delta, real sigma) {
    real sqrt_2pi = sqrt(2 * pi());
    return (sqrt_2pi * sigma / 2) + delta;
  }

  // For y=0 class: uniform from L up to mu0, Gaussian tail above mu0
  real uniform_left_gaussian_log_density_unnormalised(real y, real mu, real sigma, real L) {
    real log_density_unnormalised;
    if (y <= mu) {
      log_density_unnormalised = 0;
    } else {
      log_density_unnormalised = -0.5 * square((y - mu) / sigma);
    }
    return log_density_unnormalised;
  }

  // For y=1 class: Gaussian tail below mu1, uniform from mu1 up to U
  real uniform_right_gaussian_log_density_unnormalised(real y, real mu, real sigma, real U) {
    real log_density_unnormalised;
    if (y >= mu) {
      log_density_unnormalised = 0;
    } else {
      log_density_unnormalised = -0.5 * square((y - mu) / sigma);
    }
    return log_density_unnormalised;
  }
}

data {
  int<lower=1> N;                         // number of observations
  int<lower=1> J;                         // number of populations (groups)
  real L;                                 // lower truncation bound
  real U;                                 // upper truncation bound

  array[N] real<lower=L, upper=U> x;      // observed predictor (length)
  array[N] int<lower=0, upper=1> y;       // response: 0 = immature, 1 = mature
  array[N] int<lower=1, upper=J> pop;     // population index

  // prior hyperparameters (passed from R)
  real prior_mu_m50_mu;
  real<lower=0> prior_mu_m50_tau;
  real prior_d_mu;
  real<lower=0> prior_d_tau;
  real<lower=0> prior_sigma_x_tau;
  real<lower=0> prior_sigma_alpha_tau;
}

parameters {
  real<lower=L, upper=U> mu_m50;   // global mean transition point
  real<lower=0> d;                 // distance between mu0 and mu1

  vector[J] z;                     // non-centred group effects on m50

  real<lower=0> sigma_x;           // SD of x around group transition
  real<lower=0> sigma_alpha;       // SD of population effects on m50
}

transformed parameters {
  vector[J] alpha;     // population deviations
  vector[J] m50_pop;   // population-specific transition points
  vector[J] mu0_pop;   // lower (immature) Gaussian means
  vector[J] mu1_pop;   // upper (mature) Gaussian means
  vector[J] C0_pop;    // normalising constants for y=0 class
  vector[J] C1_pop;    // normalising constants for y=1 class

  for (j in 1:J) {
    alpha[j]   = z[j] * sigma_alpha;        // non-centred
    m50_pop[j] = mu_m50 + alpha[j];
    mu0_pop[j] = m50_pop[j] - d / 2;
    mu1_pop[j] = m50_pop[j] + d / 2;

    C0_pop[j] = uniform_gaussian_normalising_constant(mu0_pop[j] - L, sigma_x);
    C1_pop[j] = uniform_gaussian_normalising_constant(U - mu1_pop[j], sigma_x);
  }
}

model {
  // Priors
  mu_m50      ~ normal(prior_mu_m50_mu, prior_mu_m50_tau);
  d           ~ normal(prior_d_mu, prior_d_tau);
  sigma_x     ~ normal(0, prior_sigma_x_tau);
  sigma_alpha ~ normal(0, prior_sigma_alpha_tau);
  z           ~ normal(0, 1);

  // Likelihood
  for (i in 1:N) {
    int j = pop[i];
    if (y[i] == 0) {
      target +=
        uniform_left_gaussian_log_density_unnormalised(
          x[i], mu0_pop[j], sigma_x, L
        ) - log(C0_pop[j]);
    } else {
      target +=
        uniform_right_gaussian_log_density_unnormalised(
          x[i], mu1_pop[j], sigma_x, U
        ) - log(C1_pop[j]);
    }
  }
}

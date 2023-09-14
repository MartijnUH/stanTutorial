//
data {
  int n_schools;
  array[n_schools] real y;
  vector<lower = 0>[n_schools] sigma;
}

parameters {
  real mu;
  real<lower = 0> tau;
  // vector[n_schools] theta;
}

model {
  // Priors
  mu ~ normal(5, 3);
  tau ~ normal(0, 10);
  
  // Marginal Likelihood
  y ~ normal(mu, sqrt(square(tau) + square(sigma)));  // p(y| mu, tau)
  
}

generated quantities {
  array [n_schools] real theta;
  vector [n_schools] conjugate_variance = 1 / (1/square(sigma) + 1/square(tau));
  vector [n_schools] conjugate_mean = 
  to_vector(y) ./ square(sigma) + mu/square(tau) .* conjugate_variance;
  
  // p(theta| mu, tau, y)
  theta = normal_rng(conjugate_mean, sqrt(conjugate_variance));
}


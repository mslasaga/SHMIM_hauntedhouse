// Ordered-logit model. Latent mean carries slopes only; cutpoints carry the location.

data {
  int<lower=1> N;
  int<lower=1> N_subj;
  array[N] int<lower=1, upper=N_subj> subject_id;
  int<lower=2> K;
  array[N] int<lower=1, upper=K> y;

  int<lower=1> K_mu;
  matrix[N, K_mu] X_mu;

  int<lower=1> K_random_mu;
  matrix[N, K_random_mu] Z_mu;
}

parameters {
  vector[K_mu] beta_mu;

  matrix[N_subj, K_random_mu] u_raw_mu;
  vector<lower=0>[K_random_mu] sigma_u_mu;
  cholesky_factor_corr[K_random_mu] L_u_mu;

  ordered[K-1] cutpoints;
}

transformed parameters {
  matrix[N_subj, K_random_mu] u_mu = (diag_pre_multiply(sigma_u_mu, L_u_mu) * u_raw_mu')';

  vector[N] eta_mu = X_mu * beta_mu;
  for (n in 1:N) {
    eta_mu[n] += dot_product(Z_mu[n], u_mu[subject_id[n]]);
  }
}

model {
  beta_mu ~ normal(0, 1);

  sigma_u_mu[1] ~ normal(0, 1);
  for (k in 2:K_random_mu) sigma_u_mu[k] ~ normal(0, 0.5);
  L_u_mu ~ lkj_corr_cholesky(2);
  to_vector(u_raw_mu) ~ std_normal();

  cutpoints ~ normal(0, 2);

  for (n in 1:N) {
    target += ordered_logistic_lpmf(y[n] | eta_mu[n], cutpoints);
  }
}

generated quantities {
  vector[K_mu] beta_mu_gq = beta_mu;

  vector[K_random_mu] sigma_u_gq_mu = sigma_u_mu;
  matrix[K_random_mu, K_random_mu] Omega_mu = multiply_lower_tri_self_transpose(L_u_mu);
  matrix[N_subj, K_random_mu] u_gq_mu = u_mu;

  vector[K-1] cutpoint_intercepts_gq = cutpoints;

  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = ordered_logistic_lpmf(y[n] | eta_mu[n], cutpoints);
  }
}

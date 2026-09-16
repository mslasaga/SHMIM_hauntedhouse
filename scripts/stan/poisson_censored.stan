// Right-censored Poisson mixed model. Counts equal to cens are treated as >= cens.

data {
  int<lower=1> N;
  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N> subject_id;
  array[N_obs] int<lower=0> y;
  int<lower=1> cens;

  int<lower=1> K_count;
  matrix[N_obs, K_count] X_count;
  int<lower=1> count_intercept_pos;

  int<lower=1> K_random_count;
  matrix[N_obs, K_random_count] Z_count;
}

parameters {
  vector[K_count] beta_count;

  matrix[N, K_random_count] u_raw_count;
  vector<lower=0>[K_random_count] sigma_u_count;
  cholesky_factor_corr[K_random_count] L_u_count;
}

transformed parameters {
  matrix[N, K_random_count] u_count = (diag_pre_multiply(sigma_u_count, L_u_count) * u_raw_count')';

  vector[N_obs] eta_count = X_count * beta_count;

  for (n in 1:N_obs) {
    eta_count[n] += dot_product(Z_count[n], u_count[subject_id[n]]);
  }
}

model {
  beta_count[count_intercept_pos] ~ normal(0, 2);
  for (k in 1:K_count) {
    if (k != count_intercept_pos) beta_count[k] ~ normal(0, 1);
  }

  sigma_u_count[1] ~ normal(0, 1);
  for (k in 2:K_random_count) sigma_u_count[k] ~ normal(0, 0.5);
  L_u_count ~ lkj_corr_cholesky(2);
  to_vector(u_raw_count) ~ std_normal();

  for (n in 1:N_obs) {
    real lambda = exp(eta_count[n]);
    if (y[n] == cens) {
      target += poisson_lccdf(cens - 1 | lambda);
    } else {
      target += poisson_lpmf(y[n] | lambda);
    }
  }
}

generated quantities {
  vector[K_count] beta_count_gq = beta_count;

  vector[K_random_count] sigma_u_gq_count = sigma_u_count;
  matrix[K_random_count, K_random_count] Omega_count = multiply_lower_tri_self_transpose(L_u_count);
  matrix[N, K_random_count] u_gq_count = u_count;

  vector[N_obs] log_lik;

  for (n in 1:N_obs) {
    real lambda = exp(eta_count[n]);
    if (y[n] == cens) {
      log_lik[n] = poisson_lccdf(cens - 1 | lambda);
    } else {
      log_lik[n] = poisson_lpmf(y[n] | lambda);
    }
  }
}

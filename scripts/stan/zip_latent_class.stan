// Person-level latent-class zero-inflated censored Poisson.
// pri = inv_logit(eta_zi) = P(susceptible class); brms zi is 1 - pri.

data {
  int<lower=1> N;
  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N> subject_id;
  array[N_obs] int<lower=0> y;
  int<lower=1> cens;

  int<lower=1> K_count;
  matrix[N_obs, K_count] X_count;
  int<lower=1> count_intercept_pos;

  int<lower=1> K_zi;
  matrix[N, K_zi] X_zi;
  int<lower=1> zi_intercept_pos;

  int<lower=1> K_random;
  matrix[N_obs, K_random] Z;
  int<lower=0> lik_ind;
}

parameters {
  vector[K_count] beta_count;

  vector[K_zi] beta_zi;

  matrix[N, K_random] u_raw;
  vector<lower=0>[K_random] sigma_u;
  cholesky_factor_corr[K_random] L_u;
}

transformed parameters {
  matrix[N, K_random] u = (diag_pre_multiply(sigma_u, L_u) * u_raw')';

  vector[N_obs] eta_count = X_count * beta_count;
  vector[N] eta_zi = X_zi * beta_zi;
  vector<lower=0, upper=1>[N] pri = inv_logit(eta_zi);

  for (n in 1:N_obs) {
    eta_count[n] += dot_product(Z[n], u[subject_id[n]]);
  }
}

model {
  beta_count[count_intercept_pos] ~ normal(0, 2);
  for (k in 1:K_count) {
    if (k != count_intercept_pos) beta_count[k] ~ normal(0, 1);
  }

  beta_zi[zi_intercept_pos] ~ normal(0, 2);
  for (k in 1:K_zi) {
    if (k != zi_intercept_pos) beta_zi[k] ~ normal(0, 1);
  }

  sigma_u[1] ~ normal(0, 1);
  for (k in 2:K_random) sigma_u[k] ~ normal(0, 0.5);
  L_u ~ lkj_corr_cholesky(2);
  to_vector(u_raw) ~ std_normal();

  for (i in 1:N) {
    real log_lik_A = log(pri[i]);
    real log_lik_B = log1m(pri[i]);

    for (n in 1:N_obs) {
      if (subject_id[n] == i) {
        real lambda = exp(eta_count[n]);
        if (y[n] == cens) {
          log_lik_A += poisson_lccdf(cens - 1 | lambda);
        } else {
          log_lik_A += poisson_lpmf(y[n] | lambda);
        }
        if (lik_ind == 0) {
          log_lik_B += (y[n] == 0) ? 0 : negative_infinity();
        } else {
          log_lik_B += poisson_lpmf(0 | lambda);
        }
      }
    }
    target += log_sum_exp(log_lik_A, log_lik_B);
  }
}

generated quantities {
  vector[K_count] beta_count_gq = beta_count;
  vector[K_zi] beta_zi_gq = beta_zi;

  vector[K_random] sigma_u_gq = sigma_u;
  matrix[K_random, K_random] Omega = multiply_lower_tri_self_transpose(L_u);
  matrix[N, K_random] u_gq = u;

  vector[N] log_lik;

  for (i in 1:N) {
    real log_lik_A = log(pri[i]);
    real log_lik_B = log1m(pri[i]);

    for (n in 1:N_obs) {
      if (subject_id[n] == i) {
        real lambda = exp(eta_count[n]);
        if (y[n] == cens) {
          log_lik_A += poisson_lccdf(cens - 1 | lambda);
        } else {
          log_lik_A += poisson_lpmf(y[n] | lambda);
        }
        if (lik_ind == 0) {
          log_lik_B += (y[n] == 0) ? 0 : negative_infinity();
        } else {
          log_lik_B += poisson_lpmf(0 | lambda);
        }
      }
    }
    log_lik[i] = log_sum_exp(log_lik_A, log_lik_B);
  }
}

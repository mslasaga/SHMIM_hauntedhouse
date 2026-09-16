# 2.0 Censored Poisson. Stan: poisson_censored.stan.
# Out: models/fit_pois_cens.rds, models/fit_pois_hr_cens.rds

source('scripts/0. Setup/0.0 Libraries and Functions.R')
source('scripts/1. Processing/1.0 Data Processing.R')

mod_pois <- cmdstan_model('scripts/stan/poisson_censored.stan')

# ---- Stan data ----

sdata_pois <- within(list(), {
  N = length(unique(dat_im_all$id))
  N_obs = nrow(dat_im_all)
  subject_id = as.numeric(factor(dat_im_all$id))
  y = dat_im_all$im_frequency
  cens = 10
  X_count = model.matrix(~ 1 + day_after*scare_level + total_sh,
                         data = data.frame(day_after = dat_im_all$day_after,
                                           scare_level = dat_im_all$scare_level,
                                           total_sh = dat_im_all$total_sh))
  K_count = ncol(X_count)
  count_intercept_pos = 1
  Z_count = model.matrix(~ 1 + day_after, data = data.frame(day_after = dat_im_all$day_after))
  K_random_count = ncol(Z_count)
})

sdata_pois_hr <- within(list(), {
  N = length(unique(dat_im_hr$id))
  N_obs = nrow(dat_im_hr)
  subject_id = as.numeric(factor(dat_im_hr$id))
  y = dat_im_hr$im_frequency
  cens = 10
  X_count = model.matrix(~ 1 + day_after*scare_level + total_sh + mean_hr,
                         data = data.frame(day_after = dat_im_hr$day_after,
                                           scare_level = dat_im_hr$scare_level,
                                           total_sh = dat_im_hr$total_sh,
                                           mean_hr = dat_im_hr$mean_hr))
  K_count = ncol(X_count)
  count_intercept_pos = 1
  Z_count = model.matrix(~ 1 + day_after, data = data.frame(day_after = dat_im_hr$day_after))
  K_random_count = ncol(Z_count)
})

# ---- fit ----

fit_pois <- mod_pois$sample(data = sdata_pois,
             chains = 10,
             parallel_chains = 10,
             refresh = 50,
             seed = 2026,
             iter_warmup = 4000,
             iter_sampling = 6000)

fit_pois$save_object(file = 'models/fit_pois_cens.rds')

fit_pois_hr <- mod_pois$sample(data = sdata_pois_hr,
                  chains = 10,
                  parallel_chains = 10,
                  refresh = 50,
                  seed = 2026,
                  iter_warmup = 4000,
                  iter_sampling = 6000)

fit_pois_hr$save_object(file = 'models/fit_pois_hr_cens.rds')

# ---- summaries ----

extract_and_summarise_base(fit_pois, sdata_pois) %>%
  print(n = 100)

extract_and_summarise_base(fit_pois_hr, sdata_pois_hr) %>%
  print(n = 100)

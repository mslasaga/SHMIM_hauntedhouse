# 2.3 Ordinal distress: Day 1 (primary) + longitudinal. Stan: ordinal_probit.stan.
# Out: models/fit_ord_d1.rds, models/fit_ord_d1_hr.rds, models/fit_ord_long*.rds

source('scripts/0. Setup/0.0 Libraries and Functions.R')
source('scripts/1. Processing/1.0 Data Processing.R')

mod_ord <- cmdstan_model('scripts/stan/ordinal_probit.stan')

# ---- Stan data ----

sdata_ord_d1 <- within(list(), {
  N = nrow(dat_distress_d1_all)
  subject_id = as.numeric(factor(dat_distress_d1_all$id))
  N_subj = length(unique(subject_id))
  y = dat_distress_d1_all$distress
  K = length(unique(dat_distress_d1_all$distress))
  X_mu = model.matrix(~ 1 + scare_level + total_sh,
                      data = data.frame(scare_level = dat_distress_d1_all$scare_level,
                                        total_sh = dat_distress_d1_all$total_sh))[, -1, drop = FALSE]
  K_mu = ncol(X_mu)
  Z_mu = model.matrix(~ 1, data = dat_distress_d1_all)
  K_random_mu = ncol(Z_mu)
})

sdata_ord_d1_hr <- within(list(), {
  N = nrow(dat_distress_d1_hr)
  subject_id = as.numeric(factor(dat_distress_d1_hr$id))
  N_subj = length(unique(subject_id))
  y = dat_distress_d1_hr$distress
  K = length(unique(dat_distress_d1_hr$distress))
  X_mu = model.matrix(~ 1 + scare_level + total_sh + mean_hr,
                      data = data.frame(scare_level = dat_distress_d1_hr$scare_level,
                                        total_sh = dat_distress_d1_hr$total_sh,
                                        mean_hr = dat_distress_d1_hr$mean_hr))[, -1, drop = FALSE]
  K_mu = ncol(X_mu)
  Z_mu = model.matrix(~ 1, data = dat_distress_d1_hr)
  K_random_mu = ncol(Z_mu)
})

# ---- fit ----

fit_ord_d1 <- mod_ord$sample(data = sdata_ord_d1,
               chains = 10,
               parallel_chains = 10,
               refresh = 50,
               seed = 2026,
               iter_warmup = 4000,
               iter_sampling = 6000)

fit_ord_d1$save_object(file = 'models/fit_ord_d1.rds')

fit_ord_d1_hr <- mod_ord$sample(data = sdata_ord_d1_hr,
                    chains = 10,
                    parallel_chains = 10,
                    refresh = 50,
                    seed = 2026,
                    iter_warmup = 4000,
                    iter_sampling = 6000)

fit_ord_d1_hr$save_object(file = 'models/fit_ord_d1_hr.rds')

# ---- summaries ----

extract_and_summarise_ord(fit_ord_d1, sdata_ord_d1) %>%
  print(n = 100)

extract_and_summarise_ord(fit_ord_d1_hr, sdata_ord_d1_hr) %>%
  print(n = 100)

# ---- Stan data ----

sdata_ord_long <- within(list(), {
  N = nrow(dat_distress_all)
  subject_id = as.numeric(factor(dat_distress_all$id))
  N_subj = length(unique(subject_id))
  y = dat_distress_all$distress
  K = length(unique(dat_distress_all$distress))
  X_mu = model.matrix(~ 1 + day_after*scare_level + total_sh,
                      data = data.frame(day_after = dat_distress_all$day_after,
                                        scare_level = dat_distress_all$scare_level,
                                        total_sh = dat_distress_all$total_sh))[, -1, drop = FALSE]
  K_mu = ncol(X_mu)
  Z_mu = model.matrix(~ 1 + day_after, data = data.frame(day_after = dat_distress_all$day_after))
  K_random_mu = ncol(Z_mu)
})

sdata_ord_long_hr <- within(list(), {
  N = nrow(dat_distress_hr)
  subject_id = as.numeric(factor(dat_distress_hr$id))
  N_subj = length(unique(subject_id))
  y = dat_distress_hr$distress
  K = length(unique(dat_distress_hr$distress))
  X_mu = model.matrix(~ 1 + day_after*scare_level + total_sh + mean_hr,
                      data = data.frame(day_after = dat_distress_hr$day_after,
                                        scare_level = dat_distress_hr$scare_level,
                                        total_sh = dat_distress_hr$total_sh,
                                        mean_hr = dat_distress_hr$mean_hr))[, -1, drop = FALSE]
  K_mu = ncol(X_mu)
  Z_mu = model.matrix(~ 1 + day_after, data = data.frame(day_after = dat_distress_hr$day_after))
  K_random_mu = ncol(Z_mu)
})

# ---- fit ----

fit_ord_long <- mod_ord$sample(data = sdata_ord_long,
                   chains = 10,
                   parallel_chains = 10,
                   refresh = 50,
                   seed = 2026,
                   adapt_delta = 0.95,
                   iter_warmup = 4000,
                   iter_sampling = 6000)

fit_ord_long$save_object(file = 'models/fit_ord_long.rds')

fit_ord_long_hr <- mod_ord$sample(data = sdata_ord_long_hr,
                      chains = 10,
                      parallel_chains = 10,
                      refresh = 50,
                      seed = 2026,
                      iter_warmup = 4000,
                      iter_sampling = 6000)

fit_ord_long_hr$save_object(file = 'models/fit_ord_long_hr.rds')

# ---- summaries ----

extract_and_summarise_ord(fit_ord_long, sdata_ord_long) %>%
  print(n = 100)

extract_and_summarise_ord(fit_ord_long_hr, sdata_ord_long_hr) %>%
  print(n = 100)

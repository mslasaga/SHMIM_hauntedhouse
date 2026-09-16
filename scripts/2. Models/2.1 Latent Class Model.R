# 2.1 Latent-class ZIP. Stan: zip_latent_class.stan.
# Out: models/fit_zip_lc.rds, models/fit_zip_lc_hr.rds

source('scripts/0. Setup/0.0 Libraries and Functions.R')
source('scripts/1. Processing/1.0 Data Processing.R')

mod_zip_lc <- cmdstan_model('scripts/stan/zip_latent_class.stan')

# ---- Stan data ----

sdata_lc <- within(list(), {
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
  X_zi = model.matrix(~ 1 + total_sh,
                      data = data.frame(total_sh = dat_im_all[!duplicated(dat_im_all$id), ]$total_sh))
  K_zi = ncol(X_zi)
  zi_intercept_pos = 1
  Z = model.matrix(~ 1 + day_after, data = data.frame(day_after = dat_im_all$day_after))
  K_random = ncol(Z)
  lik_ind = 0
})

sdata_lc_hr <- within(list(), {
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
  X_zi = model.matrix(~ 1 + total_sh,
                      data = data.frame(total_sh = dat_im_hr[!duplicated(dat_im_hr$id), ]$total_sh))
  K_zi = ncol(X_zi)
  zi_intercept_pos = 1
  Z = model.matrix(~ 1 + day_after, data = data.frame(day_after = dat_im_hr$day_after))
  K_random = ncol(Z)
  lik_ind = 0
})

# ---- fit ----

fit_zip_lc <- mod_zip_lc$sample(data = sdata_lc,
               chains = 10,
               parallel_chains = 10,
               refresh = 50,
               seed = 2026,
               iter_warmup = 4000,
               iter_sampling = 6000)

fit_zip_lc$save_object(file = 'models/fit_zip_lc.rds')

fit_zip_lc_hr <- mod_zip_lc$sample(data = sdata_lc_hr,
                    chains = 10,
                    parallel_chains = 10,
                    refresh = 50,
                    seed = 2026,
                    iter_warmup = 4000,
                    iter_sampling = 6000)

fit_zip_lc_hr$save_object(file = 'models/fit_zip_lc_hr.rds')

# ---- summaries ----

extract_and_summarise_lc(fit_zip_lc, sdata_lc) %>%
  print(n = 100)

extract_and_summarise_lc(fit_zip_lc_hr, sdata_lc_hr) %>%
  print(n = 100)

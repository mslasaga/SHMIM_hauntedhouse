# 1.0 Raw CSVs -> analysis data sets.

source('scripts/0. Setup/0.0 Libraries and Functions.R')

# ---- raw questionnaire data ----

dat_diary <- read.csv('data/post_7_days_IM.csv') %>%
  janitor::clean_names()

dat_quest <- read.csv('data/Pre_post_questionnaire.csv') %>%
  janitor::clean_names() %>%
  mutate(across(c(cog_dis1_distress:control3_cantcontrol), ~ recode(.,
                                                                    "Strongly disagree"            = 1,
                                                                    "Disagree"                     = 2,
                                                                    "Somewhat disagree"            = 3,
                                                                    "Neither agree nor disagree"   = 4,
                                                                    "Somewhat agree"               = 5,
                                                                    "Agree"                        = 6,
                                                                    "Strongly agree"               = 7
  ))) %>%
  mutate(general_im = ordered(general_im,
                              levels = c('<10%', '10-20%', '21-40%', '41-60%', '61-80%', '>80%')),

         sh_cogdis = rowMeans(select(., cog_dis1_distress:cog_dis3_uncomfortable), na.rm = TRUE),
         sh_viv    = rowMeans(select(., vivid1_vivid:vivid3_lifelike), na.rm = TRUE),
         sh_phys   = rowMeans(select(., phys1_temperature:phys3_breathing), na.rm = TRUE),
         sh_int    = rowMeans(select(., intent1_unintentional:intent3_notthinking), na.rm = TRUE),
         sh_con    = rowMeans(select(., control1_cantstop:control3_cantcontrol), na.rm = TRUE))

dat_all <- dat_diary %>%
  left_join(dat_quest) %>%
  mutate(day_after_zero = day_after - 1,
         total_sh = rowMeans(select(., sh_cogdis:sh_con)))

# ---- heart-rate data ----

dat_hr_raw <- read.csv('data/Processed HR Data.csv', skip = 0, header = TRUE) %>%
  mutate(id = as.integer(sub("^0*([0-9]+)_.*", "\\1", FileName)),
         mean_hr = as.numeric(S1_Mean.HR..bpm.)) %>%
  filter(!is.na(mean_hr)) %>%
  select(id, mean_hr)

dup_hr_ids <- unique(dat_hr_raw$id[duplicated(dat_hr_raw$id)])
if (length(dup_hr_ids) > 0) {
  message("Duplicate HR ids detected (first row kept): ",
          paste(dup_hr_ids, collapse = ", "))
}

dat_hr <- dat_hr_raw %>%
  distinct(id, .keep_all = TRUE)

# ---- analysis data sets ----

dat_im <- dat_all %>%
  select(id, im_frequency, day_after_zero, scare_level, sh_cogdis:sh_con, total_sh) %>%
  filter(complete.cases(.)) %>%
  mutate(cens_type = ifelse(im_frequency == 10, "right", "none")) %>%
  mutate_at(vars(sh_cogdis:sh_con, total_sh, scare_level), standardize)

dat_im_all <- dat_im %>%
  mutate(day_after = day_after_zero + 1) %>%
  arrange(id)

dat_im_hr <- dat_im_all %>%
  left_join(dat_hr, by = 'id') %>%
  filter(!is.na(mean_hr)) %>%
  mutate(mean_hr = standardize(mean_hr))

dat_distress_all <- dat_all %>%
  select(id, im_frequency, day_after_zero, scare_level, total_sh, ends_with('distress')) %>%
  select(-cog_dis1_distress) %>%
  mutate(cens_type = ifelse(im_frequency == 10, "right", "none")) %>%
  mutate_at(vars(total_sh, scare_level), standardize) %>%
  mutate(day_after = day_after_zero + 1) %>%
  pivot_longer(ends_with('distress'), names_to = 'intrusion_num', values_to = 'distress') %>%
  filter(!is.na(distress)) %>%
  filter(complete.cases(.)) %>%
  mutate(distress = distress + 1) %>%
  arrange(id)

dat_distress_hr <- dat_distress_all %>%
  left_join(dat_hr, by = 'id') %>%
  filter(!is.na(mean_hr)) %>%
  mutate(mean_hr = standardize(mean_hr))

dat_distress_d1_all <- dat_distress_all %>% filter(day_after == 1)
dat_distress_d1_hr  <- dat_distress_hr  %>% filter(day_after == 1)


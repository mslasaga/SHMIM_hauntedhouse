# 3.1 LCA figure. Prereq: 2.1. Output: figures/fig_lca.jpg

source('scripts/0. Setup/0.0 Libraries and Functions.R')
source('scripts/1. Processing/1.0 Data Processing.R')

if (!exists('fit_zip_lc')) fit_zip_lc <- readRDS('models/fit_zip_lc.rds')

# ---- fixed-effect draws ----
bmat <- fit_zip_lc$draws("beta_count", format = "draws_matrix")
zmat <- fit_zip_lc$draws("beta_zi", format = "draws_matrix")
bd <- bind_cols(
  as_tibble(matrix(as.numeric(bmat), nrow = nrow(bmat),
                   dimnames = list(NULL, c("b_int", "b_day", "b_scare",
                                           "b_sh", "b_dayXscare")))),
  as_tibble(matrix(as.numeric(zmat), nrow = nrow(zmat),
                   dimnames = list(NULL, c("z_int", "z_sh")))))

# ---- participants, in the order Stan indexed them (sorted id) ----
persons <- dat_im_all %>%
  distinct(id, scare_level, total_sh) %>%
  arrange(id) %>%
  mutate(scare_group = factor(ntile(scare_level, 3), levels = 1:3,
                              labels = c("Low", "Mid", "High")))
stopifnot(identical(persons$id, sort(unique(dat_im_all$id))))

U     <- re_draws(fit_zip_lc, "u", nrow(persons))
Mmat  <- re_M(fit_zip_lc, "u", days = 1:7)
dat_g <- dat_im_all %>% left_join(select(persons, id, scare_group), by = "id")

pri_i <- plogis(bd$z_int + outer(bd$z_sh, persons$total_sh))

# ---- P1: day_after x scare interaction (marginal count) ----

line_int <- gcomp_group_days(bd, U, persons, dat_g, days = 1:7,
                             observed_only = TRUE, pri = pri_i)

raw_int <- dat_g %>%
  group_by(scare_group, day_after) %>%
  summarise(mean_y = mean(im_frequency), se = sd(im_frequency) / sqrt(n()),
            n = n(), .groups = "drop")

dodge <- position_dodge(width = 0.45)

p_int <- ggplot(line_int, aes(day_after, .value, colour = scare_group, fill = scare_group)) +
  geom_ribbon(aes(ymin = .lower95, ymax = .upper95), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_linerange(data = raw_int, position = dodge, linewidth = 0.4,
                 aes(x = day_after, ymin = mean_y - se, ymax = mean_y + se,
                     colour = scare_group), inherit.aes = FALSE) +
  geom_point(data = raw_int, position = dodge, shape = 16,
             aes(day_after, mean_y, colour = scare_group, size = n),
             inherit.aes = FALSE) +
  scale_colour_manual(values = haunt_pal$scare, name = "Scare Level") +
  scale_fill_manual(values = haunt_pal$scare, guide = "none") +
  scale_size_area(max_size = 3, guide = "none") +
  scale_x_continuous(breaks = 1:7, expand = c(0.02, 0)) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_haunt() +
  labs(x = "Day After Exposure", y = "Intrusions",
       subtitle = "Day After Exposure × Scare Level")

# ---- P2: total_sh main effect (marginal count) ----

Eta0 <- outer(bd$b_int, rep(1, 7)) + outer(bd$b_day, 1:7)
A  <- rowMeans(exp(Eta0) * Mmat)
S0 <- rowSums(exp(Eta0))
sh_grid <- seq(min(dat_im_all$total_sh), max(dat_im_all$total_sh), length.out = 200)

line_sh <- tibble(A = A, b_sh = bd$b_sh, z_int = bd$z_int, z_sh = bd$z_sh) %>%
  crossing(total_sh = sh_grid) %>%
  mutate(epred = plogis(z_int + z_sh * total_sh) * A * exp(b_sh * total_sh)) %>%
  group_by(total_sh) %>%
  reframe(post_band(epred))

raw_sh <- binned_mean(dat_im_all, total_sh, im_frequency, nbin = 6)

p_sh <- ggplot(line_sh, aes(total_sh, .value)) +
  geom_ribbon(aes(ymin = .lower95, ymax = .upper95), fill = haunt_pal$fill95) +
  geom_line(colour = haunt_pal$line, linewidth = 0.9) +
  geom_linerange(data = raw_sh, aes(x = x, ymin = mean_y - se, ymax = mean_y + se),
                 colour = haunt_pal$line, linewidth = 0.4, inherit.aes = FALSE) +
  geom_point(data = raw_sh, aes(x, mean_y, size = n), colour = haunt_pal$line,
             shape = 16, inherit.aes = FALSE) +
  scale_size_area(max_size = 3, guide = "none") +
  scale_x_continuous(expand = c(0.02, 0)) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_haunt() +
  labs(x = "Total SHMIM (Scaled)", y = "Intrusions",
       subtitle = "Total SHMIM")

# ---- P3: total_sh -> P(any intrusion) + latent P(susceptible) ----

lz <- tibble(z_int = bd$z_int, z_sh = bd$z_sh, b_sh = bd$b_sh, S0 = S0)

line_any <- lz %>% crossing(total_sh = sh_grid) %>%
  mutate(epred = plogis(z_int + z_sh * total_sh) *
           (1 - exp(-exp(b_sh * total_sh) * S0))) %>%
  group_by(total_sh) %>% reframe(post_band(epred))

line_pri <- lz %>% crossing(total_sh = sh_grid) %>%
  mutate(epred = plogis(z_int + z_sh * total_sh)) %>%
  group_by(total_sh) %>% reframe(post_band(epred))

raw_zi <- dat_im_all %>%
  group_by(id, total_sh) %>%
  summarise(any_im = as.integer(any(im_frequency > 0)), .groups = "drop")

p_zi <- ggplot(line_pri, aes(total_sh, .value)) +
  geom_ribbon(aes(ymin = .lower95, ymax = .upper95), fill = haunt_pal$fill95) +
  geom_line(colour = haunt_pal$line, linewidth = 0.9) +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.02, 0)) +
  theme_haunt() +
  labs(x = "Total SHMIM (Scaled)", y = "P(Susceptible Class)",
       subtitle = "Class Membership")

# ---- compose ----
fig_lca <- p_int + p_sh + p_zi

ggsave('figures/fig_lca.jpg', fig_lca, width = 13.5, height = 4.4, dpi = 700)

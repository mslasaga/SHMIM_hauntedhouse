# 3.0 Figure 1. Prereq: 2.0. Output: figures/fig_base.jpg

source('scripts/0. Setup/0.0 Libraries and Functions.R')
source('scripts/1. Processing/1.0 Data Processing.R')

if (!exists('fit_pois')) fit_pois <- readRDS('models/fit_pois_cens.rds')

# ---- fixed-effect draws ----
bmat <- fit_pois$draws("beta_count", format = "draws_matrix")
bd <- as_tibble(matrix(as.numeric(bmat), nrow = nrow(bmat),
                       dimnames = list(NULL, c("b_int", "b_day", "b_scare",
                                               "b_sh", "b_dayXscare"))))

# ---- participants, in the order Stan indexed them (sorted id) ----
persons <- dat_im_all %>%
  distinct(id, scare_level, total_sh) %>%
  arrange(id) %>%
  mutate(scare_group = factor(ntile(scare_level, 3), levels = 1:3,
                              labels = c("Low", "Mid", "High")))
stopifnot(identical(persons$id, sort(unique(dat_im_all$id))))

U      <- re_draws(fit_pois, "u_count", nrow(persons))
Mmat   <- re_M(fit_pois, "u_count", days = 1:7)
dat_g  <- dat_im_all %>% left_join(select(persons, id, scare_group), by = "id")

# ---- P1: day_after x scare interaction ----

line_int <- gcomp_group_days(bd, U, persons, dat_g, days = 1:7, observed_only = TRUE)

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

# ---- P2: total_sh main effect ----

Eta0 <- outer(bd$b_int, rep(1, 7)) + outer(bd$b_day, 1:7)
A <- rowMeans(exp(Eta0) * Mmat)
sh_grid <- seq(min(dat_im_all$total_sh), max(dat_im_all$total_sh), length.out = 200)

line_sh <- tibble(A = A, b_sh = bd$b_sh) %>%
  crossing(total_sh = sh_grid) %>%
  mutate(epred = A * exp(b_sh * total_sh)) %>%
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

# ---- compose ----
fig_base <- p_int + p_sh

ggsave('figures/fig_base.jpg', fig_base, width = 9.5, height = 4.2, dpi = 700)

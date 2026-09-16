# 3.2 Figure 2. Prereq: 2.3. Output: figures/fig_distress_d1.jpg

source('scripts/0. Setup/0.0 Libraries and Functions.R')
source('scripts/1. Processing/1.0 Data Processing.R')

if (!exists('fit_ord_d1')) fit_ord_d1 <- readRDS('models/fit_ord_d1.rds')

K <- length(unique(dat_distress_d1_all$distress))

# ---- draws ----
as_plain <- function(dm) matrix(as.numeric(dm), nrow = nrow(dm), dimnames = list(NULL, colnames(dm)))
bm      <- as_plain(fit_ord_d1$draws("beta_mu",   format = "draws_matrix"))
cp      <- as_plain(fit_ord_d1$draws("cutpoints", format = "draws_matrix"))
b_sh    <- bm[, "beta_mu[2]"]

er_curve <- function(gvec, coef) {
  map_dfr(gvec, function(g) {
    eta <- coef * g
    er  <- 1 + rowSums(plogis(eta - cp))
    bind_cols(tibble(x = g), post_band(er))
  })
}

er_panel <- function(line_df, raw_bin, xlab) {
  ggplot(line_df, aes(x, .value)) +
    geom_ribbon(aes(ymin = .lower95, ymax = .upper95), fill = haunt_pal$fill95) +
    geom_line(colour = haunt_pal$line, linewidth = 0.9) +
    geom_linerange(data = raw_bin, aes(x = x, ymin = mean_y - se, ymax = mean_y + se),
                   colour = haunt_pal$line, linewidth = 0.4, inherit.aes = FALSE) +
    geom_point(data = raw_bin, aes(x, mean_y, size = n), colour = haunt_pal$line,
               shape = 16, inherit.aes = FALSE) +
    scale_size_area(max_size = 3, guide = "none") +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(breaks = 1:K, labels = 0:(K - 1), expand = c(0.02, 0)) + 
    coord_cartesian(ylim = c(1, K)) +
    theme_haunt() +
    labs(x = xlab, y = "Distress Rating")
}

# ---- total_sh ----
sh_grid <- seq(min(dat_distress_d1_all$total_sh), max(dat_distress_d1_all$total_sh),
               length.out = 200)
p_sh <- er_panel(
  er_curve(sh_grid, b_sh),
  binned_mean(dat_distress_d1_all, total_sh, distress, nbin = 6),
  "Total SHMIM (Scaled)")

ggsave('figures/fig_distress_d1.jpg', p_sh, width = 4.2, height = 4.2, dpi = 700)

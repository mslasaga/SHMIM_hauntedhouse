# 0.0 Packages and shared helpers. Sourced by every 1.x-4.x script.

# ---- libraries ----

library(cmdstanr)
library(brms)
library(tidybayes)
library(posterior)
library(modelr)
library(tidyverse)
library(patchwork)

# ---- general helpers ----

standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

check_overdispersion <- function(model) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Please install the 'brms' package to use this function.")
  }

  ll_mat <- log_lik(model, combine_chains = FALSE)

  deviance_draws <- -2 * rowSums(ll_mat)

  mean_deviance <- mean(deviance_draws)

  n_data <- nrow(model$data)

  num_fixed <- length(fixef(model))

  ranef_list <- VarCorr(model)
  num_ranpars <- sum(sapply(ranef_list, function(x) length(x[[1]])))

  param_count <- num_fixed + num_ranpars

  df_approx <- max(1, n_data - param_count)

  ratio_dev <- mean_deviance / df_approx

  list(
    mean_deviance = mean_deviance,
    df_approx = df_approx,
    ratio_dev = ratio_dev,
    interpretation = ifelse(ratio_dev > 1.5,
                            "Overdispersion likely present. Consider a negative binomial or zero-inflated model.",
                            "No strong evidence of overdispersion under this heuristic.")
  )
}

# ---- figure helpers ----

haunt_pal <- list(
  line   = "#31688E",
  scare  = c(Low = "#414487", Mid = "#25848E", High = "#43BF71"),
  fill95 = "#E4E4E4"
)

theme_haunt <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line     = element_line(linewidth = 0.4, colour = "grey25"),
      axis.ticks    = element_line(linewidth = 0.4, colour = "grey25"),
      axis.title    = element_text(colour = "grey10"),
      axis.text     = element_text(colour = "grey30"),
      plot.subtitle = element_text(face = "bold", size = rel(0.95), colour = "grey10"),
      plot.tag      = element_text(face = "bold", size = rel(1.1)),
      legend.position   = "bottom",
      legend.title      = element_text(face = "bold", size = rel(0.9)),
      legend.key.height = unit(10, "pt"),
      plot.margin       = margin(6, 10, 6, 6)
    )
}

theme_haunt_min <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
      axis.title    = element_text(colour = "grey10"),
      axis.text     = element_text(colour = "grey30"),
      axis.ticks    = element_line(linewidth = 0.3, colour = "grey60"),
      plot.subtitle = element_text(face = "bold", size = rel(0.95), colour = "grey10"),
      plot.tag      = element_text(face = "bold", size = rel(1.1)),
      legend.position   = "bottom",
      legend.title      = element_text(face = "bold", size = rel(0.9)),
      legend.key.height = unit(10, "pt"),
      plot.margin       = margin(6, 10, 6, 6)
    )
}

post_band <- function(draws) {
  q <- unname(quantile(draws, c(0.025, 0.25, 0.5, 0.75, 0.975)))
  tibble(.value = q[3], .lower50 = q[2], .upper50 = q[4],
         .lower95 = q[1], .upper95 = q[5])
}

re_M <- function(fit, uvar, days = 1:7) {
  m  <- fit$draws(uvar, format = "draws_matrix")
  cn <- colnames(m)
  pick <- function(jj) {
    cols <- grep(paste0(",", jj, "\\]$"), cn, value = TRUE)
    cols <- cols[order(as.integer(sub(".*\\[(\\d+),.*", "\\1", cols)))]
    matrix(as.numeric(m[, cols]), nrow = nrow(m))
  }
  Uint <- pick(1); Uday <- pick(2)
  sapply(days, function(d) rowMeans(exp(Uint + d * Uday)))
}

binned_mean <- function(df, xvar, yvar, nbin = 6) {
  df %>%
    mutate(.bin = ntile({{ xvar }}, nbin)) %>%
    group_by(.bin) %>%
    summarise(x = mean({{ xvar }}), mean_y = mean({{ yvar }}),
              se = sd({{ yvar }}) / sqrt(n()), n = n(), .groups = "drop")
}

# ---- extraction helpers for the raw Stan models ----

extract_and_summarise_lc <- function(fit, sdata) {

  draws_df <- fit$draws(c('beta_count_gq', 'beta_zi_gq', 'sigma_u_gq', 'Omega')) %>%
    as_draws_df()

  if ("X_count" %in% names(sdata) && !is.null(colnames(sdata$X_count))) {
    bc_names <- paste0("count_", colnames(sdata$X_count))
    names_map <- setNames(grep("^beta_count_gq\\[", names(draws_df), value = TRUE), bc_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)
  }

  if ("X_zi" %in% names(sdata) && !is.null(colnames(sdata$X_zi))) {
    zi_names <- paste0("zi_", colnames(sdata$X_zi))
    names_map <- setNames(grep("^beta_zi_gq\\[", names(draws_df), value = TRUE), zi_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)
  }

  if ("Z" %in% names(sdata) && !is.null(colnames(sdata$Z))) {
    ru_names <- paste0("sd_", colnames(sdata$Z))
    names_map <- setNames(grep("^sigma_u_gq\\[", names(draws_df), value = TRUE), ru_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)
    K <- ncol(sdata$Z)
    omega_cols <- grep("^Omega\\[", names(draws_df), value = TRUE)

    omega_cols_keep <- sapply(omega_cols, function(col) {
      indices <- gsub("Omega\\[|\\]", "", col) %>% strsplit(",") %>% unlist() %>% as.integer()
      i <- indices[1]
      j <- indices[2]
      return(i > j)
    })

    omega_cols <- omega_cols[omega_cols_keep]

    omega_names <- sapply(omega_cols, function(col) {
      indices <- gsub("Omega\\[|\\]", "", col) %>% strsplit(",") %>% unlist() %>% as.integer()
      i <- indices[1]
      j <- indices[2]
      paste0("cor(", colnames(sdata$Z)[j], ",", colnames(sdata$Z)[i], ")")
    })

    names_map <- setNames(omega_cols, omega_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)

    all_omega_cols <- grep("^Omega\\[", names(draws_df), value = TRUE)
    cols_to_remove <- setdiff(all_omega_cols, omega_cols)
    draws_df <- draws_df[, !names(draws_df) %in% cols_to_remove]
  }

  draws_df %>%
    summarise_draws(mean, median, sd,
                    ~quantile2(.x, probs = c(0.025, 0.975)),
                    rhat, ess_bulk, ess_tail)
}

extract_and_summarise_ord <- function(fit, sdata) {

  draws_df <- fit$draws(c('beta_mu_gq', 'sigma_u_gq_mu', 'Omega_mu',
                          'cutpoint_intercepts_gq')) %>%
    as_draws_df()

  if (!is.null(sdata$X_mu) && !is.null(colnames(sdata$X_mu))) {
    bc_names <- paste0("mu_", colnames(sdata$X_mu))
    names_map <- setNames(grep("^beta_mu_gq\\[", names(draws_df), value = TRUE), bc_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)
  }

  if (!is.null(sdata$Z_mu) && !is.null(colnames(sdata$Z_mu))) {
    ru_names <- paste0("sd_", colnames(sdata$Z_mu))
    names_map <- setNames(grep("^sigma_u_gq_mu\\[", names(draws_df), value = TRUE), ru_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)

    K <- ncol(sdata$Z_mu)
    omega_cols <- grep("^Omega_mu\\[", names(draws_df), value = TRUE)
    omega_cols_keep <- sapply(omega_cols, function(col) {
      indices <- gsub("Omega_mu\\[|\\]", "", col) %>% strsplit(",") %>% unlist() %>% as.integer()
      i <- indices[1]; j <- indices[2]
      return(i > j)
    })
    omega_cols <- omega_cols[omega_cols_keep]

    omega_names <- sapply(omega_cols, function(col) {
      indices <- gsub("Omega_mu\\[|\\]", "", col) %>% strsplit(",") %>% unlist() %>% as.integer()
      i <- indices[1]; j <- indices[2]
      paste0("cor(", colnames(sdata$Z_mu)[j], ",", colnames(sdata$Z_mu)[i], ")")
    })
    names_map <- setNames(omega_cols, omega_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)

    all_omega_cols <- grep("^Omega_mu\\[", names(draws_df), value = TRUE)
    cols_to_remove <- setdiff(all_omega_cols, omega_cols)
    draws_df <- draws_df[, !names(draws_df) %in% cols_to_remove]
  }

  cp_cols <- grep("^cutpoint_intercepts_gq\\[", names(draws_df), value = TRUE)
  if (length(cp_cols)) {
    cp_names <- paste0("cutpoint_", seq_along(cp_cols))
    draws_df <- dplyr::rename(draws_df, !!!setNames(cp_cols, cp_names))
  }

  draws_df %>%
    summarise_draws(
      mean, median, sd,
      ~quantile2(.x, probs = c(0.025, 0.975)),
      rhat, ess_bulk, ess_tail
    )
}

extract_and_summarise_base <- function(fit, sdata) {

  draws_df <- fit$draws(c('beta_count_gq', 'sigma_u_gq_count', 'Omega_count')) %>%
    as_draws_df()

  if ("X_count" %in% names(sdata) && !is.null(colnames(sdata$X_count))) {
    bc_names <- paste0("count_", colnames(sdata$X_count))
    names_map <- setNames(grep("^beta_count_gq\\[", names(draws_df), value = TRUE), bc_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)
  }

  if ("X_zi" %in% names(sdata) && !is.null(colnames(sdata$X_zi))) {
    zi_names <- paste0("zi_", colnames(sdata$X_zi))
    names_map <- setNames(grep("^beta_zi_gq\\[", names(draws_df), value = TRUE), zi_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)
  }

  if ("Z_count" %in% names(sdata) && !is.null(colnames(sdata$Z_count))) {
    ru_names <- paste0("sd_count_", colnames(sdata$Z_count))
    names_map <- setNames(grep("^sigma_u_gq_count\\[", names(draws_df), value = TRUE), ru_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)
    K <- ncol(sdata$Z_count)
    omega_cols <- grep("^Omega_count\\[", names(draws_df), value = TRUE)

    omega_cols_keep <- sapply(omega_cols, function(col) {
      indices <- gsub("Omega_count\\[|\\]", "", col) %>% strsplit(",") %>% unlist() %>% as.integer()
      i <- indices[1]
      j <- indices[2]
      return(i > j)
    })

    omega_cols <- omega_cols[omega_cols_keep]

    omega_names <- sapply(omega_cols, function(col) {
      indices <- gsub("Omega_count\\[|\\]", "", col) %>% strsplit(",") %>% unlist() %>% as.integer()
      i <- indices[1]
      j <- indices[2]
      paste0("cor_count(", colnames(sdata$Z_count)[j], ",", colnames(sdata$Z_count)[i], ")")
    })

    names_map <- setNames(omega_cols, omega_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)

    all_omega_cols <- grep("^Omega_count\\[", names(draws_df), value = TRUE)
    cols_to_remove <- setdiff(all_omega_cols, omega_cols)
    draws_df <- draws_df[, !names(draws_df) %in% cols_to_remove]
  }

  if ("Z_zi" %in% names(sdata) && !is.null(colnames(sdata$Z_zi))) {
    ru_names <- paste0("sd_zi_", colnames(sdata$Z_zi))
    names_map <- setNames(grep("^sigma_u_gq_zi\\[", names(draws_df), value = TRUE), ru_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)
    K <- ncol(sdata$Z_zi)
    omega_cols <- grep("^Omega_zi\\[", names(draws_df), value = TRUE)

    omega_cols_keep <- sapply(omega_cols, function(col) {
      indices <- gsub("Omega_zi\\[|\\]", "", col) %>% strsplit(",") %>% unlist() %>% as.integer()
      i <- indices[1]
      j <- indices[2]
      return(i > j)
    })

    omega_cols <- omega_cols[omega_cols_keep]

    omega_names <- sapply(omega_cols, function(col) {
      indices <- gsub("Omega_zi\\[|\\]", "", col) %>% strsplit(",") %>% unlist() %>% as.integer()
      i <- indices[1]
      j <- indices[2]
      paste0("cor_zi(", colnames(sdata$Z_zi)[j], ",", colnames(sdata$Z_zi)[i], ")")
    })

    names_map <- setNames(omega_cols, omega_names)
    draws_df <- dplyr::rename(draws_df, !!!names_map)

    all_omega_cols <- grep("^Omega_zi\\[", names(draws_df), value = TRUE)
    cols_to_remove <- setdiff(all_omega_cols, omega_cols)
    draws_df <- draws_df[, !names(draws_df) %in% cols_to_remove]
  }

  draws_df %>%
    summarise_draws(mean, median, sd,
                    ~quantile2(.x, probs = c(0.025, 0.975)),
                    rhat, ess_bulk, ess_tail)
}

# ---- G-computation over the fitted participants ----

re_draws <- function(fit, uvar, np) {
  m <- fit$draws(uvar, format = "draws_matrix")
  list(int = matrix(as.numeric(m[, paste0(uvar, "[", seq_len(np), ",1]")]), nrow = nrow(m)),
       day = matrix(as.numeric(m[, paste0(uvar, "[", seq_len(np), ",2]")]), nrow = nrow(m)))
}

epred_persons <- function(bd, U, persons, d, k, pri = 1) {
  p <- if (is.matrix(pri)) pri[, k, drop = FALSE] else pri
  p * exp(outer(bd$b_int, rep(1, length(k))) + bd$b_day * d +
            outer(bd$b_scare + d * bd$b_dayXscare, persons$scare_level[k]) +
            outer(bd$b_sh, persons$total_sh[k]) +
            U$int[, k, drop = FALSE] + d * U$day[, k, drop = FALSE])
}

gcomp_group_days <- function(bd, U, persons, dat_g, days = 1:7,
                             observed_only = TRUE, pri = 1) {
  NP  <- nrow(persons)
  idx <- setNames(seq_len(NP), persons$id)
  grp <- persons$scare_group
  map_dfr(days, function(d) {
    map_dfr(levels(grp), function(g) {
      k <- if (observed_only) {
        unname(idx[as.character(dat_g$id[dat_g$scare_group == g & dat_g$day_after == d])])
      } else which(grp == g)
      bind_cols(tibble(day_after = d, scare_group = g),
                post_band(rowMeans(epred_persons(bd, U, persons, d, k, pri))))
    })
  }) %>% mutate(scare_group = factor(scare_group, levels = levels(grp)))
}

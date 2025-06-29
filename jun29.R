# ---------------------------------------------------------------------------
# Entry/Exit Times using Negative Log-Likelihood under GAL Noise Model
# ---------------------------------------------------------------------------
# Author: Said Maanan | Updated: 2025-06-21
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(gt)

# ---------------------- GAL Log Density -----------------------------------
neg_log_gal_density <- function(e, tau, z, b) {
  log_term <- ifelse(e < 0,
                     -log(tau) - log(1 + z) + abs((1 + z) * e) / b,
                     -log(1 - tau) - log(1 - z) + abs((1 - z) * e) / b)
  log_term + log(b)
}

gal_loss <- function(I, a, l, s, tau, z, b) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  e <- I - fitted
  mean(neg_log_gal_density(e, tau, z, b))
}

# --------------------- Estimation Functions -------------------------------
med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  o <- order(x)
  cum_w <- cumsum(w[o])
  cutoff <- tau * sum(w)
  x[o][min(which(cum_w >= cutoff))]
}

compute_a <- function(I, l, s, tau) med_weighted_tau(I, l, s, tau)

search_best_s <- function(I, l, tau) {
  f <- function(s) {
    a <- compute_a(I, l, s, tau)
    fitted <- a * exp(-((1:length(I)) - l)^2 / s^2)
    u <- I - fitted
    mean(u * (tau - (u < 0)))
  }
  optimize(f, c(1, 20), tol = 1e-4)$minimum
}

estimate_params <- function(I, tau, h) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau)
    a <- compute_a(I, l, s, tau)
    fitted <- a * exp(-((1:length(I)) - l)^2 / s^2)
    u <- I - fitted
    mean(u * (tau - (u < 0)))
  })
  l_hat <- which.min(losses)
  s_hat <- search_best_s(I, l_hat, tau)
  a_hat <- compute_a(I, l_hat, s_hat, tau)
  list(l_hat = l_hat, s_hat = s_hat, a_hat = a_hat)
}

# ---------------------- Noise Generator -----------------------------------
rgen_asy_laplace <- function(n, tau, z, b) {
  u <- runif(n)
  ifelse(u < tau,
         -rexp(n, rate = (1 + z) / b),
         rexp(n, rate = (1 - z) / b))
}

asy_laplace_variance <- function(tau, z) {
  mu <- tau / (1 + z) - (1 - tau) / (1 - z)
  ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
  ex2 - mu^2
}

compute_noise_scale <- function(A, snr, tau, z) {
  sqrt((A / snr)^2 / asy_laplace_variance(tau, z))
}

# --------------------- Entry/Exit Detection -------------------------------
find_stable_time <- function(l_path, T_vals, tol) {
  for (i in seq_along(T_vals)) {
    future_vals <- l_path[i:length(l_path)]
    if (max(future_vals) - min(future_vals) <= tol) {
      return(T_vals[i])
    }
  }
  return(NA)
}

# ------------------ Main Wrapper Function ---------------------------------
compute_entry_exit_table <- function(snr_value, h = 35, seed = 280) {
  A_true <- 1
  L_true <- 22
  S_true <- 10
  T_vals <- 4:h
  z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
  taus <- c(0.25, 0.5, 0.75, 0.95)
  threshold <- 4
  threshold_prime <- 0
  results <- data.frame()
  
  for (tau in taus) {
    for (z in z_vals) {
      set.seed(seed)
      signal <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
      b <- compute_noise_scale(A_true, snr_value, tau, z)
      noise <- rgen_asy_laplace(h, tau, z, b)
      I <- signal + noise
      
      l_path <- sapply(T_vals, function(T) estimate_params(I[1:T], tau, h)$l_hat)
      
      T_entry <- find_stable_time(l_path, T_vals, threshold)
      T_exit <- if (!is.na(T_entry)) {
        find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)
      } else NA
      
      results <- rbind(results, data.frame(SNR = snr_value, tau = tau, z = z, T_entry = T_entry, T_exit = T_exit))
    }
  }
  
  # Pivot and format table
  results %>%
    mutate(tau = paste0("τ = ", tau)) %>%
    arrange(z, tau) %>%
    pivot_wider(names_from = tau, values_from = c(T_entry, T_exit)) %>%
    gt(rowname_col = "z") %>%
    tab_header(
      title = md(glue::glue("**Table: Entry/Exit Times using -log(f<sub>τ,z,b</sub>) Loss (SNR = {snr_value})**"))
    ) %>%
    fmt_number(columns = everything(), decimals = 1) %>%
    tab_options(
      table.font.size = "small",
      heading.title.font.size = 14,
      data_row.padding = px(3)
    )
}

# ------------------ Generate and Save Both Tables --------------------------
table_snr100 <- compute_entry_exit_table(snr_value = 100)
table_snr40  <- compute_entry_exit_table(snr_value = 40)

gtsave(table_snr100, "table_entry_exit_snr100.png")
gtsave(table_snr40,  "table_entry_exit_snr40.png")

writeLines(as.character(as_latex(table_snr100)), "table_entry_exit_snr100.tex")
writeLines(as.character(as_latex(table_snr40)),  "table_entry_exit_snr40.tex")

cat("Saved both entry/exit tables for SNR = 100 and SNR = 40.\n")

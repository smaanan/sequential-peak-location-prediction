
# ---------------------------------------------------------------------------
# Estimation Paths + Entry/Exit Times under Low SNR = 40
# Author: Said Maanan | Date: 2025-06-21
# ---------------------------------------------------------------------------

# --- Required Libraries ----------------------------------------------------
library(dplyr)
library(gt)
library(tidyr)

# --- Noise Generator and Helpers -------------------------------------------

rgen_asy_laplace <- function(n, tau, z, b = 1) {
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

quantile_loss <- function(u, tau) u * (tau - (u < 0))

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

objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  mean(quantile_loss(I - a * exp(-((t - l)^2) / s^2), tau))
}

search_best_s <- function(I, l, tau) {
  f <- function(s) objective_function(I, compute_a(I, l, s, tau), l, s, tau)
  optimize(f, c(1, 20), tol = 1e-4)$minimum
}

estimate_params <- function(I, tau, h) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau)
    a <- compute_a(I, l, s, tau)
    objective_function(I, a, l, s, tau)
  })
  l_hat <- which.min(losses)
  s_hat <- search_best_s(I, l_hat, tau)
  a_hat <- compute_a(I, l_hat, s_hat, tau)
  list(l_hat = l_hat, s_hat = s_hat, a_hat = a_hat)
}

# --- Stability Function ----------------------------------------------------
find_stable_time <- function(l_path, T_vals, tol) {
  for (i in seq_along(T_vals)) {
    future_vals <- l_path[i:length(l_path)]
    if (max(future_vals) - min(future_vals) <= tol) {
      return(T_vals[i])
    }
  }
  return(NA)
}

# --- Main Settings ---------------------------------------------------------
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 40
seed <- 280
h <- 35
T_vals <- 4:h
taus <- c(0.5, 0.75, 0.95)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
colors <- c("orange", "red", "black", "blue", "darkgreen")
threshold <- 4
threshold_prime <- 0

results <- data.frame()

# --- Plotting Setup --------------------------------------------------------
png("fig_lhat_paths_snr40_tau_z_overlay.png", width = 1200, height = 400, res = 120)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

for (tau in taus) {
  lhat_matrix <- matrix(NA, nrow = length(T_vals), ncol = length(z_vals))

  for (i in seq_along(z_vals)) {
    z <- z_vals[i]
    set.seed(seed)
    signal <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
    b <- compute_noise_scale(A_true, snr, tau, z)
    noise <- rgen_asy_laplace(h, tau, z, b)
    I <- signal + noise

    l_path <- sapply(T_vals, function(T) estimate_params(I[1:T], tau, h)$l_hat)
    lhat_matrix[, i] <- l_path

    T_entry <- find_stable_time(l_path, T_vals, threshold)
    T_exit <- find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)

    results <- rbind(results, data.frame(tau = tau, z = z, T_entry = T_entry, T_exit = T_exit))
  }

  ylim_range <- range(lhat_matrix, na.rm = TRUE)

  plot(NA, xlim = range(T_vals), ylim = ylim_range,
       xlab = "T", ylab = expression(hat(l)(T)),
       main = bquote(hat(l)(T) ~ "for" ~ tau == .(tau)))
  abline(h = L_true, col = "gray40", lty = 2)
  grid()

  for (i in seq_along(z_vals)) {
    lines(T_vals, lhat_matrix[, i], type = "b", col = colors[i], pch = 19)
  }

  legend("bottomright", legend = paste("z =", z_vals),
         col = colors, lwd = 2, pch = 19, cex = 0.8, bty = "o")
}

dev.off()
cat("Saved plot as 'fig_lhat_paths_snr40_tau_z_overlay.png'\n")

# --- Create Table of Entry/Exit Times --------------------------------------
table_out <- results %>%
  mutate(tau = paste0("tau = ", tau)) %>%
  arrange(z, tau) %>%
  pivot_wider(names_from = tau, values_from = c(T_entry, T_exit))

gt_table <- table_out %>%
  gt(rowname_col = "z") %>%
  tab_header(
    title = md("**Table: T_entry and T_exit for each (τ, z)**"),
    subtitle = "Computed from paths of l̂(T) under SNR = 40"
  ) %>%
  fmt_number(columns = everything(), decimals = 1) %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )

gtsave(gt_table, "T_entry_exit_table.png")
cat("Saved table as 'T_entry_exit_table.png'\n")

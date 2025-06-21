# ---------------------------------------------------------------------------
# Compute T(tau, z): minimum T such that l_hat(T') == 22 for all T' >= T
# Author: Said Maanan | Date: 2025-06-21
# ---------------------------------------------------------------------------

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

# --- Main Configuration ----------------------------------------------------

A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
seed <- 280
h <- 30
T_vals <- 4:h
# taus <- c(0.5, 0.75, 0.95)
taus <- c(0.5, 0.75, 0.95, 0.99)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)

# --- Core Loop: Compute T(tau, z) ------------------------------------------

results <- expand.grid(tau = taus, z = z_vals)
results$T_star <- NA

for (i in seq_len(nrow(results))) {
  tau <- results$tau[i]
  z <- results$z[i]
  
  set.seed(seed)
  signal <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
  b <- compute_noise_scale(A_true, snr, tau, z)
  noise <- rgen_asy_laplace(h, tau, z, b)
  I <- signal + noise
  
  l_path <- sapply(T_vals, function(T) estimate_params(I[1:T], tau, h)$l_hat)
  
  # Find smallest T such that all later estimates equal 22
  for (j in seq_along(T_vals)) {
    if (all(l_path[j:length(l_path)] == L_true)) {
      results$T_star[i] <- T_vals[j]
      break
    }
  }
}

# --- Display Output --------------------------------------------------------

print(results)

# Optional: Save to CSV
# write.csv(results, file = "T_star_table.csv", row.names = FALSE)
# cat("Saved table as 'T_star_table.csv'\n")























# install.packages("gt")
# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("webshot2")
# webshot2::install_phantomjs()  # only once, required for saving as image




library(gt)
library(dplyr)
library(tidyr)

# Load the results if not already in the environment
# results <- read.csv("T_star_table.csv")

# Pivot to wide format: rows = z, columns = tau
wide_table <- results %>%
  mutate(tau = paste0("tau = ", tau)) %>%
  pivot_wider(names_from = tau, values_from = T_star) %>%
  arrange(z)

# Generate gt table
gt_table <- wide_table %>%
  gt(rowname_col = "z") %>%
  tab_header(
    title = md("**Table: Minimal T(τ, z) such that l̂(T') = 22 ∀ T' ≥ T(τ,z)**"),
    subtitle = "Simulated under SNR = 1000, using asymmetric Laplace noise"
  ) %>%
  cols_label(
    `tau = 0.5` = md("**τ = 0.5**"),
    `tau = 0.75` = md("**τ = 0.75**"),
    `tau = 0.95` = md("**τ = 0.95**"),
    `tau = 0.99` = md("**τ = 0.99**")
  ) %>%
  fmt_number(columns = everything(), decimals = 2) %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )

# Save as PNG
gtsave(gt_table, "T_star_table.png")
cat("Saved table as T_star_table.png\n")

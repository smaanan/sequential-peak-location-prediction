# ---------------------------------------------------------------------------
# Entry/Exit Time Estimation under Varying tau, b, z using GAL Noise
# Author: Said Maanan (rewritten version)
# ---------------------------------------------------------------------------

# --- Libraries ---
library(dplyr)
library(tidyr)
library(gt)
library(ggplot2)
library(purrr)
library(tibble)
library(glue)

# --- Noise Generator and Helpers ---
rgen_gal <- function(n, tau, z, b) {
  u <- runif(n)
  ifelse(u < tau,
         -rexp(n, rate = (1 + z) / b),
         rexp(n, rate = (1 - z) / b))
}

gal_variance <- function(tau, z) {
  mu <- tau / (1 + z) - (1 - tau) / (1 - z)
  ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
  ex2 - mu^2
}

compute_noise_scale <- function(A, snr, tau, z) {
  sqrt((A / snr)^2 / gal_variance(tau, z))
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

find_stable_time <- function(l_path, T_vals, tol) {
  for (i in seq_along(T_vals)) {
    future_vals <- l_path[i:length(l_path)]
    if (max(future_vals) - min(future_vals) <= tol) {
      return(T_vals[i])
    }
  }
  return(NA)
}

# --- Main Grid Evaluation --------------------------------------------------
run_sim_grid <- function(A = 1, L = 22, S = 10, snr = 40,
                         taus = c(0.5, 0.75, 0.95),
                         z_vals = c(-0.9, -0.5, 0, 0.5, 0.9),
                         b_vals = c(0.25, 0.5, 1, 2),
                         T_range = 4:35,
                         threshold = 4, threshold_prime = 0,
                         h = 35, seed = 42) {

  results <- expand.grid(tau = taus, z = z_vals, b = b_vals) %>%
    mutate(data = pmap(list(tau, z, b), function(tau, z, b) {
      set.seed(seed)
      t_full <- 1:h
      signal <- A * exp(-((t_full - L)^2) / S^2)
      noise <- rgen_gal(h, tau, z, b)
      I <- signal + noise

      l_path <- sapply(T_range, function(T) estimate_params(I[1:T], tau, h)$l_hat)
      T_entry <- find_stable_time(l_path, T_range, threshold)
      T_exit  <- find_stable_time(l_path[T_range >= T_entry], T_range[T_range >= T_entry], threshold_prime)
      list(l_path = l_path, T_vals = T_range, T_entry = T_entry, T_exit = T_exit)
    }))
  
  return(results)
}



# ---- Prepare Data for Plotting ----
entry_exit_long <- run_sim_grid() %>%
  unnest_wider(data) %>%
  mutate(b = as.factor(b)) %>%
  pivot_longer(cols = c(T_entry, T_exit),
               names_to = "type", values_to = "T_val")

# ---- Heatmap Plot ----
ggplot(entry_exit_long, aes(x = factor(z), y = b, fill = T_val)) +
  geom_tile(color = "white") +
  facet_grid(type ~ tau, labeller = label_both) +
  scale_fill_viridis_c(option = "B", direction = -1, na.value = "grey90") +
  labs(
    title = "Entry and Exit Times under General Asymmetric Laplace Noise",
    x = "Asymmetry (z)", y = "Scale (b)", fill = "T",
    caption = "Each facet shows Entry (T_entry) or Exit (T_exit) for a given τ"
  ) +
  theme_minimal(base_size = 12)


# ---- Table Summary ----
summary_table <- run_sim_grid() %>%
  unnest_wider(data) %>%
  select(tau, z, b, T_entry, T_exit) %>%
  arrange(tau, z, b)

# ---- Generate gt Table ----
gt_table <- summary_table %>%
  gt() %>%
  tab_header(
    title = md("**Entry and Exit Times (T_entry, T_exit)**"),
    subtitle = "Across all combinations of τ, z, and b (SNR = 40)"
  ) %>%
  fmt_number(columns = c(T_entry, T_exit), decimals = 1) %>%
  cols_label(
    tau = md("**τ**"),
    z = md("**z**"),
    b = md("**b**"),
    T_entry = md("**T<sub>entry</sub>**"),
    T_exit = md("**T<sub>exit</sub>**")
  ) %>%
  tab_options(table.font.size = "small")

# ---- Save and Export ----
gtsave(gt_table, "table_entry_exit_tau_z_b.png")
cat("Saved 'table_entry_exit_tau_z_b.png'\n")

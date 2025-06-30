# ---------------------------------------------------------------------------
# Entry/Exit Times for Real Data using -log(f_{τ,z,b}) Objective
# ---------------------------------------------------------------------------

library(COVID19)
library(dplyr)
library(tidyr)
library(gt)
library(glue)

# --- Parameters ---
h <- 35
T_vals <- 4:h
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
taus <- c(0.25, 0.5, 0.75, 0.95)
threshold <- 3
threshold_prime <- 0
b <- 1  # fixed for real data

# --- Load China COVID-19 data ---
china <- covid19("China", level = 1, start = "2020-01-22", end = "2020-03-22", verbose = FALSE) %>%
  arrange(date) %>%
  mutate(new_cases = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_cases)) %>%
  mutate(t = row_number()) %>%
  filter(t <= 60)

I_full <- china$new_cases

# --- Asymmetric Laplace log-density ---
log_f_asy_laplace <- function(x, tau, z, b = 1) {
  log((tau * (1 - z) + (1 - tau) * (1 + z)) / 2) -
    (x < 0) * (1 + z) * abs(x) / b -
    (x >= 0) * (1 - z) * abs(x) / b - log(b)
}

# --- Weighted quantile ---
weighted_quantile <- function(x, w, tau) {
  o <- order(x)
  cum_w <- cumsum(w[o])
  cutoff <- tau * sum(w)
  x[o][min(which(cum_w >= cutoff))]
}

# --- Compute a ---
compute_a <- function(I, l, s, tau) {
  t <- seq_along(I)
  w <- exp(-((t - l)^2) / s^2)
  weighted_quantile(I, w, tau)
}

# --- Objective function ---
loglik_objective <- function(I, a, l, s, tau, z, b) {
  t <- seq_along(I)
  u <- I - a * exp(-((t - l)^2) / s^2)
  -mean(log_f_asy_laplace(u, tau, z, b))
}

# --- Find best s given l ---
search_best_s <- function(I, l, tau, z, b) {
  f <- function(s) {
    a <- compute_a(I, l, s, tau)
    loglik_objective(I, a, l, s, tau, z, b)
  }
  optimize(f, c(1, 20), tol = 1e-4)$minimum
}

# --- Estimate l_hat(T) ---
estimate_lhat <- function(I, tau, z, b, h) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau, z, b)
    a <- compute_a(I, l, s, tau)
    loglik_objective(I, a, l, s, tau, z, b)
  })
  which.min(losses)
}

# --- Stability functions ---
find_stable_time <- function(l_path, T_vals, tol) {
  for (i in seq_along(T_vals)) {
    if (max(l_path[i:length(l_path)]) - min(l_path[i:length(l_path)]) <= tol) {
      return(T_vals[i])
    }
  }
  return(NA)
}

# --- Main estimation loop ---
entry_exit_results <- data.frame()

for (tau in taus) {
  for (z in z_vals) {
    l_path <- sapply(T_vals, function(T) {
      I <- I_full[1:T]
      estimate_lhat(I, tau, z, b, h)
    })
    T_entry <- find_stable_time(l_path, T_vals, threshold)
    T_exit <- if (!is.na(T_entry)) {
      find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)
    } else NA
    
    entry_exit_results <- bind_rows(entry_exit_results,
                                    data.frame(tau = tau, z = z, T_entry = T_entry, T_exit = T_exit))
  }
}

# --- Pivot wide for display ---
table_wide <- entry_exit_results %>%
  mutate(tau = paste0("τ = ", tau)) %>%
  pivot_wider(names_from = tau, values_from = c(T_entry, T_exit)) %>%
  arrange(z)

# --- Create gt table ---
gt_table <- table_wide %>%
  gt(rowname_col = "z") %>%
  tab_header(
    title = md("**Table: Entry/Exit Times from -log(f<sub>τ,z,b</sub>) Objective (Real Data)**"),
    subtitle = md("China, January–March 2020 | Thresholds: ε = 3, ε' = 0")
  ) %>%
  fmt_number(columns = everything(), decimals = 1) %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )

# --- Export PNG and LaTeX ---
gtsave(gt_table, "table_entry_exit_realdata.png")
writeLines(as.character(as_latex(gt_table)), "table_entry_exit_realdata.tex")
cat("✓ Entry/Exit table saved as PNG and LaTeX.\n")

# ---------------------------------------------------------------------------
# Entry/Exit Times for Real Data using GAL Loss
# ---------------------------------------------------------------------------

library(COVID19)
library(dplyr)
library(tidyr)
library(gt)
library(glue)

# --- Parameters ---
h <- 35
T_vals <- 4:h
taus <- c(0.25, 0.5, 0.75, 0.95)
lambdas <- c(1, 0.5, 0.1)
threshold <- 3
threshold_prime <- 0
b <- 1  # fixed scale

# --- Load China COVID-19 data ---
china <- covid19("China", level = 1, start = "2020-01-22", end = "2020-03-22", verbose = FALSE) %>%
  arrange(date) %>%
  mutate(new_cases = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_cases)) %>%
  mutate(t = row_number()) %>%
  filter(t <= 60)

I_full <- china$new_cases

# --- GAL loss function ---
# gal_loss <- function(u, z, lambda, b = 1) {
#   ((z - lambda * u / b)^2) / (2 * lambda) + abs(u) / lambda
# }

gal_loss <- function(u, z, lambda, b = 1) {
  ((z - (lambda * u / b))^2) / (2 * lambda) + abs(u) / b
}

# --- Weighted quantile ---
weighted_quantile <- function(x, w, tau) {
  o <- order(x)
  cum_w <- cumsum(w[o])
  cutoff <- tau * sum(w)
  x[o][min(which(cum_w >= cutoff))]
}

# --- Compute a(l,s) ---
compute_a <- function(I, l, s, tau) {
  t <- seq_along(I)
  w <- exp(-((t - l)^2) / s^2)
  weighted_quantile(I, w, tau)
}

# --- GAL-based objective ---
gal_objective <- function(I, a, l, s, z, lambda, b = 1) {
  t <- seq_along(I)
  u <- I - a * exp(-((t - l)^2) / s^2)
  mean(gal_loss(u, z, lambda, b))
}

# --- Find best s given l ---
search_best_s_gal <- function(I, l, tau, lambda, b) {
  z <- 1 - 2 * tau
  f <- function(s) {
    a <- compute_a(I, l, s, tau)
    gal_objective(I, a, l, s, z, lambda, b)
  }
  optimize(f, c(1, 20), tol = 1e-4)$minimum
}

# --- Estimate l_hat(T) ---
estimate_lhat_gal <- function(I, tau, lambda, b, h) {
  z <- 1 - 2 * tau
  losses <- sapply(1:h, function(l) {
    s <- search_best_s_gal(I, l, tau, lambda, b)
    a <- compute_a(I, l, s, tau)
    gal_objective(I, a, l, s, z, lambda, b)
  })
  which.min(losses)
}

# --- Stability check ---
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
  for (lambda in lambdas) {
    l_path <- sapply(T_vals, function(T) {
      I <- I_full[1:T]
      estimate_lhat_gal(I, tau, lambda, b, h)
    })
    
    T_entry <- find_stable_time(l_path, T_vals, threshold)
    T_exit <- if (!is.na(T_entry)) {
      find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)
    } else NA
    
    entry_exit_results <- bind_rows(entry_exit_results,
                                    data.frame(tau = tau, lambda = lambda, T_entry = T_entry, T_exit = T_exit)
    )
  }
}

# --- Format table for display ---
table_wide <- entry_exit_results %>%
  mutate(tau = paste0("τ = ", tau)) %>%
  pivot_wider(names_from = tau, values_from = c(T_entry, T_exit)) %>%
  arrange(lambda)

# --- Create gt table ---
gt_table <- table_wide %>%
  gt(rowname_col = "lambda") %>%
  tab_header(
    title = md("**Table: Entry/Exit Times from GAL Loss Objective (Real Data)**"),
    subtitle = md("China, January–March 2020 | z = 1 – 2τ | Thresholds: ε = 3, ε' = 0")
  ) %>%
  fmt_number(columns = everything(), decimals = 1) %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )

# --- Save PNG and LaTeX ---
gtsave(gt_table, "china_entry_exit_gal.png")
writeLines(as.character(as_latex(gt_table)), "china_entry_exit_gal.tex")
cat("✓ China GAL-based entry/exit table saved as PNG and LaTeX.\n")

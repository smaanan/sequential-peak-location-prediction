# ---------------------------------------------------------------------------
# Entry/Exit Times for Germany using GAL Loss
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(gt)
library(readr)
library(glue)
library(zoo)

# --- Parameters ---
threshold <- 15
threshold_prime <- 2
h <- 160
T_vals <- 40:75
taus <- c(0.25, 0.5, 0.75, 0.95)
lambdas <- c(1, 0.5, 0.1)
b <- 1  # fixed scale

# --- Load Germany Data from ECDC ---
get_de_cases_ecdc <- function(start_date = "2020-01-28", end_date = "2020-07-05") {
  url <- "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
  read_csv(url, col_types = cols()) %>%
    filter(countriesAndTerritories == "Germany") %>%
    mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
    filter(date >= as.Date(start_date) & date <= as.Date(end_date)) %>%
    arrange(date) %>%
    select(date, daily_cases = cases) %>%
    mutate(daily_cases = ifelse(daily_cases < 0, NA, daily_cases),
           daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE))
}

de <- get_de_cases_ecdc()
I_full <- de$daily_cases

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

# --- Optimize over s given l ---
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

# --- Main loop ---
results <- data.frame()

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
    
    results <- rbind(results, data.frame(tau = tau, lambda = lambda, T_entry = T_entry, T_exit = T_exit))
  }
}

# --- Format and display table ---
table_germany <- results %>%
  mutate(tau = paste0("τ = ", tau)) %>%
  arrange(lambda, tau) %>%
  pivot_wider(names_from = tau, values_from = c(T_entry, T_exit)) %>%
  gt(rowname_col = "lambda") %>%
  tab_header(
    title = md("**Table: Entry/Exit Times using GAL Loss (Germany, Real Data)**"),
    subtitle = md("ECDC dataset | z = 1 – 2τ | Thresholds: ε = 5, ε' = 1")
  ) %>%
  fmt_number(columns = everything(), decimals = 1) %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )

# --- Save PNG and LaTeX ---
gtsave(table_germany, "germany_entry_exit_gal.png")
writeLines(as.character(as_latex(table_germany)), "germany_entry_exit_gal.tex")
cat("✓ Germany GAL-based entry/exit table saved.\n")

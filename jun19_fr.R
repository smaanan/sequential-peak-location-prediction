# --- Libraries ---
library(dplyr)
library(tidyr)
library(gt)
library(readr)
library(glue)
library(zoo)

# --- France Parameters ---
threshold <- 5
threshold_prime <- 1
h <- 120
T_vals <- 50:89
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
taus <- c(0.25, 0.5, 0.75, 0.95)

# --- Load France Data from ECDC ---
get_france_cases_ecdc <- function(start_date = "2020-01-15", end_date = "2020-05-31") {
  url <- "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
  read_csv(url, col_types = cols()) %>%
    filter(countriesAndTerritories == "France") %>%
    mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
    filter(date >= as.Date(start_date) & date <= as.Date(end_date)) %>%
    arrange(date) %>%
    select(date, daily_cases = cases) %>%
    mutate(daily_cases = ifelse(daily_cases < 0, NA, daily_cases),
           daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE))
}


france <- get_france_cases_ecdc()
I_full <- france$daily_cases

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

# --- Compute a(l,s) ---
compute_a <- function(I, l, s, tau) {
  t <- seq_along(I)
  w <- exp(-((t - l)^2) / s^2)
  weighted_quantile(I, w, tau)
}

# --- Log-likelihood-based loss ---
loglik_objective <- function(I, a, l, s, tau, z, b) {
  t <- seq_along(I)
  u <- I - a * exp(-((t - l)^2) / s^2)
  -mean(log_f_asy_laplace(u, tau, z, b))
}

# --- Optimize over s given l ---
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

# --- Stability check ---
find_stable_time <- function(l_path, T_vals, tol) {
  for (i in seq_along(T_vals)) {
    if (max(l_path[i:length(l_path)]) - min(l_path[i:length(l_path)]) <= tol) {
      return(T_vals[i])
    }
  }
  return(NA)
}

# --- Loop over tau and z combinations ---
results <- data.frame()

for (tau in taus) {
  for (z in z_vals) {
    b <- 1  # fixed scale as in previous script
    l_path <- sapply(T_vals, function(T) {
      I <- I_full[1:T]
      estimate_lhat(I, tau, z, b, h)
    })
    
    T_entry <- find_stable_time(l_path, T_vals, threshold)
    T_exit <- if (!is.na(T_entry)) {
      find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)
    } else NA
    
    results <- rbind(results, data.frame(tau = tau, z = z, T_entry = T_entry, T_exit = T_exit))
  }
}

# --- Format table ---
table_france <- results %>%
  mutate(tau = paste0("τ = ", tau)) %>%
  arrange(z, tau) %>%
  pivot_wider(names_from = tau, values_from = c(T_entry, T_exit)) %>%
  gt(rowname_col = "z") %>%
  tab_header(
    title = md("**Table: Entry/Exit Times using -log(f<sub>τ,z,b</sub>) Loss (France, Real Data)**")
  ) %>%
  fmt_number(columns = everything(), decimals = 1) %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = 14,
    data_row.padding = px(3)
  )

# --- Save PNG and LaTeX ---
gtsave(table_france, "table_entry_exit_france.png")
writeLines(as.character(as_latex(table_france)), "table_entry_exit_france.tex")

cat("✓ Entry/Exit table saved for France.\n")

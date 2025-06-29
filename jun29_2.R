# ---------------------------------------------------------------------------
# Section 7.2 — Stability Analysis Using Negative Log-Likelihood (Real Data)
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(gt)
library(tibble)
library(glue)
library(readr)
library(zoo)

# ------------------- GAL Log-Likelihood Loss -------------------------
neg_log_gal_density <- function(e, tau, z, b) {
  ifelse(e < 0,
         -log(tau) - log(1 + z) + abs((1 + z) * e) / b,
         -log(1 - tau) - log(1 - z) + abs((1 - z) * e) / b) + log(b)
}

gal_loss <- function(I, a, l, s, tau, z, b) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  e <- I - fitted
  mean(neg_log_gal_density(e, tau, z, b))
}

# ------------- Estimation Functions (Quantile-based Peak) --------------
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

search_best_s <- function(I, l, tau) {
  f <- function(s) {
    a <- compute_a(I, l, s, tau)
    fitted <- a * exp(-((1:length(I)) - l)^2 / s^2)
    u <- I - fitted
    mean(u * (tau - (u < 0)))
  }
  optimize(f, c(1, 40), tol = 1e-4)$minimum
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

asy_laplace_variance <- function(tau, z) {
  mu <- tau / (1 + z) - (1 - tau) / (1 - z)
  ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
  ex2 - mu^2
}

compute_noise_scale <- function(A, snr, tau, z) {
  sqrt((A / snr)^2 / asy_laplace_variance(tau, z))
}

# ---------------------- Entry/Exit Detection --------------------------
find_stable_time <- function(l_path, T_vals, tol) {
  for (i in seq_along(T_vals)) {
    future_vals <- l_path[i:length(l_path)]
    if (max(future_vals) - min(future_vals) <= tol) {
      return(T_vals[i])
    }
  }
  return(NA)
}

# --------------------- Main Analysis Function --------------------------
compute_real_entry_exit_table <- function(I_full, country_name, T_vals, snr = 1000,
                                          threshold, threshold_prime,
                                          taus = c(0.25, 0.5, 0.75, 0.95),
                                          z_vals = c(-0.9, -0.5, 0, 0.5, 0.9)) {
  results <- data.frame()
  A <- max(I_full)  # approximate scale
  h <- length(I_full)
  
  for (tau in taus) {
    for (z in z_vals) {
      b <- compute_noise_scale(A, snr, tau, z)
      l_path <- sapply(T_vals, function(T) {
        estimate_params(I_full[1:T], tau, h)$l_hat
      })
      T_entry <- find_stable_time(l_path, T_vals, threshold)
      T_exit <- if (!is.na(T_entry)) {
        find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)
      } else NA
      results <- rbind(results, data.frame(country = country_name, tau, z, T_entry, T_exit))
    }
  }
  
  # Format Table
  results %>%
    mutate(tau = paste0("τ = ", tau)) %>%
    pivot_wider(names_from = tau, values_from = c(T_entry, T_exit)) %>%
    arrange(z) %>%
    gt(rowname_col = "z") %>%
    tab_header(title = md(glue("**Table: Entry/Exit using -log(f<sub>τ,z,b</sub>) Loss ({country_name}, SNR = {snr})**"))) %>%
    fmt_number(columns = everything(), decimals = 1) %>%
    tab_options(table.font.size = "small", heading.title.font.size = 14, data_row.padding = px(3))
}


# --- China Parameters ---
threshold <- 3
threshold_prime <- 0
h <- 35
T_vals <- 4:h

# --- Load Data: JHU Source (60 days) ---
china <- covid19("China", level = 1, start = "2020-01-22", end = "2020-03-22", verbose = FALSE) %>%
  arrange(date) %>%
  mutate(new_cases = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_cases)) %>%
  mutate(t = row_number()) %>%
  filter(t <= 60)

I_full <- china$new_cases

# --- Compute and Save Entry/Exit Table ---
china_table <- compute_real_entry_exit_table(
  I_full = I_full,
  country_name = "China",
  T_vals = T_vals,
  snr = 40,
  threshold = threshold,
  threshold_prime = threshold_prime
)

gtsave(china_table, "china_table7_entry_exit_snr40.png")
writeLines(as.character(as_latex(china_table)), "china_table7_entry_exit_snr40.tex")





# --- France Parameters ---
threshold <- 5
threshold_prime <- 1
h <- 120
T_vals <- 50:89

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

# --- Compute and Save Entry/Exit Table ---
france_table <- compute_real_entry_exit_table(
  I_full = I_full,
  country_name = "France",
  T_vals = T_vals,
  snr = 40,
  threshold = threshold,
  threshold_prime = threshold_prime
)

gtsave(france_table, "france_table7_entry_exit_snr40.png")
writeLines(as.character(as_latex(france_table)), "france_table7_entry_exit_snr40.tex")




# --- Germany Parameters ---
threshold <- 15
threshold_prime <- 2
h <- 160
T_vals <- 40:75

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

# --- Compute and Save Entry/Exit Table ---
germany_table <- compute_real_entry_exit_table(
  I_full = I_full,
  country_name = "Germany",
  T_vals = T_vals,
  snr = 40,
  threshold = threshold,
  threshold_prime = threshold_prime
)

gtsave(germany_table, "germany_table7_entry_exit_snr40.png")
writeLines(as.character(as_latex(germany_table)), "germany_table7_entry_exit_snr40.tex")

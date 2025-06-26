# ---------------------------------------------------------------------
# Table 8-style Output Using Real COVID-19 Data from France (ECDC)
# ---------------------------------------------------------------------

# --- Libraries ---
library(dplyr)
library(tidyr)
library(gt)
library(tibble)
library(glue)
library(readr)

# --- Parameters ---
taus <- c(0.5, 0.75, 0.95)
threshold <- 5
threshold_prime <- 1
h <- 120
T_vals <- 50:89  # From Table 10 in the paper

# --- Load ECDC France COVID-19 Data (Jan 15 to May 31) ---
get_france_cases_ecdc <- function(start_date = "2020-01-15", end_date = "2020-05-31") {
  url <- "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
  
  data <- read_csv(url, col_types = cols()) %>%
    filter(countriesAndTerritories == "France") %>%
    mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
    filter(date >= as.Date(start_date) & date <= as.Date(end_date)) %>%
    arrange(date) %>%
    select(date, daily_cases = cases) %>%
    mutate(daily_cases = ifelse(daily_cases < 0, NA, daily_cases),
           daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE),
           log_cases = log(pmax(daily_cases, 1)))
  
  return(data)
}

france <- get_france_cases_ecdc()
# I_full <- france$log_cases  # Log-transformed daily new cases
I_full <- france$daily_cases

# --- Estimation Functions ---
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
  optimize(f, c(1, 40), tol = 1e-4)$minimum
}

estimate_params <- function(I, tau, h) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau)
    a <- compute_a(I, l, s, tau)
    objective_function(I, a, l, s, tau)
  })
  l_hat <- which.min(losses)
  list(l_hat = l_hat)
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

# --- Main Loop Over τ Values ---
for (tau in taus) {
  l_path <- sapply(T_vals, function(T) {
    I <- I_full[1:T]
    estimate_params(I, tau, h)$l_hat
  })
  
  T_entry <- find_stable_time(l_path, T_vals, threshold)
  T_exit <- find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)
  
  # Table structure
  table_data <- tibble(
    T_1 = T_vals[1:20],
    lhat_1 = l_path[1:20],
    T_2 = T_vals[21:40],
    lhat_2 = l_path[21:40]
  )
  
  transposed_matrix <- t(table_data)
  colnames(transposed_matrix) <- paste0("V", seq_len(ncol(transposed_matrix)))
  transposed_data <- as.data.frame(transposed_matrix) %>%
    rownames_to_column(var = "Label") %>%
    as_tibble()
  
  transposed_data$Label <- recode(transposed_data$Label,
                                  "T_1" = "T", "T_2" = "T",
                                  "lhat_1" = "l̂(T)", "lhat_2" = "l̂(T)")
  
  subtitle_text <- glue(
    "COVID-19 in France, τ = {tau}. ",
    "Entry at **T_entry = {T_entry}** (ε = {threshold}), ",
    "exit at **T_exit = {T_exit}** (ε' = {threshold_prime})."
  )
  
  gt_table <- transposed_data %>%
    gt() %>%
    tab_header(
      title = md("**Table: Sequence T → l̂(T)**"),
      subtitle = md(subtitle_text)
    ) %>%
    cols_label(.list = setNames(rep("", ncol(transposed_data)), names(transposed_data))) %>%
    tab_options(
      column_labels.hidden = TRUE,
      table.font.size = "small",
      heading.title.font.size = 14,
      heading.subtitle.font.size = 11,
      data_row.padding = px(3)
    )
  
  filename <- glue("france_table8_lhat_path_tau{tau}.png")
  gtsave(gt_table, filename)
  cat(glue("Saved: {filename}\n"))
}

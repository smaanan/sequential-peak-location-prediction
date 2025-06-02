# -------------------------------
# 0. Libraries
# -------------------------------
library(dplyr)
library(readr)
library(lubridate)
library(zoo)
library(quantreg)
library(pracma)  # for hilbert()

# -------------------------------
# 1. Load and prepare data (ECDC)
# -------------------------------
data <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv",
                 na.strings = "", fileEncoding = "UTF-8-BOM")

france_cases <- data %>%
  filter(countriesAndTerritories == "France") %>%
  mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
  filter(date >= as.Date("2020-01-15") & date <= as.Date("2020-05-31")) %>%
  arrange(date) %>%
  select(date, daily_cases = cases) %>%
  mutate(
    daily_cases = ifelse(daily_cases < 0, NA, daily_cases),
    daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE),
    log_cases = log(pmax(daily_cases, 1)),
    t = row_number(),
    t_sq = t^2
  )

# Define data sets
data_T70 <- france_cases %>% filter(t <= 70)
data_T77 <- france_cases %>% filter(t <= 77)

# -------------------------------
# 2. Estimation function using Hilbert-based covariance
# -------------------------------
estimate_poly_hilbert <- function(df, t_peak = 78, tau = 0.95, f0 = 1/2) {
  T_len <- nrow(df)
  
  # Step 1: Quantile regression
  model <- rq(log_cases ~ t + t_sq, tau = tau, data = df, method = "fn")
  beta <- coef(model)
  
  # Step 2: Hilbert-based covariance matrix
  scaling <- c(T_len^(1/2), T_len^(3/2), T_len^(5/2))
  H_inv <- solve(hilbert(3))
  omega2 <- (tau * (1 - tau)) / f0^2
  Sigma_beta <- omega2 * diag(1 / scaling) %*% H_inv %*% diag(1 / scaling)
  
  # Step 3: Transform beta to (log a, l, s²)
  beta0 <- beta[1]
  beta1 <- beta[2]
  beta2 <- beta[3]
  l_hat <- -beta1 / (2 * beta2)
  s2_hat <- -1 / beta2
  log_a_hat <- beta0 + l_hat^2 / s2_hat
  
  # Step 4: Delta method for log Î(t_peak)
  X_peak <- matrix(c(1, t_peak, t_peak^2), nrow = 1)
  log_peak_pred <- as.numeric(X_peak %*% beta)
  se_peak <- sqrt(as.numeric(X_peak %*% Sigma_beta %*% t(X_peak)))
  ci <- c(log_peak_pred - 1.96 * se_peak,
          log_peak_pred + 1.96 * se_peak)
  
  return(list(
    beta = beta,
    log_a = log_a_hat,
    l = l_hat,
    s2 = s2_hat,
    log_peak_pred = log_peak_pred,
    ci = ci,
    se = se_peak
  ))
}

# -------------------------------
# 3. Run and print results
# -------------------------------
res_T70 <- estimate_poly_hilbert(data_T70, t_peak = 78)
res_T77 <- estimate_poly_hilbert(data_T77, t_peak = 78)

cat("=== Résultats pour T = 70 ===\n")
cat(sprintf("log(a) = %.3f → a ≈ %.0f\n", res_T70$log_a, exp(res_T70$log_a)))
cat(sprintf("l = %.2f\n", res_T70$l))
cat(sprintf("s² = %.2f → s ≈ %.1f\n", res_T70$s2, sqrt(res_T70$s2)))
cat(sprintf("log(Î(78)) = %.3f → Î(78) ≈ %.0f\n", res_T70$log_peak_pred, exp(res_T70$log_peak_pred)))
cat(sprintf("95%% CI (log-scale): [%.3f, %.3f]\n\n", res_T70$ci[1], res_T70$ci[2]))

cat("=== Résultats pour T = 77 ===\n")
cat(sprintf("log(a) = %.3f → a ≈ %.0f\n", res_T77$log_a, exp(res_T77$log_a)))
cat(sprintf("l = %.2f\n", res_T77$l))
cat(sprintf("s² = %.2f → s ≈ %.1f\n", res_T77$s2, sqrt(res_T77$s2)))
cat(sprintf("log(Î(78)) = %.3f → Î(78) ≈ %.0f\n", res_T77$log_peak_pred, exp(res_T77$log_peak_pred)))
cat(sprintf("95%% CI (log-scale): [%.3f, %.3f]\n", res_T77$ci[1], res_T77$ci[2]))

res_T70

res_T77

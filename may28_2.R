# ----------------------------------------
# Hybrid Method: LAD + Hilbert-based CI
# ----------------------------------------
library(dplyr)
library(lubridate)
library(readr)
library(zoo)
library(pracma)  # for hilbert()

# ----------------------------------------
# 1. Load and prepare data (ECDC)
# ----------------------------------------
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
    t = row_number()
  )

# Restrict to T = 70 and T = 77
data_T70 <- france_cases %>% filter(t <= 70)
data_T77 <- france_cases %>% filter(t <= 77)

# ----------------------------------------
# 2. LAD Estimation in (log_a, l, s2)
# ----------------------------------------
lad_objective <- function(theta, t_vec, log_I, tau = 0.95) {
  loga <- theta[1]
  l <- theta[2]
  s2 <- theta[3]
  if (s2 <= 0) return(1e10)
  y_pred <- loga - (t_vec - l)^2 / s2
  rho <- ifelse(log_I >= y_pred,
                tau * (log_I - y_pred),
                (tau - 1) * (log_I - y_pred))
  return(sum(rho))
}

estimate_hybrid <- function(df, t_peak = 78, tau = 0.95, f0 = 1/2) {
  t_vec <- df$t
  log_I <- df$log_cases
  T_len <- length(t_vec)
  
  # Step 1: LAD estimation (log_a, l, s2)
  theta0 <- c(log(max(df$daily_cases)), median(t_vec), 50)
  fit <- optim(par = theta0,
               fn = lad_objective,
               method = "L-BFGS-B",
               lower = c(log(1), 1, 1e-2),
               upper = c(log(1e6), max(t_vec)+20, 1e4),
               t_vec = t_vec,
               log_I = log_I,
               tau = tau,
               control = list(maxit = 1000))
  theta <- fit$par
  names(theta) <- c("log_a", "l", "s2")
  
  # Step 2: Convert to beta0, beta1, beta2
  beta2 <- -1 / theta["s2"]
  beta1 <- 2 * theta["l"] / theta["s2"]
  beta0 <- theta["log_a"] - theta["l"]^2 / theta["s2"]
  beta <- c(beta0, beta1, beta2)
  
  # Step 3: Hilbert-based covariance matrix for beta
  scaling <- c(T_len^(1/2), T_len^(3/2), T_len^(5/2))
  H_inv <- solve(hilbert(3))
  omega2 <- (tau * (1 - tau)) / f0^2
  Sigma_beta <- omega2 * diag(1 / scaling) %*% H_inv %*% diag(1 / scaling)
  
  # Step 4: CI for log Î(78)
  t_peak_sq <- t_peak^2
  X <- c(1, t_peak, t_peak_sq)
  log_I_pred <- sum(beta * X)
  var_log_I <- X %*% Sigma_beta %*% X
  se <- sqrt(as.numeric(var_log_I))
  ci <- c(log_I_pred - 1.96 * se, log_I_pred + 1.96 * se)
  
  return(list(
    theta = theta,
    beta = beta,
    log_peak_pred = log_I_pred,
    ci = ci,
    se = se
  ))
}

# ----------------------------------------
# 3. Run and print results
# ----------------------------------------
res_T70 <- estimate_hybrid(data_T70, t_peak = 78)
res_T77 <- estimate_hybrid(data_T77, t_peak = 78)

cat("=== Résultats pour T = 70 ===\n")
print(res_T70$theta)
cat(sprintf("log(Î(78)) = %.3f → Î(78) ≈ %.0f\n", res_T70$log_peak_pred, exp(res_T70$log_peak_pred)))
cat(sprintf("95%% CI (log): [%.3f, %.3f]\n\n", res_T70$ci[1], res_T70$ci[2]))

cat("=== Résultats pour T = 77 ===\n")
print(res_T77$theta)
cat(sprintf("log(Î(78)) = %.3f → Î(78) ≈ %.0f\n", res_T77$log_peak_pred, exp(res_T77$log_peak_pred)))
cat(sprintf("95%% CI (log): [%.3f, %.3f]\n", res_T77$ci[1], res_T77$ci[2]))

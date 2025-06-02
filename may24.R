# -------------------------------
# 1. Load and prepare data (ECDC)
# -------------------------------
library(dplyr)
library(lubridate)
library(readr)

# Load data from ECDC
data <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", 
                 na.strings = "", fileEncoding = "UTF-8-BOM")

# Extract daily new cases in France from Jan 15 to May 31, 2020
france_cases <- data %>%
  filter(countriesAndTerritories == "France") %>%
  mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
  filter(date >= as.Date("2020-01-15") & date <= as.Date("2020-05-31")) %>%
  arrange(date) %>%
  select(date, daily_cases = cases) %>%
  mutate(
    daily_cases = ifelse(daily_cases < 0, NA, daily_cases),  # remove negative artifacts
    daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE),  # fill NA downward
    log_cases = log(pmax(daily_cases, 1)),                   # avoid log(0)
    t = row_number()
  )

# Restrict to T = 70 and T = 77
data_T70 <- france_cases %>% filter(t <= 70)
data_T77 <- france_cases %>% filter(t <= 77)






# -------------------------------
# 2. Define LAD objective function
# -------------------------------
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

# -------------------------------
# 3. Estimation function
# -------------------------------
estimate_peak_model <- function(df, t_peak = 78, tau = 0.95) {
  t_vec <- df$t
  log_I <- df$log_cases
  
  # Initial values
  theta0 <- c(log(max(df$daily_cases)), median(t_vec), 50)
  
  # Optimization
  fit <- optim(par = theta0,
               fn = lad_objective,
               method = "L-BFGS-B",
               lower = c(log(1), 1, 1e-2),
               upper = c(log(1e6), max(t_vec)+20, 1e4),
               t_vec = t_vec,
               log_I = log_I,
               tau = tau,
               control = list(maxit = 1000))
  
  est <- fit$par
  names(est) <- c("log_a", "l", "s2")
  
  # Peak prediction
  log_peak_pred <- est["log_a"] - (t_peak - est["l"])^2 / est["s2"]
  
  # Numerical gradient for delta method
  grad_g <- c(1,
              2 * (t_peak - est["l"]) / est["s2"],
              (t_peak - est["l"])^2 / est["s2"]^2)
  
  # Approximate Hessian of the loss (to estimate parameter variance)
  library(numDeriv)
  hess <- hessian(func = lad_objective,
                  x = est,
                  t_vec = t_vec,
                  log_I = log_I,
                  tau = tau)
  
  cov_matrix <- tryCatch({
    solve(hess)
  }, error = function(e) {
    message("Hessian not invertible, returning NA.")
    matrix(NA, 3, 3)
  })
  
  # Delta method standard error
  se_peak <- sqrt(as.numeric(t(grad_g) %*% cov_matrix %*% grad_g))
  ci_lower <- log_peak_pred - 1.96 * se_peak
  ci_upper <- log_peak_pred + 1.96 * se_peak
  
  return(list(
    est = est,
    log_peak_pred = log_peak_pred,
    ci = c(lower = ci_lower, upper = ci_upper),
    se = se_peak
  ))
}

# -------------------------------
# 4. Run estimation for T=70 and T=77
# -------------------------------
set.seed(42)
result_T70 <- estimate_peak_model(data_T70, t_peak = 78)
result_T77 <- estimate_peak_model(data_T77, t_peak = 78)

# -------------------------------
# 5. Print results
# -------------------------------
cat("Results for T = 70:\n")
print(result_T70$est)
cat(sprintf("Predicted log-peak at day 78: %.3f\n", result_T70$log_peak_pred))
cat(sprintf("95%% CI: [%.3f, %.3f]\n\n", result_T70$ci[1], result_T70$ci[2]))

cat("Results for T = 77:\n")
print(result_T77$est)
cat(sprintf("Predicted log-peak at day 78: %.3f\n", result_T77$log_peak_pred))
cat(sprintf("95%% CI: [%.3f, %.3f]\n", result_T77$ci[1], result_T77$ci[2]))

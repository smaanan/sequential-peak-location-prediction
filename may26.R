# -------------------------------
# 0. Libraries
# -------------------------------
library(dplyr)
library(readr)
library(lubridate)
library(zoo)
library(quantreg)

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
# 2. Estimation function
# -------------------------------
estimate_poly_model <- function(df, t_peak = 78, tau = 0.95) {
  model <- rq(log_cases ~ t + t_sq, tau = tau, data = df, method = "fn")
  model_sum <- summary(model, se = "ker", cov = TRUE)  # get VCV matrix
  
  beta <- coef(model)
  cov_beta <- model_sum$cov  # directly use this
  
  # Recover parameters from beta
  beta0 <- beta[1]
  beta1 <- beta[2]
  beta2 <- beta[3]
  
  l_hat <- -beta1 / (2 * beta2)
  s2_hat <- -1 / beta2
  log_a_hat <- beta0 + l_hat^2 / s2_hat
  
  # Predict log I(78)
  X_peak <- matrix(c(1, t_peak, t_peak^2), nrow = 1)
  pred_log_peak <- X_peak %*% beta
  pred_se <- sqrt(X_peak %*% cov_beta %*% t(X_peak))
  ci <- c(pred_log_peak - 1.96 * pred_se,
          pred_log_peak + 1.96 * pred_se)
  
  return(list(
    beta = beta,
    log_a = log_a_hat,
    l = l_hat,
    s2 = s2_hat,
    log_peak_pred = as.numeric(pred_log_peak),
    ci = as.numeric(ci),
    se = as.numeric(pred_se)
  ))
}

# -------------------------------
# 3. Run and print results
# -------------------------------
res_T70 <- estimate_poly_model(data_T70, t_peak = 78)
res_T77 <- estimate_poly_model(data_T77, t_peak = 78)

cat("Results for T = 70:\n")
print(res_T70[c("log_a", "l", "s2")])
cat(sprintf("Predicted log I(78): %.3f\n", res_T70$log_peak_pred))
cat(sprintf("95%% CI: [%.3f, %.3f]\n\n", res_T70$ci[1], res_T70$ci[2]))

cat("Results for T = 77:\n")
print(res_T77[c("log_a", "l", "s2")])
cat(sprintf("Predicted log I(78): %.3f\n", res_T77$log_peak_pred))
cat(sprintf("95%% CI: [%.3f, %.3f]\n", res_T77$ci[1], res_T77$ci[2]))

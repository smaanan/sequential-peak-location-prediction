# ----------------------------------------
# Revised Covariance Comparison Script
# ----------------------------------------
library(dplyr)
library(lubridate)
library(readr)
library(zoo)
library(numDeriv)

# Load data
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

data_T70 <- france_cases %>% filter(t <= 70)
t_vec <- data_T70$t
log_I <- data_T70$log_cases
T_len <- nrow(data_T70)

# LAD objective
lad_objective <- function(theta, t_vec, log_I, tau = 0.95) {
  loga <- theta[1]; l <- theta[2]; s2 <- theta[3]
  if (s2 <= 0) return(1e10)
  y_pred <- loga - (t_vec - l)^2 / s2
  rho <- ifelse(log_I >= y_pred,
                tau * (log_I - y_pred),
                (tau - 1) * (log_I - y_pred))
  return(sum(rho))
}

# Step 1: Estimate theta and compute Hessian
theta0 <- c(log(max(data_T70$daily_cases)), median(t_vec), 50)

fit <- optim(par = theta0,
             fn = lad_objective,
             method = "L-BFGS-B",
             lower = c(log(1), 1, 1e-2),
             upper = c(log(1e6), max(t_vec)+20, 1e4),
             t_vec = t_vec,
             log_I = log_I,
             tau = 0.95,
             control = list(maxit = 1000))

theta_hat <- fit$par
names(theta_hat) <- c("log_a", "l", "s2")

# Compute Hessian-based covariance for theta
hess <- hessian(f = lad_objective,
                x = theta_hat,
                t_vec = t_vec,
                log_I = log_I,
                tau = 0.95)
cov_theta <- solve(hess)

# Step 2: Transform to polynomial basis (beta)
log_a <- theta_hat["log_a"]; l_val <- theta_hat["l"]; s2_val <- theta_hat["s2"]

# Polynomial coefficients
beta0 <- log_a - l_val^2 / s2_val
beta1 <- 2 * l_val / s2_val
beta2 <- -1 / s2_val
beta <- c(beta0, beta1, beta2)

# Jacobian matrix (dbeta/dtheta)
J <- matrix(c(
  1,                 -2*l_val/s2_val,        l_val^2/s2_val^2,
  0,                  2/s2_val,             -2*l_val/s2_val^2,
  0,                  0,                     1/s2_val^2
), nrow = 3, byrow = TRUE)

# Delta method covariance for beta
cov_beta_delta <- J %*% cov_theta %*% t(J)

# Step 3: Hilbert-based covariance (using paper's method)
# (a) Construct exact scaling matrix
scaling <- c(
  sqrt(T_len),          # sum(1) = T
  sqrt(sum(t_vec^2)),   # sum(t^2)
  sqrt(sum(t_vec^4))    # sum(t^4)
)

# (b) Build Hilbert matrix H[i,j] = 1/(i+j+1) for i,j=0,1,2
H3 <- hilbert(3)
# H3 <- matrix(c(
#   1,   1/2, 1/3,
#   1/2, 1/3, 1/4,
#   1/3, 1/4, 1/5
# ), nrow = 3, byrow = TRUE)
H_inv <- solve(H3)

# (c) Estimate f(0) from residuals
y_pred <- theta_hat[1] - (t_vec - theta_hat[2])^2 / theta_hat[3]
residuals <- log_I - y_pred
f0 <- density(residuals, from = 0, to = 0, n = 1, kernel = "gaussian")$y

# (d) Compute Hilbert covariance
tau <- 0.95
omega2 <- (tau * (1 - tau)) / (f0^2)
cov_beta_hilbert <- omega2 * diag(1/scaling) %*% H_inv %*% diag(1/scaling)

# Step 4: Compare matrices
cat("=== Parameter Estimates ===\n")
print(round(theta_hat, 4))
cat(sprintf("\nPolynomial coefficients: beta0=%.4f, beta1=%.4f, beta2=%.6f\n", 
            beta0, beta1, beta2))
cat(sprintf("f(0) estimate: %.4f\n", f0))

cat("\n=== Delta Method Covariance (Hessian-based) ===\n")
print(round(cov_beta_delta, 8))

cat("\n=== Hilbert-Based Covariance (Paper) ===\n")
print(round(cov_beta_hilbert, 8))

cat("\n=== Relative Difference (Delta - Hilbert)/Hilbert ===\n")
print(round((cov_beta_delta - cov_beta_hilbert)/cov_beta_hilbert, 4))

cat("\n=== Diagonal Elements Comparison ===\n")
diag_comparison <- data.frame(
  Parameter = c("beta0", "beta1", "beta2"),
  Delta_Var = diag(cov_beta_delta),
  Hilbert_Var = diag(cov_beta_hilbert),
  Ratio = diag(cov_beta_delta)/diag(cov_beta_hilbert)
)
print(diag_comparison)



# Find log(I(78)) from actual data
obs_log_I_78 <- france_cases %>% filter(t == 78) %>% pull(log_cases)
cat(sprintf("Observed log(I(78)) = %.3f\n", obs_log_I_78))

# Compare with Hilbert-based CI from previous result
ci_lower <- res_T70$ci[1]
ci_upper <- res_T70$ci[2]
in_ci <- obs_log_I_78 >= ci_lower & obs_log_I_78 <= ci_upper
cat(sprintf("Is observed value within CI? %s\n", ifelse(in_ci, "YES", "NO")))


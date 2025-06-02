
# ----------------------------------------
# Compare Delta and Hilbert methods on China data
# ----------------------------------------
library(dplyr)
library(zoo)
library(COVID19)
library(numDeriv)
library(pracma)

# 1. Load and prepare China data
x <- covid19(country = "China", level = 1, start = "2020-01-22", end = "2020-04-01", verbose = FALSE)

x <- x %>%
  arrange(date) %>%
  mutate(new_confirmed = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_confirmed)) %>%
  mutate(
    log_cases = log(pmax(new_confirmed, 1)),
    t = row_number()
  )

data_T18 <- x %>% filter(t <= 18)
data_T21 <- x %>% filter(t <= 21)

# 2. LAD objective function
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

# 3. Delta method estimator
estimate_peak_model <- function(df, t_peak = 22, tau = 0.95) {
  t_vec <- df$t
  log_I <- df$log_cases
  theta0 <- c(log(max(exp(log_I))), median(t_vec), 50)
  fit <- optim(par = theta0,
               fn = lad_objective,
               method = "L-BFGS-B",
               lower = c(log(1), 1, 1e-2),
               upper = c(log(1e6), max(t_vec)+20, 1e4),
               t_vec = t_vec,
               log_I = log_I,
               tau = tau)
  est <- fit$par
  names(est) <- c("log_a", "l", "s2")
  log_peak_pred <- est["log_a"] - (t_peak - est["l"])^2 / est["s2"]
  grad_g <- c(1, 2 * (t_peak - est["l"]) / est["s2"], (t_peak - est["l"])^2 / est["s2"]^2)
  hess <- hessian(f = lad_objective,
                  x = est,
                  t_vec = t_vec,
                  log_I = log_I,
                  tau = tau)
  cov_matrix <- tryCatch(solve(hess), error = function(e) matrix(NA, 3, 3))
  se_peak <- sqrt(as.numeric(t(grad_g) %*% cov_matrix %*% grad_g))
  ci <- c(log_peak_pred - 1.96 * se_peak, log_peak_pred + 1.96 * se_peak)
  return(list(log_peak_pred = log_peak_pred, ci = ci, se = se_peak))
}

# 4. Hilbert-based hybrid method
estimate_hybrid <- function(df, t_peak = 22, tau = 0.95, f0 = 0.5) {
  t_vec <- df$t
  log_I <- df$log_cases
  T_len <- length(t_vec)
  theta0 <- c(log(max(exp(log_I))), median(t_vec), 50)
  fit <- optim(par = theta0,
               fn = lad_objective,
               method = "L-BFGS-B",
               lower = c(log(1), 1, 1e-2),
               upper = c(log(1e6), max(t_vec)+20, 1e4),
               t_vec = t_vec,
               log_I = log_I,
               tau = tau)
  theta <- fit$par
  names(theta) <- c("log_a", "l", "s2")
  beta2 <- -1 / theta["s2"]
  beta1 <- 2 * theta["l"] / theta["s2"]
  beta0 <- theta["log_a"] - theta["l"]^2 / theta["s2"]
  beta <- c(beta0, beta1, beta2)
  scaling <- c(T_len^(1/2), T_len^(3/2), T_len^(5/2))
  H_inv <- solve(hilbert(3))
  omega2 <- (tau * (1 - tau)) / f0^2
  Sigma_beta <- omega2 * diag(1 / scaling) %*% H_inv %*% diag(1 / scaling)
  X <- c(1, t_peak, t_peak^2)
  log_I_pred <- sum(beta * X)
  var_log_I <- X %*% Sigma_beta %*% X
  se <- sqrt(as.numeric(var_log_I))
  ci <- c(log_I_pred - 1.96 * se, log_I_pred + 1.96 * se)
  return(list(log_peak_pred = log_I_pred, ci = ci, se = se))
}

# 5. Run comparisons for T = 18 and T = 21
res_delta_18 <- estimate_peak_model(data_T18, t_peak = 22)
res_delta_21 <- estimate_peak_model(data_T21, t_peak = 22)
res_hilbert_18 <- estimate_hybrid(data_T18, t_peak = 22)
res_hilbert_21 <- estimate_hybrid(data_T21, t_peak = 22)

# 6. True value for day 22
true_logI22 <- log(x$new_confirmed[x$t == 22])

# 7. Print summary
cat("=== Day 22 true log(I) ===\\n")
cat(sprintf("log(I[22]) = %.3f\\n\\n", true_logI22))

cat("=== Delta Method ===\\n")
cat(sprintf("T = 18: Pred = %.3f, CI = [%.3f, %.3f]\\n", res_delta_18$log_peak_pred, res_delta_18$ci[1], res_delta_18$ci[2]))
cat(sprintf("T = 21: Pred = %.3f, CI = [%.3f, %.3f]\\n\\n", res_delta_21$log_peak_pred, res_delta_21$ci[1], res_delta_21$ci[2]))

cat("=== Hilbert-Based Method ===\\n")
cat(sprintf("T = 18: Pred = %.3f, CI = [%.3f, %.3f]\\n", res_hilbert_18$log_peak_pred, res_hilbert_18$ci[1], res_hilbert_18$ci[2]))
cat(sprintf("T = 21: Pred = %.3f, CI = [%.3f, %.3f]\\n", res_hilbert_21$log_peak_pred, res_hilbert_21$ci[1], res_hilbert_21$ci[2]))

library(dplyr)
library(zoo)
library(pracma)
library(quantreg)
library(COVID19)

# Step 1: Prepare data (China COVID-19)
china <- covid19("China", start = "2020-01-22", end = "2020-04-01", level = 1, verbose = FALSE)
china <- china %>%
  arrange(date) %>%
  mutate(new_confirmed = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_confirmed)) %>%
  mutate(log_cases = log(pmax(new_confirmed, 1)),
         t = row_number(),
         t_sq = t^2)

data_T21 <- china %>% filter(t <= 21)

# Step 2: Quantile regression
model <- rq(log_cases ~ t + t_sq, tau = 0.95, data = data_T21)
beta <- coef(model)
names(beta) <- c("b0", "b1", "b2")

# Step 3: Hilbert-based covariance matrix
t_vec <- data_T21$t
T_len <- length(t_vec)
scaling <- c(sqrt(T_len), sqrt(sum(t_vec^2)), sqrt(sum(t_vec^4)))
H_inv <- solve(hilbert(3))
tau <- 0.95
# f0 <- 0.5  # can be refined

residuals <- model$residuals
f0 <- density(residuals, from = 0, to = 0, n = 1, kernel = "gaussian")$y

omega2 <- (tau * (1 - tau)) / (f0^2)
Sigma_beta <- omega2 * diag(1 / scaling) %*% H_inv %*% diag(1 / scaling)

# Step 4: Compute L and a
L_hat <- -beta["b1"] / (2 * beta["b2"])
log_a_hat <- beta["b0"] - beta["b1"]^2 / (4 * beta["b2"])
a_hat <- exp(log_a_hat)

# Step 5: Jacobian of (a, L) with respect to (b0, b1, b2)
J <- matrix(c(
  a_hat,                            -a_hat * beta["b1"] / (2 * beta["b2"]),   a_hat * beta["b1"]^2 / (4 * beta["b2"]^2),
  0,                                -1 / (2 * beta["b2"]),                    beta["b1"] / (2 * beta["b2"]^2)
), nrow = 2, byrow = TRUE)

# Step 6: Delta method for covariance of (a, L)
Sigma_aL <- J %*% Sigma_beta %*% t(J)
se_a <- sqrt(Sigma_aL[1,1])
se_L <- sqrt(Sigma_aL[2,2])

# Confidence intervals (approximate, on original scale)
ci_a <- c(a_hat - 1.96 * se_a, a_hat + 1.96 * se_a)
ci_L <- c(L_hat - 1.96 * se_L, L_hat + 1.96 * se_L)

# Output
cat("=== Estimated Amplitude and Location of Peak ===\n")
cat(sprintf("â = %.2f, 95%% CI = [%.2f, %.2f]\n", a_hat, ci_a[1], ci_a[2]))
cat(sprintf("L̂ = %.2f, 95%% CI = [%.2f, %.2f]\n", L_hat, ci_L[1], ci_L[2]))

cat(sprintf("Standard Error for â: %.2f\n", se_a))
cat(sprintf("Standard Error for L̂: %.2f\n", se_L))

true_peak <- china %>% filter(t == 22) %>% pull(new_confirmed)
cat(sprintf("Observed peak I(22) = %d → log(I) = %.2f\n", true_peak, log(true_peak)))

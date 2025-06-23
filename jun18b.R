# ---------------------------------------------------------------------------
# Generate Figures 3, 4a, and 4b for tau = 0.5, 0.75, 0.95
# ---------------------------------------------------------------------------

# Laplace Noise Generator
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# Quantile Loss
quantile_loss <- function(u, tau) {
  u * ( tau - (u < 0) )
}

# Weighted Quantile
med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  o <- order(x)
  x_sorted <- x[o]
  w_sorted <- w[o]
  cum_w <- cumsum(w_sorted)
  cutoff <- tau * sum(w_sorted)
  k <- min(which(cum_w >= cutoff))
  return(x_sorted[k]) 
}

# Estimators
compute_a <- function(I, l, s, tau) {
  med_weighted_tau(I, l, s, tau)
}

objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  u <- I - fitted
  mean(quantile_loss(u, tau))
}

search_best_s <- function(I, l, tau) {
  f_obj <- function(s) {
    a_hat <- compute_a(I, l, s, tau)
    objective_function(I, a_hat, l, s, tau)
  }
  optimize(f_obj, interval = c(1, 20), tol = 1e-4)$minimum
}

estimate_params <- function(I, tau, h) {
  losses <- sapply(1:h, function(l) {
    s_hat <- search_best_s(I, l, tau)
    a_hat <- compute_a(I, l, s_hat, tau)
    objective_function(I, a_hat, l, s_hat, tau)
  })
  l_hat <- which.min(losses)
  s_hat <- search_best_s(I, l_hat, tau)
  a_hat <- compute_a(I, l_hat, s_hat, tau)
  list(l_hat = l_hat, s_hat = s_hat, a_hat = a_hat)
}

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
seed <- 280
h <- 30
T_vals <- 4:h
taus <- c(0.5, 0.75, 0.95)

# ----------------------------------------
# Combined 3x3 Plot for tau = 0.5, 0.75, 0.95
# ----------------------------------------

png("fig3_4a_4b_combined.png", width = 1200, height = 800, res = 120)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))  # 3 rows (tau), 3 columns (l, a, s)

for (tau in taus) {
  set.seed(seed)
  signal_full <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
  b <- A_true / (sqrt(2) * snr)
  noise_full <- rlaplace(h, scale = b)
  I_full <- signal_full + noise_full
  
  l_hats <- numeric(length(T_vals))
  a_hats <- numeric(length(T_vals))
  s_hats <- numeric(length(T_vals))
  
  for (i in seq_along(T_vals)) {
    T <- T_vals[i]
    I_T <- I_full[1:T]
    est <- estimate_params(I_T, tau, h)
    l_hats[i] <- est$l_hat
    a_hats[i] <- est$a_hat
    s_hats[i] <- est$s_hat
  }
  
  # Row = tau, Column 1 = l_hat(T)
  plot(T_vals, l_hats, type = "b", col = "blue", pch = 19,
       xlab = "T", ylab = expression(hat(l)),
       main = bquote(hat(l)(T) ~ "for" ~ tau == .(tau)))
  abline(h = L_true, col = "red", lty = 2)
  grid()
  
  # Column 2 = a_hat(T)
  plot(T_vals, a_hats, type = "b", col = "darkgreen", pch = 19,
       xlab = "T", ylab = expression(hat(a)),
       main = bquote(hat(a)(T) ~ "for" ~ tau == .(tau)))
  abline(h = A_true, col = "red", lty = 2)
  grid()
  
  # Column 3 = s_hat(T)
  plot(T_vals, s_hats, type = "b", col = "purple", pch = 19,
       xlab = "T", ylab = expression(hat(s)),
       main = bquote(hat(s)(T) ~ "for" ~ tau == .(tau)))
  abline(h = S_true, col = "red", lty = 2)
  grid()
}

dev.off()
cat("Combined 3x3 plot saved as 'fig3_4a_4b_combined.png'\n")

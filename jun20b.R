# ---------------------------------------------------------------------------
# Convergence check: T ↦ l̂(T) under different SNRs (tau = 0.5, z = 0)
# ---------------------------------------------------------------------------

# Asymmetric Laplace variance
asy_laplace_variance <- function(tau, z) {
  mu <- tau / (1 + z) - (1 - tau) / (1 - z)
  ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
  ex2 - mu^2
}

# Compute noise scale b from SNR
compute_noise_scale <- function(A, snr, tau, z) {
  sigma2 <- asy_laplace_variance(tau, z)
  sqrt((A / snr)^2 / sigma2)
}

# Asymmetric Laplace noise generator
r_asy_laplace <- function(n, tau, scale = 1) {
  u <- runif(n)
  ifelse(u < tau, -rexp(n, rate = 1 / scale), rexp(n, rate = 1 / scale))
}

# Quantile loss and weighted quantile
quantile_loss <- function(u, tau) u * (tau - (u < 0))

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
  x_sorted[k]
}

# Estimation routines
compute_a <- function(I, l, s, tau) med_weighted_tau(I, l, s, tau)

objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  mean(quantile_loss(I - fitted, tau))
}

search_best_s <- function(I, l, tau) {
  f <- function(s) objective_function(I, compute_a(I, l, s, tau), l, s, tau)
  optimize(f, c(1, 20), tol = 1e-4)$minimum
}

estimate_params <- function(I, tau, h) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau)
    a <- compute_a(I, l, s, tau)
    objective_function(I, a, l, s, tau)
  })
  l_hat <- which.min(losses)
  s_hat <- search_best_s(I, l_hat, tau)
  a_hat <- compute_a(I, l_hat, s_hat, tau)
  list(l_hat = l_hat, s_hat = s_hat, a_hat = a_hat)
}

# Settings
A_true <- 1
L_true <- 22
S_true <- 10
tau <- 0.5
z <- 0
h <- 30
T_vals <- 4:h
snr_vals <- c(0.9, 5, 40, 1000)
colors <- c("black", "red", "blue", "darkgreen")

# Fixed seed for reproducibility
set.seed(280)

# Plotting
png("lhat_convergence_by_snr.png", width = 900, height = 600, res = 120)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (i in seq_along(snr_vals)) {
  snr <- snr_vals[i]
  b <- compute_noise_scale(A_true, snr, tau, z)
  signal <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
  noise <- r_asy_laplace(h, tau, scale = b)
  I <- signal + noise
  
  l_hats <- sapply(T_vals, function(T) {
    estimate_params(I[1:T], tau, h)$l_hat
  })
  
  plot(T_vals, l_hats, type = "b", col = colors[i], pch = 19,
       xlab = "T", ylab = expression(hat(l)(T)),
       main = paste0("SNR = ", snr))
  abline(h = L_true, col = "gray40", lty = 2)
  grid()
}

dev.off()
cat("Saved convergence plot as 'lhat_convergence_by_snr.png'\\n")

# Quantile-Based Peak Estimation for Gaussian Trend with Additive Laplace Noise
# Cleaned & Modular Version
# Author: Said Maanan
# Date: 2025-06-12

# ---- Libraries ----
# None needed beyond base R

# ---- Laplace Noise Generator ----
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# ---- Quantile Loss Function ----
quantile_loss <- function(u, tau) {
  u * (tau - (u < 0))
}

# ---- Weighted Quantile Function ----
weighted_quantile <- function(x, w, tau) {
  stopifnot(all(w >= 0), length(x) == length(w), tau > 0, tau < 1)
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  cum_w <- cumsum(w_sorted)
  threshold <- tau * sum(w_sorted)
  k <- max(which(cum_w <= threshold), 0)
  if (k == 0) return(x_sorted[1])
  if (cum_w[k] < threshold && k < length(x_sorted)) return(x_sorted[k + 1])
  return(x_sorted[k])
}

# ---- Weighted Median (tau = 0.5) ----
med_weighted <- function(I, l, s) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  o <- order(x)
  w <- exp(-((1:T) - l)^2 / s^2)
  wo <- w[o]
  W <- sum(w)
  y <- rep(0, T)
  for (t in 1:T) {
    if (sum(wo[1:t]) <= W / 2) y[t] <- t
  }
  return(x[o[max(y) + 1]])
}

# ---- Compute a(T, l, s) ----
compute_a <- function(I, l, s, tau) {
  if (tau == 0.5) {
    return(med_weighted(I, l, s))
  }
  t <- seq_along(I)
  w <- exp(-((t - l)^2) / s^2)
  x <- I * exp(((t - l)^2) / s^2)
  weighted_quantile(x, w, tau)
}

# ---- Objective Function f_tau ----
objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  u <- I - fitted
  mean(quantile_loss(u, tau))
}

# ---- Optimize s for fixed l ----
search_best_s <- function(I, l, tau) {
  f_obj <- function(s) {
    a_hat <- compute_a(I, l, s, tau)
    objective_function(I, a_hat, l, s, tau)
  }
  optimize(f_obj, interval = c(1, 20), tol = 1e-4)$minimum
}

# ---- Simulation Function ----
simulate_quantile_objective <- function(
    T, A, L, S, h, tau,
    snr = Inf, seed = NULL,
    l_range = NULL, force_l_gt_T = FALSE
) {
  t <- 1:T
  signal <- A * exp(-((t - L)^2) / S^2)
  
  if (!is.null(seed)) set.seed(seed)
  if (is.finite(snr)) {
    b <- A / (sqrt(2) * snr)
    noise <- rlaplace(T, scale = b)
    I <- signal + noise
  } else {
    I <- signal
  }
  
  if (is.null(l_range)) {
    l_grid <- if (force_l_gt_T) seq(T + 1, h) else if (is.finite(snr)) seq(L - 4, L + 3) else seq(1, h)
  } else {
    l_grid <- l_range
  }
  
  losses <- sapply(l_grid, function(l) {
    s_hat <- search_best_s(I, l, tau)
    a_hat <- compute_a(I, l, s_hat, tau)
    objective_function(I, a_hat, l, s_hat, tau)
  })
  
  data.frame(l = l_grid, loss = losses, T = T)
}

# ---- Plotting Function (Base R) ----
plot_quantile_objectives <- function(results, tau, filename) {
  Ts <- unique(results$T)
  nT <- length(Ts)
  nrow <- ceiling(nT / 2)
  ncol <- if (nT == 1) 1 else 2
  
  png(filename, width = 900, height = 600, res = 120)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
  on.exit({ par(old_par); dev.off() })
  
  for (T_val in Ts) {
    subset_data <- subset(results, T == T_val)
    plot(subset_data$l, subset_data$loss, type = "l", lwd = 2, col = "blue",
         xlab = "Peak Location Candidate l",
         ylab = expression("Objective " * f[tau]),
         main = paste0("T = ", T_val))
    abline(v = 22, col = "red", lty = 2)
    grid()
  }
}

# ---- Parameters ----
A_true <- 15000; L_true <- 22; S_true <- 10; h <- 30
Ts_free <- c(2, 3, 10, 20)

# ---- Run: Noise-Free, tau = 0.5 (Median Case) ----
tau_val <- 0.5
results_free_median <- do.call(rbind, lapply(Ts_free, function(T) {
  simulate_quantile_objective(T, A = A_true, L = L_true, S = S_true, h = h,
                              tau = tau_val, snr = Inf, force_l_gt_T = TRUE)
}))
plot_quantile_objectives(results_free_median, tau = tau_val, filename = "plot_noisefree_tau05.png")

# ---- Run: Noise-Free, tau = 0.95 (Quantile Case) ----
tau_val <- 0.95
results_free <- do.call(rbind, lapply(Ts_free, function(T) {
  simulate_quantile_objective(T, A = A_true, L = L_true, S = S_true, h = h,
                              tau = tau_val, snr = Inf)
}))
plot_quantile_objectives(results_free, tau = tau_val, filename = "plot_noisefree_tau095.png")
# Quantile-Based Peak Estimation for Gaussian Trend with Additive Laplace Noise
#
# This script implements quantile regression for estimating the amplitude, location,
# and width of a Gaussian peak in the presence of Laplace noise, as described in the
# provided LaTeX manuscript. The approach generalizes LAD (Least Absolute Deviations)
# to any quantile level tau in (0,1), using the check loss and weighted quantile estimator.
#
# Author: [Said Maanan]
# Date: [12-06-2025]
#
# ------------------------------------------------------------------------------

# ---- Laplace Noise Generator ----
# Returns n samples from Laplace(location, scale)
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# ---- Quantile Loss (Check Function) ----
# Implements rho_tau(u) = u * (tau - 1_{u<0})
quantile_loss <- function(u, tau) {
  u * (tau - (u < 0))
}

# ---- Weighted Quantile Function ----
# Given values x, positive weights w, and quantile level tau in (0,1),
# returns the weighted tau-quantile as described in the manuscript/proposition.
weighted_quantile <- function(x, w, tau) {
  # Ensure no negative weights
  stopifnot(all(w >= 0), length(x) == length(w), tau > 0, tau < 1)
  # Sort x and reorder corresponding weights
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  cum_w <- cumsum(w_sorted)
  total_w <- sum(w_sorted)
  threshold <- tau * total_w
  
  # Find largest k such that cum_w[k] <= threshold
  k <- max(which(cum_w <= threshold), 0)
  # If threshold is never exceeded, return min(x)
  if (k == 0) return(x_sorted[1])
  # If cum_w[k] < threshold, quantile is at x[k+1], handle edge case
  if (cum_w[k] < threshold && k < length(x_sorted)) {
    return(x_sorted[k + 1])
  } else {
    return(x_sorted[k])
  }
}

# ---- Compute a(T, l, s) for Given I(t), l, s, tau ----
# Returns the weighted tau-quantile of transformed data
compute_a <- function(I, l, s, tau) {
  t <- seq_along(I)
  w <- exp(-((t - l)^2) / s^2)
  x <- I * exp(((t - l)^2) / s^2)
  weighted_quantile(x, w, tau)
}

# ---- Objective Function f_tau(T, a, l, s) ----
# Mean quantile loss over the data
objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  u <- I - fitted
  mean(quantile_loss(u, tau))
}

# ---- Search Best s for Given l ----
# Finds s in s_grid that minimizes the quantile objective for fixed l
search_best_s <- function(I, l, tau, s_grid) {
  losses <- sapply(s_grid, function(s) {
    a_hat <- compute_a(I, l, s, tau)
    objective_function(I, a_hat, l, s, tau)
  })
  s_grid[which.min(losses)]
}

# ---- Main Simulation Function ----
# Simulates the quantile objective curve l ↦ f_tau(T, a(T,l,s(T,l)), l, s(T,l))
# for a given T (number of observations), optionally with Laplace noise.
simulate_quantile_objective <- function(
    T,
    A = 1,
    L = 22,
    S = 10,
    h = 30,
    tau = 0.95,
    snr = Inf,     # Use Inf for noise-free; otherwise, SNR used for Laplace noise
    seed = NULL,
    l_range = NULL # Optional override for l grid (for Laplace case)
) {
  t <- 1:T
  signal <- A * exp(-((t - L)^2) / S^2)
  
  # Generate noise
  if (!is.null(seed)) set.seed(seed)
  if (is.finite(snr)) {
    b <- A / (sqrt(2) * snr)
    noise <- rlaplace(T, scale = b)
    I <- signal + noise
  } else {
    I <- signal
  }
  
  # l_grid: which l to scan (default: full for noise-free, local for noisy)
  if (is.null(l_range)) {
    l_grid <- if (is.finite(snr)) seq(L-4, L+3) else seq(1, h)
  } else {
    l_grid <- l_range
  }
  # s_grid: reasonable range for width parameter search
  s_grid <- seq(3, 20, length.out = 50)
  
  # For each l, search best s and compute minimized objective
  losses <- sapply(l_grid, function(l) {
    s_hat <- search_best_s(I, l, tau, s_grid)
    a_hat <- compute_a(I, l, s_hat, tau)
    objective_function(I, a_hat, l, s_hat, tau)
  })
  
  data.frame(l = l_grid, loss = losses, T = T)
}

# ---- Example 1: Noise-Free Case ----
# Visualize the quantile objective curve for various T, no noise
# Use a moderate quantile (e.g., tau = 0.95)
tau_val <- 0.95
A_true <- 15000
L_true <- 22
S_true <- 10
horizon <- 30
Ts_free <- c(2, 3, 10, 20) # Number of observations to test

results_free <- do.call(rbind, lapply(Ts_free, function(T) {
  simulate_quantile_objective(
    T, A = A_true, L = L_true, S = S_true, h = horizon, tau = tau_val, snr = Inf
  )
}))

# ---- Example 2: Laplace Noise Case ----
# Simulate with Laplace noise, SNR = 1000, for various T
SNR_val <- 1000
seed_val <- 280
Ts_laplace <- c(4, 6, 8, 9, 10)

results_laplace <- do.call(rbind, lapply(Ts_laplace, function(T) {
  simulate_quantile_objective(
    T, A = 1, L = L_true, S = S_true, h = horizon, tau = tau_val,
    snr = SNR_val, seed = seed_val, l_range = seq(L_true-4, L_true+3)
  )
}))

# ---- Base R Plotting Function ----
# Plots the objective function versus l for various T using base R

# Noise-Free Case Plot — No Legend, No Main Title, Saved as PNG
plot_quantile_objectives_noisefree <- function(results, tau, filename = "plot_noisefree.png") {
  Ts <- unique(results$T)
  nT <- length(Ts)
  
  # Layout: auto rows and columns
  nrow <- ceiling(nT / 2)
  ncol <- if (nT == 1) 1 else 2
  
  png(filename, width = 900, height = 600, res = 120)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1), oma = c(0, 0, 0, 0))
  on.exit({
    par(old_par)
    dev.off()
  })
  
  for (T_val in Ts) {
    subset_data <- subset(results, T == T_val)
    plot(subset_data$l, subset_data$loss, type = "l", lwd = 2, col = "blue",
         xlab = "Peak Location Candidate l",
         ylab = expression("Objective " ~ f[tau]),
         main = paste0("T = ", T_val))
    
    # idx_min <- which.min(subset_data$loss)
    # points(subset_data$l[idx_min], subset_data$loss[idx_min], col = "darkgreen", pch = 19)
    # 
    abline(v = 22, col = "red", lty = 2)
    grid()
  }
}


# Laplace Noise Case Plot — 2x3 Layout, No Legend, No Main Title, Saved as PNG
plot_quantile_objectives_laplace <- function(results, tau, filename = "plot_laplace.png") {
  Ts <- unique(results$T)
  
  png(filename, width = 900, height = 600, res = 120)
  old_par <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 0, 0))
  on.exit({
    par(old_par)
    dev.off()
  })
  
  for (T_val in Ts) {
    subset_data <- subset(results, T == T_val)
    plot(subset_data$l, subset_data$loss, type = "l", lwd = 2, col = "blue",
         xlab = "Peak Location Candidate l",
         ylab = expression("Objective " ~ f[tau]),
         main = paste0("T = ", T_val))
    
    # idx_min <- which.min(subset_data$loss)
    # points(subset_data$l[idx_min], subset_data$loss[idx_min], col = "darkgreen", pch = 19)
    # 
    abline(v = 22, col = "red", lty = 2)
    grid()
  }
}


plot_quantile_objectives_noisefree(results_free, tau = tau_val)
plot_quantile_objectives_laplace(results_laplace, tau = tau_val)

# Tau=0.5 Case

tau_val <- 0.5

# Run the noise-free simulation at tau = 0.5
results_free_median <- do.call(rbind, lapply(Ts_free, function(T) {
  simulate_quantile_objective(
    T, A = A_true, L = L_true, S = S_true, h = horizon, tau = tau_val, snr = Inf
  )
}))

# Plot
plot_quantile_objectives_noisefree(results_free_median, tau = tau_val, filename = "plot_noisefree_tau05.png")

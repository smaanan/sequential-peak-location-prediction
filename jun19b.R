# ---------------------------------------------------------------------------
# Quantile-Based Peak Estimation with shifted Laplace to match tau (updated SNR)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Laplace Quantile (inverse CDF) 
# ---------------------------------------------------------------------------
qlaplace <- function(tau, scale = 1) {
  if (tau < 0.5) scale * log(2 * tau) else -scale * log(2 * (1 - tau))
}

# ---------------------------------------------------------------------------
# Laplace with shifted center to match P(e < 0) = tau
# ---------------------------------------------------------------------------
shifted_laplace <- function(n, tau, scale = 1) {
  rlap <- rlaplace(n, scale = scale)
  rlap - qlaplace(tau, scale = scale)
}

# ---------------------------------------------------------------------------
# Laplace Noise Generator (original, centered at 0) 
# ---------------------------------------------------------------------------
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# ---------------------------------------------------------------------------
# Asymmetric Laplace Variance
# ---------------------------------------------------------------------------
asy_laplace_variance <- function(tau, z) {
  stopifnot(tau > 0, tau < 1, abs(z) < 1)
  mu <- tau / (1 + z) - (1 - tau) / (1 - z)
  E_X2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
  return(E_X2 - mu^2)
}

# ---------------------------------------------------------------------------
# Quantile Loss
# ---------------------------------------------------------------------------
quantile_loss <- function(u, tau) {
  u * (tau - (u < 0))
}

# ---------------------------------------------------------------------------
# Weighted Quantile
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Quantile-Based Amplitude
# ---------------------------------------------------------------------------
compute_a <- function(I, l, s, tau) {
  med_weighted_tau(I, l, s, tau)
}

# ---------------------------------------------------------------------------
# Quantile Loss objective
# ---------------------------------------------------------------------------
objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  u <- I - fitted
  mean(quantile_loss(u, tau))
}

# ---------------------------------------------------------------------------
# Optimize s for a given l
# ---------------------------------------------------------------------------
search_best_s <- function(I, l, tau) {
  f_obj <- function(s) {
    a_hat <- compute_a(I, l, s, tau)
    objective_function(I, a_hat, l, s, tau)
  }
  optimize(f_obj, interval = c(1, 20), tol = 1e-4)$minimum
}

# ---------------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------------
simulate_quantile_objective <- function(T, A, L, S, tau, z, snr = Inf, seed = NULL) {
  t <- 1:T
  signal <- A * exp(-((t - L)^2) / S^2)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.finite(snr)) {
    sigma2 <- asy_laplace_variance(tau, z)
    b <- sqrt((A / snr)^2 / sigma2)
    noise <- shifted_laplace(T, tau, scale = b)
    I <- signal + noise
  } else {
    I <- signal
  }
  
  l_grid <- 18:25
  losses <- sapply(l_grid, function(l) {
    s_hat <- search_best_s(I, l, tau)
    a_hat <- compute_a(I, l, s_hat, tau)
    objective_function(I, a_hat, l, s_hat, tau)
  })
  
  data.frame(l = l_grid, loss = losses, T = T)
}

# ---------------------------------------------------------------------------
# Plotting utility
# ---------------------------------------------------------------------------
plot_quantile_objectives <- function(results, tau, A_true, L_true, S_true) {
  Ts <- unique(results$T)
  nrow <- 2
  ncol <- 3
  
  png(filename = paste0("shifted_obj_fun_tau_", tau, ".png"),
      width = 900, height = 600, res = 120)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
  on.exit({ par(old_par); dev.off() })
  
  for (T_val in Ts) {
    subset_data <- subset(results, T == T_val)
    plot(subset_data$l, subset_data$loss, type = "l", lwd = 2, col = "blue",
         xlab = "Peak Location Candidate l",
         ylab = expression("Objective " * f),
         main = paste0("T = ", T_val, ", tau = ", tau))
    abline(v = L_true, col = "red", lty = 2)
    grid()
  }
}

# ---------------------------------------------------------------------------
# Run simulation and plot
# ---------------------------------------------------------------------------
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
z <- 0
seed <- 280
T_vals <- c(4, 6, 8, 9, 10)

for (tau_val in c(0.5, 0.75, 0.95)) {
  results <- do.call(rbind, lapply(T_vals, function(T) {
    simulate_quantile_objective(T, A = A_true, L = L_true, S = S_true,
                                tau = tau_val, z = z, snr = snr, seed = seed)
  }))
  
  plot_quantile_objectives(results, tau_val, A_true, L_true, S_true)
  cat("Plot for tau =", tau_val, "generated and saved.\n")
}






# ------------------------------------------
# (1) Updated simulate_quantile_objective
# ------------------------------------------
simulate_quantile_objective <- function(T, A, L, S, tau, z, snr = Inf, seed = NULL) {
  t <- 1:T
  signal <- A * exp(-((t - L)^2) / S^2)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.finite(snr)) {
    sigma2 <- asy_laplace_variance(tau, z)
    b <- sqrt((A / snr)^2 / sigma2)
    noise <- shifted_laplace(T, tau, scale = b)
    I <- signal + noise
  } else {
    I <- signal
  }
  
  l_grid <- 18:25
  losses <- sapply(l_grid, function(l) {
    s_hat <- search_best_s(I, l, tau)
    a_hat <- compute_a(I, l, s_hat, tau)
    objective_function(I, a_hat, l, s_hat, tau)
  })
  
  data.frame(l = l_grid, loss = losses, T = T, z = z)
}

# ------------------------------------------
# (2) Plotting utility with z overlay
# ------------------------------------------
plot_quantile_objectives_superposed <- function(results, tau, A_true, L_true, S_true) {
  Ts <- sort(unique(results$T))
  z_vals <- sort(unique(results$z))
  colors <- c("black", "red", "blue", "darkgreen", "orange")
  stopifnot(length(z_vals) <= length(colors))
  
  png(filename = paste0("shifted_obj_fun_tau_", tau, "_zcompare.png"),
      width = 900, height = 600, res = 120)
  old_par <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  on.exit({ par(old_par); dev.off() })
  
  for (T_val in Ts) {
    df_T <- subset(results, T == T_val)
    x_range <- range(df_T$l)
    y_range <- range(df_T$loss)
    
    plot(NA, xlim = x_range, ylim = y_range,
         xlab = "Peak Location Candidate l",
         ylab = expression("Objective " * f),
         main = paste0("T = ", T_val, ", tau = ", tau))
    abline(v = L_true, col = "gray40", lty = 2)
    grid()
    
    for (i in seq_along(z_vals)) {
      df_z <- subset(df_T, z == z_vals[i])
      lines(df_z$l, df_z$loss, col = colors[i], lwd = 2)
    }
    
    legend("topright", legend = paste("z =", z_vals),
           col = colors[1:length(z_vals)], lwd = 2, bty = "o", cex = 0.8)
  }
}


# ------------------------------------------
# (3) Main execution loop for superposed Figure 2
# ------------------------------------------
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
seed <- 280
T_vals <- c(4, 6, 8, 9, 10)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)

for (tau_val in c(0.5, 0.75, 0.95)) {
  results_all <- do.call(rbind, lapply(z_vals, function(z) {
    res_z <- do.call(rbind, lapply(T_vals, function(T) {
      simulate_quantile_objective(T, A = A_true, L = L_true, S = S_true,
                                  tau = tau_val, z = z, snr = snr, seed = seed)
    }))
    res_z$z <- z
    res_z
  }))
  
  plot_quantile_objectives_superposed(results_all, tau_val, A_true, L_true, S_true)
  cat("Figure 2 for tau =", tau_val, "generated and saved with all z values.\n")
}

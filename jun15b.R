# ---------------------------------------------------------------------------
# Quantile-Based Peak Estimation with shifted Laplace to match tau
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Laplace Quantile (inverse CDF) 
# ---------------------------------------------------------------------------
qlaplace <- function(tau, scale = 1) {
  if ( tau < 0.5 )
    return( scale * log( 2 * tau ) )
  else
    return( - scale * log( 2 * (1 - tau) ) )
}

# ---------------------------------------------------------------------------
# Laplace with shifted center to match P(e < 0)= tau
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
# Quantile Loss
# ---------------------------------------------------------------------------
quantile_loss <- function(u, tau) {
  u * ( tau - (u < 0) )
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
  
  return(x_sorted[k]) 
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
simulate_quantile_objective <- function(T, A, L, S, tau, snr = Inf, seed = NULL) {
  t <- 1:T
  signal <- A * exp(-((t - L)^2) / S^2)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.finite(snr)) {
    b <- A / (sqrt(2) * snr)
    noise <- shifted_laplace(T, tau, scale = b )
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
# Combined Plotting Utility for All Tau and T Values (Shifted Laplace Case)
# ---------------------------------------------------------------------------
plot_all_shifted_objectives <- function(results_all, L_true) {
  Ts <- sort(unique(results_all$T))
  taus <- sort(unique(results_all$tau))
  
  png(filename = "plot_shifted_all.png", width = 900, height = 1200, res = 120)
  old_par <- par(mfrow = c(length(Ts), length(taus)), mar = c(4, 4, 3, 1))
  on.exit({ par(old_par); dev.off() })
  
  for (T_val in Ts) {
    for (tau_val in taus) {
      subset_data <- subset(results_all, T == T_val & tau == tau_val)
      plot(subset_data$l, subset_data$loss, type = "l", lwd = 2, col = "blue",
           xlab = "Peak Location Candidate l",
           ylab = expression("Objective " * f),
           main = bquote("T = " * .(T_val) * ",  " * tau * " = " * .(tau_val)))
      abline(v = L_true, col = "red", lty = 2)
      grid()
    }
  }
}

# ---------------------------------------------------------------------------
# Run simulation and produce combined plot
# ---------------------------------------------------------------------------
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
seed <- 280
T_vals <- c(4, 6, 8, 9, 10)
tau_values <- c(0.5, 0.75, 0.95)

results_all <- do.call(rbind, lapply(tau_values, function(tau_val) {
  do.call(rbind, lapply(T_vals, function(T) {
    df <- simulate_quantile_objective(
      T, A = A_true, L = L_true, S = S_true,
      tau = tau_val, snr = snr, seed = seed
    )
    df$tau <- tau_val
    df
  }))
}))

plot_all_shifted_objectives(results_all, L_true)
cat("Combined plot saved to plot_shifted_all.png\n")

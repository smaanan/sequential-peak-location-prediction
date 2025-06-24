# ---------------------------------------------------------------------------
# Quantile-Based Peak Estimation with Asymmetric Laplace to match tau
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Asymmetric Laplace with P(e < 0)= tau
# ---------------------------------------------------------------------------
rasy_laplace <- function(n, tau, z) {
  u <- runif(n)
  sample <- numeric(n)
  
  for (i in 1:n) {
    if (u[i] < tau) { 
      sample[i] <- -rexp(1, rate = 1 + z) 
    } else { 
      sample[i] <- rexp(1, rate = 1 - z) 
    }
  }
  
  return(sample)
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
simulate_quantile_objective_asy <- function(T, A, L, S, tau, z, snr = Inf, seed = NULL) {
  t <- 1:T
  signal <- A * exp(-((t - L)^2) / S^2)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.finite(snr)) {
    b <- A / (sqrt(2) * snr)
    noise <- rasy_laplace(T, tau, z) * b
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
# Combined Plotting Utility for All Tau and T Values (Asymmetric Laplace)
# ---------------------------------------------------------------------------
plot_all_asymmetric_objectives <- function(results_all, L_true) {
  Ts <- sort(unique(results_all$T))
  taus <- sort(unique(results_all$tau))
  pdf(file = "plot_asymmetric_all.pdf", width = 7.5, height = 10)
  # png(filename = "plot_asymmetric_all.png", width = 900, height = 1200, res = 120)
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
# Run simulation and generate merged plot
# ---------------------------------------------------------------------------
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
seed <- 280
T_vals <- c(4, 6, 8, 9, 10)
tau_values <- c(0.5, 0.75, 0.95)
z <- 0  # Controls asymmetry; z = 0 gives symmetric Laplace

results_all <- do.call(rbind, lapply(tau_values, function(tau_val) {
  do.call(rbind, lapply(T_vals, function(T) {
    df <- simulate_quantile_objective_asy(
      T, A = A_true, L = L_true, S = S_true,
      tau = tau_val, z = z, snr = snr, seed = seed
    )
    df$tau <- tau_val
    df
  }))
}))

plot_all_asymmetric_objectives(results_all, L_true)
cat("Combined plot saved to plot_asymmetric_all.png\n")

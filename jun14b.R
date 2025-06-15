
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

quantile_loss <- function(u, tau) {
  u * ( tau - (u < 0) )
}


med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I)
  
  # Transform signals and weights first
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


simulate_quantile_objective <- function(T, A, L, S, h, tau, snr = Inf, seed = NULL, l_range = NULL, force_l_gt_T = FALSE) {
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

# ---------------------------------------------------------------------------
# Plotting utility
# ---------------------------------------------------------------------------
plot_quantile_objectives <- function(results, tau, A_true, L_true, S_true, h) {
  Ts <- unique(results$T)
  nT <- length(Ts)
  nrow <- ceiling(nT / 2)
  ncol <- if (nT == 1) 1 else 2
  
  png(filename = paste0("plot_quantile_tau", tau, ".png"),
      width = 900, height = 600, res = 120)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
  on.exit({ par(old_par); dev.off() })
  
  for (T_val in Ts) {
    subset_data <- subset(results, T == T_val)
    plot(subset_data$l, subset_data$loss, type = "l", lwd = 2, col = "blue",
         xlab = "Peak Location Candidate l",
         ylab = expression("Objective " * f[tau]),
         main = paste0("T = ", T_val, ", tau = ", tau))
    abline(v = L_true, col = "red", lty = 2)
    grid()
  }
}


A_true <- 15000
L_true <- 22
S_true <- 10
h <- 30
Ts_free <- c(2, 3, 10, 20)

for ( tau_val in c(0.5, 0.75, 0.95) ) {
  
  results <- do.call(rbind, lapply(Ts_free, function(T) {
    simulate_quantile_objective(T, A = A_true, L = L_true, S = S_true, h = h,
                                tau = tau_val, snr = Inf, force_l_gt_T = ( tau_val == 0.5 )) 
  }))
  
  plot_file <- paste0("plot_quantile_tau", tau_val, ".png")
  plot_quantile_objectives(results, tau_val, A_true, L_true, S_true, h)
  
  cat("Plot for tau =", tau_val, "generated and saved to", plot_file, "\n")
}


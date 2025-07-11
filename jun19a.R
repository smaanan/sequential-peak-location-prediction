# ---------------------------------------------------------------------------
# Quantile-Based Peak Estimation with Laplace Noise (Multiple z, tau)
# ---------------------------------------------------------------------------
# Author: Said Maanan
# Refactored: 2025-06-19
# ---------------------------------------------------------------------------

# --- Core Functions ---
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# asy_laplace_variance <- function(tau, z) {
#   stopifnot(tau > 0, tau < 1, abs(z) < 1)
#   mu <- tau / (1 + z) - (1 - tau) / (1 - z)
#   ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
#   ex2 - mu^2
# }

asy_laplace_variance <- function(tau, z) {
  t1 <- tau / (1 + z)^2
  t2 <- (1 - tau) / (1 - z)^2
  mu <- (1 - tau) / (1 - z) - tau / (1 + z)
  t1 + t2 - 0.25 * mu^2
}

compute_noise_scale <- function(A, snr, tau, z) {
  sqrt((A / snr)^2 / asy_laplace_variance(tau, z))
}

quantile_loss <- function(u, tau) u * (tau - (u < 0))

med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  o <- order(x)
  x_sorted <- x[o]; w_sorted <- w[o]
  cum_w <- cumsum(w_sorted)
  k <- min(which(cum_w >= tau * sum(w_sorted)))
  x_sorted[k]
}

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

# --- Simulation ---
simulate_quantile_objective <- function(T, A, L, S, h, tau, snr, z, seed = NULL, force_l_gt_T = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  t <- 1:T
  signal <- A * exp(-((t - L)^2) / S^2)
  b <- compute_noise_scale(A, snr, tau, z)
  noise <- rlaplace(T, scale = b)
  I <- signal + noise
  l_grid <- if (force_l_gt_T) seq(T + 1, h) else seq(L - 4, L + 3)
  
  data.frame(
    l = l_grid,
    loss = sapply(l_grid, function(l) {
      s <- search_best_s(I, l, tau)
      a <- compute_a(I, l, s, tau)
      objective_function(I, a, l, s, tau)
    }),
    T = T, z = z
  )
}

# --- Config ---
A_true <- 15000
L_true <- 22
S_true <- 10
h <- 30
Ts_free <- c(2, 3, 10, 20, 21)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
snr <- 1000
seed <- 280
tau_vals <- c(0.25, 0.5, 0.75, 0.95)

# --- Combined Simulation and Plotting ---
results_all <- do.call(rbind, lapply(tau_vals, function(tau) {
  do.call(rbind, lapply(z_vals, function(z) {
    res_z <- do.call(rbind, lapply(Ts_free, function(T) {
      simulate_quantile_objective(T, A_true, L_true, S_true, h,
                                  tau, snr, z, seed, force_l_gt_T = (tau == 0.5))
    }))
    res_z$z <- z
    res_z$tau <- tau
    res_z
  }))
}))

plot_all_quantile_objectives_zcompare <- function(results_all, L_true) {
  Ts <- sort(unique(results_all$T))
  taus <- sort(unique(results_all$tau))
  z_vals <- sort(unique(results_all$z))
  colors <- c("blue", "red", "darkgreen", "black", "grey60")
  
  nrow <- length(Ts)
  ncol <- length(taus)
  
  pdf(file = "plot_all_tau_zcompare.pdf", width = 2.5 * ncol, height = 2.5 * nrow)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(2, 2, 1.5, 1), oma = c(4, 4, 3, 0), bg = "white")
  on.exit({ par(old_par); dev.off() })
  
  for (i in seq_along(Ts)) {
    for (j in seq_along(taus)) {
      T_val <- Ts[i]
      tau_val <- taus[j]
      df_sub <- subset(results_all, T == T_val & tau == tau_val)
      
      x_range <- range(df_sub$l)
      y_range <- range(df_sub$loss)
      
      plot(NA, xlim = x_range, ylim = y_range,
           xlab = "", ylab = "", main = "")
      abline(v = L_true, col = "gray40", lty = 2)
      grid()
      
      for (k in seq_along(z_vals)) {
        df_z <- subset(df_sub, z == z_vals[k])
        lines(df_z$l, df_z$loss, col = colors[k], lwd = 2)
      }
      
      # Add tau label on top
      if (i == 1) mtext(bquote(tau == .(tau_val)), side = 3, line = 0.5, cex = 1.1)
      # Add T label on the left
      if (j == 1) mtext(bquote(T == .(T_val)), side = 2, line = 3, cex = 1.1)
      
      # Legend only in top row
      if (i == 1) {
        legend("topright", legend = paste("z =", z_vals),
               col = colors, lwd = 2, cex = 0.8)
      }
    }
  }
  
  # mtext("Peak Location Candidate l", side = 1, outer = TRUE, line = 2.5, cex = 1.2)
  # mtext(expression("Objective Function " * f[tau]), side = 2, outer = TRUE, line = 2.5, cex = 1.2)
  # mtext("Quantile-Based Peak Estimation — Laplace Noise", side = 3, outer = TRUE, line = 1.5, cex = 1.4)
}

plot_all_quantile_objectives_zcompare(results_all, L_true)
cat("Merged plot saved to plot_all_tau_zcompare.png\n")

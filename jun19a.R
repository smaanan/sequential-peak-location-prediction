# # ---------------------------------------------------------------------------
# # Quantile-Based Peak Estimation with Laplace Noise
# # ---------------------------------------------------------------------------
# # Author: Said Maanan
# # Date: 2025-06-12
# # ---------------------------------------------------------------------------
# 
# # ---------------------------------------------------------------------------
# # Laplace Noise Generator
# # ---------------------------------------------------------------------------
# rlaplace <- function(n, location = 0, scale = 1) {
#   u <- runif(n, -0.5, 0.5)
#   location - scale * sign(u) * log(1 - 2 * abs(u))
# }
# 
# # ---------------------------------------------------------------------------
# # Asymmetric Laplace Variance
# # ---------------------------------------------------------------------------
# asy_laplace_variance <- function(tau, z) {
#   stopifnot(tau > 0, tau < 1, abs(z) < 1)
#   mu <- tau / (1 + z) - (1 - tau) / (1 - z)
#   E_X2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
#   return(E_X2 - mu^2)
# }
# 
# # ---------------------------------------------------------------------------
# # Compute scale b for desired SNR
# # ---------------------------------------------------------------------------
# compute_noise_scale <- function(A, SNR_target, tau, z) {
#   sigma2 <- asy_laplace_variance(tau, z)
#   sqrt((A / SNR_target)^2 / sigma2)
# }
# 
# # ---------------------------------------------------------------------------
# # Quantile Loss Function
# # ---------------------------------------------------------------------------
# quantile_loss <- function(u, tau) {
#   u * ( tau - (u < 0) )
# }
# 
# # ---------------------------------------------------------------------------
# # Weighted Quantile (Generalized Algorithm, tau in (0,1))
# # ---------------------------------------------------------------------------
# med_weighted_tau <- function(I, l, s, tau) {
#   T <- length(I)
#   x <- I * exp(((1:T) - l)^2 / s^2)
#   w <- exp(-((1:T) - l)^2 / s^2)
#   o <- order(x)
#   x_sorted <- x[o]
#   w_sorted <- w[o]
#   cum_w <- cumsum(w_sorted)
#   cutoff <- tau * sum(w_sorted)
#   k <- min(which(cum_w >= cutoff))
#   return(x_sorted[k]) 
# }
# 
# # ---------------------------------------------------------------------------
# # Quantile-Based Amplitude Estimation
# # ---------------------------------------------------------------------------
# compute_a <- function(I, l, s, tau) {
#   med_weighted_tau(I, l, s, tau)
# }
# 
# # ---------------------------------------------------------------------------
# # Quantile Loss objective
# # ---------------------------------------------------------------------------
# objective_function <- function(I, a, l, s, tau) {
#   t <- seq_along(I)
#   fitted <- a * exp(-((t - l)^2) / s^2)
#   u <- I - fitted
#   mean(quantile_loss(u, tau))
# }
# 
# # ---------------------------------------------------------------------------
# # Optimize s for a given l
# # ---------------------------------------------------------------------------
# search_best_s <- function(I, l, tau) {
#   f_obj <- function(s) {
#     a_hat <- compute_a(I, l, s, tau)
#     objective_function(I, a_hat, l, s, tau)
#   }
#   optimize(f_obj, interval = c(1, 20), tol = 1e-4)$minimum
# }
# 
# # ---------------------------------------------------------------------------
# # Main simulation function
# # ---------------------------------------------------------------------------
# simulate_quantile_objective <- function(T, A, L, S, h, tau, snr = Inf, z = 0, seed = NULL, l_range = NULL, force_l_gt_T = FALSE) {
#   t <- 1:T
#   signal <- A * exp(-((t - L)^2) / S^2)
#   
#   if (!is.null(seed)) set.seed(seed)
#   
#   if (is.finite(snr)) {
#     b <- compute_noise_scale(A, snr, tau, z)
#     noise <- rlaplace(T, scale = b)
#     I <- signal + noise
#   } else {
#     I <- signal
#   }
#   
#   if (is.null(l_range)) {
#     l_grid <- if (force_l_gt_T) seq(T + 1, h) else if (is.finite(snr)) seq(L - 4, L + 3) else seq(1, h)
#   } else {
#     l_grid <- l_range
#   }
#   
#   losses <- sapply(l_grid, function(l) {
#     s_hat <- search_best_s(I, l, tau)
#     a_hat <- compute_a(I, l, s_hat, tau)
#     objective_function(I, a_hat, l, s_hat, tau)
#   })
#   
#   data.frame(l = l_grid, loss = losses, T = T)
# }
# 
# # ---------------------------------------------------------------------------
# # Plotting utility
# # ---------------------------------------------------------------------------
# plot_quantile_objectives <- function(results, tau, A_true, L_true, S_true, h) {
#   Ts <- unique(results$T)
#   nT <- length(Ts)
#   nrow <- 2
#   ncol <- 3
#   
#   png(filename = paste0("plot_quantile_tau", tau, ".png"),
#       width = 900, height = 600, res = 120)
#   old_par <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
#   on.exit({ par(old_par); dev.off() })
#   
#   for (T_val in Ts) {
#     subset_data <- subset(results, T == T_val)
#     plot(subset_data$l, subset_data$loss, type = "l", lwd = 2, col = "blue",
#          xlab = "Peak Location Candidate l",
#          ylab = expression("Objective " * f[tau]),
#          main = paste0("T = ", T_val, ", tau = ", tau))
#     abline(v = L_true, col = "red", lty = 2)
#     grid()
#   }
# }
# 
# # ---------------------------------------------------------------------------
# # Run simulation and produce plots for tau = 0.5, 0.75, 0.95
# # ---------------------------------------------------------------------------
# A_true <- 15000
# L_true <- 22
# S_true <- 10
# h <- 30
# # Ts_free <- c(2, 3, 10, 20, 21)
# Ts_free <- c(2, 3, 10, 15, 20, 21)
# z <- 0
# snr <- 1000
# 
# for ( tau_val in c(0.5, 0.75, 0.95) ) {
#   results <- do.call(rbind, lapply(Ts_free, function(T) {
#     simulate_quantile_objective(T, A = A_true, L = L_true, S = S_true, h = h,
#                                 tau = tau_val, snr = snr, z = z, force_l_gt_T = ( tau_val == 0.5 )) 
#   }))
#   
#   plot_file <- paste0("plot_quantile_tau", tau_val, ".png")
#   plot_quantile_objectives(results, tau_val, A_true, L_true, S_true, h)
#   
#   cat("Plot for tau =", tau_val, "generated and saved to", plot_file, "\n")
# }
# 
# 
# 
# 
# 
# 
# 
# 
# plot_quantile_objectives_superposed <- function(results, tau, A_true, L_true, S_true, h) {
#   Ts <- sort(unique(results$T))
#   z_vals <- sort(unique(results$z))
#   colors <- c("blue", "red", "darkgreen", "orange", "purple") # match z_vals order
#   nrow <- 2
#   ncol <- 3
#   
#   png(filename = paste0("plot_quantile_tau", tau, "_zcompare.png"),
#       width = 900, height = 600, res = 120)
#   old_par <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
#   on.exit({ par(old_par); dev.off() })
#   
#   for (T_val in Ts) {
#     df_T <- subset(results, T == T_val)
#     x_range <- range(df_T$l)
#     y_range <- range(df_T$loss)
#     
#     plot(NA, xlim = x_range, ylim = y_range,
#          xlab = "Peak Location Candidate l",
#          ylab = expression("Objective " * f[tau]),
#          main = paste0("T = ", T_val, ", tau = ", tau))
#     abline(v = L_true, col = "gray40", lty = 2)
#     grid()
#     
#     for (i in seq_along(z_vals)) {
#       df_z <- subset(df_T, z == z_vals[i])
#       lines(df_z$l, df_z$loss, col = colors[i], lwd = 2)
#     }
#     
#     legend("topright", legend = paste("z =", z_vals),
#            col = colors[1:length(z_vals)], lwd = 2, bty = "o", cex = 0.8)
#   }
# }
# 
# 
# z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
# 
# for (tau_val in c(0.5, 0.75, 0.95)) {
#   results_all <- do.call(rbind, lapply(z_vals, function(z) {
#     res_z <- do.call(rbind, lapply(Ts_free, function(T) {
#       simulate_quantile_objective(T, A = A_true, L = L_true, S = S_true, h = h,
#                                   tau = tau_val, snr = snr, z = z, 
#                                   force_l_gt_T = (tau_val == 0.5))
#     }))
#     res_z$z <- z
#     res_z
#   }))
#   
#   plot_file <- paste0("plot_quantile_tau", tau_val, "_zcompare.png")
#   plot_quantile_objectives_superposed(results_all, tau_val, A_true, L_true, S_true, h)
#   cat("Plot for tau =", tau_val, "generated and saved to", plot_file, "\n")
# }


# ---------------------------------------------------------------------------
# Quantile-Based Peak Estimation with Laplace Noise (Multiple z, tau)
# ---------------------------------------------------------------------------
# Author: Said Maanan
# Refactored: 2025-06-13
# ---------------------------------------------------------------------------

# --- Core Functions ---
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

asy_laplace_variance <- function(tau, z) {
  stopifnot(tau > 0, tau < 1, abs(z) < 1)
  mu <- tau / (1 + z) - (1 - tau) / (1 - z)
  ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
  ex2 - mu^2
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

# --- Plotting ---
plot_quantile_objectives <- function(results, tau, A_true, L_true, S_true, h, superposed = FALSE) {
  Ts <- sort(unique(results$T))
  z_vals <- sort(unique(results$z))
  colors <- c("blue", "red", "darkgreen", "orange", "purple")
  nrow <- 2; ncol <- 3
  
  fname <- paste0("plot_quantile_tau", tau, if (superposed) "_zcompare" else "", ".png")
  png(filename = fname, width = 900, height = 600, res = 120)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
  on.exit({ par(old_par); dev.off() })
  
  for (T_val in Ts) {
    df_T <- subset(results, T == T_val)
    x_range <- range(df_T$l)
    y_range <- range(df_T$loss)
    
    plot(NA, xlim = x_range, ylim = y_range,
         xlab = "Peak Location Candidate l",
         ylab = expression("Objective " * f[tau]),
         main = paste0("T = ", T_val, ", tau = ", tau))
    abline(v = L_true, col = "gray40", lty = 2)
    grid()
    
    if (superposed) {
      for (i in seq_along(z_vals)) {
        df_z <- subset(df_T, z == z_vals[i])
        lines(df_z$l, df_z$loss, col = colors[i], lwd = 2)
      }
      legend("topright", legend = paste("z =", z_vals), col = colors, lwd = 2, bty = "o", cex = 0.8)
    } else {
      lines(df_T$l, df_T$loss, col = colors[1], lwd = 2)
    }
  }
}

# --- Config ---
A_true <- 15000
L_true <- 22
S_true <- 10
h <- 30
Ts_free <- c(2, 3, 10, 15, 20, 21)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
snr <- 1000
seed <- 280
tau_vals <- c(0.5, 0.75, 0.95)

# --- Main Loop ---
for (tau in tau_vals) {
  results_all <- do.call(rbind, lapply(z_vals, function(z) {
    res_z <- do.call(rbind, lapply(Ts_free, function(T) {
      simulate_quantile_objective(T, A_true, L_true, S_true, h,
                                  tau, snr, z, seed, force_l_gt_T = (tau == 0.5))
    }))
    res_z$z <- z
    res_z
  }))
  plot_quantile_objectives(results_all, tau, A_true, L_true, S_true, h, superposed = TRUE)
  cat("Plot for tau =", tau, "generated with z comparison.\n")
}

# # ---------------------------------------------------------------------------
# # Figure 1 — Quantile Objective under Corrected Asymmetric Laplace Noise
# # ---------------------------------------------------------------------------
# 
# rgen_corrected_asy_laplace <- function(n, z, b = 1) {
#   stopifnot(abs(z) < 1)
#   p_neg <- (1 + z) / 2
#   u <- runif(n)
#   x <- numeric(n)
#   idx_neg <- which(u < p_neg)
#   idx_pos <- which(u >= p_neg)
#   x[idx_neg] <- -rexp(length(idx_neg), rate = (1 + z) / b)
#   x[idx_pos] <-  rexp(length(idx_pos), rate = (1 - z) / b)
#   return(x)
# }
# 
# asy_laplace_variance <- function(z) {
#   stopifnot(abs(z) < 1)
#   b <- 1  # since the variance will be scaled later
#   p_neg <- (1 + z) / 2
#   p_pos <- (1 - z) / 2
#   rate_neg <- (1 + z) / b
#   rate_pos <- (1 - z) / b
#   var_neg <- 2 / rate_neg^2
#   var_pos <- 2 / rate_pos^2
#   mu_neg <- -1 / rate_neg
#   mu_pos <- 1 / rate_pos
#   mu <- p_neg * mu_neg + p_pos * mu_pos
#   ex2 <- p_neg * (var_neg + mu_neg^2) + p_pos * (var_pos + mu_pos^2)
#   return(ex2 - mu^2)
# }
# 
# compute_noise_scale <- function(A, snr, z) {
#   sqrt((A / snr)^2 / asy_laplace_variance(z))
# }
# 
# quantile_loss <- function(u, tau) u * (tau - (u < 0))
# 
# med_weighted_tau <- function(I, l, s, tau) {
#   T <- length(I)
#   x <- I * exp(((1:T) - l)^2 / s^2)
#   w <- exp(-((1:T) - l)^2 / s^2)
#   o <- order(x); x_sorted <- x[o]; w_sorted <- w[o]
#   cum_w <- cumsum(w_sorted)
#   k <- min(which(cum_w >= tau * sum(w_sorted)))
#   x_sorted[k]
# }
# 
# compute_a <- function(I, l, s, tau) med_weighted_tau(I, l, s, tau)
# 
# objective_function <- function(I, a, l, s, tau) {
#   t <- seq_along(I)
#   mean(quantile_loss(I - a * exp(-((t - l)^2) / s^2), tau))
# }
# 
# search_best_s <- function(I, l, tau) {
#   optimize(function(s) objective_function(I, compute_a(I, l, s, tau), l, s, tau),
#            interval = c(1, 20), tol = 1e-4)$minimum
# }
# 
# simulate_quantile_objective <- function(T, A, L, S, h, tau, snr, z, seed = NULL, force_l_gt_T = FALSE) {
#   if (!is.null(seed)) set.seed(seed)
#   signal <- A * exp(-((1:T) - L)^2 / S^2)
#   b <- compute_noise_scale(A, snr, z)
#   noise <- rgen_corrected_asy_laplace(T, z, b)
#   I <- signal + noise
#   l_grid <- if (force_l_gt_T) seq(T + 1, h) else seq(L - 4, L + 3)
#   
#   data.frame(
#     l = l_grid,
#     loss = sapply(l_grid, function(l) {
#       s <- search_best_s(I, l, tau)
#       a <- compute_a(I, l, s, tau)
#       objective_function(I, a, l, s, tau)
#     }),
#     T = T, z = z
#   )
# }
# 
# # --- Config ---
# A_true <- 15000; L_true <- 22; S_true <- 10; h <- 30
# Ts_free <- c(2, 3, 10, 20, 21)
# z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
# snr <- 1000; seed <- 280; tau_vals <- c(0.25, 0.5, 0.75, 0.95)
# 
# # --- Updated Combined Panel Plot for Figure 1 ---
# plot_quantile_all_tau_T <- function(results_all, L_true) {
#   Ts <- sort(unique(results_all$T))
#   taus <- sort(unique(results_all$tau))
#   z_vals <- sort(unique(results_all$z))
#   colors <- c("blue", "red", "darkgreen", "black", "grey40")
#   
#   nrow <- length(Ts)
#   ncol <- length(taus)
#   
#   pdf(file = "figure1_new_asy_laplace_all_tau_T.pdf", width = 2.5 * ncol, height = 2.5 * nrow)
#   old_par <- par(mfrow = c(nrow, ncol), mar = c(2, 2, 1.5, 1), oma = c(4, 4, 3, 1), bg = "white")
#   on.exit({ par(old_par); dev.off() })
#   
#   for (i in seq_along(Ts)) {
#     for (j in seq_along(taus)) {
#       T_val <- Ts[i]
#       tau_val <- taus[j]
#       df_sub <- subset(results_all, T == T_val & tau == tau_val)
#       
#       plot(NA, xlim = range(df_sub$l), ylim = range(df_sub$loss),
#            xlab = "", ylab = "", main = "", axes = FALSE)
#       axis(1); axis(2)
#       box()
#       
#       for (k in seq_along(z_vals)) {
#         df_z <- subset(df_sub, z == z_vals[k])
#         lines(df_z$l, df_z$loss, col = colors[k], lwd = 2)
#       }
#       
#       abline(v = L_true, col = "gray40", lty = 2)
#       grid()
#       
#       if (j == 1) mtext(bquote(T == .(T_val)), side = 2, line = 3, cex = 1.1)
#       if (i == 1) mtext(bquote(tau == .(tau_val)), side = 3, line = 1, cex = 1.1)
#       if (i == 1) {
#         legend("topright", legend = paste("z =", z_vals),
#                col = colors, lwd = 2, cex = 0.8)
#       }
#     }
#   }
# }
# 
# # --- Rerun simulation including tau = 0.25 ---
# 
# results_all <- do.call(rbind, lapply(tau_vals, function(tau) {
#   do.call(rbind, lapply(z_vals, function(z) {
#     do.call(rbind, lapply(Ts_free, function(T) {
#       simulate_quantile_objective(T, A_true, L_true, S_true, h,
#                                   tau, snr, z, seed, force_l_gt_T = (tau == 0.5))
#     })) -> res_z
#     res_z$z <- z
#     res_z$tau <- tau
#     res_z
#   }))
# }))
# 
# plot_quantile_all_tau_T(results_all, L_true)
# cat("Figure 1 saved to figure1_new_asy_laplace_all_tau_T.pdf\n")


# ---------------------------------------------------------------------------
# Figure 1 — Quantile Objective under Generalized Asymmetric Laplace Noise
# ---------------------------------------------------------------------------

# Generate noise from the GL_{tau, z, b} law for arbitrary tau, z, b
rgen_gal <- function(n, tau, z, b = 1) {
  stopifnot(abs(z) < 1, tau > 0, tau < 1, b > 0)
  # Probabilities for negative and positive sides
  p_neg <- tau
  p_pos <- 1 - tau
  # Draw uniform
  u <- runif(n)
  x <- numeric(n)
  # Negative side: x < 0
  idx_neg <- which(u < p_neg)
  # Positive side: x >= 0
  idx_pos <- which(u >= p_neg)
  # Negative exponential: rate = (1+z)/b, shifted to left
  if(length(idx_neg) > 0) {
    x[idx_neg] <- -rexp(length(idx_neg), rate = (1 + z) / b)
  }
  # Positive exponential: rate = (1-z)/b, shifted to right
  if(length(idx_pos) > 0) {
    x[idx_pos] <-  rexp(length(idx_pos), rate = (1 - z) / b)
  }
  return(x)
}

# Theoretical variance as in your section
gal_variance <- function(tau, z, b = 1) {
  stopifnot(abs(z) < 1, tau > 0, tau < 1, b > 0)
  v <- tau / (1 + z)^2 + (1 - tau) / (1 - z)^2 -
    0.25 * ( (1 - tau)/(1 - z) - tau/(1 + z) )^2
  return(b^2 * v)
}

# Compute b so that SNR is as specified (matches your formula)
compute_noise_scale <- function(A, snr, tau, z) {
  # Ensure SNR = A / sqrt(Var(noise))
  v <- gal_variance(tau, z, b=1)
  b <- sqrt( (A / snr)^2 / v )
  return(b)
}

# Quantile loss function (check loss)
quantile_loss <- function(u, tau) u * (tau - (u < 0))

# Weighted tau-median estimator for a
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
  mean(quantile_loss(I - a * exp(-((t - l)^2) / s^2), tau))
}

search_best_s <- function(I, l, tau) {
  optimize(function(s) objective_function(I, compute_a(I, l, s, tau), l, s, tau),
           interval = c(1, 20), tol = 1e-4)$minimum
}

simulate_quantile_objective <- function(T, A, L, S, h, tau, snr, z, seed = NULL, force_l_gt_T = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  signal <- A * exp(-((1:T) - L)^2 / S^2)
  b <- compute_noise_scale(A, snr, tau, z)
  noise <- rgen_gal(T, tau, z, b)
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

# --- Run simulation ---
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

# --- Plotting ---
plot_quantile_all_tau_T <- function(results_all, L_true) {
  Ts <- sort(unique(results_all$T))
  taus <- sort(unique(results_all$tau))
  z_vals <- sort(unique(results_all$z))
  colors <- c("blue", "red", "darkgreen", "black", "grey40")
  
  nrow <- length(Ts)
  ncol <- length(taus)
  
  pdf(file = "figure1_new_asy_laplace_all_tau_T.pdf", width = 2.5 * ncol, height = 2.5 * nrow)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(2, 2, 1.5, 1), oma = c(4, 4, 3, 1), bg = "white")
  on.exit({ par(old_par); dev.off() })
  
  for (i in seq_along(Ts)) {
    for (j in seq_along(taus)) {
      T_val <- Ts[i]
      tau_val <- taus[j]
      df_sub <- subset(results_all, T == T_val & tau == tau_val)
      
      plot(NA, xlim = range(df_sub$l), ylim = range(df_sub$loss),
           xlab = "", ylab = "", main = "", axes = FALSE)
      axis(1); axis(2)
      box()
      
      for (k in seq_along(z_vals)) {
        df_z <- subset(df_sub, z == z_vals[k])
        lines(df_z$l, df_z$loss, col = colors[k], lwd = 2)
      }
      
      abline(v = L_true, col = "gray40", lty = 2)
      grid()
      
      if (j == 1) mtext(bquote(T == .(T_val)), side = 2, line = 3, cex = 1.1)
      if (i == 1) mtext(bquote(tau == .(tau_val)), side = 3, line = 1, cex = 1.1)
      if (i == 1) {
        legend("topright", legend = paste("z =", z_vals),
               col = colors, lwd = 2, cex = 0.8)
      }
    }
  }
}

plot_quantile_all_tau_T(results_all, L_true)
cat("Figure 1 saved to figure1_new_asy_laplace_all_tau_T.pdf\n")
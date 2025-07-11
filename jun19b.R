# ---------------------------------------------------------------------------
# Figure 2 — Quantile-Based Peak Estimation with Shifted Laplace Noise
# Author: Said Maanan | Date: 2025-06-19
# ---------------------------------------------------------------------------

# ----------------- Distribution & Noise Utilities --------------------------
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

qlaplace <- function(tau, scale = 1) {
  if (tau < 0.5) scale * log(2 * tau) else -scale * log(2 * (1 - tau))
}

shifted_laplace <- function(n, tau, scale = 1) {
  rlaplace(n, scale = scale) - qlaplace(tau, scale)
}

# asy_laplace_variance <- function(tau, z) {
#   mu <- tau / (1 + z) - (1 - tau) / (1 - z)
#   ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
#   ex2 - mu^2
# }

asy_laplace_variance <- function(tau, z) {
  mean_U <- (1 - tau) / (1 - z) - tau / (1 + z)
  ex2 <- (tau / (1 + z)^2) + ((1 - tau) / (1 - z)^2)
  ex2 - 0.25 * mean_U^2
}


# ----------------- Model Estimation Core --------------------------
quantile_loss <- function(u, tau) u * (tau - (u < 0))

med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  o <- order(x)
  k <- min(which(cumsum(w[o]) >= tau * sum(w)))
  x[o][k]
}

compute_a <- function(I, l, s, tau) med_weighted_tau(I, l, s, tau)

objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  mean(quantile_loss(I - fitted, tau))
}

search_best_s <- function(I, l, tau) {
  optimize(function(s) objective_function(I, compute_a(I, l, s, tau), l, s, tau),
           interval = c(1, 20), tol = 1e-4)$minimum
}

# ----------------- Simulation Wrapper --------------------------
simulate_shifted_laplace_objective <- function(T, A, L, S, tau, z, snr, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  signal <- A * exp(-((1:T) - L)^2 / S^2)
  b <- sqrt((A / snr)^2 / asy_laplace_variance(tau, z))
  noise <- shifted_laplace(T, tau, scale = b)
  I <- signal + noise
  
  l_grid <- 18:25
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

# ----------------- Main Execution Loop --------------------------
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
seed <- 280
T_vals <- c(4, 6, 8, 9, 10)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)

results_all <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 0.95), function(tau_val) {
  do.call(rbind, lapply(z_vals, function(z) {
    do.call(rbind, lapply(T_vals, function(T) {
      simulate_shifted_laplace_objective(T, A_true, L_true, S_true, tau_val, z, snr, seed)
    })) -> res_z
    res_z$z <- z
    res_z$tau <- tau_val
    res_z
  }))
}))


# ----------------- Combined Plotting for All Tau and T --------------------------
plot_shifted_all_tau_T <- function(results_all, L_true) {
  Ts <- sort(unique(results_all$T))
  taus <- sort(unique(results_all$tau))
  z_vals <- sort(unique(results_all$z))
  colors <- c("black", "red", "blue", "darkgreen", "grey")
  nrow <- length(Ts)
  ncol <- length(taus)
  
  pdf(file = "figure2_shifted_laplace_all_tau_T.pdf", width = 2.5 * ncol, height = 2.5 * nrow)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(2, 2, 1.5, 1), oma = c(4, 4, 3, 0))
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
      
      if (j == 1) mtext(bquote(T == .(T_val)), side = 2, line = 3, cex = 1.1)
      if (i == 1) mtext(bquote(tau == .(tau_val)), side = 3, line = 0.5, cex = 1.1)
      
      if (i == 1) {
        legend("topright", legend = paste("z =", z_vals),
               col = colors, lwd = 2, cex = 0.8, bty = "o", bg = "white")
      }
      
    }
  }
}


plot_shifted_all_tau_T(results_all, L_true)
cat("Figure 2 saved to figure2_shifted_laplace_all_tau_T.png\n")

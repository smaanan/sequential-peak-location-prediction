# ---------------------------------------------------------------------------
# Figure 1 — Quantile Objective under New Asymmetric Laplace Noise
# ---------------------------------------------------------------------------

rgen_asy_laplace <- function(n, tau, z, b = 1) {
  u <- runif(n)
  ifelse(u < tau,
         -rexp(n, rate = (1 + z) / b),
         rexp(n, rate = (1 - z) / b))
}

asy_laplace_variance <- function(tau, z) {
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
  o <- order(x); x_sorted <- x[o]; w_sorted <- w[o]
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
  noise <- rgen_asy_laplace(T, tau, z, b)
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
A_true <- 15000; L_true <- 22; S_true <- 10; h <- 30
Ts_free <- c(2, 3, 10, 20, 21)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
snr <- 1000; seed <- 280; tau_vals <- c(0.5, 0.75, 0.95)

# Combined superposed panel plot for all T and tau values
plot_quantile_all_tau_T <- function(results_all, L_true) {
  Ts <- sort(unique(results_all$T))
  taus <- sort(unique(results_all$tau))
  z_vals <- sort(unique(results_all$z))
  colors <- c("blue", "red", "darkgreen", "black", "grey")
  pdf(file = "figure1_new_asy_laplace_all_tau_T.pdf", width = 7.5, height = 10)
  # png("figure1_new_asy_laplace_all_tau_T.png", width = 900, height = 1200, res = 120)
  old_par <- par(mfrow = c(length(Ts), length(taus)), mar = c(4, 4, 3, 1), bg = "white")
  on.exit({ par(old_par); dev.off() })
  
  for (T_val in Ts) {
    for (tau_val in taus) {
      df_sub <- subset(results_all, T == T_val & tau == tau_val)
      plot(NA, xlim = range(df_sub$l), ylim = range(df_sub$loss),
           xlab = "Peak Location Candidate l", ylab = expression("Objective " * f),
           main = bquote("T = " * .(T_val) * ",  " * tau * " = " * .(tau_val)))
      abline(v = L_true, col = "gray40", lty = 2)
      grid()
      for (i in seq_along(z_vals)) {
        df_z <- subset(df_sub, z == z_vals[i])
        lines(df_z$l, df_z$loss, col = colors[i], lwd = 2)
      }
      legend("topright", legend = paste("z =", z_vals), col = colors, lwd = 2, cex = 0.8)
    }
  }
}

results_all <- do.call(rbind, lapply(tau_vals, function(tau) {
  do.call(rbind, lapply(z_vals, function(z) {
    do.call(rbind, lapply(Ts_free, function(T) {
      simulate_quantile_objective(T, A_true, L_true, S_true, h,
                                  tau, snr, z, seed, force_l_gt_T = (tau == 0.5))
    })) -> res_z
    res_z$z <- z
    res_z$tau <- tau
    res_z
  }))
}))

plot_quantile_all_tau_T(results_all, L_true)
cat("Figure 1 saved to figure1_new_asy_laplace_all_tau_T.png\n")

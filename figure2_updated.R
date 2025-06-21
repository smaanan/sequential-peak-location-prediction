# ---------------------------------------------------------------------------
# Figure 2 — Shifted Asymmetric Laplace Noise, Quantile Objective
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
  o <- order(x)
  k <- min(which(cumsum(w[o]) >= tau * sum(w)))
  x[o][k]
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

simulate_objective <- function(T, A, L, S, tau, z, snr, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  signal <- A * exp(-((1:T) - L)^2 / S^2)
  b <- compute_noise_scale(A, snr, tau, z)
  noise <- rgen_asy_laplace(T, tau, z, b)
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

plot_shifted_obj_superposed <- function(results, tau, L_true) {
  Ts <- sort(unique(results$T))
  z_vals <- sort(unique(results$z))
  colors <- c("black", "red", "blue", "darkgreen", "orange")
  
  png(paste0("shifted_obj_fun_tau_", tau, "_zcompare.png"), width = 900, height = 600, res = 120)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  on.exit(dev.off())
  
  for (T_val in Ts) {
    df_T <- subset(results, T == T_val)
    plot(NA, xlim = range(df_T$l), ylim = range(df_T$loss),
         xlab = "Peak Location Candidate l", ylab = expression("Objective " * f),
         main = paste0("T = ", T_val, ", tau = ", tau))
    abline(v = L_true, col = "gray40", lty = 2)
    grid()
    
    for (i in seq_along(z_vals)) {
      lines(subset(df_T, z == z_vals[i])[, c("l", "loss")], col = colors[i], lwd = 2)
    }
    
    legend("topright", legend = paste("z =", z_vals),
           col = colors[1:length(z_vals)], lwd = 2, cex = 0.8)
  }
}

# --- Main Parameters ---
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
seed <- 280
T_vals <- c(4, 6, 8, 9, 10)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)

for (tau_val in c(0.5, 0.75, 0.95)) {
  results_all <- do.call(rbind, lapply(z_vals, function(z) {
    do.call(rbind, lapply(T_vals, function(T) {
      simulate_objective(T, A_true, L_true, S_true, tau_val, z, snr, seed)
    }))
  }))
  
  plot_shifted_obj_superposed(results_all, tau_val, L_true)
  cat("Figure 2 for tau =", tau_val, "saved.\n")
}

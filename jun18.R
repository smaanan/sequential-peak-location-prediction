simulate_quantile_objective_fixed_noise <- function(T, A, L, S, tau, I_fixed) {
  l_grid <- 18:25
  losses <- sapply(l_grid, function(l) {
    s_hat <- search_best_s(I_fixed, l, tau)
    a_hat <- compute_a(I_fixed, l, s_hat, tau)
    objective_function(I_fixed, a_hat, l, s_hat, tau)
  })
  data.frame(l = l_grid, loss = losses, T = T)
}

find_min_T_fixed_noise <- function(tau, A, L_true, S, I_fixed, T_vals, tolerance = 1) {
  for (T in T_vals) {
    I_sub <- I_fixed[1:T]  # truncate signal+noise to first T
    results <- simulate_quantile_objective_fixed_noise(T, A, L_true, S, tau, I_sub)
    min_l <- results$l[which.min(results$loss)]
    
    if (abs(min_l - L_true) <= tolerance) return(T)
  }
  return(NA)
}

# Settings
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
z <- 0
seed <- 280
T_max <- 30
T_vals <- 4:T_max
tolerance <- 1
taus <- 50:99/100

# Step 1: Generate fixed noise
set.seed(seed)
signal <- A_true * exp(-((1:T_max) - L_true)^2 / S_true^2)
b <- A_true / (sqrt(2) * snr)
noise_fixed <- rasy_laplace(T_max, tau = 0.5, z = z) * b
I_fixed <- signal + noise_fixed

# Step 2: Loop over taus
min_Ts <- sapply(taus, function(tau) {
  find_min_T_fixed_noise(tau, A_true, L_true, S_true, I_fixed, T_vals, tolerance)
})

# Step 3: Store and plot
results_fixed <- data.frame(tau = taus, min_T = min_Ts)

png(filename = "min_T_vs_tau_fixed_noise.png", width = 900, height = 600, res = 120)
plot(results_fixed$tau, results_fixed$min_T,
     type = "b", pch = 19, col = "blue",
     xlab = expression(tau),
     ylab = "Minimum T (same noise trajectory)",
     main = "Minimum T vs tau (fixed noise trajectory)")
grid()
dev.off()

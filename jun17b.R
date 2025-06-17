# ------------------------------------------------------------------------------
# (using previously defined rasy_laplace, quantile_loss, compute_a, etc.)
# ------------------------------------------------------------------------------
find_min_T <- function(tau, A, L_true, S, z, snr, seed, T_vals, tolerance = 1) {
  for (T in T_vals) {
    results <- simulate_quantile_objective_asy(T, A, L_true, S, tau, z, snr, seed)
    min_l <- results$l[which.min(results$loss)]
    
    if (abs(min_l - L_true) <= tolerance) {
      return(T)  # first T where minimum is close enough
    }
  }
  return(NA)  # if true location not recovered
}

find_best_tau <- function(taus, A, L_true, S, z, snr, seed_vals, T_vals, tolerance = 1) {
  min_vals_matrix <- sapply(taus, function(tau) {
    sapply(seed_vals, function(seed) {
      find_min_T(tau, A, L_true, S, z, snr, seed, T_vals, tolerance)
    })
  })
  
  avg_vals <- apply(min_vals_matrix, 2, function(x) mean(na.omit(x)))
  sd_vals <- apply(min_vals_matrix, 2, function(x) sd(na.omit(x)))
  
  data.frame(tau = taus, average_T = avg_vals, sd_T = sd_vals)
}

# ------------------------------------------------------------------------------
# Run pipeline
# ------------------------------------------------------------------------------

A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
z <- 0
seed_vals <- 280:289  # 10 different simulations
T_vals <- 4:30
tolerance <- 1

taus <- c(0.5, 0.75, 0.95)
taus <- 50:99/100
results <- find_best_tau(taus, A_true, L_true, S_true, z, snr, seed_vals, T_vals, tolerance)

print("Summary of minimum T required for each tau:")
print(results)

# ------------------------------------------------------------------------------
# Plotting minimum T against tau
# ---------------------------------------------------------------------------
png(filename = "min_T_vs_tau2.png", width = 900, height = 600, res = 120)

plot(results$tau, results$average_T, 
     type = "b", pch = 19, col = "blue",
     xlab = expression(tau), 
     ylab = "Average minimum T",
     main = "Average minimum T vs tau")
# arrows(results$tau, 
#        results$average_T - results$sd_T, 
#        results$tau, 
#        results$average_T + results$sd_T, 
#        angle = 90, code = 3, length = 0.05, col = "darkgrey")
grid()
dev.off()

cat("Plot 'min_T_vs_tau.png' generated in directory.\n")

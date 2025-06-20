# ---------------------------------------------------------
# SNR as a function of z (for each fixed tau)
# and as a function of tau (for each fixed z)
# ---------------------------------------------------------

# Asymmetric Laplace variance
asy_laplace_variance <- function(tau, z) {
  mu <- tau / (1 + z) - (1 - tau) / (1 - z)
  ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
  ex2 - mu^2
}

# SNR if b = 1
compute_snr <- function(A, tau, z, b = 1) {
  A / (b * sqrt(asy_laplace_variance(tau, z)))
}

# Settings
A <- 1
z_seq <- seq(-0.99, 0.99, length.out = 200)
tau_seq <- seq(0.01, 0.99, length.out = 200)
tau_vals <- c(0.5, 0.75, 0.95)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
colors <- c("black", "red", "blue", "darkgreen", "orange")

# --- Plot SNR(z) for each fixed tau
png("snr_vs_z.png", 800, 600, res = 120)
plot(NA, xlim = c(-1, 1), ylim = c(0, 2),
     xlab = "z", ylab = "SNR(z, tau)", main = "SNR(z, tau) for fixed tau")
grid()
for (i in seq_along(tau_vals)) {
  tau <- tau_vals[i]
  snr_vals <- sapply(z_seq, function(z) compute_snr(A, tau, z))
  lines(z_seq, snr_vals, col = colors[i], lwd = 2)
}
legend("topright", legend = paste("tau =", tau_vals), col = colors, lwd = 2, cex = 0.8)
dev.off()

# --- Plot SNR(tau) for each fixed z
png("snr_vs_tau.png", 800, 600, res = 120)
plot(NA, xlim = c(0, 1), ylim = c(0, 2),
     xlab = "tau", ylab = "SNR(z, tau)", main = "SNR(z, tau) for fixed z")
grid()
for (i in seq_along(z_vals)) {
  z <- z_vals[i]
  snr_vals <- sapply(tau_seq, function(tau) compute_snr(A, tau, z))
  lines(tau_seq, snr_vals, col = colors[i], lwd = 2)
}
legend("topright", legend = paste("z =", z_vals), col = colors, lwd = 2, cex = 0.8)
dev.off()

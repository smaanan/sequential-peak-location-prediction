# ---------------------------------------------------------------------------
# Generate Figures 3, 4a, and 4b for tau = 0.25, 0.5, 0.75, 0.95
# ---------------------------------------------------------------------------

# Correct GL_{tau,b} noise sampler
rGL_tau_b <- function(n, tau, b) {
  u <- runif(n)
  q <- tau
  sgn <- ifelse(u < q, -1, 1)
  abs_x <- ifelse(sgn == -1,
                  rexp(n, rate = (1 - tau) / b),
                  rexp(n, rate = tau / b))
  sgn * abs_x
}

# Correct GL variance from document
# gl_variance <- function(tau) {
#   z <- 1 - 2 * tau
#   t1 <- tau / (1 + z)^2
#   t2 <- (1 - tau) / (1 - z)^2
#   mu <- (1 - tau) / (1 - z) - tau / (1 + z)
#   t1 + t2 - 0.25 * mu^2
# }

# Variance coefficient for the tau-aligned asymmetric Laplace distribution:
#
#   Var(e) =
#     b^2 * [tau^2 + (1-tau)^2]
#         / [tau^2 (1-tau)^2]

al_variance <- function(tau) {
  (tau^2 + (1 - tau)^2) /
    (tau^2 * (1 - tau)^2)
}

# Quantile loss
quantile_loss <- function(u, tau) u * (tau - (u < 0))

# Weighted quantile amplitude
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

# Estimators
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

estimate_params <- function(I, tau, h) {
  losses <- sapply(1:h, function(l) {
    s_hat <- search_best_s(I, l, tau)
    a_hat <- compute_a(I, l, s_hat, tau)
    objective_function(I, a_hat, l, s_hat, tau)
  })
  l_hat <- which.min(losses)
  s_hat <- search_best_s(I, l_hat, tau)
  a_hat <- compute_a(I, l_hat, s_hat, tau)
  list(l_hat = l_hat, s_hat = s_hat, a_hat = a_hat)
}

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
seed <- 280
h <- 30
T_vals <- 4:h
taus <- c(0.25, 0.5, 0.75, 0.95)

pdf(file = "fig3_4a_4b_combined2.pdf", width = 8, height = 8)
par(mfrow = c(4, 3), mar = c(2, 2, 1, 1), oma = c(3, 3, 2.5, 1))  # reduce inner margins, use outer margins

# for (tau in taus) {
#   set.seed(seed)
#   signal_full <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
#   b <- sqrt((A_true / snr)^2 / al_variance(tau))
#   noise_full <- rGL_tau_b(h, tau, b)
#   I_full <- signal_full + noise_full
#   
#   l_hats <- numeric(length(T_vals))
#   a_hats <- numeric(length(T_vals))
#   s_hats <- numeric(length(T_vals))
#   
#   for (i in seq_along(T_vals)) {
#     T <- T_vals[i]
#     I_T <- I_full[1:T]
#     est <- estimate_params(I_T, tau, h)
#     l_hats[i] <- est$l_hat
#     a_hats[i] <- est$a_hat
#     s_hats[i] <- est$s_hat
#   }
#   
#   # Column 1: l_hat(T)
#   plot(T_vals, l_hats, type = "b", col = "blue", pch = 19,
#        xlab = "", ylab = "", main = "", axes = FALSE)
#   abline(h = L_true, col = "red", lty = 2)
#   grid(); box(); axis(1); axis(2)
#   if (tau == taus[1]) mtext(expression(hat(l)), side = 3, line = 1)
#   mtext(bquote(tau == .(tau)), side = 2, line = 3)
#   
#   # Column 2: a_hat(T)
#   plot(T_vals, a_hats, type = "b", col = "darkgreen", pch = 19,
#        xlab = "", ylab = "", main = "", axes = FALSE)
#   abline(h = A_true, col = "red", lty = 2)
#   grid(); box(); axis(1); axis(2)
#   if (tau == taus[1]) mtext(expression(hat(a)), side = 3, line = 1)
#   
#   # Column 3: s_hat(T)
#   plot(T_vals, s_hats, type = "b", col = "purple", pch = 19,
#        xlab = "", ylab = "", main = "", axes = FALSE)
#   abline(h = S_true, col = "red", lty = 2)
#   grid(); box(); axis(1); axis(2)
#   if (tau == taus[1]) mtext(expression(hat(s)), side = 3, line = 1)
# }

for (tau in taus) {
  set.seed(seed)
  signal_full <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
  b <- sqrt((A_true / snr)^2 / al_variance(tau))
  noise_full <- rGL_tau_b(h, tau, b)
  I_full <- signal_full + noise_full
  
  l_hats <- numeric(length(T_vals))
  a_hats <- numeric(length(T_vals))
  s_hats <- numeric(length(T_vals))
  
  for (i in seq_along(T_vals)) {
    T <- T_vals[i]
    I_T <- I_full[1:T]
    est <- estimate_params(I_T, tau, h)
    l_hats[i] <- est$l_hat
    a_hats[i] <- est$a_hat
    s_hats[i] <- est$s_hat
  }
  
  # Column 1: l_hat(T)
  plot(T_vals, l_hats, type = "o", pch = 21, col = "gray35", bg = "darkgoldenrod1",
       xlab = "", ylab = "", main = "", axes = FALSE)
  abline(h = L_true, col = "gray", lty = 2)
  grid(); box(); axis(1); axis(2)
  if (tau == taus[1]) mtext(expression(hat(l)), side = 3, line = 1)
  mtext(bquote(tau == .(tau)), side = 2, line = 3)
  
  # Column 2: a_hat(T)
  plot(T_vals, a_hats, type = "o", pch = 21, col = "gray35", bg = "darkmagenta",
       xlab = "", ylab = "", main = "", axes = FALSE)
  abline(h = A_true, col = "gray", lty = 2)
  grid(); box(); axis(1); axis(2)
  if (tau == taus[1]) mtext(expression(hat(a)), side = 3, line = 1)
  
  # Column 3: s_hat(T)
  plot(T_vals, s_hats, type = "o", pch = 21, col = "gray35", bg = "darkslategray1",
       xlab = "", ylab = "", main = "", axes = FALSE)
  abline(h = S_true, col = "gray", lty = 2)
  grid(); box(); axis(1); axis(2)
  if (tau == taus[1]) mtext(expression(hat(s)), side = 3, line = 1)
}

dev.off()
cat("Updated combined plot saved as 'fig3_4a_4b_combined.pdf'\n")

# ---------------------------------------------------------------------------
# Quantile-Based Peak Estimation with Laplace Noise
# ---------------------------------------------------------------------------
# Author: Said Maanan
# Date: 2025-06-14
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Laplace Noise Generator
# ---------------------------------------------------------------------------
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# ---------------------------------------------------------------------------
# Quantile Loss Function
# ---------------------------------------------------------------------------
quantile_loss <- function(u, tau) {
  u * ( tau - (u < 0) )
}

# ---------------------------------------------------------------------------
# Weighted Quantile (Generalized Algorithm, tau in (0,1))
# ---------------------------------------------------------------------------
med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I)
  
  # Transform signals and weights first
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  
  o <- order(x)
  x_sorted <- x[o]
  w_sorted <- w[o]
  
  cum_w <- cumsum(w_sorted)
  cutoff <- tau * sum(w_sorted)
  
  k <- min(which(cum_w >= cutoff))
  
  return(x_sorted[k]) 
}

# ---------------------------------------------------------------------------
# Quantile-Based Amplitude Estimation
# ---------------------------------------------------------------------------
compute_a <- function(I, l, s, tau) {
  med_weighted_tau(I, l, s, tau)
}

# ---------------------------------------------------------------------------
# Quantile Loss objective
# ---------------------------------------------------------------------------
objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  u <- I - fitted
  mean(quantile_loss(u, tau))
}

# ---------------------------------------------------------------------------
# Optimize s for a given l
# ---------------------------------------------------------------------------
search_best_s <- function(I, l, tau) {
  f_obj <- function(s) {
    a_hat <- compute_a(I, l, s, tau)
    objective_function(I, a_hat, l, s, tau)
  }
  optimize(f_obj, interval = c(1, 20), tol = 1e-4)$minimum
}

# ---------------------------------------------------------------------------
# Main simulation function
# ---------------------------------------------------------------------------
simulate_quantile_objective <- function(T, A, L, S, h, tau, snr = Inf, seed = NULL, l_range = NULL, force_l_gt_T = FALSE) {
  t <- 1:T
  signal <- A * exp(-((t - L)^2) / S^2)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.finite(snr)) {
    b <- A / (sqrt(2) * snr)
    noise <- rlaplace(T, scale = b)
    I <- signal + noise
  } else {
    I <- signal
  }
  
  if (is.null(l_range)) {
    l_grid <- if (force_l_gt_T) seq(T + 1, h) else if (is.finite(snr)) seq(L - 4, L + 3) else seq(1, h)
  } else {
    l_grid <- l_range
  }
  
  losses <- sapply(l_grid, function(l) {
    s_hat <- search_best_s(I, l, tau)
    a_hat <- compute_a(I, l, s_hat, tau)
    objective_function(I, a_hat, l, s_hat, tau)
  })
  
  data.frame(l = l_grid, loss = losses, T = T)
}

# ---------------------------------------------------------------------------
# Run simulation for all tau values and combine plots into a single figure
# ---------------------------------------------------------------------------
A_true <- 15000
L_true <- 22
S_true <- 10
h <- 30
Ts_free <- c(2, 3, 10, 20, 21)
tau_values <- c(0.25, 0.5, 0.75, 0.95)

# Simulate objective losses
results_all <- do.call(rbind, lapply(tau_values, function(tau_val) {
  do.call(rbind, lapply(Ts_free, function(T) {
    df <- simulate_quantile_objective(
      T, A = A_true, L = L_true, S = S_true, h = h,
      tau = tau_val, snr = Inf, force_l_gt_T = TRUE
    )
    df$tau <- tau_val
    df
  }))
}))

plot_all_quantile_objectives <- function(results_all, A_true, L_true) {
  Ts <- sort(unique(results_all$T))
  taus <- sort(unique(results_all$tau))
  nrow <- length(Ts)
  ncol <- length(taus)
  
  pdf(file = "plot_quantile_all.pdf", width = 2.5 * ncol, height = 2.5 * nrow)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(2, 2, 1.5, 1), oma = c(4, 4, 3, 0))
  on.exit({ par(old_par); dev.off() })
  
  for (i in seq_along(Ts)) {
    for (j in seq_along(taus)) {
      T_val <- Ts[i]
      tau_val <- taus[j]
      subset_data <- subset(results_all, T == T_val & tau == tau_val)
      
      plot(subset_data$l, subset_data$loss, type = "l", lwd = 2, col = "blue",
           xlab = "", ylab = "", main = "")  # <- only these are blank, ticks still show
      
      abline(v = L_true, col = "red", lty = 2)
      grid()
      
      # Add T label on the left (row titles)
      if (j == 1) mtext(bquote(T == .(T_val)), side = 2, line = 3, cex = 1.1)
      # Add tau label on top (column titles)
      if (i == 1) mtext(bquote(tau == .(tau_val)), side = 3, line = 0.5, cex = 1.1)
    }
  }
}

# Run the new plotting function
plot_all_quantile_objectives(results_all, A_true, L_true)
cat("Saved: plot_quantile_all.pdf\n")


# --- Identify tau with lowest objective value at each T ---
library(dplyr)

# For each T, find the minimal loss and corresponding tau
summary_by_tau <- results_all %>%
  group_by(T, tau) %>%
  summarise(min_loss = min(loss), .groups = "drop") %>%
  group_by(T) %>%
  slice_min(min_loss, with_ties = FALSE) %>%
  arrange(T)
#
print(summary_by_tau)

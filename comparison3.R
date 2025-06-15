library(isotone)
library(laeken)
isotone::weighted.median

laeken::weightedMedian
laeken::weightedQuantile
matrixStats::weightedMedian

weighted_quantile_v2 <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  o <- order(x)
  w <- exp(-((1:T) - l)^2 / s^2)
  wo <- w[o]
  W <- sum(w)
  
  y <- rep(0, T)
  for (t in 1:T) {
    if (sum(wo[1:t]) <= tau * W) {
      y[t] <- t
    }
  }
  return(x[o[max(y) + 1]])
}

weighted_quantile_v2 <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  W <- sum(w)
  
  cum <- 0
  idx <- 0
  for (t in 1:T) {
    cum <- cum + w[t]
    if (cum > tau * W) {
      idx <- t
      break
    }
  }
  return(x[idx])  # first observation to pass 50%
}


compute_a <- function(I, l, s, tau) {
  t <- seq_along(I)
  w <- exp(-((t - l)^2) / s^2)
  x <- I * exp(((t - l)^2) / s^2)
  return(weighted_quantile_v2(I, l, s, tau))
}

tau_vals <- c(0.5, 0.75, 0.95)

results_tau <- lapply(tau_vals, function(tau) {
  do.call(rbind, lapply(Ts_free, function(T) {
    simulate_quantile_objective(T, A = A_true, L = L_true, S = S_true, h = horizon, tau = tau, snr = Inf)
  }))
})

for (i in seq_along(tau_vals)) {
  tau <- tau_vals[i]
  results <- results_tau[[i]]
  plot_file <- paste0("plot_noisefree_tau", tau*100, ".png")
  plot_quantile_objectives(results, tau = tau, filename = plot_file)
}

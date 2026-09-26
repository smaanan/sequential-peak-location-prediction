# ---------------------------------------------------------------------------
# Figure 2 — Quantile-Based Peak Estimation with Asymmetric Laplace Noise
# Author: Said Maanan
# ---------------------------------------------------------------------------

# ----------------- Asymmetric Laplace Sampler --------------------------

# Standard asymmetric Laplace distribution aligned with the tau-quantile.
#
# Density:
#
#   f_tau,b(x) = tau(1-tau)/b * exp(-rho_tau(x)/b)
#
# with
#
#   rho_tau(x) = x * (tau - I{x < 0})
#
# The sampler is centered at the tau-quantile, which is zero.

rasymmetric_laplace <- function(n, tau, b) {
  u <- runif(n)
  
  sgn <- ifelse(u < tau, -1, 1)
  
  abs_x <- ifelse(
    sgn == -1,
    rexp(n, rate = (1 - tau) / b),
    rexp(n, rate = tau / b)
  )
  
  sgn * abs_x
}


# ----------------- Quantile Loss and Estimation --------------------------

quantile_loss <- function(u, tau) {
  u * (tau - (u < 0))
}


med_weighted_tau <- function(I, l, s, tau) {
  
  T <- length(I)
  
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  
  o <- order(x)
  
  k <- min(which(cumsum(w[o]) >= tau * sum(w)))
  
  x[o][k]
}


compute_a <- function(I, l, s, tau) {
  med_weighted_tau(I, l, s, tau)
}


objective_function <- function(I, a, l, s, tau) {
  
  t <- seq_along(I)
  
  fitted <- a * exp(-((t - l)^2) / s^2)
  
  mean(quantile_loss(I - fitted, tau))
}


search_best_s <- function(I, l, tau) {
  
  optimize(
    function(s)
      objective_function(
        I,
        compute_a(I, l, s, tau),
        l,
        s,
        tau
      ),
    interval = c(1, 20),
    tol = 1e-4
  )$minimum
}


# ----------------- Simulation Wrapper --------------------------

simulate_asymmetric_laplace_objective <- function(
    T, A, L, S, tau, snr, seed = NULL) {
  
  if (!is.null(seed))
    set.seed(seed)
  
  # True Gaussian peak signal
  signal <- A * exp(-((1:T) - L)^2 / S^2)
  
  # ----------------------------------------------------------------
  # Variance calibration
  #
  # For the tau-aligned asymmetric Laplace distribution,
  #
  # Var(e) =
  #   b^2 * [tau^2 + (1-tau)^2]
  #       / [tau^2 (1-tau)^2].
  #
  # We impose:
  #
  #   Var(e) = (A / SNR)^2.
  # ----------------------------------------------------------------
  
  vtarget <- (A / snr)^2
  
  b <- sqrt(
    vtarget *
      tau^2 * (1 - tau)^2 /
      (tau^2 + (1 - tau)^2)
  )
  
  # Equivalent expression:
  #
  # b <- (A / snr) *
  #      tau * (1 - tau) /
  #      sqrt(tau^2 + (1 - tau)^2)
  
  # Generate asymmetric Laplace noise
  noise <- rasymmetric_laplace(T, tau, b)
  
  # Observed signal
  I <- signal + noise
  
  # Candidate peak locations
  l_grid <- 18:25
  
  data.frame(
    l = l_grid,
    loss = sapply(l_grid, function(l) {
      
      s <- search_best_s(I, l, tau)
      
      a <- compute_a(I, l, s, tau)
      
      objective_function(
        I,
        a,
        l,
        s,
        tau
      )
    }),
    
    T = T,
    tau = tau
  )
}

# ----------------- Plotting Function --------------------------

plot_asymmetric_all_tau_T <- function(results_all, L_true) {
  
  Ts <- sort(unique(results_all$T))
  taus <- sort(unique(results_all$tau))
  
  nrow <- length(Ts)
  ncol <- length(taus)
  
  pdf(
    file = "figure2_asymmetric_laplace_all_tau_T.pdf",
    width = 2.5 * ncol,
    height = 2.5 * nrow
  )
  
  old_par <- par(
    mfrow = c(nrow, ncol),
    mar = c(2, 2, 1.5, 1),
    oma = c(4, 4, 3, 0)
  )
  
  on.exit({
    par(old_par)
    dev.off()
  })
  
  for (i in seq_along(Ts)) {
    for (j in seq_along(taus)) {
      T_val <- Ts[i]
      tau_val <- taus[j]
      
      subset_data <- subset(results_all, T == T_val & tau == tau_val)
      
      plot(
        subset_data$l,
        subset_data$loss,
        type = "l",
        lwd = 2,
        col = "blue",
        xlab = "",
        ylab = "",
        main = ""
      )
      
      # Grid
      grid()
      
      # True peak location
      abline(
        v = L_true,
        col = "red",
        lty = 2
      )
      
      # Add T label on the left (row titles)
      if (j == 1) {
        mtext(
          bquote(T == .(T_val)),
          side = 2,
          line = 3,
          cex = 1.1
        )
      }
      
      # Add tau label on top (column titles)
      if (i == 1) {
        mtext(
          bquote(tau == .(tau_val)),
          side = 3,
          line = 0.5,
          cex = 1.1
        )
      }
    }
  }
}


# ----------------- Main Execution --------------------------

A_true <- 1
L_true <- 22
S_true <- 10

snr <- 1000

seed <- 280

T_vals <- c(4, 6, 8, 9, 10)

tau_vals <- c(0.25, 0.5, 0.75, 0.95)


# One simulation for each (tau, T) pair
results_all <- do.call(
  rbind,
  lapply(tau_vals, function(tau_val) {
    
    do.call(
      rbind,
      lapply(T_vals, function(T) {
        
        simulate_asymmetric_laplace_objective(
          T = T,
          A = A_true,
          L = L_true,
          S = S_true,
          tau = tau_val,
          snr = snr,
          seed = seed
        )
        
      })
    )
  })
)

# Generate Figure 2
plot_asymmetric_all_tau_T(
  results_all,
  L_true
)


cat(
  "Figure 2 (asymmetric Laplace) saved to ",
  "figure2_asymmetric_laplace_all_tau_T.pdf\n",
  sep = ""
)
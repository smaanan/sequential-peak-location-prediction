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
  
  # Grayscale colors from light to dark
  # colors <- c(
  #   "gray90",
  #   "gray70",
  #   "gray50",
  #   "gray30",
  #   "black"
  # )
  
  colors <- c(
    "blue",
    "lightblue",
    "forestgreen",
    "orange",
    "red"
  )
  
  pdf(
    file = "figure2_asymmetric_laplace_all_tau_T.pdf",
    width = 8,
    height = 6
  )
  
  # old_par <- par(
  #   mfrow = c(1, 4),
  #   mar = c(4.5, 4.5, 2.5, 1),
  #   oma = c(1, 1, 3, 1)
  # )
  
  old_par <- par(
    mfrow = c(2, 2),
    mar = c(4.5, 4.5, 2.5, 1),
    oma = c(1, 1, 3, 1)
  )
  
  on.exit({
    par(old_par)
    dev.off()
  })
  
  for (tau_val in taus) {
    
    df_tau <- subset(
      results_all,
      tau == tau_val
    )
    
    x_range <- range(df_tau$l)
    y_range <- range(df_tau$loss)
    
    plot(
      NA,
      xlim = x_range,
      ylim = y_range,
      xlab = expression(l),
      ylab = expression(f[tau](T,l)),
      main = bquote(tau == .(tau_val)),
      cex.lab = 1.1,
      cex.main = 1.2
    )
    
    # Grid
    grid()
    
    # True peak location
    abline(
      v = L_true,
      col = "gray40",
      lty = 2,
      lwd = 1.5
    )
    
    # Plot one curve for each T
    for (k in seq_along(Ts)) {
      
      T_val <- Ts[k]
      
      df_T <- subset(
        df_tau,
        T == T_val
      )
      
      lines(
        df_T$l,
        df_T$loss,
        col = colors[k],
        lwd = 2
      )
    }
    
    # Legend in the first panel only
    # if (tau_val == taus[1]) {
      
      legend(
        "topright",
        legend = paste("T =", Ts),
        col = colors,
        lwd = 2,
        cex = 0.8,
        bty = "o",
        bg = "white"
      )
    # }
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
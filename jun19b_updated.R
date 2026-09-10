# ---------------------------------------------------------------------------
# Figure 2 — Quantile-Based Peak Estimation with Shifted Laplace Noise
# Author: Said Maanan
# ---------------------------------------------------------------------------


# ----------------- Distribution & Noise Utilities --------------------------

rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

qlaplace <- function(tau, scale = 1) {
  if (tau < 0.5) {
    scale * log(2 * tau)
  } else {
    -scale * log(2 * (1 - tau))
  }
}

shifted_laplace <- function(n, tau, scale = 1) {
  rlaplace(n, scale = scale) - qlaplace(tau, scale)
}


# ----------------- Model Estimation Core --------------------------

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

simulate_shifted_laplace_objective <- function(
    T, A, L, S, tau, snr, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Signal
  signal <- A * exp(-((1:T) - L)^2 / S^2)
  
  # Laplace variance = 2 b^2
  #
  # SNR = A / sqrt(Var(e))
  #
  # Therefore:
  # b = A / (sqrt(2) * SNR)
  
  b <- A / (sqrt(2) * snr)
  
  # Shifted symmetric Laplace noise
  noise <- shifted_laplace(
    T,
    tau = tau,
    scale = b
  )
  
  # Observations
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


# ----------------- Main Execution Loop --------------------------

A_true <- 1
L_true <- 22
S_true <- 10

snr <- 1000

seed <- 280

T_vals <- c(4, 6, 8, 9, 10)

tau_vals <- c(0.25, 0.5, 0.75, 0.95)


results_all <- do.call(
  rbind,
  lapply(tau_vals, function(tau_val) {
    
    do.call(
      rbind,
      lapply(T_vals, function(T) {
        
        simulate_shifted_laplace_objective(
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

# ----------------- Plotting: Four Panels, One per tau --------------------------

plot_shifted_all_tau <- function(results_all, L_true) {
  
  Ts <- sort(unique(results_all$T))
  taus <- sort(unique(results_all$tau))
  
  # Grey shades for the different T values, lightest to darkest
  greys <- grey(seq(0.9, 0.1, length.out = length(Ts)))
  
  # 2 x 2 layout: tau subplots side by side
  pdf(
    file = "figure2_shifted_laplace_all_tau_T.pdf",
    width = 8,
    height = 6
  )

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
    
    # Common x-range
    x_range <- range(df_tau$l)
    
    # Common y-range within this tau panel
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
    
    # Plot one curve for each T, lightest grey -> darkest grey
    for (k in seq_along(Ts)) {
      
      T_val <- Ts[k]
      
      df_T <- subset(
        df_tau,
        T == T_val
      )
      
      lines(
        df_T$l,
        df_T$loss,
        col = greys[k],
        lwd = 2
      )
    }
    
    # Legend in the first panel only
    if (tau_val == taus[1]) {
      
      legend(
        "topright",
        legend = paste("T =", Ts),
        col = greys,
        lwd = 2,
        cex = 0.8,
        bty = "o",
        bg = "white"
      )
    }
  }
}


# ----------------- Generate Figure --------------------------

plot_shifted_all_tau(
  results_all,
  L_true
)

cat(
  "Figure 2 saved to figure2_shifted_laplace_all_tau_T.pdf\n"
)
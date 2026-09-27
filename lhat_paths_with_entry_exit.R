# ---------------------------------------------------------------------------
# Estimation Paths + Stabilization Times under Low SNR = 40
# tau-aligned Asymmetric Laplace Noise
#
# Author: Said Maanan
# ---------------------------------------------------------------------------

# --- Required Libraries ----------------------------------------------------
library(dplyr)
library(gt)
library(tidyr)

# --- Noise Generator and Helpers -------------------------------------------

# tau-aligned asymmetric Laplace random generator
#
# P(e < 0) = tau
#
# Density:
#   f_tau,b(x) = tau(1-tau)/b * exp(-rho_tau(x)/b)
#
# where
#   rho_tau(x) = x * (tau - 1{x < 0})
#
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


# Variance factor of the tau-aligned asymmetric Laplace distribution
#
# Var(e) =
#   b^2 [tau^2 + (1-tau)^2] /
#       [tau^2 (1-tau)^2]
#
al_variance_factor <- function(tau) {
  (tau^2 + (1 - tau)^2) /
    (tau^2 * (1 - tau)^2)
}


# Compute b from the desired SNR:
#
# Var(e) = (A / SNR)^2
#
compute_noise_scale <- function(A, snr, tau) {
  
  sqrt(
    (A / snr)^2 /
      al_variance_factor(tau)
  )
}


# --- Quantile Loss ---------------------------------------------------------

quantile_loss <- function(u, tau) {
  u * (tau - (u < 0))
}


# --- Weighted tau-Quantile -------------------------------------------------

med_weighted_tau <- function(I, l, s, tau) {
  
  T <- length(I)
  
  x <- I * exp(((1:T) - l)^2 / s^2)
  
  w <- exp(-((1:T) - l)^2 / s^2)
  
  o <- order(x)
  
  cum_w <- cumsum(w[o])
  
  cutoff <- tau * sum(w)
  
  x[o][min(which(cum_w >= cutoff))]
}


compute_a <- function(I, l, s, tau) {
  med_weighted_tau(I, l, s, tau)
}


# --- Objective Function ----------------------------------------------------

objective_function <- function(I, a, l, s, tau) {
  
  t <- seq_along(I)
  
  fitted <- a * exp(-((t - l)^2) / s^2)
  
  mean(
    quantile_loss(I - fitted, tau)
  )
}


# --- Search over s ---------------------------------------------------------

search_best_s <- function(I, l, tau) {
  
  f <- function(s) {
    
    a <- compute_a(I, l, s, tau)
    
    objective_function(
      I = I,
      a = a,
      l = l,
      s = s,
      tau = tau
    )
  }
  
  optimize(
    f,
    interval = c(1, 20),
    tol = 1e-4
  )$minimum
}


# --- Parameter Estimation --------------------------------------------------

estimate_params <- function(I, tau, h) {
  
  losses <- sapply(1:h, function(l) {
    
    s <- search_best_s(I, l, tau)
    
    a <- compute_a(I, l, s, tau)
    
    objective_function(
      I = I,
      a = a,
      l = l,
      s = s,
      tau = tau
    )
  })
  
  l_hat <- which.min(losses)
  
  s_hat <- search_best_s(
    I = I,
    l = l_hat,
    tau = tau
  )
  
  a_hat <- compute_a(
    I = I,
    l = l_hat,
    s = s_hat,
    tau = tau
  )
  
  list(
    l_hat = l_hat,
    s_hat = s_hat,
    a_hat = a_hat
  )
}


# --- Stabilization-Time Function -------------------------------------------
#
# Earliest T from which all subsequent estimates remain within
# a range of width <= tol.
#
# This is a stabilization/entry time, not an exit time.

find_stable_time <- function(l_path, T_vals, tol) {
  
  for (i in seq_along(T_vals)) {
    
    future_vals <- l_path[i:length(l_path)]
    
    if (
      max(future_vals) - min(future_vals) <= tol
    ) {
      return(T_vals[i])
    }
  }
  
  return(NA)
}


# --- Main Settings ---------------------------------------------------------

A_true <- 1
L_true <- 22
S_true <- 10

snr <- 40
seed <- 280

h <- 35
T_vals <- 4:h

taus <- c(
  0.25,
  0.50,
  0.75,
  0.95
)

threshold <- 4


# --- Results Container -----------------------------------------------------

results <- data.frame()


# --- Plotting Setup --------------------------------------------------------

pdf(
  file = "fig_lhat_paths_snr40_tau.pdf",
  width = 8.5,
  height = 6.5
)

par(
  mfrow = c(2, 2),
  mar = c(4, 4.5, 3, 1),
  bg = "white"
)


# --- Main Loop -------------------------------------------------------------

for (tau in taus) {
  
  # Reproducible simulation for each tau
  set.seed(seed)
  
  signal <- A_true *
    exp(
      -((1:h) - L_true)^2 / S_true^2
    )
  
  b <- compute_noise_scale(
    A = A_true,
    snr = snr,
    tau = tau
  )
  
  noise <- rasymmetric_laplace(
    n = h,
    tau = tau,
    b = b
  )
  
  I <- signal + noise
  
  
  # Estimation path
  l_path <- sapply(
    T_vals,
    function(T) {
      estimate_params(
        I = I[1:T],
        tau = tau,
        h = h
      )$l_hat
    }
  )
  
  
  # Stabilization time
  T_entry <- find_stable_time(
    l_path = l_path,
    T_vals = T_vals,
    tol = threshold
  )
  
  
  # Store results
  results <- rbind(
    results,
    data.frame(
      tau = tau,
      T_entry = T_entry
    )
  )
  
  
  # --- Plot ---------------------------------------------------------------
  
  ylim_range <- range(
    l_path,
    L_true,
    na.rm = TRUE
  )
  
  plot(
    NA,
    xlim = range(T_vals),
    ylim = ylim_range,
    xlab = "T",
    ylab = expression(hat(l)(T)),
    main = bquote(tau == .(tau))
  )
  
  abline(
    h = L_true,
    col = "gray40",
    lty = 2
  )
  
  grid()
  
  lines(
    T_vals,
    l_path,
    type = "b",
    pch = 19,
    lwd = 2
  )
}


dev.off()

cat(
  "Saved plot as 'fig_lhat_paths_snr40_tau.pdf'\n"
)


# --- Create Entry-Time Table -----------------------------------------------

table_out <- results %>%
  mutate(
    tau = paste0("tau = ", tau)
  ) %>%
  arrange(tau) %>%
  pivot_wider(
    names_from = tau,
    values_from = T_entry
  )


gt_table <- table_out %>%
  gt() %>%
  tab_header(
    title = md("**Stabilization Times $T_{\\mathrm{entry}}$**"),
    subtitle = "Estimated peak location under tau-aligned asymmetric Laplace noise, SNR = 40"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 0
  ) %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )


# --- Save Table ------------------------------------------------------------

gtsave(
  gt_table,
  "table_entry_times_snr40_tau.png"
)

cat(
  "Saved table as 'table_entry_times_snr40_tau.png'\n"
)


latex_code <- as_latex(gt_table)

cat(
  as.character(latex_code),
  file = "table_entry_times_snr40_tau.tex"
)

cat(
  as.character(latex_code)
)
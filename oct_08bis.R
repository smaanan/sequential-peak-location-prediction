###############################################################################
# Study of epsilon'_1(T) and epsilon_1(T)
#
# Purpose:
#   Implement the construction requested by the colleague for SNR = 1000:
#
#       T_exit^(0)(T) = T
#       T_entry^(0)(T) = t(T)
#       epsilon'_0(T) = epsilon_0(T) = 0
#
#   and then study the first non-trivial candidate:
#
#       T_exit^(1)(T)
#       T_entry^(1)(T)
#       epsilon'_1(T)
#       epsilon_1(T)
#
#   The script verifies numerically:
#
#       T_exit^(1)(T) = T_entry^(0)(T) - 1
#
#   and
#
#       epsilon_1(T) >= epsilon'_1(T) > 0.
#
#   It produces:
#
#       1. epsilon1_structure_SNR1000.csv
#       2. epsilon1_structure_summary_SNR1000.csv
#       3. fig_epsilon1_prime_vs_epsilon1_SNR1000.pdf
#
#   One panel is produced for each tau level.
#
###############################################################################

###############################################################################
# 0. Packages
###############################################################################

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)

###############################################################################
# 1. Check that the required simulation functions exist
###############################################################################

required_functions <- c(
  "estimate_peak_curves",
  "ald_scale_from_snr",
  "rald_tau",
  "ald_variance",
  "weighted_quantile",
  "rho_tau"
)

missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]

if (length(missing_functions) > 0) {
  stop(
    paste0(
      "\nThe following required functions are not available:\n",
      paste(missing_functions, collapse = "\n"),
      "\n\nPlease source your existing simulation script first.\n"
    )
  )
}

###############################################################################
# 2. Simulation parameters
###############################################################################

A <- 1
L <- 22
S <- 10
h <- 30

SNR <- 1000

tau_values <- c(
  0.25,
  0.50,
  0.75,
  0.95
)

T_values <- 8:30

seed <- 280

s.interval <- c(1, 20)
tol <- 1e-4

###############################################################################
# 3. Helper function:
#    reconstruct omega_l(t,T) from the estimated location trajectory
#
#    omega_l(t,T)
#       = max_{t <= t1,t2 <= T}
#           |hat l(t1) - hat l(t2)|
#
#    Since the maximum absolute pairwise difference is simply
#
#       max(hat l) - min(hat l),
#
#    we use this equivalent expression.
###############################################################################

compute_omega_vector <- function(
    hat_l,
    T_grid,
    t_min,
    T_current
) {
  
  # Estimated locations are available only on T_grid.
  available <- T_grid[T_grid <= T_current]
  
  if (length(available) == 0) {
    stop("No estimated location values are available up to T_current.")
  }
  
  t_values <- t_min:T_current
  
  omega_values <- numeric(length(t_values))
  
  for (i in seq_along(t_values)) {
    
    t_start <- t_values[i]
    
    idx <- which(
      T_grid >= t_start &
        T_grid <= T_current
    )
    
    if (length(idx) == 0) {
      omega_values[i] <- NA_real_
    } else {
      values <- hat_l[idx]
      
      omega_values[i] <-
        max(values, na.rm = TRUE) -
        min(values, na.rm = TRUE)
    }
  }
  
  data.frame(
    t = t_values,
    omega = omega_values
  )
}

###############################################################################
# 4. Extract T_entry^(0)(T)
#
#    T_exit^(0)(T) = T
#
#    T_entry^(0)(T) is the first t in the terminal zero plateau:
#
#       T_entry^(0)(T)
#       = min{t < T : omega_l(t,T) = 0}.
#
###############################################################################

compute_entry0 <- function(
    omega_df,
    T_current
) {
  
  zero_times <- omega_df$t[
    !is.na(omega_df$omega) &
      omega_df$omega == 0
  ]
  
  zero_times <- zero_times[
    zero_times < T_current
  ]
  
  if (length(zero_times) == 0) {
    return(NA_integer_)
  }
  
  min(zero_times)
}

###############################################################################
# 5. Extract T_exit^(1)(T)
#
#    The strict-decrease points are t for which
#
#       omega_l(t,T) > omega_l(t+1,T).
#
#    T_exit^(1)(T) is the largest such t.
#
###############################################################################

compute_exit1 <- function(
    omega_df,
    T_current
) {
  
  omega_df <- omega_df %>%
    filter(t < T_current)
  
  if (nrow(omega_df) == 0) {
    return(NA_integer_)
  }
  
  strict_decrease <- which(
    omega_df$omega > omega_df$omega[
      match(omega_df$t + 1, omega_df$t)
    ]
  )
  
  if (length(strict_decrease) == 0) {
    return(NA_integer_)
  }
  
  max(
    omega_df$t[strict_decrease],
    na.rm = TRUE
  )
}

###############################################################################
# 6. Compute epsilon'_1(T)
#
#    epsilon'_1(T)
#       = omega_l(T_exit^(1)(T), T)
#
###############################################################################

compute_epsilon1_prime <- function(
    omega_df,
    T_exit1
) {
  
  if (is.na(T_exit1)) {
    return(NA_real_)
  }
  
  value <- omega_df$omega[
    omega_df$t == T_exit1
  ]
  
  if (length(value) == 0) {
    return(NA_real_)
  }
  
  value[1]
}

###############################################################################
# 7. Compute Stage-2 quantities for k = 1
#
#    E_1(T)
#      = { t < T_exit^(1)(T):
#          omega_l(t,T_exit^(1)(T))
#          >= epsilon'_1(T) }
#
#    T_1^*(T)
#      = max E_1(T)
#
#    epsilon_1(T)
#      = omega_l(T_1^*(T), T_exit^(1)(T))
#
#    T_entry^(1)(T)
#      = min { t < T_exit^(1)(T):
#          omega_l(t,T_exit^(1)(T))
#          = epsilon_1(T) }
###############################################################################

compute_stage2_k1 <- function(
    hat_l,
    T_grid,
    T_exit1,
    epsilon1_prime,
    t_min
) {
  
  if (
    is.na(T_exit1) ||
    is.na(epsilon1_prime)
  ) {
    
    return(
      list(
        T_star1 = NA_integer_,
        epsilon1 = NA_real_,
        T_entry1 = NA_integer_
      )
    )
  }
  
  # Need omega_l(t, T_exit1).
  omega_exit1 <- compute_omega_vector(
    hat_l = hat_l,
    T_grid = T_grid,
    t_min = t_min,
    T_current = T_exit1
  )
  
  # Only t < T_exit1 are allowed.
  omega_exit1 <- omega_exit1 %>%
    filter(t < T_exit1)
  
  if (nrow(omega_exit1) == 0) {
    
    return(
      list(
        T_star1 = NA_integer_,
        epsilon1 = NA_real_,
        T_entry1 = NA_integer_
      )
    )
  }
  
  ###########################################################################
  # E_1(T)
  ###########################################################################
  
  E1 <- omega_exit1$t[
    omega_exit1$omega >= epsilon1_prime
  ]
  
  if (length(E1) == 0) {
    
    return(
      list(
        T_star1 = NA_integer_,
        epsilon1 = NA_real_,
        T_entry1 = NA_integer_
      )
    )
  }
  
  ###########################################################################
  # T_1^*(T)
  ###########################################################################
  
  T_star1 <- max(E1)
  
  ###########################################################################
  # epsilon_1(T)
  ###########################################################################
  
  epsilon1 <- omega_exit1$omega[
    omega_exit1$t == T_star1
  ]
  
  epsilon1 <- epsilon1[1]
  
  ###########################################################################
  # T_entry^(1)(T)
  #
  # We take the earliest t at which omega reaches epsilon_1(T).
  ###########################################################################
  
  entry_candidates <- omega_exit1$t[
    omega_exit1$omega == epsilon1
  ]
  
  if (length(entry_candidates) == 0) {
    T_entry1 <- NA_integer_
  } else {
    T_entry1 <- min(entry_candidates)
  }
  
  list(
    T_star1 = T_star1,
    epsilon1 = epsilon1,
    T_entry1 = T_entry1
  )
}

###############################################################################
# 8. Main calculation for one tau
###############################################################################

analyse_one_tau <- function(
    tau,
    A,
    L,
    S,
    h,
    SNR,
    seed,
    T_values,
    s.interval,
    tol
) {
  
  cat("\n")
  cat("============================================================\n")
  cat("tau =", tau, "\n")
  cat("SNR =", SNR, "\n")
  cat("============================================================\n")
  
  ###########################################################################
  # Run the existing peak-estimation procedure
  ###########################################################################
  
  res <- estimate_peak_curves(
    A = A,
    L = L,
    S = S,
    h = h,
    SNR = SNR,
    tau = tau,
    seed = seed,
    s.interval = s.interval,
    tol = tol
  )
  
  ###########################################################################
  # The returned trajectory corresponds to T = 4,...,h
  ###########################################################################
  
  T_grid <- res$T
  hat_l <- res$hat.l
  
  if (length(T_grid) != length(hat_l)) {
    stop(
      "T and hat.l have different lengths."
    )
  }
  
  ###########################################################################
  # Minimum available starting point
  #
  # The current estimate_peak_curves() function starts at T = 4.
  ###########################################################################
  
  t_min <- min(T_grid)
  
  ###########################################################################
  # Store results
  ###########################################################################
  
  results_tau <- vector(
    "list",
    length(T_values)
  )
  
  ###########################################################################
  # Loop over current observation horizons
  ###########################################################################
  
  for (T_current in T_values) {
    
    #########################################################################
    # omega_l(t,T)
    #########################################################################
    
    omega_df <- compute_omega_vector(
      hat_l = hat_l,
      T_grid = T_grid,
      t_min = t_min,
      T_current = T_current
    )
    
    #########################################################################
    # Stage 0
    #########################################################################
    
    T_exit0 <- T_current
    
    T_entry0 <- compute_entry0(
      omega_df = omega_df,
      T_current = T_current
    )
    
    epsilon0_prime <- 0
    epsilon0 <- 0
    
    #########################################################################
    # Stage 1
    #########################################################################
    
    T_exit1 <- compute_exit1(
      omega_df = omega_df,
      T_current = T_current
    )
    
    epsilon1_prime <- compute_epsilon1_prime(
      omega_df = omega_df,
      T_exit1 = T_exit1
    )
    
    #########################################################################
    # Stage 2 for k = 1
    #########################################################################
    
    stage2 <- compute_stage2_k1(
      hat_l = hat_l,
      T_grid = T_grid,
      T_exit1 = T_exit1,
      epsilon1_prime = epsilon1_prime,
      t_min = t_min
    )
    
    #########################################################################
    # Store
    #########################################################################
    
    results_tau[[which(T_values == T_current)]] <-
      data.frame(
        tau = tau,
        SNR = SNR,
        T = T_current,
        
        T_exit0 = T_exit0,
        T_entry0 = T_entry0,
        
        epsilon0_prime = epsilon0_prime,
        epsilon0 = epsilon0,
        
        T_exit1 = T_exit1,
        epsilon1_prime = epsilon1_prime,
        
        T_star1 = stage2$T_star1,
        epsilon1 = stage2$epsilon1,
        T_entry1 = stage2$T_entry1,
        
        exit1_minus_entry0_minus1 =
          T_exit1 - (T_entry0 - 1),
        
        epsilon1_minus_epsilon1_prime =
          stage2$epsilon1 - epsilon1_prime,
        
        exit1_relation_holds =
          !is.na(T_exit1) &&
          !is.na(T_entry0) &&
          T_exit1 == T_entry0 - 1,
        
        epsilon_positive =
          !is.na(epsilon1_prime) &&
          epsilon1_prime > 0,
        
        epsilon_ordering_holds =
          !is.na(stage2$epsilon1) &&
          !is.na(epsilon1_prime) &&
          stage2$epsilon1 >= epsilon1_prime
      )
  }
  
  ###########################################################################
  # Combine
  ###########################################################################
  
  bind_rows(results_tau)
}

###############################################################################
# 9. Run the analysis for all four quantile levels
###############################################################################

all_results <- bind_rows(
  lapply(
    tau_values,
    function(tau) {
      
      analyse_one_tau(
        tau = tau,
        A = A,
        L = L,
        S = S,
        h = h,
        SNR = SNR,
        seed = seed,
        T_values = T_values,
        s.interval = s.interval,
        tol = tol
      )
    }
  )
)

###############################################################################
# 10. Save complete numerical results
###############################################################################

write.csv(
  all_results,
  file = "epsilon1_structure_SNR1000.csv",
  row.names = FALSE
)

###############################################################################
# 11. Verification of the colleague's first claim
#
#     T_exit^(1)(T) = T_entry^(0)(T) - 1
###############################################################################

valid_exit_relation <- all_results %>%
  filter(
    !is.na(T_exit1),
    !is.na(T_entry0)
  )

cat("\n")
cat("============================================================\n")
cat("VERIFICATION 1\n")
cat("T_exit^(1)(T) = T_entry^(0)(T) - 1\n")
cat("============================================================\n")

if (nrow(valid_exit_relation) == 0) {
  
  cat("No valid observations available for this verification.\n")
  
} else {
  
  n_valid <- nrow(valid_exit_relation)
  
  n_holds <- sum(
    valid_exit_relation$exit1_relation_holds
  )
  
  n_violations <- n_valid - n_holds
  
  proportion <- n_holds / n_valid
  
  cat(
    "Number of valid cases :", n_valid, "\n"
  )
  
  cat(
    "Number satisfying relation :", n_holds, "\n"
  )
  
  cat(
    "Number of violations :", n_violations, "\n"
  )
  
  cat(
    "Proportion satisfying relation :",
    round(proportion, 4),
    "\n"
  )
  
  cat("\n")
  
  if (n_violations > 0) {
    
    cat(
      "WARNING: violations of the proposed identity were found.\n"
    )
    
    print(
      valid_exit_relation %>%
        filter(!exit1_relation_holds) %>%
        select(
          tau,
          T,
          T_entry0,
          T_exit1,
          exit1_minus_entry0_minus1
        )
    )
    
  } else {
    
    cat(
      "The identity holds for all valid cases.\n"
    )
  }
}

###############################################################################
# 12. Verification of epsilon'_1(T) > 0
###############################################################################

valid_epsilon <- all_results %>%
  filter(
    !is.na(epsilon1_prime)
  )

cat("\n")
cat("============================================================\n")
cat("VERIFICATION 2a\n")
cat("epsilon'_1(T) > 0\n")
cat("============================================================\n")

if (nrow(valid_epsilon) == 0) {
  
  cat("No valid epsilon'_1 values available.\n")
  
} else {
  
  n_valid_eps <- nrow(valid_epsilon)
  
  n_positive <- sum(
    valid_epsilon$epsilon_positive
  )
  
  n_nonpositive <- n_valid_eps - n_positive
  
  cat(
    "Number of valid cases :", n_valid_eps, "\n"
  )
  
  cat(
    "Number with epsilon'_1(T) > 0 :",
    n_positive,
    "\n"
  )
  
  cat(
    "Number with epsilon'_1(T) <= 0 :",
    n_nonpositive,
    "\n"
  )
  
  if (n_nonpositive > 0) {
    
    cat("\nCases where epsilon'_1(T) <= 0:\n")
    
    print(
      valid_epsilon %>%
        filter(!epsilon_positive) %>%
        select(
          tau,
          T,
          T_exit1,
          epsilon1_prime
        )
    )
  }
}

###############################################################################
# 13. Verification of epsilon_1(T) >= epsilon'_1(T)
###############################################################################

valid_epsilon_order <- all_results %>%
  filter(
    !is.na(epsilon1),
    !is.na(epsilon1_prime)
  )

cat("\n")
cat("============================================================\n")
cat("VERIFICATION 2b\n")
cat("epsilon_1(T) >= epsilon'_1(T)\n")
cat("============================================================\n")

if (nrow(valid_epsilon_order) == 0) {
  
  cat("No valid epsilon pairs available.\n")
  
} else {
  
  n_valid_order <- nrow(valid_epsilon_order)
  
  n_holds_order <- sum(
    valid_epsilon_order$epsilon_ordering_holds
  )
  
  n_violations_order <- n_valid_order - n_holds_order
  
  proportion_order <- n_holds_order / n_valid_order
  
  cat(
    "Number of valid cases :",
    n_valid_order,
    "\n"
  )
  
  cat(
    "Number satisfying inequality :",
    n_holds_order,
    "\n"
  )
  
  cat(
    "Number of violations :",
    n_violations_order,
    "\n"
  )
  
  cat(
    "Proportion satisfying inequality :",
    round(proportion_order, 4),
    "\n"
  )
  
  if (n_violations_order > 0) {
    
    cat(
      "\nWARNING: violations of epsilon_1 >= epsilon'_1 were found.\n"
    )
    
    print(
      valid_epsilon_order %>%
        filter(!epsilon_ordering_holds) %>%
        select(
          tau,
          T,
          T_exit1,
          T_entry1,
          epsilon1_prime,
          epsilon1,
          epsilon1_minus_epsilon1_prime
        )
    )
    
  } else {
    
    cat(
      "The inequality holds for all valid cases.\n"
    )
  }
}

###############################################################################
# 14. Compact summary by tau
###############################################################################

summary_results <- all_results %>%
  group_by(tau) %>%
  summarise(
    
    n_T = n(),
    
    n_valid_exit_relation =
      sum(
        !is.na(T_exit1) &
          !is.na(T_entry0)
      ),
    
    n_exit_relation_holds =
      sum(
        exit1_relation_holds,
        na.rm = TRUE
      ),
    
    proportion_exit_relation =
      ifelse(
        n_valid_exit_relation > 0,
        n_exit_relation_holds /
          n_valid_exit_relation,
        NA_real_
      ),
    
    n_positive_epsilon_prime =
      sum(
        epsilon_positive,
        na.rm = TRUE
      ),
    
    n_valid_epsilon_order =
      sum(
        !is.na(epsilon1) &
          !is.na(epsilon1_prime)
      ),
    
    n_epsilon_order_holds =
      sum(
        epsilon_ordering_holds,
        na.rm = TRUE
      ),
    
    proportion_epsilon_order =
      ifelse(
        n_valid_epsilon_order > 0,
        n_epsilon_order_holds /
          n_valid_epsilon_order,
        NA_real_
      ),
    
    min_epsilon_prime =
      ifelse(
        all(is.na(epsilon1_prime)),
        NA_real_,
        min(epsilon1_prime, na.rm = TRUE)
      ),
    
    max_epsilon_prime =
      ifelse(
        all(is.na(epsilon1_prime)),
        NA_real_,
        max(epsilon1_prime, na.rm = TRUE)
      ),
    
    min_epsilon =
      ifelse(
        all(is.na(epsilon1)),
        NA_real_,
        min(epsilon1, na.rm = TRUE)
      ),
    
    max_epsilon =
      ifelse(
        all(is.na(epsilon1)),
        NA_real_,
        max(epsilon1, na.rm = TRUE)
      ),
    
    .groups = "drop"
  )

###############################################################################
# 15. Save summary
###############################################################################

write.csv(
  summary_results,
  file = "epsilon1_structure_summary_SNR1000.csv",
  row.names = FALSE
)

###############################################################################
# 16. Print summary
###############################################################################

cat("\n")
cat("============================================================\n")
cat("SUMMARY BY TAU\n")
cat("============================================================\n")

print(summary_results)

###############################################################################
# 17. Prepare data for plotting
###############################################################################

plot_data_prime <- all_results %>%
  select(
    tau,
    T,
    epsilon1_prime
  ) %>%
  arrange(tau, T)

plot_data_eps <- all_results %>%
  select(
    tau,
    T,
    epsilon1
  ) %>%
  arrange(tau, T)

###############################################################################
# 18. Plot using base R
#
#    Four panels (2x2):
#
#       tau = 0.25
#       tau = 0.50
#       tau = 0.75
#       tau = 0.95
#
#    In each panel:
#
#       solid line  = epsilon'_1(T) [red]
#       dashed line = epsilon_1(T)  [blue]
#
###############################################################################

pdf(
  file = "fig_epsilon1_prime_vs_epsilon1_SNR1000.pdf",
  width = 10,
  height = 8
)

old.par <- par(
  mfrow = c(2, 2),
  mar = c(4, 4.5, 3, 1)
)

for (tau in tau_values) {
  
  data_tau_prime <- plot_data_prime %>%
    filter(tau == !!tau) %>%
    arrange(T)
  
  data_tau_eps <- plot_data_eps %>%
    filter(tau == !!tau) %>%
    arrange(T)
  
  y_min <- min(
    c(
      data_tau_prime$epsilon1_prime,
      data_tau_eps$epsilon1
    ),
    na.rm = TRUE
  )
  
  y_max <- max(
    c(
      data_tau_prime$epsilon1_prime,
      data_tau_eps$epsilon1
    ),
    na.rm = TRUE
  )
  
  y_range <- c(y_min, y_max)
  
  plot(
    NA,
    xlim = c(
      min(T_values),
      max(T_values)
    ),
    ylim = y_range,
    xlab = "T",
    ylab = expression(epsilon),
    main = bquote(tau == .(tau))
  )
  
  lines(
    data_tau_prime$T,
    data_tau_prime$epsilon1_prime,
    type = "o",
    pch = 16,
    lty = 1,
    col = "red",
    lwd = 1.5
  )
  
  lines(
    data_tau_eps$T,
    data_tau_eps$epsilon1,
    type = "o",
    pch = 16,
    lty = 2,
    col = "blue",
    lwd = 1.5
  )
  
  legend(
    "topright",
    legend = c(
      expression(epsilon[1] * "'(T)"),
      expression(epsilon[1] * "(T)")
    ),
    lty = c(1, 2),
    pch = c(16, 16),
    col = c("red", "blue"),
    bty = "o",
    bg = "white",
    cex = 0.85
  )
}

par(old.par)

dev.off()

###############################################################################
# 19. Also save PNG for quick inspection
###############################################################################

png(
  filename = "fig_epsilon1_prime_vs_epsilon1_SNR1000.png",
  width = 1000,
  height = 600,
  # units = "in",
  res = 120
)

old.par <- par(
  mfrow = c(2, 2),
  mar = c(4, 4.5, 3, 1)
)

for (tau in tau_values) {
  
  data_tau_prime <- plot_data_prime %>%
    filter(tau == !!tau) %>%
    arrange(T)
  
  data_tau_eps <- plot_data_eps %>%
    filter(tau == !!tau) %>%
    arrange(T)
  
  y_min <- min(
    c(
      data_tau_prime$epsilon1_prime,
      data_tau_eps$epsilon1
    ),
    na.rm = TRUE
  )
  
  y_max <- max(
    c(
      data_tau_prime$epsilon1_prime,
      data_tau_eps$epsilon1
    ),
    na.rm = TRUE
  )
  
  y_range <- c(y_min, y_max)
  
  plot(
    NA,
    xlim = c(
      min(T_values),
      max(T_values)
    ),
    ylim = y_range,
    xlab = "T",
    ylab = expression(epsilon),
    main = bquote(tau == .(tau))
  )
  
  lines(
    data_tau_prime$T,
    data_tau_prime$epsilon1_prime,
    type = "o",
    pch = 16,
    lty = 1,
    col = "red",
    lwd = 1.5
  )
  
  lines(
    data_tau_eps$T,
    data_tau_eps$epsilon1,
    type = "o",
    pch = 16,
    lty = 2,
    col = "blue",
    lwd = 1.5
  )
  
  legend(
    "topright",
    legend = c(
      expression(epsilon[1] * "'(T)"),
      expression(epsilon[1] * "(T)")
    ),
    lty = c(1, 2),
    pch = c(16, 16),
    col = c("red", "blue"),
    bty = "o",
    bg = "white",
    cex = 0.85
  )
}

par(old.par)

dev.off()

###############################################################################
# 20. Print confirmation
###############################################################################

cat(
  "\nPlots saved:\n",
  "  - fig_epsilon1_prime_vs_epsilon1_SNR1000.pdf\n",
  "  - fig_epsilon1_prime_vs_epsilon1_SNR1000.png\n"
)

###############################################################################
# 21. Additional diagnostic:
#     show T_entry^(0), T_exit^(1), epsilon'_1, epsilon_1
#     for T = 8,...,30
###############################################################################

diagnostic_table <- all_results %>%
  select(
    tau,
    T,
    T_entry0,
    T_exit1,
    epsilon1_prime,
    epsilon1,
    T_entry1
  )

cat("\n")
cat("============================================================\n")
cat("DETAILED DIAGNOSTIC TABLE\n")
cat("============================================================\n")

print(
  diagnostic_table,
  # n = Inf
)

###############################################################################
# 22. Final diagnostic messages
###############################################################################

cat("\n")
cat("============================================================\n")
cat("FILES CREATED\n")
cat("============================================================\n")

cat(
  "1. epsilon1_structure_SNR1000.csv\n"
)

cat(
  "2. epsilon1_structure_summary_SNR1000.csv\n"
)

cat(
  "3. fig_epsilon1_prime_vs_epsilon1_SNR1000.pdf\n"
)

cat(
  "4. fig_epsilon1_prime_vs_epsilon1_SNR1000.png\n"
)

cat("\n")
cat("Analysis completed.\n")
cat("============================================================\n")
############################################################
# COVID-19 FRANCE — QUANTILE PEAK-LOCATION ESTIMATION
# AND ENTRY / EXIT ANALYSIS
#
# Data:
#   ECDC daily COVID-19 cases
#   January 15 to May 31, 2020
#
# Quantile levels:
#   tau = 0.25, 0.50, 0.75, 0.95
#
# Observation horizons:
#   T = 50,...,89
#
# The sequential estimator produces
#
#   T -> hat{l}(T)
#
# The entry / exit analysis follows the operational
# construction of Section 7:
#
#   1. Future-tail oscillation omega_l(T,h)
#   2. Recursive Stage 1 construction of
#      (tau'_k, T_exit^(k))
#   3. Stage-2 construction of E_k
#   4. T_k^*
#   5. tau_k
#   6. T_entry^(k)
#   7. w_k
#
# IMPORTANT:
# The true peak location is NOT used in the operational
# entry / exit procedure.
#
# For the real-data application, the observed peak location
# is used only retrospectively to select, among the Stage-2
# candidates whose windows strictly contain the observed peak,
# the candidate with the smallest window width.
############################################################


############################################################
# 1. LIBRARIES
############################################################

library(dplyr)
library(tibble)
library(readr)
library(zoo)


############################################################
# 2. PARAMETERS
############################################################

taus <- c(
  0.25,
  0.50,
  0.75,
  0.95
)

# Observation horizons used in the sequential analysis
T_vals <- 50:89

# Analysis horizon: the search range for l is {1,...,h},
# with h = max(T_vals) = 89, matching the Section 7
# convention and the China and Germany scripts.
h <- max(T_vals)

# Numerical tolerances used in Section 7
STAGE1_TOLERANCE <- 1e-10
# E_TOLERANCE      <- 1e-10
# STAGE2_TOLERANCE <- 1e-10
E_TOLERANCE      <- 0
STAGE2_TOLERANCE <- 0

############################################################
# 3. DOWNLOAD AND PREPARE FRANCE COVID-19 DATA
############################################################

get_france_cases_ecdc <- function(
    file = "6YbQxT6Z",
    start_date = "2020-01-15",
    end_date   = "2020-05-31"
) {
  
  dat <- readr::read_csv(
    file,
    col_types = readr::cols(),
    progress = FALSE
  )
  
  dat %>%
    dplyr::filter(
      countriesAndTerritories == "France"
    ) %>%
    dplyr::mutate(
      date = as.Date(
        sprintf(
          "%04d-%02d-%02d",
          year,
          month,
          day
        )
      )
    ) %>%
    dplyr::filter(
      date >= as.Date(start_date),
      date <= as.Date(end_date)
    ) %>%
    dplyr::arrange(date) %>%
    dplyr::select(
      date,
      daily_cases = cases
    )
}

france <- get_france_cases_ecdc()

# Raw daily cases. Only the first h = 89 days are ever used
# by the estimator, since T ranges over T_vals = 50:89.
I_full <- france$daily_cases

stopifnot(length(I_full) >= h)


############################################################
# 4. QUANTILE LOSS
############################################################

quantile_loss <- function(
    u,
    tau
) {
  
  u * (
    tau -
      (u < 0)
  )
}


############################################################
# 5. WEIGHTED QUANTILE FOR THE AMPLITUDE
############################################################

med_weighted_tau <- function(
    I,
    l,
    s,
    tau
) {
  
  T <- length(I)
  
  x <- I *
    exp(
      ((1:T) - l)^2 /
        s^2
    )
  
  w <- exp(
    -((1:T) - l)^2 /
      s^2
  )
  
  # When l is far outside {1,...,T} and s is small,
  # x may overflow to Inf. Those entries carry
  # numerically-zero weight and do not affect the
  # weighted quantile. Drop them before ordering.
  finite <- is.finite(x) & is.finite(w)
  
  if (!any(finite)) {
    return(NA_real_)
  }
  
  x <- x[finite]
  w <- w[finite]
  
  if (any(w < 0) || sum(w) <= 0) {
    return(NA_real_)
  }
  
  o <- order(x)
  
  cum_w <- cumsum(
    w[o]
  )
  
  cutoff <- tau *
    sum(w)
  
  x[o][
    min(
      which(
        cum_w >= cutoff
      )
    )
  ]
}


compute_a <- function(
    I,
    l,
    s,
    tau
) {
  
  med_weighted_tau(
    I = I,
    l = l,
    s = s,
    tau = tau
  )
}


############################################################
# 6. QUANTILE-LOSS OBJECTIVE
############################################################

objective_function <- function(
    I,
    a,
    l,
    s,
    tau
) {
  
  t <- seq_along(I)
  
  fitted <- a *
    exp(
      -((t - l)^2) /
        s^2
    )
  
  mean(
    quantile_loss(
      I - fitted,
      tau
    )
  )
}


############################################################
# 7. OPTIMIZATION OVER THE WIDTH
############################################################

search_best_s <- function(
    I,
    l,
    tau
) {
  
  f <- function(s) {
    
    a <- compute_a(
      I = I,
      l = l,
      s = s,
      tau = tau
    )
    
    objective_function(
      I = I,
      a = a,
      l = l,
      s = s,
      tau = tau
    )
  }
  
  # France uses s in [1,40].
  optimize(
    f,
    c(1, 40),
    tol = 1e-4
  )$minimum
}


############################################################
# 8. ESTIMATE PEAK LOCATION
#
# At observation horizon T:
#
#     l in {1,...,h}
#
# where h = max(T_vals) = 89.
############################################################

estimate_params <- function(
    I,
    tau,
    h
) {
  
  T <- length(I)
  
  losses <- sapply(
    
    1:h,
    
    function(l) {
      
      s <- search_best_s(
        I = I,
        l = l,
        tau = tau
      )
      
      a <- compute_a(
        I = I,
        l = l,
        s = s,
        tau = tau
      )
      
      objective_function(
        I = I,
        a = a,
        l = l,
        s = s,
        tau = tau
      )
    }
  )
  
  l_hat <- which.min(
    losses
  )
  
  list(
    l_hat = l_hat,
    loss = min(losses)
  )
}


############################################################
# 9. SEQUENTIAL ESTIMATION PATH
############################################################

estimate_peak_path <- function(
    I_full,
    T_values,
    tau,
    h
) {
  
  hat.l <- sapply(
    
    T_values,
    
    function(T) {
      
      I <- I_full[
        1:T
      ]
      
      estimate_params(
        I = I,
        tau = tau,
        h = h
      )$l_hat
    }
  )
  
  
  data.frame(
    T = T_values,
    hat.l = hat.l
  )
}


############################################################
# 10. INTERVAL OSCILLATION
############################################################

interval_oscillation <- function(
    T.values,
    hat.l,
    T1,
    T2
) {
  
  if (
    length(T.values) !=
    length(hat.l)
  ) {
    
    stop(
      "T.values and hat.l must have the same length."
    )
  }
  
  
  if (
    length(T.values) == 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  idx <- which(
    T.values >= T1 &
      T.values <= T2
  )
  
  
  if (
    length(idx) == 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  if (
    length(idx) == 1
  ) {
    
    return(
      0
    )
  }
  
  
  max(
    abs(
      outer(
        hat.l[idx],
        hat.l[idx],
        "-"
      )
    )
  )
}


############################################################
# 11. COMPLETE FUTURE-TAIL OSCILLATION TABLE
############################################################

compute_tail_oscillation <- function(
    T.values,
    hat.l
) {
  
  h <- max(
    T.values
  )
  
  T.grid <- sort(
    unique(T.values)
  )
  
  
  omega <- sapply(
    
    T.grid,
    
    function(T) {
      
      interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = T,
        T2 = h
      )
    }
  )
  
  
  data.frame(
    T = T.grid,
    omega = omega
  )
}


############################################################
# 12. STAGE 1
#
# Recursive construction of
#
#   (tau'_k, T_exit^(k))
############################################################

stage1_tau_prime <- function(
    T.values,
    hat.l,
    tolerance = STAGE1_TOLERANCE
) {
  
  h <- max(
    T.values
  )
  
  T.grid <- sort(
    unique(T.values)
  )
  
  tail.table <- compute_tail_oscillation(
    T.values = T.values,
    hat.l = hat.l
  )
  
  
  candidates <- data.frame(
    step = 0,
    tau_prime = 0,
    T_exit = h,
    T_boundary = h
  )
  
  
  tau.current <- 0
  
  boundary.current <- h
  
  max.iterations <-
    length(T.grid) + 1
  
  
  for (
    step in seq_len(max.iterations)
  ) {
    
    earlier.T <- T.grid[
      T.grid < boundary.current
    ]
    
    
    if (
      length(earlier.T) == 0
    ) {
      
      break
    }
    
    
    possible.T <- tail.table$T[
      
      !is.na(
        tail.table$omega
      ) &
        
        abs(
          tail.table$omega -
            tau.current
        ) <=
        tolerance &
        
        tail.table$T <
        boundary.current
    ]
    
    
    ########################################################
    # CASE 1:
    #
    # The current tau' level is attained earlier.
    ########################################################
    
    if (
      length(possible.T) > 0
    ) {
      
      T.prime <- min(
        possible.T
      )
      
      previous.T <- T.grid[
        T.grid < T.prime
      ]
      
      if (
        length(previous.T) == 0
      ) {
        
        break
      }
      
      T.previous <- max(
        previous.T
      )
      
      tau.next <- interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = T.previous,
        T2 = h
      )
      
      T.exit.next <- T.previous
      
      T.boundary.next <- T.prime
      
    } else {
      
      ######################################################
      # CASE 2:
      #
      # No earlier T' satisfies
      #
      # omega_l(T',h) = tau'_current.
      ######################################################
      
      T.previous <- max(
        earlier.T
      )
      
      T.exit.next <- T.previous
      
      tau.next <- interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = T.previous,
        T2 = h
      )
      
      T.boundary.next <- T.previous
    }
    
    
    if (
      is.na(tau.next)
    ) {
      
      break
    }
    
    
    if (
      tau.next <=
      tau.current +
      tolerance
    ) {
      
      break
    }
    
    
    candidates <- rbind(
      
      candidates,
      
      data.frame(
        step = step,
        tau_prime = tau.next,
        T_exit = T.exit.next,
        T_boundary = T.boundary.next
      )
    )
    
    
    tau.current <- tau.next
    
    boundary.current <-
      T.boundary.next
  }
  
  
  ##########################################################
  # Definition check
  ##########################################################
  
  candidates$omega_at_exit <-
    
    sapply(
      
      candidates$T_exit,
      
      function(T.exit) {
        
        interval_oscillation(
          T.values = T.values,
          hat.l = hat.l,
          T1 = T.exit,
          T2 = h
        )
      }
    )
  
  
  candidates$definition_check <-
    
    abs(
      candidates$tau_prime -
        candidates$omega_at_exit
    ) <=
    tolerance
  
  
  list(
    h = h,
    tail.table = tail.table,
    tau.prime.table = candidates
  )
}


############################################################
# 13. STAGE 2
#
# For each Stage-1 candidate:
#
#   E_k
#   T_k^*
#   tau_k
#   T_entry^(k)
#   w_k
############################################################

stage2_one_candidate <- function(
    T.values,
    hat.l,
    tau.prime,
    T.exit,
    tolerance = STAGE2_TOLERANCE
) {
  
  candidate.T <- sort(
    unique(
      T.values[
        T.values < T.exit
      ]
    )
  )
  
  
  if (
    length(candidate.T) == 0
  ) {
    
    return(
      list(
        E.k = numeric(0),
        T.star = NA_real_,
        tau = NA_real_,
        T.entry = NA_real_,
        w.star = NA_real_
      )
    )
  }
  
  
  omega.values <- sapply(
    
    candidate.T,
    
    function(T.entry) {
      
      interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = T.entry,
        T2 = T.exit
      )
    }
  )
  
  
  valid <- !is.na(
    omega.values
  )
  
  T.valid <- candidate.T[
    valid
  ]
  
  omega.valid <- omega.values[
    valid
  ]
  
  
  if (
    length(T.valid) == 0
  ) {
    
    return(
      list(
        E.k = numeric(0),
        T.star = NA_real_,
        tau = NA_real_,
        T.entry = NA_real_,
        w.star = NA_real_
      )
    )
  }
  
  
  ##########################################################
  # E_k
  ##########################################################
  
  E.k <- T.valid[
    
    omega.valid >=
      tau.prime -
      tolerance
  ]
  
  
  if (
    length(E.k) == 0
  ) {
    
    return(
      list(
        E.k = numeric(0),
        T.star = NA_real_,
        tau = NA_real_,
        T.entry = NA_real_,
        w.star = NA_real_
      )
    )
  }
  
  
  ##########################################################
  # T_k^* = max E_k
  ##########################################################
  
  T.star <- max(
    E.k
  )
  
  
  ##########################################################
  # tau_k = omega(T_k^*, T_exit)
  ##########################################################
  
  index.star <- which(
    T.valid == T.star
  )[1]
  
  tau <- omega.valid[
    index.star
  ]
  
  
  ##########################################################
  # Find first occurrence of tau_k
  ##########################################################
  
  T.equal <- T.valid[
    
    abs(
      omega.valid -
        tau
    ) <=
      tolerance
  ]
  
  
  if (
    length(T.equal) == 0
  ) {
    
    return(
      list(
        E.k = E.k,
        T.star = T.star,
        tau = tau,
        T.entry = NA_real_,
        w.star = NA_real_
      )
    )
  }
  
  
  T.entry <- min(
    T.equal
  )
  
  
  ##########################################################
  # Candidate width
  ##########################################################
  
  w.star <-
    T.exit -
    T.entry
  
  
  list(
    E.k = E.k,
    T.star = T.star,
    tau = tau,
    T.entry = T.entry,
    w.star = w.star
  )
}


############################################################
# 14. APPLY STAGE 2 TO ALL STAGE-1 CANDIDATES
############################################################

stage2_all_candidates <- function(
    T.values,
    hat.l,
    stage1.result,
    tolerance = STAGE2_TOLERANCE
) {
  
  tau.prime.table <-
    stage1.result$tau.prime.table
  
  n.candidates <-
    nrow(
      tau.prime.table
    )
  
  results <- vector(
    "list",
    n.candidates
  )
  
  
  for (
    i in seq_len(n.candidates)
  ) {
    
    step <-
      tau.prime.table$step[i]
    
    tau.prime <-
      tau.prime.table$tau_prime[i]
    
    T.exit <-
      tau.prime.table$T_exit[i]
    
    
    S2 <- stage2_one_candidate(
      T.values = T.values,
      hat.l = hat.l,
      tau.prime = tau.prime,
      T.exit = T.exit,
      tolerance = tolerance
    )
    
    
    results[[i]] <- data.frame(
      
      step = step,
      
      tau_prime = tau.prime,
      
      T_exit = T.exit,
      
      E_k = paste(
        S2$E.k,
        collapse = ", "
      ),
      
      E_k_min =
        if (
          length(S2$E.k) > 0
        )
          min(S2$E.k)
      else
        NA_real_,
      
      E_k_max =
        if (
          length(S2$E.k) > 0
        )
          max(S2$E.k)
      else
        NA_real_,
      
      n_E_k =
        length(S2$E.k),
      
      E_empty =
        length(S2$E.k) == 0,
      
      T_star =
        S2$T.star,
      
      tau =
        S2$tau,
      
      T_entry =
        S2$T.entry,
      
      w_star =
        S2$w.star,
      
      stringsAsFactors = FALSE
    )
  }
  
  
  do.call(
    rbind,
    results
  )
}


############################################################
# 15. COMPLETE ENTRY / EXIT PROCEDURE
############################################################

run_entry_exit <- function(
    T.values,
    hat.l,
    tolerance.stage1 = STAGE1_TOLERANCE,
    tolerance.stage2 = STAGE2_TOLERANCE
) {
  
  S1 <- stage1_tau_prime(
    T.values = T.values,
    hat.l = hat.l,
    tolerance = tolerance.stage1
  )
  
  
  S2 <- stage2_all_candidates(
    T.values = T.values,
    hat.l = hat.l,
    stage1.result = S1,
    tolerance = tolerance.stage2
  )
  
  
  list(
    stage1 = S1,
    stage2 = S2
  )
}


############################################################
# 16. RUN FRANCE ANALYSIS FOR ALL TAU
############################################################

all_paths <- list()

all_stage1 <- list()

all_stage2 <- list()


for (
  tau in taus
) {
  
  cat(
    "\n\n========================================\n"
  )
  
  cat(
    "TAU =",
    tau,
    "\n"
  )
  
  cat(
    "========================================\n"
  )
  
  
  ##########################################################
  # Sequential peak-location path
  ##########################################################
  
  path <- estimate_peak_path(
    I_full = I_full,
    T_values = T_vals,
    tau = tau,
    h = h
  )
  
  
  all_paths[[as.character(tau)]] <- path
  
  
  ##########################################################
  # Entry / exit analysis
  ##########################################################
  
  result <- run_entry_exit(
    T.values = path$T,
    hat.l = path$hat.l
  )
  
  
  all_stage1[[as.character(tau)]] <-
    result$stage1$tau.prime.table
  
  all_stage2[[as.character(tau)]] <-
    result$stage2
  
  
  ##########################################################
  # Display estimation path
  ##########################################################
  
  cat(
    "\nSequential peak-location path:\n"
  )
  
  print(
    path
  )
  
  
  ##########################################################
  # Display Stage 1
  ##########################################################
  
  cat(
    "\nStage-1 candidates:\n"
  )
  
  print(
    result$stage1$tau.prime.table
  )
  
  
  ##########################################################
  # Display Stage 2
  ##########################################################
  
  cat(
    "\nStage-2 candidates:\n"
  )
  
  print(
    result$stage2
  )
}


############################################################
# 17. COMBINE STAGE-2 RESULTS
############################################################

france_stage2_table <- bind_rows(
  
  lapply(
    
    names(all_stage2),
    
    function(tau.name) {
      
      data <- all_stage2[[tau.name]]
      
      data$tau <- as.numeric(
        tau.name
      )
      
      data
    }
  )
) %>%
  select(
    tau,
    everything()
  )


############################################################
# 18. COMBINE STAGE-1 RESULTS
############################################################

france_stage1_table <- bind_rows(
  
  lapply(
    
    names(all_stage1),
    
    function(tau.name) {
      
      data <- all_stage1[[tau.name]]
      
      data$tau <- as.numeric(
        tau.name
      )
      
      data
    }
  )
) %>%
  select(
    tau,
    everything()
  )


############################################################
# 19. SAVE STAGE-1 AND STAGE-2 RESULTS
############################################################

write.csv(
  france_stage1_table,
  "france_section8_stage1_candidates.csv",
  row.names = FALSE
)

write.csv(
  france_stage2_table,
  "france_section8_stage2_candidates.csv",
  row.names = FALSE
)


############################################################
# 20. FINAL STAGE-2 OUTPUT
############################################################

cat(
  "\n\n========================================\n"
)

cat(
  "FRANCE — STAGE 2 ENTRY / EXIT CANDIDATES\n"
)

cat(
  "========================================\n\n"
)

print(
  france_stage2_table
)


############################################################
# 21. RETROSPECTIVE EX-POST SELECTION
#
# In the real-data application, the observed peak location
# is known from the complete infection trajectory and is
# therefore used only retrospectively.
#
# We use the strict containment criterion
#
#   T_entry^(k) < L < T_exit^(k)
#
# Among the Stage-2 candidates satisfying this criterion,
# we select the candidate with the smallest window width.
#
# This selection is NOT part of the operational procedure
# when the true peak location is unknown.
############################################################

# Observed peak location for the complete France trajectory
L_observed <- 78


select_expost_candidate <- function(
    stage2_table,
    L
) {
  
  # Keep only Stage-2 candidates with a valid entry time
  # and exit time.
  
  candidates <- stage2_table %>%
    filter(
      !is.na(T_entry),
      !is.na(T_exit),
      !is.na(w_star)
    )
  
  
  ##########################################################
  # Strict containment:
  #
  #     T_entry < L < T_exit
  ##########################################################
  
  peak_containing <- candidates %>%
    filter(
      T_entry < L,
      L < T_exit
    )
  
  
  ##########################################################
  # No candidate contains the observed peak
  ##########################################################
  
  if (
    nrow(peak_containing) == 0
  ) {
    
    return(
      data.frame(
        step = NA_integer_,
        tau_prime = NA_real_,
        T_exit = NA_real_,
        T_star = NA_real_,
        tau = NA_real_,
        T_entry = NA_real_,
        w_star = NA_real_,
        L = L,
        peak_contained = FALSE
      )
    )
  }
  
  
  ##########################################################
  # Select the narrowest peak-containing candidate.
  #
  # In case of a tie, the candidate with the smallest
  # Stage-2 step is selected.
  ##########################################################
  
  selected <- peak_containing %>%
    arrange(
      w_star,
      step
    ) %>%
    slice(1)
  
  
  data.frame(
    step = selected$step,
    tau_prime = selected$tau_prime,
    T_exit = selected$T_exit,
    T_star = selected$T_star,
    tau = selected$tau,
    T_entry = selected$T_entry,
    w_star = selected$w_star,
    L = L,
    peak_contained = TRUE
  )
}


############################################################
# 22. APPLY EX-POST SELECTION FOR ALL FOUR TAU VALUES
############################################################

france_selected_candidates <- bind_rows(
  
  lapply(
    
    names(all_stage2),
    
    function(tau.name) {
      
      tau.value <- as.numeric(
        tau.name
      )
      
      stage2.data <- all_stage2[[tau.name]]
      
      selected <- select_expost_candidate(
        stage2_table = stage2.data,
        L = L_observed
      )
      
      selected$tau_level <- tau.value
      
      selected
    }
  )
) %>%
  select(
    tau_level,
    step,
    tau_prime,
    T_exit,
    T_star,
    tau,
    T_entry,
    w_star,
    L,
    peak_contained
  ) %>%
  arrange(
    tau_level
  )


############################################################
# 23. SAVE RETROSPECTIVE SELECTION
############################################################

write.csv(
  france_selected_candidates,
  "france_section8_selected_candidates.csv",
  row.names = FALSE
)


############################################################
# 24. FINAL RETROSPECTIVE OUTPUT
############################################################

cat(
  "\n\n========================================\n"
)

cat(
  "FRANCE — RETROSPECTIVE SELECTION\n"
)

cat(
  "========================================\n\n"
)

cat(
  "Observed peak location L =",
  L_observed,
  "\n\n"
)

print(
  france_selected_candidates
)


############################################################
# 25. COMPACT MANUSCRIPT TABLE
############################################################

france_manuscript_table <- france_selected_candidates %>%
  transmute(
    tau = tau_level,
    T_entry = T_entry,
    T_exit = T_exit,
    width = w_star,
    L = L,
    peak_contained = peak_contained
  )


cat(
  "\n\n========================================\n"
)

cat(
  "FRANCE — COMPACT MANUSCRIPT TABLE\n"
)

cat(
  "========================================\n\n"
)

print(
  france_manuscript_table
)
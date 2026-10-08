############################################################
#
# ENTRY / EXIT ANALYSIS OF THE PEAK REGION
#
# CURRENT MODEL:
#
#   I(t) = A exp(-(t-L)^2/S^2) + e(t)
#
# The estimator is the current tau-quantile estimator with
# tau-aligned asymmetric Laplace noise.
#
# OPERATIONAL ENTRY / EXIT PROCEDURE:
#
# For each CURRENT observation horizon T < h:
#
#   1. Use only the sequential estimates available up to T:
#
#        hat.l(t),   t <= T.
#
#   2. Compute the retrospective oscillation
#
#        omega_l(t,T)
#        = max_{t <= t1,t2 <= T}
#          |hat.l(t1)-hat.l(t2)|.
#
#   3. Construct the candidate exit times:
#
#        T_exit^(0)(T) = T,
#
#        T_exit^(1)(T) > T_exit^(2)(T) > ...
#
#      where k >= 1 correspond to the strict-decrease points
#      of t -> omega_l(t,T), ordered in decreasing order.
#
#   4. For each candidate exit time, construct:
#
#        epsilon'_k(T)
#        E_k(T)
#        T_k^*(T)
#        epsilon_k(T)
#        T_entry^(k)(T)
#        w_k(T).
#
# The true peak location L is NOT used in the operational
# construction. It is used only EX POST to determine whether
# a candidate interval contains the true peak.
#
############################################################


############################################################
# TOLERANCES
############################################################

STAGE1_TOLERANCE <- 1e-10
E_TOLERANCE      <- 0
STAGE2_TOLERANCE <- 0


############################################################
# 1. INTERVAL OSCILLATION
#
# Computes
#
#   omega_l(T1,T2)
#   = max_{T1 <= t1,t2 <= T2}
#     |hat.l(t1)-hat.l(t2)|.
#
# This function is generic and is used both for:
#
#   - omega_l(t,T.current) in Stage 1
#   - omega_l(t,T.exit)    in Stage 2.
#
############################################################

interval_oscillation <- function(
    T.values,
    hat.l,
    T1,
    T2
) {
  
  if (length(T.values) != length(hat.l)) {
    stop("T.values and hat.l must have the same length.")
  }
  
  if (length(T.values) == 0) {
    return(NA_real_)
  }
  
  if (T1 > T2) {
    stop("T1 must be less than or equal to T2.")
  }
  
  idx <- which(
    T.values >= T1 &
      T.values <= T2
  )
  
  if (length(idx) == 0) {
    return(NA_real_)
  }
  
  if (length(idx) == 1) {
    return(0)
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
# 2. RETROSPECTIVE OSCILLATION TABLE FOR A FIXED CURRENT T
#
# For a fixed current observation horizon T.current,
# compute
#
#   t -> omega_l(t,T.current)
#
# for all available t <= T.current.
#
############################################################

compute_retrospective_oscillation <- function(
    T.values,
    hat.l,
    T.current
) {
  
  T.grid <- sort(
    unique(
      T.values[
        T.values <= T.current
      ]
    )
  )
  
  if (length(T.grid) == 0) {
    stop("No observation horizons are available up to T.current.")
  }
  
  omega <- sapply(
    T.grid,
    function(t) {
      
      interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = t,
        T2 = T.current
      )
    }
  )
  
  data.frame(
    T_current = T.current,
    t = T.grid,
    omega = omega,
    stringsAsFactors = FALSE
  )
}


############################################################
# 3. STAGE 1
#
# For a FIXED current observation horizon T.current:
#
#   omega_l(t,T.current)
#
# is nonincreasing in t.
#
# Strict-decrease point:
#
#   omega_l(t,T.current)
#       >
#   omega_l(t+1,T.current).
#
# The initialization candidate is
#
#   T_exit^(0)(T.current) = T.current.
#
# The nontrivial candidates are the strict-decrease points,
# ordered in decreasing order.
#
############################################################

stage1_exit_candidates <- function(
    T.values,
    hat.l,
    T.current,
    tolerance = STAGE1_TOLERANCE
) {
  
  T.grid <- sort(
    unique(
      T.values[
        T.values <= T.current
      ]
    )
  )
  
  if (length(T.grid) == 0) {
    stop("No observation horizons are available up to T.current.")
  }
  
  ##########################################################
  # Retrospective oscillation table
  ##########################################################
  
  retrospective.table <- compute_retrospective_oscillation(
    T.values = T.values,
    hat.l = hat.l,
    T.current = T.current
  )
  
  ##########################################################
  # Initialization candidate:
  #
  # T_exit^(0)(T.current) = T.current
  #
  ##########################################################
  
  candidates <- data.frame(
    T_current = T.current,
    step = 0,
    T_exit = T.current,
    epsilon_prime = 0,
    strict_decrease = FALSE,
    definition_check = TRUE,
    stringsAsFactors = FALSE
  )
  
  ##########################################################
  # Find strict-decrease points
  #
  # We only consider t for which t+1 is also available.
  ##########################################################
  
  earlier.t <- T.grid[T.grid < T.current]
  
  if (length(earlier.t) == 0) {
    
    strict_points <- numeric(0)
    
  } else {
    
    is.strict <- vapply(
      earlier.t,
      function(t) {
        
        omega.t <- interval_oscillation(
          T.values = T.values,
          hat.l    = hat.l,
          T1       = t,
          T2       = T.current
        )
        
        omega.t1 <- interval_oscillation(
          T.values = T.values,
          hat.l    = hat.l,
          T1       = t + 1,
          T2       = T.current
        )
        
        if (is.na(omega.t) || is.na(omega.t1)) {
          return(FALSE)
        }
        
        omega.t > omega.t1 + tolerance
      },
      logical(1)
    )
    
    strict_points <- earlier.t[is.strict]
  }
  
  ##########################################################
  # Order strict-decrease points in decreasing order
  #
  # Example:
  #
  #   t3 < t2 < t1 < T.current
  #
  # gives
  #
  #   T_exit^(1)(T.current) = t1
  #   T_exit^(2)(T.current) = t2
  #   T_exit^(3)(T.current) = t3
  #
  ##########################################################
  
  strict_points <- sort(
    strict_points,
    decreasing = TRUE
  )
  
  ##########################################################
  # Add nontrivial candidates
  ##########################################################
  
  if (length(strict_points) > 0) {
    
    strict_rows <- lapply(
      seq_along(strict_points),
      function(k) {
        
        T.exit <- strict_points[k]
        
        epsilon.prime <- interval_oscillation(
          T.values = T.values,
          hat.l = hat.l,
          T1 = T.exit,
          T2 = T.current
        )
        
        data.frame(
          T_current = T.current,
          step = k,
          T_exit = T.exit,
          epsilon_prime = epsilon.prime,
          strict_decrease = TRUE,
          definition_check = !is.na(epsilon.prime),
          stringsAsFactors = FALSE
        )
      }
    )
    
    candidates <- rbind(
      candidates,
      do.call(rbind, strict_rows)
    )
  }
  
  ##########################################################
  # Return
  ##########################################################
  
  list(
    T_current = T.current,
    retrospective.table = retrospective.table,
    strict_decrease_points = strict_points,
    exit.table = candidates
  )
}


############################################################
# 4. COMPUTE E_k(T) FOR ONE STAGE-1 CANDIDATE
#
# For a fixed current horizon T.current and candidate
# exit time T.exit:
#
#   E_k(T.current)
#   =
#   {
#       t < T.exit :
#       omega_l(t,T.exit)
#       >= epsilon'_k(T.current)
#   }.
#
############################################################

compute_E_for_candidate <- function(
    T.values,
    hat.l,
    T.current,
    epsilon.prime,
    T.exit,
    tolerance = E_TOLERANCE
) {
  
  candidate.t <- sort(
    unique(
      T.values[
        T.values < T.exit &
          T.values <= T.current
      ]
    )
  )
  
  if (length(candidate.t) == 0) {
    return(
      list(
        T_values_in_set = numeric(0),
        T_min = NA_real_,
        T_max = NA_real_,
        n_T = 0,
        empty_set = TRUE,
        omega.values = numeric(0)
      )
    )
  }
  
  ##########################################################
  # omega_l(t,T.exit)
  ##########################################################
  
  omega.values <- sapply(
    candidate.t,
    function(t) {
      
      interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = t,
        T2 = T.exit
      )
    }
  )
  
  valid <- !is.na(omega.values)
  
  if (!any(valid)) {
    return(
      list(
        T_values_in_set = numeric(0),
        T_min = NA_real_,
        T_max = NA_real_,
        n_T = 0,
        empty_set = TRUE,
        omega.values = omega.values
      )
    )
  }
  
  ##########################################################
  # E_k(T.current)
  ##########################################################
  
  in.set <-
    valid &
    omega.values >= epsilon.prime - tolerance
  
  T.values.in.set <- candidate.t[in.set]
  
  empty.set <- length(T.values.in.set) == 0
  
  if (empty.set) {
    T.min <- NA_real_
    T.max <- NA_real_
    n.T <- 0
  } else {
    T.min <- min(T.values.in.set)
    T.max <- max(T.values.in.set)
    n.T <- length(T.values.in.set)
  }
  
  list(
    T_values_in_set = T.values.in.set,
    T_min = T.min,
    T_max = T.max,
    n_T = n.T,
    empty_set = empty.set,
    omega.values = omega.values
  )
}


############################################################
# 5. APPLY E_k(T) CONSTRUCTION TO ALL STAGE-1 CANDIDATES
#
# This is retained as a separate diagnostic function.
#
############################################################

stage1_E_sets <- function(
    T.values,
    hat.l,
    stage1.result,
    Data,
    L,
    tolerance = E_TOLERANCE
) {
  
  exit.table <- stage1.result$exit.table
  T.current <- stage1.result$T_current
  
  n.candidates <- nrow(exit.table)
  results <- vector("list", n.candidates)
  
  for (i in seq_len(n.candidates)) {
    
    step <- exit.table$step[i]
    epsilon.prime <- exit.table$epsilon_prime[i]
    T.exit <- exit.table$T_exit[i]
    
    E.result <- compute_E_for_candidate(
      T.values = T.values,
      hat.l = hat.l,
      T.current = T.current,
      epsilon.prime = epsilon.prime,
      T.exit = T.exit,
      tolerance = tolerance
    )
    
    peak.in.set <-
      !is.na(L) &
      L %in% E.result$T_values_in_set
    
    results[[i]] <- data.frame(
      Data = Data,
      T_current = T.current,
      step = step,
      epsilon_prime = epsilon.prime,
      T_exit = T.exit,
      T_values_in_set = paste(
        E.result$T_values_in_set,
        collapse = ", "
      ),
      T_min = E.result$T_min,
      T_max = E.result$T_max,
      n_T = E.result$n_T,
      empty_set = E.result$empty_set,
      L = L,
      peak_in_set = peak.in.set,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, results)
}


############################################################
# 6. STAGE 2 FOR ONE CANDIDATE
#
# For a fixed current observation horizon T.current and
# candidate exit time T.exit:
#
#   STEP A:
#       E_k(T.current)
#
#   STEP B:
#       T_k^*(T.current) = max E_k(T.current)
#
#   STEP C:
#       epsilon_k(T.current)
#       =
#       omega_l(T_k^*(T.current), T.exit)
#
#   STEP D/E:
#       T_entry^(k)(T.current)
#       =
#       first t such that
#       omega_l(t,T.exit) = epsilon_k(T.current)
#
#   STEP F:
#       w_k(T.current)
#       =
#       T_exit - T_entry.
#
############################################################

stage2_one_candidate <- function(
    T.values,
    hat.l,
    T.current,
    step,
    epsilon.prime,
    T.exit,
    L = NA_real_,
    tolerance = STAGE2_TOLERANCE
) {
  
  ##########################################################
  # Candidate times strictly before T.exit
  ##########################################################
  
  candidate.t <- sort(
    unique(
      T.values[
        T.values < T.exit &
          T.values <= T.current
      ]
    )
  )
  
  if (length(candidate.t) == 0) {
    return(
      list(
        T_current = T.current,
        step = step,
        epsilon.prime = epsilon.prime,
        T.exit = T.exit,
        E.k = numeric(0),
        T.star = NA_real_,
        epsilon = NA_real_,
        T.entry = NA_real_,
        w.star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        candidate.t = numeric(0),
        omega.values = numeric(0)
      )
    )
  }
  
  ##########################################################
  # Compute omega_l(t,T.exit)
  ##########################################################
  
  omega.values <- sapply(
    candidate.t,
    function(t) {
      
      interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = t,
        T2 = T.exit
      )
    }
  )
  
  valid <- !is.na(omega.values)
  
  T.valid <- candidate.t[valid]
  omega.valid <- omega.values[valid]
  
  if (length(T.valid) == 0) {
    return(
      list(
        T_current = T.current,
        step = step,
        epsilon.prime = epsilon.prime,
        T.exit = T.exit,
        E.k = numeric(0),
        T.star = NA_real_,
        epsilon = NA_real_,
        T.entry = NA_real_,
        w.star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        candidate.t = numeric(0),
        omega.values = numeric(0)
      )
    )
  }
  
  ##########################################################
  # STEP A:
  #
  # E_k(T.current)
  #
  ##########################################################
  
  in.E <-
    omega.valid >= epsilon.prime - tolerance
  
  E.k <- T.valid[in.E]
  
  if (length(E.k) == 0) {
    return(
      list(
        T_current = T.current,
        step = step,
        epsilon.prime = epsilon.prime,
        T.exit = T.exit,
        E.k = numeric(0),
        T.star = NA_real_,
        epsilon = NA_real_,
        T.entry = NA_real_,
        w.star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        candidate.t = T.valid,
        omega.values = omega.valid
      )
    )
  }
  
  ##########################################################
  # STEP B:
  #
  # T_k^*(T.current) = max E_k(T.current)
  #
  ##########################################################
  
  T.star <- max(E.k)
  
  ##########################################################
  # STEP C:
  #
  # epsilon_k(T.current)
  # =
  # omega_l(T_k^*(T.current), T.exit)
  #
  ##########################################################
  
  index.star <- which(T.valid == T.star)[1]
  epsilon <- omega.valid[index.star]
  
  ##########################################################
  # STEP D:
  #
  # Find all t such that
  #
  # omega_l(t,T.exit) = epsilon_k(T.current).
  #
  ##########################################################
  
  equal.epsilon <-
    abs(omega.valid - epsilon) <= tolerance
  
  T.equal <- T.valid[equal.epsilon]
  
  if (length(T.equal) == 0) {
    return(
      list(
        T_current = T.current,
        step = step,
        epsilon.prime = epsilon.prime,
        T.exit = T.exit,
        E.k = E.k,
        T.star = T.star,
        epsilon = epsilon,
        T.entry = NA_real_,
        w.star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        candidate.t = T.valid,
        omega.values = omega.valid
      )
    )
  }
  
  ##########################################################
  # STEP E:
  #
  # T_entry^(k)(T.current)
  # =
  # min T.equal
  #
  ##########################################################
  
  T.entry <- min(T.equal)
  
  ##########################################################
  # STEP F:
  #
  # w_k(T.current)
  #
  ##########################################################
  
  w.star <- T.exit - T.entry
  
  ##########################################################
  # STEP G:
  #
  # EX-POST peak containment.
  #
  # The manuscript criterion is:
  #
  #   T_entry^(k)(T) < L < T_exit^(k)(T).
  #
  ##########################################################
  
  peak.in.window <-
    !is.na(L) &
    T.entry < L &
    L < T.exit
  
  ##########################################################
  # STEP H:
  #
  # Return
  #
  ##########################################################
  
  list(
    T_current = T.current,
    step = step,
    epsilon.prime = epsilon.prime,
    T.exit = T.exit,
    E.k = E.k,
    T.star = T.star,
    epsilon = epsilon,
    T.entry = T.entry,
    w.star = w.star,
    L = L,
    peak.in.window = peak.in.window,
    candidate.t = T.valid,
    omega.values = omega.valid
  )
}


############################################################
# 7. APPLY STAGE 2 TO ALL STAGE-1 CANDIDATES
#
# For a fixed CURRENT observation horizon T.current.
#
############################################################

stage2_all_candidates <- function(
    T.values,
    hat.l,
    stage1.result,
    Data,
    L,
    tolerance = STAGE2_TOLERANCE
) {
  
  exit.table <- stage1.result$exit.table
  T.current <- stage1.result$T_current
  
  n.candidates <- nrow(exit.table)
  results <- vector("list", n.candidates)
  
  for (i in seq_len(n.candidates)) {
    
    step <- exit.table$step[i]
    epsilon.prime <- exit.table$epsilon_prime[i]
    T.exit <- exit.table$T_exit[i]
    
    S2 <- stage2_one_candidate(
      T.values = T.values,
      hat.l = hat.l,
      T.current = T.current,
      step = step,
      epsilon.prime = epsilon.prime,
      T.exit = T.exit,
      L = L,
      tolerance = tolerance
    )
    
    results[[i]] <- data.frame(
      Data = Data,
      T_current = T.current,
      k = step,
      epsilon_prime = epsilon.prime,
      T_exit = T.exit,
      E_k = paste(
        S2$E.k,
        collapse = ", "
      ),
      E_k_min =
        if (length(S2$E.k) > 0)
          min(S2$E.k)
      else
        NA_real_,
      E_k_max =
        if (length(S2$E.k) > 0)
          max(S2$E.k)
      else
        NA_real_,
      n_E_k = length(S2$E.k),
      E_empty = length(S2$E.k) == 0,
      T_star = S2$T.star,
      epsilon = S2$epsilon,
      T_entry = S2$T.entry,
      w_k = S2$w.star,
      L = L,
      peak.in.window = S2$peak.in.window,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, results)
}


############################################################
# 8. COMPLETE PROCEDURE FOR ONE CURRENT T
#
# This function performs:
#
#   Stage 1:
#       t -> omega_l(t,T.current)
#
#   Stage 2:
#       E_k(T.current), T_entry^(k)(T.current),
#       w_k(T.current)
#
############################################################

run_entry_exit_at_T <- function(
    T.values,
    hat.l,
    T.current,
    Data,
    L,
    tolerance.stage1 = STAGE1_TOLERANCE,
    tolerance.E = E_TOLERANCE,
    tolerance.stage2 = STAGE2_TOLERANCE
) {
  
  ##########################################################
  # Stage 1
  ##########################################################
  
  S1 <- stage1_exit_candidates(
    T.values = T.values,
    hat.l = hat.l,
    T.current = T.current,
    tolerance = tolerance.stage1
  )
  
  ##########################################################
  # E_k diagnostics
  ##########################################################
  
  E.table <- stage1_E_sets(
    T.values = T.values,
    hat.l = hat.l,
    stage1.result = S1,
    Data = Data,
    L = L,
    tolerance = tolerance.E
  )
  
  ##########################################################
  # Stage 2
  ##########################################################
  
  S2 <- stage2_all_candidates(
    T.values = T.values,
    hat.l = hat.l,
    stage1.result = S1,
    Data = Data,
    L = L,
    tolerance = tolerance.stage2
  )
  
  list(
    T_current = T.current,
    stage1 = S1,
    E_sets = E.table,
    stage2 = S2
  )
}


############################################################
# 9. EX-POST SIMULATION SELECTION AT A FIXED CURRENT T
#
# Select the peak-containing candidate with the smallest
# width w_k(T).
#
# This is an EX-POST simulation diagnostic because L is used
# here only to determine whether the candidate contains the
# true peak.
#
############################################################

select_final_candidate <- function(
    stage2.table,
    L
) {
  
  out.cols <- c(
    "T_current", "k", "epsilon_prime", "T_exit", "T_star",
    "epsilon", "T_entry", "w_k", "L", "peak.in.window"
  )
  
  empty_row <- function(T.current) {
    data.frame(
      T_current     = T.current,
      k             = NA_real_,
      epsilon_prime = NA_real_,
      T_exit        = NA_real_,
      T_star        = NA_real_,
      epsilon       = NA_real_,
      T_entry       = NA_real_,
      w_k           = NA_real_,
      L             = L,
      peak.in.window = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  if (is.null(stage2.table) || nrow(stage2.table) == 0) {
    return(empty_row(NA_real_))
  }
  
  T.current <- unique(stage2.table$T_current)[1]
  
  valid <-
    !is.na(stage2.table$T_entry) &
    !is.na(stage2.table$T_exit) &
    !is.na(stage2.table$w_k) &
    !is.na(stage2.table$peak.in.window) &
    stage2.table$peak.in.window
  
  if (!any(valid)) {
    return(empty_row(T.current))
  }
  
  candidates <- stage2.table[valid, , drop = FALSE]
  
  candidates <- candidates[
    order(candidates$w_k, candidates$T_entry, candidates$T_exit),
    ,
    drop = FALSE
  ]
  
  selected <- candidates[1, out.cols, drop = FALSE]
  selected$L <- L
  rownames(selected) <- NULL
  
  selected
}

############################################################
# 10. SIMULATION PARAMETERS
############################################################

A <- 1
L <- 22
S <- 10
h <- 30
seed <- 280


############################################################
# QUANTILE LEVELS
############################################################

tau.values <- c(
  0.25,
  0.50,
  0.75,
  0.95
)


############################################################
# SNR VALUES
############################################################

SNR.values <- c(
  1000,
  100,
  1
)


############################################################
# 11. RUN COMPLETE PROCEDURE
#
# IMPORTANT:
#
# For every (tau,SNR), we first obtain the complete
# sequential trajectory
#
#     T -> hat.l(T),
#
# and then, for EVERY current observation horizon
#
#     T.current = T0,...,h,
#
# we independently construct:
#
#     t -> omega_l(t,T.current)
#
# followed by the corresponding candidate exit and entry
# times.
#
############################################################

all_selected <- vector(
  "list",
  length(tau.values) *
    length(SNR.values)
)

all_diagnostics <- vector(
  "list",
  length(tau.values) *
    length(SNR.values)
)

all_E_diagnostics <- vector(
  "list",
  length(tau.values) *
    length(SNR.values)
)

counter.selected <- 1
counter.diagnostics <- 1
counter.E <- 1


for (tau.index in seq_along(tau.values)) {
  
  tau <- tau.values[tau.index]
  
  cat("\n\n##################################################\n")
  cat("# tau =", tau, "\n")
  cat("##################################################\n")
  
  
  for (SNR.index in seq_along(SNR.values)) {
    
    SNR <- SNR.values[SNR.index]
    
    cat("\n========================================\n")
    cat("tau =", tau, "\n")
    cat("SNR =", SNR, "\n")
    cat("========================================\n")
    
    
    ########################################################
    # Estimate the complete sequential trajectory
    ########################################################
    
    res <- estimate_peak_curves(
      A = A,
      L = L,
      S = S,
      h = h,
      SNR = SNR,
      tau = tau,
      seed = seed
    )
    
    T.values <- res$T
    hat.l <- res$hat.l
    
    
    ########################################################
    # Current observation horizons
    #
    # The estimator is available for these T values.
    ########################################################
    
    current.T.values <- sort(
      unique(T.values)
    )
    
    
    ########################################################
    # Storage for this (tau,SNR) pair
    ########################################################
    
    selected.at.T <- vector(
      "list",
      length(current.T.values)
    )
    
    diagnostics.at.T <- vector(
      "list",
      length(current.T.values)
    )
    
    E.at.T <- vector(
      "list",
      length(current.T.values)
    )
    
    
    ########################################################
    # LOOP OVER CURRENT OBSERVATION HORIZONS
    ########################################################
    
    for (T.current in current.T.values) {
      
      cat(
        "\n----------------------------------------\n"
      )
      
      cat(
        "Current observation horizon T =",
        T.current,
        "\n"
      )
      
      
      ######################################################
      # Run the complete operational construction at T
      ######################################################
      
      result <- run_entry_exit_at_T(
        T.values = T.values,
        hat.l = hat.l,
        T.current = T.current,
        Data = SNR,
        L = L,
        tolerance.stage1 = STAGE1_TOLERANCE,
        tolerance.E = E_TOLERANCE,
        tolerance.stage2 = STAGE2_TOLERANCE
      )
      
      
      ######################################################
      # Stage-1 diagnostics
      ######################################################
      
      cat("\nStage-1 exit candidates:\n")
      
      print(
        result$stage1$exit.table
      )
      
      
      ######################################################
      # Retrospective oscillation table
      ######################################################
      
      cat(
        "\nRetrospective oscillation table:\n"
      )
      
      print(
        result$stage1$retrospective.table
      )
      
      
      ######################################################
      # E_k(T) diagnostics
      ######################################################
      
      cat("\nE_k(T) sets:\n")
      
      print(
        result$E_sets
      )
      
      
      ######################################################
      # Stage-2 candidates
      ######################################################
      
      cat("\nStage-2 candidates:\n")
      
      print(
        result$stage2
      )
      
      
      ######################################################
      # Store E_k diagnostics
      ######################################################
      
      E.at.T[[which(current.T.values == T.current)[1]]] <-
        result$E_sets
      
      
      ######################################################
      # Store Stage-2 diagnostics
      ######################################################
      
      stage2.table <- result$stage2
      
      
      ######################################################
      # Add explicit strict-containment diagnostic
      ######################################################
      
      stage2.table$peak.in.strict <-
        !is.na(stage2.table$T_entry) &
        !is.na(stage2.table$T_exit) &
        stage2.table$T_entry < L &
        L < stage2.table$T_exit
      
      stage2.table$tau_level <- tau
      stage2.table$SNR_level <- SNR
      
      
      ######################################################
      # Add explicit width diagnostic
      ######################################################
      
      stage2.table$w_k_check <-
        stage2.table$T_exit -
        stage2.table$T_entry
      
      
      diagnostics.at.T[[
        which(current.T.values == T.current)[1]
      ]] <- stage2.table
      
      
      ######################################################
      # EX-POST selection
      ######################################################
      
      selected <- select_final_candidate(
        stage2.table = stage2.table,
        L = L
      )
      
      
      cat("\nSelected candidate at T =", T.current, ":\n")
      
      print(
        selected
      )
      
      
      ######################################################
      # Store selected candidate
      ######################################################
      
      selected$tau <- tau
      selected$SNR <- SNR
      
      selected.at.T[[
        which(current.T.values == T.current)[1]
      ]] <- selected
    }
    
    
    ########################################################
    # Combine results over current T
    ########################################################
    
    selected.pair <- do.call(
      rbind,
      selected.at.T
    )
    
    diagnostics.pair <- do.call(
      rbind,
      diagnostics.at.T
    )
    
    E.pair <- do.call(
      rbind,
      E.at.T
    )
    
    
    ########################################################
    # Store by (tau,SNR)
    ########################################################
    
    all_selected[[counter.selected]] <-
      selected.pair
    
    all_diagnostics[[counter.diagnostics]] <-
      diagnostics.pair
    
    all_E_diagnostics[[counter.E]] <-
      E.pair
    
    counter.selected <- counter.selected + 1
    counter.diagnostics <- counter.diagnostics + 1
    counter.E <- counter.E + 1
  }
}


############################################################
# 12. COMBINE ALL RESULTS
############################################################

table.entry.interpretation <- do.call(
  rbind,
  all_selected
)

rownames(
  table.entry.interpretation
) <- NULL


combined.diagnostics <- do.call(
  rbind,
  all_diagnostics
)

rownames(
  combined.diagnostics
) <- NULL


combined.E.diagnostics <- do.call(
  rbind,
  all_E_diagnostics
)

rownames(
  combined.E.diagnostics
) <- NULL


############################################################
# 13. DISPLAY ALL STAGE-2 CANDIDATES
############################################################

cat(
  "\n\n========================================\n"
)

cat(
  "ALL STAGE-2 CANDIDATES\n"
)

cat(
  "========================================\n\n"
)


print(
  combined.diagnostics[
    ,
    c(
      "tau_level",
      "SNR_level",
      "T_current",
      "k",
      "epsilon_prime",
      "T_exit",
      "T_star",
      "epsilon",
      "T_entry",
      "w_k",
      "peak.in.strict"
    )
  ]
)


############################################################
# 14. FIRST PEAK-CONTAINING k BY (tau,SNR,T_current)
#
# Two variants:
#
#   first_k_incl0:
#       smallest k >= 0 satisfying
#
#       T_entry^(k)(T) < L < T_exit^(k)(T)
#
#   first_k_excl0:
#       smallest k >= 1 satisfying the same condition.
#
############################################################

groups <- split(
  combined.diagnostics,
  interaction(
    combined.diagnostics$tau_level,
    combined.diagnostics$SNR_level,
    combined.diagnostics$T_current,
    drop = TRUE
  )
)


first_k_table <- do.call(
  rbind,
  lapply(
    groups,
    function(df) {
      
      valid.incl0 <-
        df[
          df$peak.in.strict,
          ,
          drop = FALSE
        ]
      
      valid.excl0 <-
        df[
          df$k >= 1 &
            df$peak.in.strict,
          ,
          drop = FALSE
        ]
      
      data.frame(
        tau = unique(df$tau_level),
        SNR = unique(df$SNR_level),
        T_current = unique(df$T_current),
        
        first_k_incl0 =
          if (nrow(valid.incl0) > 0)
            min(valid.incl0$k)
        else
          NA_integer_,
        
        first_k_excl0 =
          if (nrow(valid.excl0) > 0)
            min(valid.excl0$k)
        else
          NA_integer_,
        
        stringsAsFactors = FALSE
      )
    }
  )
)


rownames(first_k_table) <- NULL


first_k_table <- first_k_table[
  order(
    first_k_table$tau,
    -first_k_table$SNR,
    first_k_table$T_current
  ),
  ,
  drop = FALSE
]


rownames(first_k_table) <- NULL


cat(
  "\n\n========================================\n"
)

cat(
  "FIRST PEAK-CONTAINING k BY (tau,SNR,T_current)\n"
)

cat(
  "========================================\n\n"
)


cat(
  paste0(
    "first_k_incl0 : smallest k >= 0 such that ",
    "T_entry^(k)(T) < L < T_exit^(k)(T)\n"
  )
)


cat(
  paste0(
    "first_k_excl0 : smallest k >= 1 such that ",
    "T_entry^(k)(T) < L < T_exit^(k)(T)\n\n"
  )
)


print(
  first_k_table
)


############################################################
# 15. RELATION WITH THE TRUE PEAK
#
# This is EX POST only.
############################################################

table.entry.interpretation$relation <-
  with(
    table.entry.interpretation,
    
    ifelse(
      is.na(T_entry) |
        is.na(T_exit),
      
      NA_character_,
      
      ifelse(
        T_exit < L,
        "T_exit < L",
        
        ifelse(
          T_entry <= L &
            L <= T_exit,
          "T_entry <= L <= T_exit",
          "L < T_entry"
        )
      )
    )
  )


############################################################
# 16. INTERPRETATION OF THE SELECTED CANDIDATE
#
# This is EX POST only.
############################################################

table.entry.interpretation$interpretation <-
  with(
    table.entry.interpretation,
    
    ifelse(
      is.na(T_entry) |
        is.na(T_exit),
      
      "No valid peak-containing candidate",
      
      ifelse(
        T_exit < L,
        "Pre-peak candidate",
        
        ifelse(
          T_entry <= L &
            L <= T_exit,
          "Peak-containing candidate",
          "Post-peak candidate"
        )
      )
    )
  )


############################################################
# 17. DISPLAY FINAL SELECTED TABLE
############################################################

cat(
  "\n\n========================================\n"
)

cat(
  "FINAL ENTRY-INTERPRETATION TABLE\n"
)

cat(
  "========================================\n\n"
)


print(
  table.entry.interpretation
)


############################################################
# 18. SAVE CSV FILES
############################################################

write.csv(
  table.entry.interpretation,
  "entry-interpretation-updated-all-tau-T.csv",
  row.names = FALSE
)


write.csv(
  combined.diagnostics,
  "simulation_all_stage2_candidates-T-dependent.csv",
  row.names = FALSE
)


write.csv(
  combined.E.diagnostics,
  "simulation_all_E_candidates-T-dependent.csv",
  row.names = FALSE
)


write.csv(
  first_k_table,
  "simulation_first_k_summary-T-dependent.csv",
  row.names = FALSE
)


############################################################
# 19. LATEX-READY SELECTED TABLE
############################################################

cat(
  "\n\n========================================\n"
)

cat(
  "LATEX TABLE — SELECTED CANDIDATES\n"
)

cat(
  "========================================\n\n"
)


print(
  knitr::kable(
    table.entry.interpretation,
    format = "latex",
    booktabs = TRUE,
    digits = 6
  )
)


############################################################
# 20. LATEX-READY FIRST-k TABLE
############################################################

cat(
  "\n\n========================================\n"
)

cat(
  "LATEX TABLE — FIRST k SUMMARY\n"
)

cat(
  "========================================\n\n"
)


print(
  knitr::kable(
    first_k_table,
    format = "latex",
    booktabs = TRUE,
    digits = 6
  )
)


############################################################
# 21. OPTIONAL: DISPLAY ONE SELECTED TABLE PER TAU
############################################################

for (tau in tau.values) {
  
  cat(
    "\n\n========================================\n"
  )
  
  cat(
    "FINAL SELECTED TABLE FOR tau =",
    tau,
    "\n"
  )
  
  cat(
    "========================================\n\n"
  )
  
  
  tau.table <- table.entry.interpretation[
    table.entry.interpretation$tau == tau,
    ,
    drop = FALSE
  ]
  
  
  print(
    tau.table
  )
  
  
  cat(
    "\nLaTeX table for tau =",
    tau,
    ":\n\n"
  )
  
  
  print(
    knitr::kable(
      tau.table,
      format = "latex",
      booktabs = TRUE,
      digits = 6
    )
  )
}
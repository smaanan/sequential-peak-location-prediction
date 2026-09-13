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
# The operational entry/exit procedure uses ONLY
#
#   T -> hat.l(T)
#
# and does NOT use the true peak location L.
#
# L is used only EX POST in the simulation to determine
# whether a candidate window contains the true peak.
#
############################################################


############################################################
# TOLERANCES
############################################################
#
# These tolerances are used for numerical equality in the
# recursive construction and in the definition of E_k.
#
############################################################

STAGE1_TOLERANCE <- 1e-10
E_TOLERANCE      <- 1e-10
STAGE2_TOLERANCE <- 1e-10


############################################################
# 1. INTERVAL OSCILLATION
############################################################
#
# Computes
#
#   omega_l(T1,T2)
#     = max_{T1 <= s,t <= T2}
#         |\hat l(s)-\hat l(t)|.
#
############################################################

interval_oscillation <- function(
    T.values,
    hat.l,
    T1,
    T2
) {
  
  ##########################################################
  # Basic validity checks
  ##########################################################
  
  if (length(T.values) != length(hat.l)) {
    
    stop(
      "T.values and hat.l must have the same length."
    )
  }
  
  
  if (length(T.values) == 0) {
    
    return(NA_real_)
  }
  
  
  ##########################################################
  # Observations belonging to [T1,T2]
  ##########################################################
  
  idx <- which(
    T.values >= T1 &
      T.values <= T2
  )
  
  
  ##########################################################
  # No observations
  ##########################################################
  
  if (length(idx) == 0) {
    
    return(NA_real_)
  }
  
  
  ##########################################################
  # One observation
  ##########################################################
  
  if (length(idx) == 1) {
    
    return(0)
  }
  
  
  ##########################################################
  # Maximum pairwise absolute difference
  ##########################################################
  
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
# 2. COMPUTE COMPLETE TAIL-OSCILLATION TABLE
############################################################
#
# Computes omega_l(T,h) for every available T.
#
############################################################

compute_tail_oscillation <- function(
    T.values,
    hat.l
) {
  
  ##########################################################
  # Final observation
  ##########################################################
  
  h <- max(T.values)
  
  
  ##########################################################
  # Available T-grid
  ##########################################################
  
  T.grid <- sort(
    unique(T.values)
  )
  
  
  ##########################################################
  # Compute omega_l(T,h)
  ##########################################################
  
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
  
  
  ##########################################################
  # Return
  ##########################################################
  
  data.frame(
    
    T = T.grid,
    
    omega = omega,
    
    stringsAsFactors = FALSE
  )
}


############################################################
# 3. STAGE 1
#
# CORRECTED RECURSIVE CONSTRUCTION OF tau'
############################################################
#
# Initial candidate:
#
#   tau'_0 = 0
#   T_exit(0) = h
#
# At each step, search for the earliest T' before the
# current boundary such that
#
#   omega_l(T',h) = tau'_current.
#
# If such a T' exists, the observation immediately preceding
# T' becomes the next exit time.
#
# If no such T' exists, the observation immediately preceding
# the current boundary becomes the next exit time.
#
############################################################

stage1_tau_prime <- function(
    T.values,
    hat.l,
    tolerance = STAGE1_TOLERANCE
) {
  
  
  ##########################################################
  # Basic quantities
  ##########################################################
  
  h <- max(T.values)
  
  T.grid <- sort(
    unique(T.values)
  )
  
  
  ##########################################################
  # Compute complete tail-oscillation table once
  ##########################################################
  
  tail.table <- compute_tail_oscillation(
    
    T.values = T.values,
    
    hat.l = hat.l
  )
  
  
  ##########################################################
  # Initial candidate
  #
  # tau'_0 = 0
  # T_exit(0) = h
  ##########################################################
  
  candidates <- data.frame(
    
    step = 0,
    
    tau_prime = 0,
    
    T_exit = h,
    
    T_boundary = h,
    
    stringsAsFactors = FALSE
  )
  
  
  ##########################################################
  # Current recursive state
  ##########################################################
  
  tau.current <- 0
  
  boundary.current <- h
  
  
  ##########################################################
  # Safety limit
  ##########################################################
  
  max.iterations <- length(T.grid) + 1
  
  
  ##########################################################
  # Recursive construction
  ##########################################################
  
  for (step in seq_len(max.iterations)) {
    
    
    ########################################################
    # Observations strictly before current boundary
    ########################################################
    
    earlier.T <- T.grid[
      T.grid < boundary.current
    ]
    
    
    ########################################################
    # If no earlier observation exists, stop.
    ########################################################
    
    if (length(earlier.T) == 0) {
      
      break
    }
    
    
    ########################################################
    # Find T' < current boundary such that
    #
    # omega_l(T',h) = tau'_current
    #
    # Equality is implemented numerically with tolerance.
    ########################################################
    
    possible.T <- tail.table$T[
      
      !is.na(tail.table$omega) &
        
        abs(
          tail.table$omega -
            tau.current
        ) <= tolerance &
        
        tail.table$T <
        boundary.current
    ]
    
    
    ########################################################
    # CASE 1:
    #
    # The current tau' level is attained earlier.
    ########################################################
    
    if (length(possible.T) > 0) {
      
      
      ######################################################
      # Choose the SMALLEST such T'
      ######################################################
      
      T.prime <- min(
        possible.T
      )
      
      
      ######################################################
      # Previous available observation before T'
      ######################################################
      
      previous.T <- T.grid[
        T.grid < T.prime
      ]
      
      
      if (length(previous.T) == 0) {
        
        break
      }
      
      
      T.previous <- max(
        previous.T
      )
      
      
      ######################################################
      # Construct next tau'
      ######################################################
      
      tau.next <- interval_oscillation(
        
        T.values = T.values,
        
        hat.l = hat.l,
        
        T1 = T.previous,
        
        T2 = h
      )
      
      
      ######################################################
      # New exit time and new boundary
      ######################################################
      
      T.exit.next <- T.previous
      
      T.boundary.next <- T.prime
      
    } else {
      
      
      ######################################################
      # CASE 2:
      #
      # No earlier T' satisfies
      #
      # omega_l(T',h) = tau'_current.
      #
      # The recursion does NOT stop.
      #
      # Instead take the observation immediately preceding
      # the current boundary.
      ######################################################
      
      T.previous <- max(
        earlier.T
      )
      
      
      ######################################################
      # This becomes the next exit time
      ######################################################
      
      T.exit.next <- T.previous
      
      
      ######################################################
      # Define next tau'
      ######################################################
      
      tau.next <- interval_oscillation(
        
        T.values = T.values,
        
        hat.l = hat.l,
        
        T1 = T.previous,
        
        T2 = h
      )
      
      
      ######################################################
      # New boundary
      ######################################################
      
      T.boundary.next <- T.previous
    }
    
    
    ########################################################
    # Safety check
    ########################################################
    
    if (is.na(tau.next)) {
      
      break
    }
    
    
    ########################################################
    # tau' should increase
    ########################################################
    
    if (
      tau.next <=
      tau.current + tolerance
    ) {
      
      warning(
        paste0(
          "Stage 1 stopped at step ",
          step,
          ": tau' did not increase. ",
          "Current = ",
          tau.current,
          ", next = ",
          tau.next
        )
      )
      
      break
    }
    
    
    ########################################################
    # Store new candidate
    ########################################################
    
    candidates <- rbind(
      
      candidates,
      
      data.frame(
        
        step = step,
        
        tau_prime = tau.next,
        
        T_exit = T.exit.next,
        
        T_boundary = T.boundary.next,
        
        stringsAsFactors = FALSE
      )
    )
    
    
    ########################################################
    # Update recursive state
    ########################################################
    
    tau.current <- tau.next
    
    boundary.current <- T.boundary.next
  }
  
  
  ##########################################################
  # Diagnostic:
  #
  # Verify
  #
  # tau' = omega_l(T_exit,h)
  ##########################################################
  
  candidates$omega_at_exit <- sapply(
    
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
  
  
  ##########################################################
  # Definition check
  ##########################################################
  
  candidates$definition_check <-
    
    abs(
      candidates$tau_prime -
        candidates$omega_at_exit
    ) <= tolerance
  
  
  ##########################################################
  # Return
  ##########################################################
  
  list(
    
    h = h,
    
    tail.table = tail.table,
    
    tau.prime.table = candidates
  )
}


############################################################
# 4. COMPUTE E_k FOR ONE STAGE-1 CANDIDATE
############################################################
#
# For a Stage-1 candidate (tau',T_exit):
#
# E_k =
# {
#     T < T_exit :
#     omega_l(T,T_exit) >= tau'
# }.
#
############################################################

compute_E_for_candidate <- function(
    T.values,
    hat.l,
    tau.prime,
    T.exit,
    tolerance = E_TOLERANCE
) {
  
  
  ##########################################################
  # Candidate T values strictly before T_exit
  ##########################################################
  
  candidate.T <- sort(
    unique(
      T.values[
        T.values < T.exit
      ]
    )
  )
  
  
  ##########################################################
  # No T < T_exit
  ##########################################################
  
  if (length(candidate.T) == 0) {
    
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
  # Compute omega_l(T,T_exit)
  ##########################################################
  
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
  
  
  ##########################################################
  # Keep valid omega values
  ##########################################################
  
  valid <- !is.na(
    omega.values
  )
  
  
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
  # Construct E_k
  ##########################################################
  
  in.set <-
    
    valid &
    
    omega.values >=
    tau.prime - tolerance
  
  
  T.values.in.set <-
    
    candidate.T[
      in.set
    ]
  
  
  ##########################################################
  # Check whether E_k is empty
  ##########################################################
  
  empty.set <-
    
    length(
      T.values.in.set
    ) == 0
  
  
  ##########################################################
  # Summary
  ##########################################################
  
  if (empty.set) {
    
    T.min <- NA_real_
    
    T.max <- NA_real_
    
    n.T <- 0
    
  } else {
    
    T.min <- min(
      T.values.in.set
    )
    
    T.max <- max(
      T.values.in.set
    )
    
    n.T <- length(
      T.values.in.set
    )
  }
  
  
  ##########################################################
  # Return
  ##########################################################
  
  list(
    
    T_values_in_set =
      T.values.in.set,
    
    T_min =
      T.min,
    
    T_max =
      T.max,
    
    n_T =
      n.T,
    
    empty_set =
      empty.set,
    
    omega.values =
      omega.values
  )
}


############################################################
# 5. APPLY E_k CONSTRUCTION TO ALL STAGE-1 CANDIDATES
############################################################

stage1_E_sets <- function(
    T.values,
    hat.l,
    stage1.result,
    Data,
    L,
    tolerance = E_TOLERANCE
) {
  
  
  ##########################################################
  # Stage-1 candidate table
  ##########################################################
  
  tau.prime.table <-
    stage1.result$tau.prime.table
  
  
  ##########################################################
  # Number of candidates
  ##########################################################
  
  n.candidates <-
    nrow(
      tau.prime.table
    )
  
  
  ##########################################################
  # Storage
  ##########################################################
  
  results <-
    vector(
      "list",
      n.candidates
    )
  
  
  ##########################################################
  # Loop over all candidates
  ##########################################################
  
  for (i in seq_len(n.candidates)) {
    
    
    ########################################################
    # Stage-1 quantities
    ########################################################
    
    step <-
      tau.prime.table$step[i]
    
    tau.prime <-
      tau.prime.table$tau_prime[i]
    
    T.exit <-
      tau.prime.table$T_exit[i]
    
    
    ########################################################
    # Compute E_k
    ########################################################
    
    E.result <-
      
      compute_E_for_candidate(
        
        T.values =
          T.values,
        
        hat.l =
          hat.l,
        
        tau.prime =
          tau.prime,
        
        T.exit =
          T.exit,
        
        tolerance =
          tolerance
      )
    
    
    ########################################################
    # Check whether true peak L belongs to E_k
    #
    # IMPORTANT:
    #
    # E_k contains ENTRY TIMES T, not the complete window.
    #
    # This is an EX-POST diagnostic only.
    ########################################################
    
    peak.in.set <-
      
      !is.na(L) &
      
      L %in%
      E.result$T_values_in_set
    
    
    ########################################################
    # Store
    ########################################################
    
    results[[i]] <-
      
      data.frame(
        
        Data =
          Data,
        
        step =
          step,
        
        tau_prime =
          tau.prime,
        
        T_exit =
          T.exit,
        
        T_values_in_set =
          paste(
            E.result$T_values_in_set,
            collapse = ", "
          ),
        
        T_min =
          E.result$T_min,
        
        T_max =
          E.result$T_max,
        
        n_T =
          E.result$n_T,
        
        empty_set =
          E.result$empty_set,
        
        L =
          L,
        
        peak_in_set =
          peak.in.set,
        
        stringsAsFactors =
          FALSE
      )
  }
  
  
  ##########################################################
  # Combine
  ##########################################################
  
  do.call(
    rbind,
    results
  )
}


############################################################
# 6. UPDATED STAGE 2 FOR ONE CANDIDATE
############################################################
#
# E_k is first constructed.
#
# If E_k != empty:
#
#   T_star = max E_k
#
#   tau = omega_l(T_star,T_exit)
#
#   T_entry = first T attaining tau
#
############################################################

stage2_one_candidate <- function(
    T.values,
    hat.l,
    tau.prime,
    T.exit,
    L = NA_real_,
    tolerance = STAGE2_TOLERANCE
) {
  
  
  ##########################################################
  # Candidate entry times
  ##########################################################
  
  candidate.T <- sort(
    unique(
      T.values[
        T.values < T.exit
      ]
    )
  )
  
  
  ##########################################################
  # CASE 1:
  # No candidate entry time
  ##########################################################
  
  if (length(candidate.T) == 0) {
    
    return(
      list(
        
        tau.prime =
          tau.prime,
        
        T.exit =
          T.exit,
        
        E.k =
          numeric(0),
        
        T.star =
          NA_real_,
        
        tau =
          NA_real_,
        
        T.entry =
          NA_real_,
        
        w.star =
          NA_real_,
        
        L =
          L,
        
        peak.in.window =
          FALSE,
        
        candidate.T =
          numeric(0),
        
        omega.values =
          numeric(0)
      )
    )
  }
  
  
  ##########################################################
  # Compute omega_l(T,T_exit)
  ##########################################################
  
  omega.values <- sapply(
    
    candidate.T,
    
    function(T.entry) {
      
      interval_oscillation(
        
        T.values =
          T.values,
        
        hat.l =
          hat.l,
        
        T1 =
          T.entry,
        
        T2 =
          T.exit
      )
    }
  )
  
  
  ##########################################################
  # Keep valid omega values
  ##########################################################
  
  valid <-
    !is.na(
      omega.values
    )
  
  
  T.valid <-
    candidate.T[
      valid
    ]
  
  
  omega.valid <-
    omega.values[
      valid
    ]
  
  
  ##########################################################
  # CASE 2:
  # No valid omega values
  ##########################################################
  
  if (length(T.valid) == 0) {
    
    return(
      list(
        
        tau.prime =
          tau.prime,
        
        T.exit =
          T.exit,
        
        E.k =
          numeric(0),
        
        T.star =
          NA_real_,
        
        tau =
          NA_real_,
        
        T.entry =
          NA_real_,
        
        w.star =
          NA_real_,
        
        L =
          L,
        
        peak.in.window =
          FALSE,
        
        candidate.T =
          numeric(0),
        
        omega.values =
          omega.values
      )
    )
  }
  
  
  ##########################################################
  # STEP A
  #
  # Construct E_k
  ##########################################################
  
  in.E <-
    
    omega.valid >=
    tau.prime - tolerance
  
  
  E.k <-
    
    T.valid[
      in.E
    ]
  
  
  ##########################################################
  # CASE 3:
  # E_k is empty
  ##########################################################
  
  if (length(E.k) == 0) {
    
    return(
      list(
        
        tau.prime =
          tau.prime,
        
        T.exit =
          T.exit,
        
        E.k =
          numeric(0),
        
        T.star =
          NA_real_,
        
        tau =
          NA_real_,
        
        T.entry =
          NA_real_,
        
        w.star =
          NA_real_,
        
        L =
          L,
        
        peak.in.window =
          FALSE,
        
        candidate.T =
          T.valid,
        
        omega.values =
          omega.valid
      )
    )
  }
  
  
  ##########################################################
  # STEP B
  #
  # T_star = max E_k
  ##########################################################
  
  T.star <-
    
    max(
      E.k
    )
  
  
  ##########################################################
  # STEP C
  #
  # tau = omega_l(T_star,T_exit)
  ##########################################################
  
  index.star <-
    
    which(
      T.valid ==
        T.star
    )[1]
  
  
  tau <-
    
    omega.valid[
      index.star
    ]
  
  
  ##########################################################
  # STEP D
  #
  # Find all T for which omega(T,T_exit) = tau
  #
  # Numerically:
  #
  # abs(omega - tau) <= tolerance
  ##########################################################
  
  equal.tau <-
    
    abs(
      omega.valid -
        tau
    ) <=
    tolerance
  
  
  T.equal <-
    
    T.valid[
      equal.tau
    ]
  
  
  ##########################################################
  # CASE 4:
  # Numerical equality not found
  ##########################################################
  
  if (length(T.equal) == 0) {
    
    return(
      list(
        
        tau.prime =
          tau.prime,
        
        T.exit =
          T.exit,
        
        E.k =
          E.k,
        
        T.star =
          T.star,
        
        tau =
          tau,
        
        T.entry =
          NA_real_,
        
        w.star =
          NA_real_,
        
        L =
          L,
        
        peak.in.window =
          FALSE,
        
        candidate.T =
          T.valid,
        
        omega.values =
          omega.valid
      )
    )
  }
  
  
  ##########################################################
  # STEP E
  #
  # T_entry = FIRST occurrence of tau
  ##########################################################
  
  T.entry <-
    
    min(
      T.equal
    )
  
  
  ##########################################################
  # STEP F
  #
  # Window width
  ##########################################################
  
  w.star <-
    
    T.exit -
    T.entry
  
  
  ##########################################################
  # STEP G
  #
  # True peak contained in the window:
  #
  # T_entry <= L <= T_exit
  #
  # This is an EX-POST diagnostic only.
  ##########################################################
  
  peak.in.window <-
    
    !is.na(L) &
    
    T.entry <= L &
    
    L <= T.exit
  
  
  ##########################################################
  # STEP H
  #
  # Return
  ##########################################################
  
  list(
    
    tau.prime =
      tau.prime,
    
    T.exit =
      T.exit,
    
    E.k =
      E.k,
    
    T.star =
      T.star,
    
    tau =
      tau,
    
    T.entry =
      T.entry,
    
    w.star =
      w.star,
    
    L =
      L,
    
    peak.in.window =
      peak.in.window,
    
    candidate.T =
      T.valid,
    
    omega.values =
      omega.valid
  )
}


############################################################
# 7. APPLY UPDATED STAGE 2 TO ALL STAGE-1 CANDIDATES
############################################################

stage2_all_candidates <- function(
    T.values,
    hat.l,
    stage1.result,
    Data,
    L,
    tolerance = STAGE2_TOLERANCE
) {
  
  
  ##########################################################
  # Stage-1 candidate table
  ##########################################################
  
  tau.prime.table <-
    stage1.result$tau.prime.table
  
  
  ##########################################################
  # Number of candidates
  ##########################################################
  
  n.candidates <-
    nrow(
      tau.prime.table
    )
  
  
  ##########################################################
  # Storage
  ##########################################################
  
  results <-
    vector(
      "list",
      n.candidates
    )
  
  
  ##########################################################
  # Loop
  ##########################################################
  
  for (i in seq_len(n.candidates)) {
    
    
    ########################################################
    # Stage-1 quantities
    ########################################################
    
    step <-
      tau.prime.table$step[i]
    
    tau.prime <-
      tau.prime.table$tau_prime[i]
    
    T.exit <-
      tau.prime.table$T_exit[i]
    
    
    ########################################################
    # Run updated Stage 2
    ########################################################
    
    S2 <-
      
      stage2_one_candidate(
        
        T.values =
          T.values,
        
        hat.l =
          hat.l,
        
        tau.prime =
          tau.prime,
        
        T.exit =
          T.exit,
        
        L =
          L,
        
        tolerance =
          tolerance
      )
    
    
    ########################################################
    # Store
    ########################################################
    
    results[[i]] <-
      
      data.frame(
        
        Data =
          Data,
        
        step =
          step,
        
        tau_prime =
          tau.prime,
        
        T_exit =
          T.exit,
        
        E_k =
          paste(
            S2$E.k,
            collapse = ", "
          ),
        
        E_k_min =
          ifelse(
            length(S2$E.k) > 0,
            min(S2$E.k),
            NA_real_
          ),
        
        E_k_max =
          ifelse(
            length(S2$E.k) > 0,
            max(S2$E.k),
            NA_real_
          ),
        
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
        
        L =
          L,
        
        peak.in.window =
          S2$peak.in.window,
        
        stringsAsFactors =
          FALSE
      )
  }
  
  
  ##########################################################
  # Combine
  ##########################################################
  
  do.call(
    rbind,
    results
  )
}


############################################################
# 8. COMPLETE STAGE 1 + STAGE 2 PROCEDURE
############################################################
#
# IMPORTANT:
#
# L is NOT used by Stage 1.
#
# L is passed to Stage 2 only so that the returned table can
# contain the EX-POST diagnostic
#
#   T_entry <= L <= T_exit.
#
# It does not affect the construction of E_k, T_star,
# tau, or T_entry.
#
############################################################

run_updated_entry_exit <- function(
    T.values,
    hat.l,
    Data,
    L,
    tolerance.stage1 = STAGE1_TOLERANCE,
    tolerance.E = E_TOLERANCE,
    tolerance.stage2 = STAGE2_TOLERANCE
) {
  
  
  ##########################################################
  # Stage 1
  ##########################################################
  
  S1 <-
    stage1_tau_prime(
      
      T.values =
        T.values,
      
      hat.l =
        hat.l,
      
      tolerance =
        tolerance.stage1
    )
  
  
  ##########################################################
  # Stage 1: E_k sets
  #
  # This is retained as a diagnostic representation of
  # the sets E_k.
  ##########################################################
  
  E.table <-
    
    stage1_E_sets(
      
      T.values =
        T.values,
      
      hat.l =
        hat.l,
      
      stage1.result =
        S1,
      
      Data =
        Data,
      
      L =
        L,
      
      tolerance =
        tolerance.E
    )
  
  
  ##########################################################
  # Stage 2
  ##########################################################
  
  S2 <-
    
    stage2_all_candidates(
      
      T.values =
        T.values,
      
      hat.l =
        hat.l,
      
      stage1.result =
        S1,
      
      Data =
        Data,
      
      L =
        L,
      
      tolerance =
        tolerance.stage2
    )
  
  
  list(
    
    stage1 =
      S1,
    
    E_sets =
      E.table,
    
    stage2 =
      S2
  )
}


############################################################
# 9. EX-POST SIMULATION SELECTION
############################################################
#
# Among valid Stage-2 candidates:
#
#   T_entry <= L <= T_exit
#
# select the candidate having the smallest width
#
#   w_star = T_exit - T_entry.
#
# This function MUST NOT be used operationally when L is
# unknown.
#
############################################################

select_final_candidate <- function(
    stage2.table,
    L
) {
  
  
  ##########################################################
  # Empty result
  ##########################################################
  
  if (
    is.null(stage2.table) ||
    nrow(stage2.table) == 0
  ) {
    
    return(
      data.frame(
        
        step =
          NA_real_,
        
        tau_prime =
          NA_real_,
        
        T_exit =
          NA_real_,
        
        T_star =
          NA_real_,
        
        tau =
          NA_real_,
        
        T_entry =
          NA_real_,
        
        w_star =
          NA_real_,
        
        L =
          L,
        
        peak.in.window =
          FALSE,
        
        stringsAsFactors =
          FALSE
      )
    )
  }
  
  
  ##########################################################
  # Valid Stage-2 candidates
  #
  # The peak-containment condition has already been computed
  # by stage2_one_candidate() and stored in
  #
  #   peak.in.window
  #
  # We use this indicator directly rather than recomputing
  # the condition here.
  ##########################################################
  
  valid <-
    
    !is.na(stage2.table$T_entry) &
    
    !is.na(stage2.table$T_exit) &
    
    !is.na(stage2.table$w_star) &
    
    !is.na(stage2.table$peak.in.window) &
    
    stage2.table$peak.in.window
  
  
  ##########################################################
  # No valid peak-containing candidate
  ##########################################################
  
  if (!any(valid)) {
    
    return(
      data.frame(
        
        step =
          NA_real_,
        
        tau_prime =
          NA_real_,
        
        T_exit =
          NA_real_,
        
        T_star =
          NA_real_,
        
        tau =
          NA_real_,
        
        T_entry =
          NA_real_,
        
        w_star =
          NA_real_,
        
        L =
          L,
        
        peak.in.window =
          FALSE,
        
        stringsAsFactors =
          FALSE
      )
    )
  }
  
  
  ##########################################################
  # Keep valid peak-containing candidates
  ##########################################################
  
  candidates <-
    
    stage2.table[
      valid,
      ,
      drop = FALSE
    ]
  
  
  ##########################################################
  # Order by window width.
  #
  # In case of a tie, use the earliest entry time and then
  # the earliest exit time as secondary criteria.
  ##########################################################
  
  candidates <-
    
    candidates[
      order(
        candidates$w_star,
        candidates$T_entry,
        candidates$T_exit
      ),
      ,
      drop = FALSE
    ]
  
  
  ##########################################################
  # Select narrowest peak-containing candidate
  ##########################################################
  
  selected <-
    
    candidates[
      1,
      ,
      drop = FALSE
    ]
  
  
  ##########################################################
  # Preserve the true simulation peak location
  ##########################################################
  
  selected$L <-
    L
  
  
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
# CURRENT QUANTILE LEVEL
############################################################
#
# tau = 0.5 reproduces the symmetric case.
#
# Change to 0.25, 0.75, or 0.95 to study the asymmetric
# quantile-specific formulation.
#
############################################################

tau <- 0.5


############################################################
# SNR VALUES
############################################################

SNR.values <-
  c(
    1000,
    100,
    1
  )


############################################################
# 11. RUN COMPLETE PROCEDURE FOR EACH SNR
############################################################

final.rows <-
  vector(
    "list",
    length(SNR.values)
  )


for (i in seq_along(SNR.values)) {
  
  
  SNR <-
    SNR.values[i]
  
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "tau =",
    tau,
    "\n"
  )
  
  cat(
    "SNR =",
    SNR,
    "\n"
  )
  
  cat(
    "========================================\n"
  )
  
  
  ##########################################################
  # Estimate sequential trajectory
  ##########################################################
  
  res <-
    estimate_peak_curves(
      
      A =
        A,
      
      L =
        L,
      
      S =
        S,
      
      h =
        h,
      
      SNR =
        SNR,
      
      tau =
        tau,
      
      seed =
        seed
    )
  
  
  ##########################################################
  # Extract T-grid and location estimates
  ##########################################################
  
  T.values <-
    res$T
  
  hat.l <-
    res$hat.l
  
  
  ##########################################################
  # Run complete operational Stage 1 + Stage 2 procedure
  #
  # L is NOT used to construct the operational entry/exit
  # quantities.
  #
  # It is used only for EX-POST diagnostics.
  ##########################################################
  
  result <-
    
    run_updated_entry_exit(
      
      T.values =
        T.values,
      
      hat.l =
        hat.l,
      
      Data =
        SNR,
      
      L =
        L,
      
      tolerance.stage1 =
        STAGE1_TOLERANCE,
      
      tolerance.E =
        E_TOLERANCE,
      
      tolerance.stage2 =
        STAGE2_TOLERANCE
    )
  
  
  ##########################################################
  # Extract Stage-2 table
  ##########################################################
  
  stage2.table <-
    result$stage2
  
  
  ##########################################################
  # EX-POST selection using the known simulation L
  ##########################################################
  
  selected <-
    
    select_final_candidate(
      
      stage2.table =
        stage2.table,
      
      L =
        L
    )
  
  
  ##########################################################
  # Print Stage-1 candidates
  ##########################################################
  
  cat(
    "\nStage-1 candidates:\n"
  )
  
  print(
    result$stage1$tau.prime.table
  )
  
  
  ##########################################################
  # Print E_k sets
  ##########################################################
  
  cat(
    "\nE_k sets:\n"
  )
  
  print(
    result$E_sets
  )
  
  
  ##########################################################
  # Print Stage-2 candidates
  ##########################################################
  
  cat(
    "\nStage-2 candidates:\n"
  )
  
  print(
    stage2.table
  )
  
  
  ##########################################################
  # Print selected candidate
  ##########################################################
  
  cat(
    "\nSelected candidate:\n"
  )
  
  print(
    selected
  )
  
  
  ##########################################################
  # Store result
  ##########################################################
  
  final.rows[[i]] <-
    
    data.frame(
      
      tau =
        tau,
      
      SNR =
        SNR,
      
      k =
        selected$step,
      
      tau_prime =
        selected$tau_prime,
      
      T_exit =
        selected$T_exit,
      
      T_star =
        selected$T_star,
      
      tau_entry =
        selected$tau,
      
      T_entry =
        selected$T_entry,
      
      w_star =
        selected$w_star,
      
      L =
        L,
      
      peak.in.window =
        selected$peak.in.window,
      
      stringsAsFactors =
        FALSE
    )
}


############################################################
# 12. COMBINE RESULTS
############################################################

table.entry.interpretation <-
  
  do.call(
    rbind,
    final.rows
  )


rownames(
  table.entry.interpretation
) <-
  NULL


############################################################
# 13. EX-POST RELATION WITH THE TRUE PEAK
############################################################
#
# Closed containment is used:
#
#   T_entry <= L <= T_exit.
#
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
# 14. INTERPRETATION
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
# 15. DISPLAY FINAL TABLE
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
# 16. SAVE CSV
############################################################

write.csv(
  
  table.entry.interpretation,
  
  file =
    "entry-interpretation-updated.csv",
  
  row.names =
    FALSE
)


############################################################
# 17. LATEX-READY TABLE
############################################################

cat(
  "\n\n========================================\n"
)

cat(
  "LATEX TABLE\n"
)

cat(
  "========================================\n\n"
)

print(
  
  knitr::kable(
    
    table.entry.interpretation,
    
    format =
      "latex",
    
    booktabs =
      TRUE,
    
    digits =
      6
  )
)
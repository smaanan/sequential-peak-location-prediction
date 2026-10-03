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
# NEW in this version:
#   The script now outputs, for every (tau, SNR) pair:
#     - all Stage-2 candidates with their tuple
#       (T_entry^(k), T_exit^(k), w_k)
#     - a strict-containment indicator for L
#     - the first peak-containing k, including and excluding
#       the initialization candidate k = 0
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
# 2. COMPLETE TAIL-OSCILLATION TABLE
############################################################

compute_tail_oscillation <- function(
    T.values,
    hat.l
) {
  
  h <- max(T.values)
  T.grid <- sort(unique(T.values))
  
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
    omega = omega,
    stringsAsFactors = FALSE
  )
}


############################################################
# 3. STAGE 1
#
# Recursive construction of (tau'_k, T_exit^(k))
############################################################

stage1_tau_prime <- function(
    T.values,
    hat.l,
    tolerance = STAGE1_TOLERANCE
) {
  
  h <- max(T.values)
  T.grid <- sort(unique(T.values))
  
  tail.table <- compute_tail_oscillation(
    T.values = T.values,
    hat.l = hat.l
  )
  
  candidates <- data.frame(
    step = 0,
    tau_prime = 0,
    T_exit = h,
    T_boundary = h,
    stringsAsFactors = FALSE
  )
  
  tau.current <- 0
  boundary.current <- h
  
  max.iterations <- length(T.grid) + 1
  
  for (step in seq_len(max.iterations)) {
    
    earlier.T <- T.grid[
      T.grid < boundary.current
    ]
    
    if (length(earlier.T) == 0) {
      break
    }
    
    possible.T <- tail.table$T[
      !is.na(tail.table$omega) &
        abs(tail.table$omega - tau.current) <= tolerance &
        tail.table$T < boundary.current
    ]
    
    ########################################################
    # CASE 1: current tau' attained earlier
    ########################################################
    
    if (length(possible.T) > 0) {
      
      T.prime <- min(possible.T)
      
      previous.T <- T.grid[
        T.grid < T.prime
      ]
      
      if (length(previous.T) == 0) {
        break
      }
      
      T.previous <- max(previous.T)
      
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
      # CASE 2: no earlier T' satisfies equality
      ######################################################
      
      T.previous <- max(earlier.T)
      
      T.exit.next <- T.previous
      
      tau.next <- interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = T.previous,
        T2 = h
      )
      
      T.boundary.next <- T.previous
    }
    
    if (is.na(tau.next)) {
      break
    }
    
    if (tau.next <= tau.current + tolerance) {
      break
    }
    
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
    
    tau.current <- tau.next
    boundary.current <- T.boundary.next
  }
  
  ##########################################################
  # Definition check
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
  
  candidates$definition_check <-
    abs(candidates$tau_prime - candidates$omega_at_exit) <= tolerance
  
  list(
    h = h,
    tail.table = tail.table,
    tau.prime.table = candidates
  )
}


############################################################
# 4. COMPUTE E_k FOR ONE STAGE-1 CANDIDATE
############################################################

compute_E_for_candidate <- function(
    T.values,
    hat.l,
    tau.prime,
    T.exit,
    tolerance = E_TOLERANCE
) {
  
  candidate.T <- sort(
    unique(
      T.values[
        T.values < T.exit
      ]
    )
  )
  
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
  
  in.set <-
    valid &
    omega.values >= tau.prime - tolerance
  
  T.values.in.set <- candidate.T[in.set]
  
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
  
  tau.prime.table <- stage1.result$tau.prime.table
  n.candidates <- nrow(tau.prime.table)
  results <- vector("list", n.candidates)
  
  for (i in seq_len(n.candidates)) {
    
    step <- tau.prime.table$step[i]
    tau.prime <- tau.prime.table$tau_prime[i]
    T.exit <- tau.prime.table$T_exit[i]
    
    E.result <- compute_E_for_candidate(
      T.values = T.values,
      hat.l = hat.l,
      tau.prime = tau.prime,
      T.exit = T.exit,
      tolerance = tolerance
    )
    
    peak.in.set <-
      !is.na(L) &
      L %in% E.result$T_values_in_set
    
    results[[i]] <- data.frame(
      Data = Data,
      step = step,
      tau_prime = tau.prime,
      T_exit = T.exit,
      T_values_in_set = paste(E.result$T_values_in_set, collapse = ", "),
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
############################################################

stage2_one_candidate <- function(
    T.values,
    hat.l,
    tau.prime,
    T.exit,
    L = NA_real_,
    tolerance = STAGE2_TOLERANCE
) {
  
  candidate.T <- sort(
    unique(
      T.values[
        T.values < T.exit
      ]
    )
  )
  
  if (length(candidate.T) == 0) {
    return(
      list(
        tau.prime = tau.prime,
        T.exit = T.exit,
        E.k = numeric(0),
        T.star = NA_real_,
        tau = NA_real_,
        T.entry = NA_real_,
        w.star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        candidate.T = numeric(0),
        omega.values = numeric(0)
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
  
  valid <- !is.na(omega.values)
  T.valid <- candidate.T[valid]
  omega.valid <- omega.values[valid]
  
  if (length(T.valid) == 0) {
    return(
      list(
        tau.prime = tau.prime,
        T.exit = T.exit,
        E.k = numeric(0),
        T.star = NA_real_,
        tau = NA_real_,
        T.entry = NA_real_,
        w.star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        candidate.T = numeric(0),
        omega.values = omega.values
      )
    )
  }
  
  ##########################################################
  # STEP A: Construct E_k
  ##########################################################
  
  in.E <- omega.valid >= tau.prime - tolerance
  E.k <- T.valid[in.E]
  
  if (length(E.k) == 0) {
    return(
      list(
        tau.prime = tau.prime,
        T.exit = T.exit,
        E.k = numeric(0),
        T.star = NA_real_,
        tau = NA_real_,
        T.entry = NA_real_,
        w.star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        candidate.T = T.valid,
        omega.values = omega.valid
      )
    )
  }
  
  ##########################################################
  # STEP B: T_k^* = max E_k
  ##########################################################
  
  T.star <- max(E.k)
  
  ##########################################################
  # STEP C: tau_k = omega_l(T_k^*, T_exit)
  ##########################################################
  
  index.star <- which(T.valid == T.star)[1]
  tau <- omega.valid[index.star]
  
  ##########################################################
  # STEP D: Find all T for which omega(T, T_exit) = tau_k
  ##########################################################
  
  equal.tau <- abs(omega.valid - tau) <= tolerance
  T.equal <- T.valid[equal.tau]
  
  if (length(T.equal) == 0) {
    return(
      list(
        tau.prime = tau.prime,
        T.exit = T.exit,
        E.k = E.k,
        T.star = T.star,
        tau = tau,
        T.entry = NA_real_,
        w.star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        candidate.T = T.valid,
        omega.values = omega.valid
      )
    )
  }
  
  ##########################################################
  # STEP E: T_entry = FIRST occurrence of tau_k
  ##########################################################
  
  T.entry <- min(T.equal)
  
  ##########################################################
  # STEP F: Window width
  ##########################################################
  
  w.star <- T.exit - T.entry
  
  ##########################################################
  # STEP G: EX-POST peak containment
  #
  # Diagnostic field "peak.in.window" keeps the strict form
  # T_entry < L < T_exit, matching the manuscript criterion.
  ##########################################################
  
  peak.in.window <-
    !is.na(L) &
    T.entry < L &
    L < T.exit
  
  ##########################################################
  # STEP H: Return
  ##########################################################
  
  list(
    tau.prime = tau.prime,
    T.exit = T.exit,
    E.k = E.k,
    T.star = T.star,
    tau = tau,
    T.entry = T.entry,
    w.star = w.star,
    L = L,
    peak.in.window = peak.in.window,
    candidate.T = T.valid,
    omega.values = omega.valid
  )
}


############################################################
# 7. APPLY STAGE 2 TO ALL STAGE-1 CANDIDATES
############################################################

stage2_all_candidates <- function(
    T.values,
    hat.l,
    stage1.result,
    Data,
    L,
    tolerance = STAGE2_TOLERANCE
) {
  
  tau.prime.table <- stage1.result$tau.prime.table
  n.candidates <- nrow(tau.prime.table)
  results <- vector("list", n.candidates)
  
  for (i in seq_len(n.candidates)) {
    
    step <- tau.prime.table$step[i]
    tau.prime <- tau.prime.table$tau_prime[i]
    T.exit <- tau.prime.table$T_exit[i]
    
    S2 <- stage2_one_candidate(
      T.values = T.values,
      hat.l = hat.l,
      tau.prime = tau.prime,
      T.exit = T.exit,
      L = L,
      tolerance = tolerance
    )
    
    results[[i]] <- data.frame(
      Data = Data,
      step = step,
      tau_prime = tau.prime,
      T_exit = T.exit,
      E_k = paste(S2$E.k, collapse = ", "),
      E_k_min = if (length(S2$E.k) > 0) min(S2$E.k) else NA_real_,
      E_k_max = if (length(S2$E.k) > 0) max(S2$E.k) else NA_real_,
      n_E_k = length(S2$E.k),
      E_empty = length(S2$E.k) == 0,
      T_star = S2$T.star,
      tau = S2$tau,
      T_entry = S2$T.entry,
      w_star = S2$w.star,
      L = L,
      peak.in.window = S2$peak.in.window,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, results)
}


############################################################
# 8. COMPLETE STAGE 1 + STAGE 2 PROCEDURE
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
  
  S1 <- stage1_tau_prime(
    T.values = T.values,
    hat.l = hat.l,
    tolerance = tolerance.stage1
  )
  
  E.table <- stage1_E_sets(
    T.values = T.values,
    hat.l = hat.l,
    stage1.result = S1,
    Data = Data,
    L = L,
    tolerance = tolerance.E
  )
  
  S2 <- stage2_all_candidates(
    T.values = T.values,
    hat.l = hat.l,
    stage1.result = S1,
    Data = Data,
    L = L,
    tolerance = tolerance.stage2
  )
  
  list(
    stage1 = S1,
    E_sets = E.table,
    stage2 = S2
  )
}


############################################################
# 9. EX-POST SIMULATION SELECTION
#
# Selects the peak-containing candidate with the smallest
# width w_k. This is the manuscript's rule.
############################################################

select_final_candidate <- function(
    stage2.table,
    L
) {
  
  if (is.null(stage2.table) || nrow(stage2.table) == 0) {
    return(
      data.frame(
        step = NA_real_,
        tau_prime = NA_real_,
        T_exit = NA_real_,
        T_star = NA_real_,
        tau = NA_real_,
        T_entry = NA_real_,
        w_star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }
  
  valid <-
    !is.na(stage2.table$T_entry) &
    !is.na(stage2.table$T_exit) &
    !is.na(stage2.table$w_star) &
    !is.na(stage2.table$peak.in.window) &
    stage2.table$peak.in.window
  
  if (!any(valid)) {
    return(
      data.frame(
        step = NA_real_,
        tau_prime = NA_real_,
        T_exit = NA_real_,
        T_star = NA_real_,
        tau = NA_real_,
        T_entry = NA_real_,
        w_star = NA_real_,
        L = L,
        peak.in.window = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }
  
  candidates <- stage2.table[valid, , drop = FALSE]
  
  candidates <- candidates[
    order(
      candidates$w_star,
      candidates$T_entry,
      candidates$T_exit
    ),
    ,
    drop = FALSE
  ]
  
  selected <- candidates[1, , drop = FALSE]
  selected$L <- L
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

tau.values <- c(0.25, 0.50, 0.75, 0.95)


############################################################
# SNR VALUES
############################################################

SNR.values <- c(1000, 100, 1)


############################################################
# 11. RUN COMPLETE PROCEDURE FOR EACH TAU AND EACH SNR
############################################################

final.rows <- vector(
  "list",
  length(tau.values) * length(SNR.values)
)

all_diagnostics <- vector(
  "list",
  length(tau.values) * length(SNR.values)
)

counter <- 1
diag_counter <- 1

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
    # Estimate sequential trajectory
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
    # Run complete operational Stage 1 + Stage 2
    ########################################################
    
    result <- run_updated_entry_exit(
      T.values = T.values,
      hat.l = hat.l,
      Data = SNR,
      L = L,
      tolerance.stage1 = STAGE1_TOLERANCE,
      tolerance.E = E_TOLERANCE,
      tolerance.stage2 = STAGE2_TOLERANCE
    )
    
    ########################################################
    # Extract Stage-2 table
    ########################################################
    
    stage2.table <- result$stage2
    
    ########################################################
    # DIAGNOSTIC: full Stage-2 candidate table
    #
    # Contains, for every candidate k (including the
    # initialization k = 0):
    #   step (k)
    #   tau_prime (epsilon'_k)
    #   T_exit^(k)
    #   T_star
    #   tau (epsilon_k)
    #   T_entry^(k)
    #   w_k = T_exit^(k) - T_entry^(k)
    #   peak.in.strict = 1{T_entry^(k) < L < T_exit^(k)}
    ########################################################
    
    diag <- stage2.table
    
    diag$w_k <- diag$T_exit - diag$T_entry
    
    diag$peak.in.strict <-
      !is.na(diag$T_entry) &
      !is.na(diag$T_exit) &
      diag$T_entry < L &
      L < diag$T_exit
    
    diag$tau_level <- tau
    diag$SNR_level <- SNR
    
    all_diagnostics[[diag_counter]] <- diag
    diag_counter <- diag_counter + 1
    
    ########################################################
    # EX-POST selection (manuscript rule: smallest width)
    ########################################################
    
    selected <- select_final_candidate(
      stage2.table = stage2.table,
      L = L
    )
    
    ########################################################
    # Print Stage-1 candidates
    ########################################################
    
    cat("\nStage-1 candidates:\n")
    print(result$stage1$tau.prime.table)
    
    ########################################################
    # Print E_k sets
    ########################################################
    
    cat("\nE_k sets:\n")
    print(result$E_sets)
    
    ########################################################
    # Print Stage-2 candidates
    ########################################################
    
    cat("\nStage-2 candidates:\n")
    print(stage2.table)
    
    ########################################################
    # Print selected candidate
    ########################################################
    
    cat("\nSelected candidate:\n")
    print(selected)
    
    ########################################################
    # Store result
    ########################################################
    
    final.rows[[counter]] <- data.frame(
      tau = tau,
      SNR = SNR,
      k = selected$step,
      tau_prime = selected$tau_prime,
      T_exit = selected$T_exit,
      T_star = selected$T_star,
      tau_entry = selected$tau,
      T_entry = selected$T_entry,
      w_star = selected$w_star,
      L = L,
      peak.in.window = selected$peak.in.window,
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}


############################################################
# 12. COMBINE SELECTED RESULTS
############################################################

table.entry.interpretation <- do.call(rbind, final.rows)
rownames(table.entry.interpretation) <- NULL


############################################################
# 12b. COMBINE DIAGNOSTICS
############################################################

combined.diagnostics <- do.call(rbind, all_diagnostics)
rownames(combined.diagnostics) <- NULL


############################################################
# 12c. DISPLAY FULL STAGE-2 CANDIDATE TABLE
############################################################

cat("\n\n========================================\n")
cat("ALL STAGE-2 CANDIDATES (diagnostic)\n")
cat("========================================\n\n")

print(
  combined.diagnostics[
    ,
    c(
      "tau_level", "SNR_level", "step",
      "tau_prime", "T_exit", "T_star", "tau",
      "T_entry", "w_k",
      "peak.in.strict"
    )
  ]
)


############################################################
# 12d. FIRST PEAK-CONTAINING k BY (tau, SNR)
#
# Two variants are reported:
#   first_k_incl0 : smallest k >= 0 (includes initialization)
#   first_k_excl0 : smallest k >= 1 (nontrivial candidates only)
############################################################

# Split by (tau, SNR) pairs and compute first admissible k

groups <- split(
  combined.diagnostics,
  interaction(
    combined.diagnostics$tau_level,
    combined.diagnostics$SNR_level,
    drop = TRUE
  )
)

first_k_table <- do.call(
  rbind,
  lapply(
    groups,
    function(df) {
      
      valid_incl0 <- df[df$peak.in.strict, ]
      valid_excl0 <- df[df$step >= 1 & df$peak.in.strict, ]
      
      data.frame(
        tau           = unique(df$tau_level),
        SNR           = unique(df$SNR_level),
        first_k_incl0 = if (nrow(valid_incl0) > 0)
          min(valid_incl0$step) else NA_integer_,
        first_k_excl0 = if (nrow(valid_excl0) > 0)
          min(valid_excl0$step) else NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(first_k_table) <- NULL
first_k_table <- first_k_table[order(first_k_table$tau, -first_k_table$SNR), ]
rownames(first_k_table) <- NULL

cat("\n\n========================================\n")
cat("FIRST PEAK-CONTAINING k BY (tau, SNR)\n")
cat("========================================\n\n")

cat(
  "first_k_incl0 : smallest k >= 0 such that T_entry^(k) < L < T_exit^(k)\n"
)

cat(
  "first_k_excl0 : smallest k >= 1 such that T_entry^(k) < L < T_exit^(k)\n\n"
)

print(first_k_table)


############################################################
# 13. RELATION WITH THE TRUE PEAK (selected candidate)
############################################################

table.entry.interpretation$relation <-
  with(
    table.entry.interpretation,
    ifelse(
      is.na(T_entry) | is.na(T_exit),
      NA_character_,
      ifelse(
        T_exit < L,
        "T_exit < L",
        ifelse(
          T_entry <= L & L <= T_exit,
          "T_entry <= L <= T_exit",
          "L < T_entry"
        )
      )
    )
  )


############################################################
# 14. INTERPRETATION (selected candidate)
############################################################

table.entry.interpretation$interpretation <-
  with(
    table.entry.interpretation,
    ifelse(
      is.na(T_entry) | is.na(T_exit),
      "No valid peak-containing candidate",
      ifelse(
        T_exit < L,
        "Pre-peak candidate",
        ifelse(
          T_entry <= L & L <= T_exit,
          "Peak-containing candidate",
          "Post-peak candidate"
        )
      )
    )
  )


############################################################
# 15. DISPLAY FINAL SELECTED TABLE
############################################################

cat("\n\n========================================\n")
cat("FINAL ENTRY-INTERPRETATION TABLE\n")
cat("========================================\n\n")

print(table.entry.interpretation)


############################################################
# 16. SAVE CSV
############################################################

write.csv(
  table.entry.interpretation,
  "entry-interpretation-updated-all-tau.csv",
  row.names = FALSE
)

write.csv(
  combined.diagnostics,
  "simulation_all_stage2_candidates.csv",
  row.names = FALSE
)

write.csv(
  first_k_table,
  "simulation_first_k_summary.csv",
  row.names = FALSE
)


############################################################
# 17. LATEX-READY COMBINED TABLE (selected candidate)
############################################################

cat("\n\n========================================\n")
cat("LATEX TABLE — ALL TAU VALUES (selected)\n")
cat("========================================\n\n")

print(
  knitr::kable(
    table.entry.interpretation,
    format = "latex",
    booktabs = TRUE,
    digits = 6
  )
)

cat("\n\n========================================\n")
cat("LATEX TABLE — FIRST K SUMMARY\n")
cat("========================================\n\n")

print(
  knitr::kable(
    first_k_table,
    format = "latex",
    booktabs = TRUE,
    digits = 6
  )
)


############################################################
# 18. OPTIONAL: DISPLAY ONE TABLE PER TAU (selected)
############################################################

for (tau in tau.values) {
  
  cat("\n\n========================================\n")
  cat("FINAL SELECTED TABLE FOR tau =", tau, "\n")
  cat("========================================\n\n")
  
  tau.table <- table.entry.interpretation[
    table.entry.interpretation$tau == tau,
    ,
    drop = FALSE
  ]
  
  print(tau.table)
  
  cat("\nLaTeX table for tau =", tau, ":\n\n")
  
  print(
    knitr::kable(
      tau.table,
      format = "latex",
      booktabs = TRUE,
      digits = 6
    )
  )
}
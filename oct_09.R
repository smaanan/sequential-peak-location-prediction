############################################################
#
# ENTRY / EXIT ANALYSIS OF THE PEAK REGION
#
# MODEL:
#
#   I(t) = A exp(-(t-L)^2/S^2) + e(t)
#
# The estimator is the current tau-quantile estimator with
# tau-aligned asymmetric Laplace noise.
#
# OPERATIONAL CONSTRUCTION:
#
#   1. Compute the sequential trajectory hat.l(T).
#   2. For each current horizon T, construct all candidate
#      intervals [T_entry^(k)(T), T_exit^(k)(T)].
#   3. Compute T_est:
#      the first T for which a candidate interval contains
#      the current point estimate hat.l(T).
#
# EX-POST SIMULATION DIAGNOSTICS:
#
#   4. Compute T_first:
#      the first T>L for which a selected candidate strictly
#      contains the true peak location L.
#   5. Compute T_det:
#      the first T>L such that a selected candidate strictly
#      contains L at every available horizon T,...,h.
#
# Selection conventions:
#
#   T_est:
#     Among candidates containing hat.l(T), select the
#     narrowest interval. Ties are broken by smaller k,
#     then smaller T_entry, then smaller T_exit.
#
#   Ex-post selection:
#     Among candidates strictly containing L, select the
#     narrowest interval. Ties are broken by smaller k,
#     then smaller T_entry, then smaller T_exit.
#
# IMPORTANT:
#
#   L is not used to construct the intervals or calculate
#   T_est. It is used only in the simulation diagnostics
#   T_first and T_det.
#
#   These intervals are data-driven candidate intervals.
#   Their containment of hat.l(T) does not itself establish
#   statistical coverage of L.
#
############################################################


############################################################
# 0. TOLERANCES
############################################################

STAGE1_TOLERANCE <- 1e-10
E_TOLERANCE      <- 0
STAGE2_TOLERANCE <- 0


############################################################
# 1. INTERVAL OSCILLATION
#
# omega_l(T1,T2)
#   = max_{T1 <= t1,t2 <= T2} |hat.l(t1)-hat.l(t2)|.
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
# Compute t -> omega_l(t,T.current) for every available
# observation horizon t <= T.current.
#
############################################################

compute_retrospective_oscillation <- function(
    T.values,
    hat.l,
    T.current
) {
  
  T.grid <- sort(
    unique(
      T.values[T.values <= T.current]
    )
  )
  
  if (length(T.grid) == 0) {
    stop("No observation horizons are available up to T.current.")
  }
  
  omega <- vapply(
    T.grid,
    function(t) {
      interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = t,
        T2 = T.current
      )
    },
    numeric(1)
  )
  
  data.frame(
    T_current = T.current,
    t = T.grid,
    omega = omega,
    stringsAsFactors = FALSE
  )
}


############################################################
# 3. STAGE 1: CONSTRUCT EXIT CANDIDATES
#
# Initialization:
#
#   T_exit^(0)(T) = T,
#   epsilon'_0(T) = 0.
#
# Nontrivial candidates:
#
#   T_exit^(k)(T) is a strict-decrease point of
#   t -> omega_l(t,T), ordered in decreasing order.
#
# Strict decrease is assessed between consecutive
# AVAILABLE observation horizons, not by assuming that
# every possible integer horizon is present.
#
############################################################

stage1_exit_candidates <- function(
    T.values,
    hat.l,
    T.current,
    tolerance = STAGE1_TOLERANCE
) {
  
  retrospective.table <- compute_retrospective_oscillation(
    T.values = T.values,
    hat.l = hat.l,
    T.current = T.current
  )
  
  T.grid <- retrospective.table$t
  omega <- retrospective.table$omega
  n.grid <- length(T.grid)
  
  candidates <- data.frame(
    T_current = T.current,
    step = 0L,
    T_exit = T.current,
    epsilon_prime = 0,
    strict_decrease = FALSE,
    definition_check = TRUE,
    stringsAsFactors = FALSE
  )
  
  strict.points <- numeric(0)
  
  if (n.grid >= 2) {
    
    # Compare each pair of consecutive AVAILABLE horizons.
    # Since omega(t,T) is nonincreasing in t, a strict
    # decrease occurs when omega at the earlier horizon
    # is greater than omega at the next available horizon.
    
    strict.index <- which(
      omega[-n.grid] > omega[-1L] + tolerance
    )
    
    if (length(strict.index) > 0) {
      strict.points <- T.grid[strict.index]
    }
  }
  
  # Largest exit time first.
  strict.points <- sort(
    strict.points,
    decreasing = TRUE
  )
  
  if (length(strict.points) > 0) {
    
    strict.rows <- lapply(
      seq_along(strict.points),
      function(k) {
        
        T.exit <- strict.points[k]
        
        epsilon.prime <- interval_oscillation(
          T.values = T.values,
          hat.l = hat.l,
          T1 = T.exit,
          T2 = T.current
        )
        
        data.frame(
          T_current = T.current,
          step = as.integer(k),
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
      do.call(rbind, strict.rows)
    )
  }
  
  rownames(candidates) <- NULL
  
  list(
    T_current = T.current,
    retrospective.table = retrospective.table,
    strict_decrease_points = strict.points,
    exit.table = candidates
  )
}


############################################################
# 4. COMPUTE E_k(T) FOR ONE EXIT CANDIDATE
#
# E_k(T) = {
#   t < T_exit^(k)(T):
#   omega_l(t,T_exit^(k)(T)) >= epsilon'_k(T)
# }.
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
        n_T = 0L,
        empty_set = TRUE,
        omega.values = numeric(0)
      )
    )
  }
  
  omega.values <- vapply(
    candidate.t,
    function(t) {
      interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = t,
        T2 = T.exit
      )
    },
    numeric(1)
  )
  
  valid <- !is.na(omega.values)
  
  in.set <- valid &
    omega.values >= epsilon.prime - tolerance
  
  T.values.in.set <- candidate.t[in.set]
  
  if (length(T.values.in.set) == 0) {
    return(
      list(
        T_values_in_set = numeric(0),
        T_min = NA_real_,
        T_max = NA_real_,
        n_T = 0L,
        empty_set = TRUE,
        omega.values = omega.values
      )
    )
  }
  
  list(
    T_values_in_set = T.values.in.set,
    T_min = min(T.values.in.set),
    T_max = max(T.values.in.set),
    n_T = length(T.values.in.set),
    empty_set = FALSE,
    omega.values = omega.values
  )
}


############################################################
# 5. APPLY E_k(T) CONSTRUCTION TO ALL EXIT CANDIDATES
#
# Retained for diagnostics.
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
  
  results <- vector("list", nrow(exit.table))
  
  for (i in seq_len(nrow(exit.table))) {
    
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
    
    peak.in.set <- !is.na(L) &&
      any(E.result$T_values_in_set == L)
    
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
# 6. STAGE 2: CONSTRUCT ONE CANDIDATE INTERVAL
#
# For a fixed current horizon T.current and candidate
# exit time T.exit:
#
#   T_star = max E_k(T.current)
#
#   epsilon_k(T.current)
#       = omega_l(T_star,T.exit)
#
#   T_entry^(k)(T.current)
#       = min{t < T.exit:
#             omega_l(t,T.exit) = epsilon_k(T.current)}
#
#   w_k(T.current) = T_exit - T_entry.
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
  
  E.result <- compute_E_for_candidate(
    T.values = T.values,
    hat.l = hat.l,
    T.current = T.current,
    epsilon.prime = epsilon.prime,
    T.exit = T.exit,
    tolerance = E_TOLERANCE
  )
  
  E.k <- E.result$T_values_in_set
  
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
        estimate.in.window = FALSE
      )
    )
  }
  
  T.star <- max(E.k)
  
  epsilon <- interval_oscillation(
    T.values = T.values,
    hat.l = hat.l,
    T1 = T.star,
    T2 = T.exit
  )
  
  candidate.t <- sort(
    unique(
      T.values[
        T.values < T.exit &
          T.values <= T.current
      ]
    )
  )
  
  omega.values <- vapply(
    candidate.t,
    function(t) {
      interval_oscillation(
        T.values = T.values,
        hat.l = hat.l,
        T1 = t,
        T2 = T.exit
      )
    },
    numeric(1)
  )
  
  valid <- !is.na(omega.values)
  
  equal.epsilon <- valid &
    abs(omega.values - epsilon) <= tolerance
  
  T.equal <- candidate.t[equal.epsilon]
  
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
        estimate.in.window = FALSE
      )
    )
  }
  
  T.entry <- min(T.equal)
  w.star <- T.exit - T.entry
  
  # Ex-post containment uses strict inequalities, matching
  # the definition in the existing Section 7.5.
  peak.in.window <- !is.na(L) &&
    T.entry < L &&
    L < T.exit
  
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
    estimate.in.window = FALSE
  )
}


############################################################
# 7. APPLY STAGE 2 TO ALL EXIT CANDIDATES
#
# Both the true-peak containment diagnostic and the
# point-estimate containment diagnostic are recorded.
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
  
  # Retrieve the current point estimate.
  idx.current <- which(T.values == T.current)
  
  if (length(idx.current) == 0) {
    stop("T.current is not present in T.values.")
  }
  
  hat.l.current <- hat.l[idx.current[1]]
  
  results <- vector("list", nrow(exit.table))
  
  for (i in seq_len(nrow(exit.table))) {
    
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
    
    valid.interval <- !is.na(S2$T.entry) && !is.na(S2$T.exit)
    
    estimate.in.window <- valid.interval &&
      S2$T.entry <= hat.l.current &&
      hat.l.current <= S2$T.exit
    
    results[[i]] <- data.frame(
      Data = Data,
      T_current = T.current,
      k = step,
      epsilon_prime = epsilon.prime,
      T_exit = T.exit,
      E_k = paste(S2$E.k, collapse = ", "),
      E_k_min = if (length(S2$E.k) > 0)
        min(S2$E.k) else NA_real_,
      E_k_max = if (length(S2$E.k) > 0)
        max(S2$E.k) else NA_real_,
      n_E_k = length(S2$E.k),
      E_empty = length(S2$E.k) == 0,
      T_star = S2$T.star,
      epsilon = S2$epsilon,
      T_entry = S2$T.entry,
      w_k = S2$w.star,
      L = L,
      hat_l_current = hat.l.current,
      peak.in.strict = S2$peak.in.window,
      estimate.in.window = estimate.in.window,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, results)
}


############################################################
# 8. COMPLETE PROCEDURE FOR ONE CURRENT HORIZON T
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
  
  S1 <- stage1_exit_candidates(
    T.values = T.values,
    hat.l = hat.l,
    T.current = T.current,
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
    T_current = T.current,
    stage1 = S1,
    E_sets = E.table,
    stage2 = S2
  )
}


############################################################
# 9. COMMON CANDIDATE SELECTION FUNCTION
#
# Select the narrowest interval among candidates satisfying
# a supplied containment condition.
#
# Selection order:
#   1. smallest width;
#   2. smallest k;
#   3. smallest T_entry;
#   4. smallest T_exit.
#
############################################################

select_narrowest_candidate <- function(
    stage2.table,
    containment.column,
    L
) {
  
  empty_row <- function(T.current, hat.l.current) {
    data.frame(
      T_current = T.current,
      k = NA_real_,
      epsilon_prime = NA_real_,
      T_exit = NA_real_,
      T_star = NA_real_,
      epsilon = NA_real_,
      T_entry = NA_real_,
      w_k = NA_real_,
      L = L,
      hat_l_current = hat.l.current,
      peak.in.strict = FALSE,
      estimate.in.window = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  if (is.null(stage2.table) || nrow(stage2.table) == 0) {
    return(empty_row(NA_real_, NA_real_))
  }
  
  T.current <- unique(stage2.table$T_current)[1]
  hat.l.current <- unique(stage2.table$hat_l_current)[1]
  
  condition <- stage2.table[[containment.column]]
  
  valid <- !is.na(condition) &
    condition &
    !is.na(stage2.table$T_entry) &
    !is.na(stage2.table$T_exit) &
    !is.na(stage2.table$w_k)
  
  if (!any(valid)) {
    return(empty_row(T.current, hat.l.current))
  }
  
  candidates <- stage2.table[valid, , drop = FALSE]
  
  candidates <- candidates[
    order(
      candidates$w_k,
      candidates$k,
      candidates$T_entry,
      candidates$T_exit
    ),
    ,
    drop = FALSE
  ]
  
  out.cols <- c(
    "T_current", "k", "epsilon_prime", "T_exit", "T_star",
    "epsilon", "T_entry", "w_k", "L", "hat_l_current",
    "peak.in.strict", "estimate.in.window"
  )
  
  selected <- candidates[1, out.cols, drop = FALSE]
  rownames(selected) <- NULL
  selected
}


############################################################
# 10. EX-POST SELECTION AT ONE CURRENT HORIZON
#
# Among candidates strictly containing L, select the
# narrowest interval.
#
############################################################

select_final_candidate <- function(
    stage2.table,
    L
) {
  select_narrowest_candidate(
    stage2.table = stage2.table,
    containment.column = "peak.in.strict",
    L = L
  )
}


############################################################
# 11. OPERATIONAL SELECTION AT ONE CURRENT HORIZON
#
# Among candidates satisfying
#
#   T_entry^(k)(T) <= hat.l(T) <= T_exit^(k)(T),
#
# select the narrowest interval.
#
# This selection does not use L.
#
############################################################

select_estimate_candidate <- function(
    stage2.table,
    L
) {
  select_narrowest_candidate(
    stage2.table = stage2.table,
    containment.column = "estimate.in.window",
    L = L
  )
}


############################################################
# 12. SIMULATION PARAMETERS
############################################################

A <- 1
L <- 22
S <- 10
h <- 30
seed <- 280

tau.values <- c(
  0.25,
  0.50,
  0.75,
  0.95
)

SNR.values <- c(
  1000,
  100,
  1
)


############################################################
# 13. RUN THE COMPLETE PROCEDURE
#
# Requires the existing function estimate_peak_curves().
#
# For each (tau,SNR):
#   - estimate the full sequential trajectory;
#   - construct candidate intervals at every available T;
#   - select the operational candidate containing hat.l(T);
#   - select the ex-post candidate containing L;
#   - retain the full candidate and E-set diagnostics.
#
############################################################

n.pairs <- length(tau.values) * length(SNR.values)

all.selected.expost <- vector("list", n.pairs)
all.selected.estimate <- vector("list", n.pairs)
all.diagnostics <- vector("list", n.pairs)
all.E.diagnostics <- vector("list", n.pairs)

counter <- 1L

for (tau in tau.values) {
  
  cat("\n\n##################################################\n")
  cat("# tau =", tau, "\n")
  cat("##################################################\n")
  
  for (SNR in SNR.values) {
    
    cat("\n========================================\n")
    cat("tau =", tau, "\n")
    cat("SNR =", SNR, "\n")
    cat("========================================\n")
    
    # Obtain the sequential estimate trajectory.
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
    
    if (length(T.values) != length(hat.l)) {
      stop("The returned T and hat.l vectors have different lengths.")
    }
    
    if (anyDuplicated(T.values)) {
      stop("T.values must not contain duplicate horizons.")
    }
    
    ord <- order(T.values)
    T.values <- T.values[ord]
    hat.l <- hat.l[ord]
    
    current.T.values <- sort(unique(T.values))
    
    selected.expost.at.T <- vector(
      "list",
      length(current.T.values)
    )
    
    selected.estimate.at.T <- vector(
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
    
    for (iT in seq_along(current.T.values)) {
      
      T.current <- current.T.values[iT]
      
      cat("\n----------------------------------------\n")
      cat("Current observation horizon T =", T.current, "\n")
      
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
      
      stage2.table <- result$stage2
      
      stage2.table$tau_level <- tau
      stage2.table$SNR_level <- SNR
      
      # Keep the width identity explicit for auditing.
      stage2.table$w_k_check <- with(
        stage2.table,
        T_exit - T_entry
      )
      
      diagnostics.at.T[[iT]] <- stage2.table
      E.at.T[[iT]] <- result$E_sets
      
      # Ex-post selection uses the true L.
      selected.expost <- select_final_candidate(
        stage2.table = stage2.table,
        L = L
      )
      
      selected.expost$tau <- tau
      selected.expost$SNR <- SNR
      selected.expost$hat_l_current <-
        unique(stage2.table$hat_l_current)[1]
      
      selected.expost.at.T[[iT]] <- selected.expost
      
      # Operational selection uses only hat.l(T) and the
      # candidate intervals. It does not use L for selection.
      selected.estimate <- select_estimate_candidate(
        stage2.table = stage2.table,
        L = L
      )
      
      selected.estimate$tau <- tau
      selected.estimate$SNR <- SNR
      
      selected.estimate.at.T[[iT]] <- selected.estimate
      
      cat("\nOperational candidate containing hat.l(T):\n")
      print(selected.estimate)
      
      cat("\nEx-post candidate containing L:\n")
      print(selected.expost)
    }
    
    all.selected.expost[[counter]] <- do.call(
      rbind,
      selected.expost.at.T
    )
    
    all.selected.estimate[[counter]] <- do.call(
      rbind,
      selected.estimate.at.T
    )
    
    all.diagnostics[[counter]] <- do.call(
      rbind,
      diagnostics.at.T
    )
    
    all.E.diagnostics[[counter]] <- do.call(
      rbind,
      E.at.T
    )
    
    counter <- counter + 1L
  }
}


############################################################
# 14. COMBINE RESULTS ACROSS (tau,SNR)
############################################################

table.selected.expost <- do.call(
  rbind,
  all.selected.expost
)

table.selected.estimate <- do.call(
  rbind,
  all.selected.estimate
)

combined.diagnostics <- do.call(
  rbind,
  all.diagnostics
)

combined.E.diagnostics <- do.call(
  rbind,
  all.E.diagnostics
)

rownames(table.selected.expost) <- NULL
rownames(table.selected.estimate) <- NULL
rownames(combined.diagnostics) <- NULL
rownames(combined.E.diagnostics) <- NULL


############################################################
# 15. FIRST HORIZON WITH ANY CANDIDATE CONTAINING L
#
# This is the first horizon T>L at which the selected
# ex-post candidate strictly contains L.
#
# It is not the persistent detection time T_det.
#
############################################################

first_true_containment <- function(
    selected.table,
    L
) {
  
  groups <- split(
    selected.table,
    interaction(
      selected.table$tau,
      selected.table$SNR,
      drop = TRUE
    )
  )
  
  out <- lapply(
    groups,
    function(df) {
      
      df <- df[
        order(df$T_current),
        ,
        drop = FALSE
      ]
      
      eligible <- df[
        df$T_current > L &
          !is.na(df$peak.in.strict) &
          df$peak.in.strict,
        ,
        drop = FALSE
      ]
      
      if (nrow(eligible) == 0) {
        return(
          data.frame(
            tau = unique(df$tau)[1],
            SNR = unique(df$SNR)[1],
            T_first = NA_real_,
            k_first = NA_real_,
            T_entry_first = NA_real_,
            T_exit_first = NA_real_,
            w_first = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
      
      first.row <- eligible[1, , drop = FALSE]
      
      data.frame(
        tau = unique(df$tau)[1],
        SNR = unique(df$SNR)[1],
        T_first = first.row$T_current,
        k_first = first.row$k,
        T_entry_first = first.row$T_entry,
        T_exit_first = first.row$T_exit,
        w_first = first.row$w_k,
        stringsAsFactors = FALSE
      )
    }
  )
  
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}


table.first.true.containment <- first_true_containment(
  selected.table = table.selected.expost,
  L = L
)


############################################################
# 16. OPERATIONAL ESTIMATION TIME T_est
#
# T_est is the first available horizon at which at least
# one candidate interval contains the current point estimate:
#
#   T_entry^(k)(T) <= hat.l(T) <= T_exit^(k)(T).
#
# The narrowest qualifying interval is selected at each T.
# T_est is the earliest horizon with such a selected interval.
#
############################################################

compute_T_est <- function(
    selected.table
) {
  
  groups <- split(
    selected.table,
    interaction(
      selected.table$tau,
      selected.table$SNR,
      drop = TRUE
    )
  )
  
  out <- lapply(
    groups,
    function(df) {
      
      df <- df[
        order(df$T_current),
        ,
        drop = FALSE
      ]
      
      eligible <- df[
        !is.na(df$estimate.in.window) &
          df$estimate.in.window &
          !is.na(df$T_entry) &
          !is.na(df$T_exit) &
          !is.na(df$w_k),
        ,
        drop = FALSE
      ]
      
      if (nrow(eligible) == 0) {
        return(
          data.frame(
            tau = unique(df$tau)[1],
            SNR = unique(df$SNR)[1],
            T_est = NA_real_,
            k_est = NA_real_,
            hat_l_at_T_est = NA_real_,
            T_entry_est = NA_real_,
            T_exit_est = NA_real_,
            w_est = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
      
      first.row <- eligible[1, , drop = FALSE]
      
      data.frame(
        tau = unique(df$tau)[1],
        SNR = unique(df$SNR)[1],
        T_est = first.row$T_current,
        k_est = first.row$k,
        hat_l_at_T_est = first.row$hat_l_current,
        T_entry_est = first.row$T_entry,
        T_exit_est = first.row$T_exit,
        w_est = first.row$w_k,
        stringsAsFactors = FALSE
      )
    }
  )
  
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}


table.T_est <- compute_T_est(
  selected.table = table.selected.estimate
)


############################################################
# 17. PERSISTENT EX-POST DETECTION TIME T_det
#
# T_det is the earliest available horizon T>L such that
# the selected ex-post candidate strictly contains L at
# EVERY available horizon T' >= T, up to the terminal
# horizon h.
#
# In other words, a single later failure of containment
# disqualifies an earlier horizon from being T_det.
#
############################################################

compute_T_det <- function(
    selected.table,
    L
) {
  
  groups <- split(
    selected.table,
    interaction(
      selected.table$tau,
      selected.table$SNR,
      drop = TRUE
    )
  )
  
  out <- lapply(
    groups,
    function(df) {
      
      df <- df[
        order(df$T_current),
        ,
        drop = FALSE
      ]
      
      # The detection definition is considered only after
      # the true peak location L.
      post.peak <- df[
        df$T_current > L,
        ,
        drop = FALSE
      ]
      
      if (nrow(post.peak) == 0) {
        return(
          data.frame(
            tau = unique(df$tau)[1],
            SNR = unique(df$SNR)[1],
            T_det = NA_real_,
            k_det = NA_real_,
            T_entry_det = NA_real_,
            T_exit_det = NA_real_,
            w_det = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
      
      contains.true.peak <- !is.na(post.peak$peak.in.strict) &
        post.peak$peak.in.strict
      
      # A horizon qualifies only if every subsequent available
      # horizon, including itself, has strict containment.
      persistent.from.horizon <- vapply(
        seq_len(nrow(post.peak)),
        function(i) {
          all(contains.true.peak[i:nrow(post.peak)])
        },
        logical(1)
      )
      
      qualifying.index <- which(persistent.from.horizon)
      
      if (length(qualifying.index) == 0) {
        return(
          data.frame(
            tau = unique(df$tau)[1],
            SNR = unique(df$SNR)[1],
            T_det = NA_real_,
            k_det = NA_real_,
            T_entry_det = NA_real_,
            T_exit_det = NA_real_,
            w_det = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
      
      i <- qualifying.index[1]
      det.row <- post.peak[i, , drop = FALSE]
      
      data.frame(
        tau = unique(df$tau)[1],
        SNR = unique(df$SNR)[1],
        T_det = det.row$T_current,
        k_det = det.row$k,
        T_entry_det = det.row$T_entry,
        T_exit_det = det.row$T_exit,
        w_det = det.row$w_k,
        stringsAsFactors = FALSE
      )
    }
  )
  
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}


table.T_det <- compute_T_det(
  selected.table = table.selected.expost,
  L = L
)


############################################################
# 18. TERMINAL-HORIZON CANDIDATES
#
# These are the selected candidates at the largest available
# horizon for each (tau,SNR) pair.
############################################################

terminal_candidate_table <- function(
    selected.table,
    prefix = ""
) {
  
  groups <- split(
    selected.table,
    interaction(
      selected.table$tau,
      selected.table$SNR,
      drop = TRUE
    )
  )
  
  out <- lapply(
    groups,
    function(df) {
      
      df <- df[
        order(df$T_current),
        ,
        drop = FALSE
      ]
      
      terminal <- df[nrow(df), , drop = FALSE]
      
      data.frame(
        tau = terminal$tau,
        SNR = terminal$SNR,
        T_terminal = terminal$T_current,
        k_terminal = terminal$k,
        T_entry_terminal = terminal$T_entry,
        T_exit_terminal = terminal$T_exit,
        w_terminal = terminal$w_k,
        hat_l_terminal = terminal$hat_l_current,
        stringsAsFactors = FALSE
      )
    }
  )
  
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}


table.terminal.expost <- terminal_candidate_table(
  selected.table = table.selected.expost
)

table.terminal.estimate <- terminal_candidate_table(
  selected.table = table.selected.estimate
)


############################################################
# 19. COMPARISON TABLE: T_est, T_first, T_det
#
# T_est:
#   Operational criterion based on hat.l(T).
#
# T_first:
#   First horizon at which the selected ex-post candidate
#   strictly contains the true L.
#
# T_det:
#   First horizon from which the selected ex-post candidate
#   continues to contain L at all subsequent available
#   horizons through the terminal horizon.
############################################################

table.comparison <- merge(
  table.T_est,
  table.first.true.containment,
  by = c("tau", "SNR"),
  all = TRUE,
  sort = FALSE
)

table.comparison <- merge(
  table.comparison,
  table.T_det,
  by = c("tau", "SNR"),
  all = TRUE,
  sort = FALSE
)

table.comparison <- merge(
  table.comparison,
  table.terminal.expost[
    ,
    c(
      "tau", "SNR", "T_terminal",
      "k_terminal", "T_entry_terminal",
      "T_exit_terminal", "w_terminal"
    )
  ],
  by = c("tau", "SNR"),
  all = TRUE,
  sort = FALSE
)

table.comparison <- table.comparison[
  order(table.comparison$tau, -table.comparison$SNR),
  ,
  drop = FALSE
]

rownames(table.comparison) <- NULL


############################################################
# 20. DISPLAY COMPARISON TABLE
############################################################

cat("\n\n========================================\n")
cat("COMPARISON: T_est, T_first, AND T_det\n")
cat("========================================\n\n")

print(table.comparison)


############################################################
# 21. DISPLAY SELECTED OPERATIONAL CANDIDATES
############################################################

cat("\n\n========================================\n")
cat("SELECTED OPERATIONAL CANDIDATES\n")
cat("========================================\n\n")

print(
  table.selected.estimate[
    ,
    c(
      "tau", "SNR", "T_current", "k",
      "hat_l_current", "T_entry", "T_exit",
      "w_k", "estimate.in.window"
    )
  ]
)


############################################################
# 22. DISPLAY SELECTED EX-POST CANDIDATES
############################################################

cat("\n\n========================================\n")
cat("SELECTED EX-POST CANDIDATES\n")
cat("========================================\n\n")

print(
  table.selected.expost[
    ,
    c(
      "tau", "SNR", "T_current", "k",
      "hat_l_current", "T_entry", "T_exit",
      "w_k", "peak.in.strict"
    )
  ]
)


############################################################
# 23. DISPLAY ALL STAGE-2 CANDIDATES
############################################################

cat("\n\n========================================\n")
cat("ALL STAGE-2 CANDIDATES\n")
cat("========================================\n\n")

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
      "hat_l_current",
      "estimate.in.window",
      "peak.in.strict"
    )
  ]
)


############################################################
# 24. FIRST PEAK-CONTAINING k BY (tau,SNR,T_current)
#
# These diagnostics preserve the original two variants:
#
#   first_k_incl0: smallest k >= 0 strictly containing L.
#   first_k_excl0: smallest k >= 1 strictly containing L.
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
      
      valid.incl0 <- df[
        df$peak.in.strict,
        ,
        drop = FALSE
      ]
      
      valid.excl0 <- df[
        df$k >= 1 &
          df$peak.in.strict,
        ,
        drop = FALSE
      ]
      
      data.frame(
        tau = unique(df$tau_level),
        SNR = unique(df$SNR_level),
        T_current = unique(df$T_current),
        
        first_k_incl0 = if (nrow(valid.incl0) > 0)
          min(valid.incl0$k) else NA_real_,
        
        first_k_excl0 = if (nrow(valid.excl0) > 0)
          min(valid.excl0$k) else NA_real_,
        
        stringsAsFactors = FALSE
      )
    }
  )
)

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

cat("\n\n========================================\n")
cat("FIRST PEAK-CONTAINING k BY (tau,SNR,T_current)\n")
cat("========================================\n\n")

print(first_k_table)


############################################################
# 25. EX-POST INTERPRETATION OF SELECTED CANDIDATES
############################################################

table.selected.expost$relation <- with(
  table.selected.expost,
  
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

table.selected.expost$interpretation <- with(
  table.selected.expost,
  
  ifelse(
    is.na(T_entry) | is.na(T_exit),
    "No valid peak-containing candidate",
    
    ifelse(
      T_exit < L,
      "Pre-peak candidate",
      
      ifelse(
        T_entry <= L & L <= T_exit,
        "Candidate contains L inclusively",
        "Post-peak candidate"
      )
    )
  )
)


############################################################
# 26. SAVE CSV FILES
############################################################

# write.csv(
#   table.selected.expost,
#   "entry-exit-selected-expost-all-tau-SNR-T.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   table.selected.estimate,
#   "entry-exit-selected-operational-all-tau-SNR-T.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   combined.diagnostics,
#   "simulation-all-stage2-candidates-T-dependent.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   combined.E.diagnostics,
#   "simulation-all-E-candidates-T-dependent.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   first_k_table,
#   "simulation-first-k-summary-T-dependent.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   table.T_est,
#   "simulation-T-est-summary.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   table.first.true.containment,
#   "simulation-T-first-true-containment-summary.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   table.T_det,
#   "simulation-T-det-summary.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   table.comparison,
#   "simulation-comparison-T-est-T-first-T-det.csv",
#   row.names = FALSE
# )


############################################################
# 27. LATEX-READY COMPARISON TABLE
############################################################

cat("\n\n========================================\n")
cat("LATEX TABLE — T_est, T_first, T_det\n")
cat("========================================\n\n")

print(
  knitr::kable(
    table.comparison,
    format = "latex",
    booktabs = TRUE,
    digits = 6
  )
)


############################################################
# 28. LATEX-READY TERMINAL CANDIDATE TABLES
############################################################

cat("\n\n========================================\n")
cat("LATEX TABLE — TERMINAL EX-POST CANDIDATES\n")
cat("========================================\n\n")

print(
  knitr::kable(
    table.terminal.expost,
    format = "latex",
    booktabs = TRUE,
    digits = 6
  )
)

cat("\n\n========================================\n")
cat("LATEX TABLE — TERMINAL OPERATIONAL CANDIDATES\n")
cat("========================================\n\n")

print(
  knitr::kable(
    table.terminal.estimate,
    format = "latex",
    booktabs = TRUE,
    digits = 6
  )
)


############################################################
# 29. OPTIONAL: SUMMARY FOR EACH TAU
############################################################

for (tau in tau.values) {
  
  cat("\n\n========================================\n")
  cat("COMPARISON TABLE FOR tau =", tau, "\n")
  cat("========================================\n\n")
  
  tau.table <- table.comparison[
    table.comparison$tau == tau,
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


############################################################
# END OF SCRIPT
############################################################


############################################################
# FIXED-HORIZON ENTRY / EXIT CONSTRUCTION (terminal horizon h)
# + RETROSPECTIVE SELECTION, strict and inclusive containment
#
# Shared by the China, France and Germany scripts.
# L is NOT used to build the candidates; only in selection.
############################################################

library(dplyr)

STAGE1_TOLERANCE <- 1e-10
E_TOLERANCE      <- 0
STAGE2_TOLERANCE <- 0


############################################################
# 1. INTERVAL OSCILLATION
############################################################

interval_oscillation <- function(T.values, hat.l, T1, T2) {
  if (length(T.values) != length(hat.l))
    stop("T.values and hat.l must have the same length.")
  if (length(T.values) == 0) return(NA_real_)
  idx <- which(T.values >= T1 & T.values <= T2)
  if (length(idx) == 0) return(NA_real_)
  if (length(idx) == 1) return(0)
  max(abs(outer(hat.l[idx], hat.l[idx], "-")))
}


############################################################
# 2. FUTURE-TAIL OSCILLATION
############################################################

compute_tail_oscillation <- function(T.values, hat.l) {
  h <- max(T.values)
  T.grid <- sort(unique(T.values))
  omega <- sapply(T.grid, function(T)
    interval_oscillation(T.values, hat.l, T, h))
  data.frame(T = T.grid, omega = omega)
}


############################################################
# 3. STAGE 1: recursive construction of (tau'_k, T_exit^(k))
############################################################

stage1_tau_prime <- function(T.values, hat.l,
                             tolerance = STAGE1_TOLERANCE) {
  h <- max(T.values)
  T.grid <- sort(unique(T.values))
  tail.table <- compute_tail_oscillation(T.values, hat.l)
  
  candidates <- data.frame(step = 0, tau_prime = 0,
                           T_exit = h, T_boundary = h)
  tau.current <- 0
  boundary.current <- h
  
  for (step in seq_len(length(T.grid) + 1)) {
    
    earlier.T <- T.grid[T.grid < boundary.current]
    if (length(earlier.T) == 0) break
    
    possible.T <- tail.table$T[
      !is.na(tail.table$omega) &
        abs(tail.table$omega - tau.current) <= tolerance &
        tail.table$T < boundary.current]
    
    if (length(possible.T) > 0) {            # CASE 1
      T.prime <- min(possible.T)
      previous.T <- T.grid[T.grid < T.prime]
      if (length(previous.T) == 0) break
      T.previous <- max(previous.T)
      T.boundary.next <- T.prime
    } else {                                 # CASE 2
      T.previous <- max(earlier.T)
      T.boundary.next <- T.previous
    }
    
    T.exit.next <- T.previous
    tau.next <- interval_oscillation(T.values, hat.l, T.previous, h)
    
    if (is.na(tau.next)) break
    if (tau.next <= tau.current + tolerance) break
    
    candidates <- rbind(candidates, data.frame(
      step = step, tau_prime = tau.next,
      T_exit = T.exit.next, T_boundary = T.boundary.next))
    
    tau.current <- tau.next
    boundary.current <- T.boundary.next
  }
  
  candidates$omega_at_exit <- sapply(candidates$T_exit, function(Te)
    interval_oscillation(T.values, hat.l, Te, h))
  candidates$definition_check <-
    abs(candidates$tau_prime - candidates$omega_at_exit) <= tolerance
  
  list(h = h, tail.table = tail.table, tau.prime.table = candidates)
}


############################################################
# 4. STAGE 2: E_k, T_k^*, epsilon_k, T_entry^(k), w_k, |Z_k|
############################################################

stage2_one_candidate <- function(T.values, hat.l, tau.prime, T.exit,
                                 tolerance = STAGE2_TOLERANCE) {
  
  empty <- list(E.k = numeric(0), T.star = NA_real_, tau = NA_real_,
                T.entry = NA_real_, w.star = NA_real_, n.Z = NA_integer_)
  
  candidate.T <- sort(unique(T.values[T.values < T.exit]))
  if (length(candidate.T) == 0) return(empty)
  
  omega.values <- sapply(candidate.T, function(Te)
    interval_oscillation(T.values, hat.l, Te, T.exit))
  valid <- !is.na(omega.values)
  T.valid <- candidate.T[valid]
  omega.valid <- omega.values[valid]
  if (length(T.valid) == 0) return(empty)
  
  E.k <- T.valid[omega.valid >= tau.prime - tolerance]
  if (length(E.k) == 0) return(empty)
  
  T.star <- max(E.k)
  tau <- omega.valid[which(T.valid == T.star)[1]]
  
  T.equal <- T.valid[abs(omega.valid - tau) <= tolerance]
  if (length(T.equal) == 0) {
    out <- empty; out$E.k <- E.k; out$T.star <- T.star; out$tau <- tau
    return(out)
  }
  
  T.entry <- min(T.equal)
  list(E.k = E.k, T.star = T.star, tau = tau, T.entry = T.entry,
       w.star = T.exit - T.entry, n.Z = length(T.equal))
}

stage2_all_candidates <- function(T.values, hat.l, stage1.result,
                                  tolerance = STAGE2_TOLERANCE) {
  tp <- stage1.result$tau.prime.table
  do.call(rbind, lapply(seq_len(nrow(tp)), function(i) {
    S2 <- stage2_one_candidate(T.values, hat.l, tp$tau_prime[i],
                               tp$T_exit[i], tolerance)
    data.frame(
      step = tp$step[i], tau_prime = tp$tau_prime[i],
      T_exit = tp$T_exit[i],
      E_k = paste(S2$E.k, collapse = ", "),
      E_k_min = if (length(S2$E.k)) min(S2$E.k) else NA_real_,
      E_k_max = if (length(S2$E.k)) max(S2$E.k) else NA_real_,
      n_E_k = length(S2$E.k), E_empty = length(S2$E.k) == 0,
      T_star = S2$T.star, tau = S2$tau, T_entry = S2$T.entry,
      w_star = S2$w.star, n_Z = S2$n.Z, stringsAsFactors = FALSE)
  }))
}

run_entry_exit <- function(T.values, hat.l,
                           tolerance.stage1 = STAGE1_TOLERANCE,
                           tolerance.stage2 = STAGE2_TOLERANCE) {
  S1 <- stage1_tau_prime(T.values, hat.l, tolerance.stage1)
  S2 <- stage2_all_candidates(T.values, hat.l, S1, tolerance.stage2)
  list(stage1 = S1, stage2 = S2)
}


############################################################
# 5. SEQUENTIAL PATH (uses estimate_params() defined in the
#    country script)
############################################################

estimate_peak_path <- function(I_full, T_values, tau, h) {
  hat.l <- sapply(T_values, function(T)
    estimate_params(I = I_full[1:T], tau = tau, h = h)$l_hat)
  data.frame(T = T_values, hat.l = hat.l)
}


############################################################
# 6. RETROSPECTIVE SELECTION
#
# Containment of the closed interval [T_entry, T_exit]:
#   inclusive = TRUE  :  T_entry <= L <= T_exit   (manuscript)
#   inclusive = FALSE :  T_entry <  L <  T_exit   (old; check only)
# Narrowest containing candidate; ties: smaller step.
# NOT part of the operational procedure (uses L).
############################################################

select_expost_candidate <- function(stage2_table, L, inclusive = TRUE) {
  
  candidates <- stage2_table %>%
    filter(!is.na(T_entry), !is.na(T_exit), !is.na(w_star))
  
  peak_containing <- if (inclusive) {
    candidates %>% filter(T_entry <= L, L <= T_exit)
  } else {
    candidates %>% filter(T_entry <  L, L <  T_exit)
  }
  
  if (nrow(peak_containing) == 0) {
    return(data.frame(
      step = NA_integer_, tau_prime = NA_real_, T_exit = NA_real_,
      T_star = NA_real_, tau = NA_real_, T_entry = NA_real_,
      w_star = NA_real_, n_Z = NA_integer_, L = L,
      peak_contained = FALSE))
  }
  
  selected <- peak_containing %>% arrange(w_star, step) %>% slice(1)
  
  data.frame(
    step = selected$step, tau_prime = selected$tau_prime,
    T_exit = selected$T_exit, T_star = selected$T_star,
    tau = selected$tau, T_entry = selected$T_entry,
    w_star = selected$w_star, n_Z = selected$n_Z, L = L,
    peak_contained = TRUE)
}

apply_selection <- function(all_stage2, L, inclusive) {
  bind_rows(lapply(names(all_stage2), function(tau.name) {
    sel <- select_expost_candidate(all_stage2[[tau.name]], L, inclusive)
    sel$tau_level <- as.numeric(tau.name)
    sel
  })) %>%
    select(tau_level, step, tau_prime, T_exit, T_star, tau,
           T_entry, w_star, n_Z, L, peak_contained) %>%
    arrange(tau_level)
}


############################################################
# 7. RUN ONE COUNTRY
############################################################

run_country_fixed_h <- function(name, I_full, T_vals, L_observed,
                                taus = c(0.25, 0.50, 0.75, 0.95),
                                use_cache = TRUE) {
  
  h <- max(T_vals)
  stopifnot(length(I_full) >= h)
  prefix <- tolower(name)
  
  all_paths <- list(); all_stage1 <- list(); all_stage2 <- list()
  
  for (tau in taus) {
    
    cat("\n========================================\n")
    cat(name, ": TAU =", tau, "\n")
    cat("========================================\n")
    
    f <- sprintf("%s_path_tau%.2f.rds", prefix, tau)
    path <- if (use_cache && file.exists(f)) {
      readRDS(f)
    } else {
      p <- estimate_peak_path(I_full, T_vals, tau, h)
      saveRDS(p, f); p
    }
    
    res <- run_entry_exit(path$T, path$hat.l)
    
    key <- as.character(tau)
    all_paths[[key]]  <- path
    all_stage1[[key]] <- res$stage1$tau.prime.table
    all_stage2[[key]] <- res$stage2
    
    cat("\nStage-1 candidates:\n"); print(res$stage1$tau.prime.table)
    cat("\nStage-2 candidates:\n"); print(res$stage2)
  }
  
  combine <- function(lst) bind_rows(lapply(names(lst), function(nm) {
    d <- lst[[nm]]; d$tau <- as.numeric(nm); d
  })) %>% select(tau, everything())
  
  stage1_table <- combine(all_stage1)
  stage2_table <- combine(all_stage2)
  
  write.csv(stage1_table, sprintf("%s_section8_stage1_candidates.csv", prefix),
            row.names = FALSE)
  write.csv(stage2_table, sprintf("%s_section8_stage2_candidates.csv", prefix),
            row.names = FALSE)
  
  selected_strict    <- apply_selection(all_stage2, L_observed, FALSE)
  selected_inclusive <- apply_selection(all_stage2, L_observed, TRUE)
  
  write.csv(selected_strict,
            sprintf("%s_section8_selected_candidates_strict.csv", prefix),
            row.names = FALSE)
  write.csv(selected_inclusive,
            sprintf("%s_section8_selected_candidates_inclusive.csv", prefix),
            row.names = FALSE)
  
  cat("\n========================================\n")
  cat(name, ": observed peak L =", L_observed, "\n")
  cat("\nSTRICT (must reproduce the current manuscript table):\n")
  print(selected_strict)
  cat("\nINCLUSIVE (new convention):\n")
  print(selected_inclusive)
  
  compact <- selected_inclusive %>%
    transmute(tau = tau_level, k = step, T_entry, T_exit,
              width = w_star, n_Z, L, peak_contained)
  cat("\nCOMPACT MANUSCRIPT TABLE (inclusive):\n")
  print(compact)
  
  invisible(list(paths = all_paths, stage1 = stage1_table,
                 stage2 = stage2_table, strict = selected_strict,
                 inclusive = selected_inclusive))
}

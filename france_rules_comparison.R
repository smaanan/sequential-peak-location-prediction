############################################################
# FRANCE: OPERATIONAL INTERVALS WITH BOTH SELECTION RULES
#
# Data:   ECDC daily cases, 2020-01-15 to 2020-05-31
# T:      50,...,89      h = 89      L_observed = 78
# s range for the width: [1, 40]
#
# Order: (1) data + estimator, (2) source the functions-only
# file, (3) rule-aware selection functions, (4) analysis.
############################################################

library(dplyr)
library(readr)
library(knitr)

############################################################
# 1. PARAMETERS
############################################################

taus   <- c(0.25, 0.50, 0.75, 0.95)
T_vals <- 4:89
h <- max(T_vals)

STAGE1_TOLERANCE <- 1e-10
E_TOLERANCE      <- 0
STAGE2_TOLERANCE <- 0

L_observed <- 78     # observed peak: used ONLY to check containment
rules <- c("narrowest", "first_k")

############################################################
# 2. DATA
############################################################

get_france_cases_ecdc <- function(
    file = "6YbQxT6Z",
    start_date = "2020-01-15",
    end_date   = "2020-05-31"
) {
  dat <- readr::read_csv(file, col_types = readr::cols(),
                         progress = FALSE)
  dat %>%
    dplyr::filter(countriesAndTerritories == "France") %>%
    dplyr::mutate(
      date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))
    ) %>%
    dplyr::filter(date >= as.Date(start_date),
                  date <= as.Date(end_date)) %>%
    dplyr::arrange(date) %>%
    dplyr::select(date, daily_cases = cases)
}

france <- get_france_cases_ecdc()
I_full <- france$daily_cases
stopifnot(length(I_full) >= h)

############################################################
# 3. ESTIMATOR (as in the France scripts, sections 4-9)
############################################################

quantile_loss <- function(u, tau) u * (tau - (u < 0))

med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  finite <- is.finite(x) & is.finite(w)
  if (!any(finite)) return(NA_real_)
  x <- x[finite]; w <- w[finite]
  if (any(w < 0) || sum(w) <= 0) return(NA_real_)
  o <- order(x)
  cum_w <- cumsum(w[o])
  cutoff <- tau * sum(w)
  x[o][min(which(cum_w >= cutoff))]
}

compute_a <- function(I, l, s, tau) med_weighted_tau(I, l, s, tau)

objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  mean(quantile_loss(I - fitted, tau))
}

search_best_s <- function(I, l, tau) {
  f <- function(s) {
    a <- compute_a(I, l, s, tau)
    objective_function(I, a, l, s, tau)
  }
  optimize(f, c(1, 40), tol = 1e-4)$minimum     # France: [1,40]
}

estimate_params <- function(I, tau, h) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau)
    a <- compute_a(I, l, s, tau)
    objective_function(I, a, l, s, tau)
  })
  list(l_hat = which.min(losses), loss = min(losses))
}

estimate_peak_path <- function(I_full, T_values, tau, h) {
  hat.l <- sapply(T_values, function(T) {
    estimate_params(I_full[1:T], tau, h)$l_hat
  })
  data.frame(T = T_values, hat.l = hat.l)
}

############################################################
# 4. CURRENT-HORIZON ENTRY/EXIT FUNCTIONS
#    (functions only: no data, no estimator)
############################################################

source("oct_09.R")
# source("entry_exit_functions.R")

############################################################
# 5. RULE-AWARE SELECTION FUNCTIONS
#    Defined AFTER source(), so they override older versions.
#
#   "narrowest": argmin_{k in K(T)} w_k(T)
#                (ties: smaller k, T_entry, T_exit)
#   "first_k"  : min K(T)
############################################################

select_narrowest_candidate <- function(
    stage2.table,
    containment.column,
    L,
    rule = c("narrowest", "first_k")
) {
  
  rule <- match.arg(rule)
  
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
  
  T.current     <- unique(stage2.table$T_current)[1]
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
  
  ord <- if (rule == "narrowest") {
    order(candidates$w_k, candidates$k,
          candidates$T_entry, candidates$T_exit)
  } else {
    order(candidates$k)
  }
  candidates <- candidates[ord, , drop = FALSE]
  
  out.cols <- c(
    "T_current", "k", "epsilon_prime", "T_exit", "T_star",
    "epsilon", "T_entry", "w_k", "L", "hat_l_current",
    "peak.in.strict", "estimate.in.window"
  )
  
  selected <- candidates[1, out.cols, drop = FALSE]
  rownames(selected) <- NULL
  selected
}

select_final_candidate <- function(stage2.table, L,
                                   rule = "narrowest") {
  select_narrowest_candidate(stage2.table, "peak.in.strict",
                             L, rule = rule)
}

select_estimate_candidate <- function(stage2.table, L,
                                      rule = "narrowest") {
  select_narrowest_candidate(stage2.table, "estimate.in.window",
                             L, rule = rule)
}

# Safety checks: stop early if an old version is still active.
stopifnot("rule" %in% names(formals(select_estimate_candidate)))
stopifnot("T.current" %in% names(formals(run_entry_exit_at_T)))

############################################################
# 6. RUN: PATHS, CANDIDATES, BOTH RULES AT EVERY T
############################################################

rows  <- list()
paths <- list()
cands <- list()

for (tau in taus) {
  
  cat("tau =", tau, "\n")
  
  # Cached: the path is the expensive part (40 horizons x 89 l).
  file <- sprintf("france_path_tau%.2f.rds", tau)
  path <- if (file.exists(file)) readRDS(file) else {
    p <- estimate_peak_path(I_full, T_vals, tau, h)
    saveRDS(p, file); p
  }
  paths[[as.character(tau)]] <- transform(path, tau = tau)
  
  T.values <- path$T
  hat.l    <- path$hat.l
  
  for (Tc in T.values) {
    
    r <- run_entry_exit_at_T(
      T.values, hat.l, Tc, Data = "France", L = L_observed,
      tolerance.stage1 = STAGE1_TOLERANCE,
      tolerance.E = E_TOLERANCE,
      tolerance.stage2 = STAGE2_TOLERANCE
    )
    s2 <- r$stage2
    s2$tau <- tau
    cands[[length(cands) + 1L]] <- s2
    
    for (rl in rules) {
      sel <- select_estimate_candidate(s2, L_observed, rule = rl)
      ok  <- !is.na(sel$T_entry) && !is.na(sel$T_exit)
      rows[[length(rows) + 1L]] <- data.frame(
        tau = tau, rule = rl, T = Tc,
        hat_l = sel$hat_l_current,
        selected = ok,
        k = if (ok) sel$k else NA_real_,
        T_entry = if (ok) sel$T_entry else NA_real_,
        T_exit  = if (ok) sel$T_exit  else NA_real_,
        width   = if (ok) sel$T_exit - sel$T_entry else NA_real_,
        contains_L_strict =
          ok && sel$T_entry < L_observed && L_observed < sel$T_exit,
        contains_L_incl =
          ok && sel$T_entry <= L_observed && L_observed <= sel$T_exit,
        stringsAsFactors = FALSE
      )
    }
  }
}

france_sel   <- do.call(rbind, rows)
# cat("Intervals selected:", sum(france_sel$selected), "\n")
france_cands <- do.call(rbind, cands)
france_paths <- do.call(rbind, paths)

write.csv(france_sel, "france_operational_both_rules_all_T.csv",
          row.names = FALSE)
write.csv(france_cands,
          "france_all_stage2_candidates_operational.csv",
          row.names = FALSE)
write.csv(france_paths, "france_paths_all_tau.csv",
          row.names = FALSE)

############################################################
# 7. T_est (does not depend on the rule)
############################################################

T_est <- do.call(rbind, lapply(taus, function(tt) {
  d <- france_sel[france_sel$rule == "narrowest" &
                    france_sel$tau == tt & france_sel$selected, ]
  data.frame(tau = tt,
             T_est = if (nrow(d)) min(d$T) else NA_real_)
}))
cat("\nFirst operational selection horizon T_est:\n")
print(T_est, row.names = FALSE)

############################################################
# 8. INTERVALS AT T_est AND AT THE TERMINAL HORIZON
############################################################

at_Test <- do.call(rbind, lapply(seq_len(nrow(T_est)), function(i) {
  if (is.na(T_est$T_est[i])) return(NULL)
  france_sel[france_sel$tau == T_est$tau[i] &
               france_sel$T == T_est$T_est[i], ]
}))

at_Tmax <- france_sel[france_sel$T == max(T_vals), ]

cols <- c("tau", "rule", "T", "hat_l", "k",
          "T_entry", "T_exit", "width", "contains_L_strict")

cat("\nAt T_est:\n")
if (!is.null(at_Test)) print(at_Test[, cols], row.names = FALSE)

cat("\nAt T = ", max(T_vals), ":\n", sep = "")
print(at_Tmax[, cols], row.names = FALSE)

print(knitr::kable(at_Tmax[, cols], format = "latex",
                   booktabs = TRUE, digits = 3))

write.csv(at_Tmax[, cols], "france_both_rules_terminal.csv",
          row.names = FALSE)

############################################################
# 9. CHECK: AVAILABILITY OF AN INTERVAL IS RULE INDEPENDENT
############################################################

chk <- tapply(france_sel$selected,
              list(france_sel$tau, france_sel$T, france_sel$rule),
              sum)
stopifnot(all(chk[, , "narrowest"] == chk[, , "first_k"]))
cat("\nCheck passed: selection availability is identical.\n")

############################################################
# 10. FIGURE: hat.l(T) AND SELECTED INTERVALS
############################################################

png("france_intervals_both_rules.png",
    width = 1600, height = 1200, res = 150)
old.par <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1))

for (tau in taus) {
  p <- france_paths[france_paths$tau == tau, ]
  plot(p$T, p$hat.l, type = "o", pch = 16, cex = 0.6,
       xlim = range(T_vals), ylim = c(0, h),
       xlab = expression(T), ylab = expression(hat(l)(T)),
       main = bquote(tau == .(tau)))
  abline(h = L_observed, lty = 3)
  abline(0, 1, col = "grey80")         # T_exit <= T boundary
  
  for (j in seq_along(rules)) {
    d <- france_sel[france_sel$tau == tau &
                      france_sel$rule == rules[j] &
                      france_sel$selected, ]
    if (nrow(d) == 0) next
    off <- c(-0.15, 0.15)[j]
    segments(d$T + off, d$T_entry, d$T + off, d$T_exit,
             col = c("blue", "red")[j], lwd = 1.5)
  }
  legend("topleft", c("narrowest", "first k"),
         col = c("blue", "red"), lwd = 1.5, bty = "n", cex = 0.8)
}

par(old.par)
dev.off()

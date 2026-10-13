############################################################
# CHINA: OPERATIONAL INTERVALS WITH BOTH SELECTION RULES
# Requires: entry_exit_functions.R, and the China data and
# estimator functions (sections 1-9 of the China script):
#   I_full, T_vals, taus, estimate_peak_path()
############################################################

source("oct_09.R")

############################################################
# RULE-AWARE SELECTION FUNCTIONS
# Defined AFTER source(), so they override any older
# versions (sections 9-11 of the original script).
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

# Safety check: stop early if the override did not work.
stopifnot("rule" %in% names(formals(select_estimate_candidate)))

STAGE1_TOLERANCE <- 1e-10
E_TOLERANCE      <- 0
STAGE2_TOLERANCE <- 0

L_observed <- 22        # observed peak: used ONLY to check containment
rules <- c("narrowest", "first_k")

rows  <- list()
paths <- list()
cands <- list()

for (tau in taus) {
  
  cat("tau =", tau, "\n")
  
  file <- sprintf("china_path_tau%.2f.rds", tau)
  path <- if (file.exists(file)) readRDS(file) else {
    p <- estimate_peak_path(I_full, T_vals, tau)
    saveRDS(p, file); p
  }
  paths[[as.character(tau)]] <- transform(path, tau = tau)
  
  T.values <- path$T
  hat.l    <- path$hat.l
  
  for (Tc in T.values) {
    
    r  <- run_entry_exit_at_T(
      T.values, hat.l, Tc, Data = "China", L = L_observed,
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

china_sel   <- do.call(rbind, rows)
china_cands <- do.call(rbind, cands)
china_paths <- do.call(rbind, paths)

write.csv(china_sel,   "china_operational_both_rules_all_T.csv",
          row.names = FALSE)
write.csv(china_cands, "china_all_stage2_candidates_operational.csv",
          row.names = FALSE)

# ---- T_est (does not depend on the rule) -----------------
T_est <- do.call(rbind, lapply(split(
  china_sel[china_sel$rule == "narrowest", ], ~ tau),
  function(d) {
    d <- d[d$selected, ]
    data.frame(tau = d$tau[1] ,
               T_est = if (nrow(d)) min(d$T) else NA_real_)
  }))
rownames(T_est) <- NULL
cat("\nFirst operational selection horizon T_est:\n")
print(T_est)

# ---- Intervals at T_est and at the terminal horizon ------
pick <- function(d, Tsel) d[d$T == Tsel, ]

at_Test <- do.call(rbind, lapply(seq_len(nrow(T_est)), function(i) {
  if (is.na(T_est$T_est[i])) return(NULL)
  d <- china_sel[china_sel$tau == T_est$tau[i], ]
  pick(d, T_est$T_est[i])
}))

at_Tmax <- china_sel[china_sel$T == max(T_vals), ]

cols <- c("tau", "rule", "T", "hat_l", "k",
          "T_entry", "T_exit", "width", "contains_L_strict")

cat("\nAt T_est:\n");  print(at_Test[, cols], row.names = FALSE)
cat("\nAt T = ", max(T_vals), ":\n", sep = "")
print(at_Tmax[, cols], row.names = FALSE)

print(knitr::kable(at_Tmax[, cols], format = "latex",
                   booktabs = TRUE, digits = 3))

# ---- Check: selection availability is rule independent ---
chk <- tapply(china_sel$selected,
              list(china_sel$tau, china_sel$T, china_sel$rule),
              sum)
stopifnot(all(chk[, , "narrowest"] == chk[, , "first_k"]))

# ---- Figure: hat.l(T) with the selected intervals --------
# pdf("china_intervals_both_rules.pdf", width = 10, height = 8)
png("china_intervals_both_rules.png", width = 1600, height = 1200, res = 150)
# ... same plotting code ...
# dev.off()
old.par <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1))

for (tau in taus) {
  p <- china_paths[china_paths$tau == tau, ]
  plot(p$T, p$hat.l, type = "o", pch = 16, cex = 0.6,
       xlim = range(T_vals), ylim = c(0, max(T_vals)),
       xlab = expression(T), ylab = expression(hat(l)(T)),
       main = bquote(tau == .(tau)))
  abline(h = L_observed, lty = 3)
  abline(0, 1, col = "grey80")         # T_exit <= T boundary
  
  for (j in seq_along(rules)) {
    d <- china_sel[china_sel$tau == tau &
                     china_sel$rule == rules[j] &
                     china_sel$selected, ]
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

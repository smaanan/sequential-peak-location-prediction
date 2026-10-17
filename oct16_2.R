tab.rows <- list()

for (tau in tau.values) {
  for (snr in SNR.values) {
    
    res <- estimate_peak_curves(A = A, L = L, S = S, h = h,
                                SNR = snr, tau = tau, seed = 280)
    ord <- order(res$T)
    T.values <- res$T[ord]; hat.l <- res$hat.l[ord]
    
    T.est <- NA_real_; s2.est <- NULL
    for (Tc in T.values) {
      s2 <- run_entry_exit_at_T(T.values, hat.l, Tc, Data = snr,
                                L = L)$stage2
      sel <- select_estimate_candidate(s2, L, rule = "narrowest")
      if (!is.na(sel$T_entry)) { T.est <- Tc; s2.est <- s2; break }
    }
    if (is.na(T.est)) next
    
    for (rl in c("narrowest", "first_k")) {
      sel <- select_estimate_candidate(s2.est, L, rule = rl)
      tab.rows[[length(tab.rows) + 1L]] <- data.frame(
        rule = rl, tau = tau, SNR = snr, Tstar = T.est,
        k = sel$k, hat_l = sel$hat_l_current,
        entry = sel$T_entry, exit = sel$T_exit,
        width = sel$T_exit - sel$T_entry)
    }
  }
}
tab.est <- do.call(rbind, tab.rows)

# 1. Check: the narrowest rows must match the old table
print(tab.est[tab.est$rule == "narrowest", -1], row.names = FALSE)

# 2. New rows (first_k)
new <- tab.est[tab.est$rule == "first_k", -1]
print(new, row.names = FALSE)

# 3. LaTeX rows
new <- new[order(new$tau, -new$SNR), ]
for (i in seq_len(nrow(new)))
  with(new[i, ], cat(sprintf(
    "$%.2f$ & %g & %g & %g & %g & %g & %g & %g \\\\\n",
    tau, SNR, Tstar, k, hat_l, entry, exit, width)))


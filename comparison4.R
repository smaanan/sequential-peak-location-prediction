for ( tau_val in c(0.5, 0.75, 0.95) ) {
  results <- do.call(rbind, lapply(Ts_free, function(T) {
    res <- simulate_quantile_objective(T, A = A_true, L = L_true, S = S_true, h = h,
                                       tau = tau_val, snr = Inf, force_l_gt_T = ( tau_val == 0.5 )) 
    idx_min <- which.min(res$loss)
    data.frame(T=T, tau=tau_val, l_star=res$l[idx_min], objective = res$loss[idx_min])
  }))
  
  cat(sprintf("Results for tau = %0.2f:\n", tau_val))
  print(results)
}

############################################################
# MONTE CARLO COVERAGE STUDY
#
# Requires, already defined in the session:
#   stage1_exit_candidates(), stage2_all_candidates(),
#   select_estimate_candidate(), the tolerances,
#   and estimate_peak_curves()
#   (and estimate_peak_curves_h_equals_T() if used)
############################################################

############################################################
# 0. SETTINGS
############################################################

M.rep      <- 50          # use 20-50 first to test
T.min.eval <- 4
estimator  <- "fixed_h"     # or "nested" (h = T version)
n.cores    <- 1             # set > 1 for parallel
seed.base  <- 280

A <- 1; L <- 22; S <- 10; h <- 30
tau.values <- c(0.25, 0.50, 0.75, 0.95)
SNR.values <- c(1000, 100, 1)

############################################################
# 1. TRAJECTORY WRAPPER
############################################################

get_trajectory <- function(tau, SNR, seed) {
  if (estimator == "fixed_h") {
    res <- estimate_peak_curves(
      A = A, L = L, S = S, h = h,
      SNR = SNR, tau = tau, seed = seed
    )
  } else {
    res <- estimate_peak_curves_h_equals_T(
      A = A, L = L, S = S, h.max = h,
      SNR = SNR, tau = tau, seed = seed,
      T.min = T.min.eval
    )
  }
  ord <- order(res$T)
  list(T = res$T[ord], hat.l = res$hat.l[ord])
}

############################################################
# 2. ONE REPLICATE
############################################################

one_replicate <- function(m, tau, SNR) {
  
  traj <- try(get_trajectory(tau, SNR, seed.base + m), silent = TRUE)
  if (inherits(traj, "try-error")) return(NULL)
  
  T.values <- traj$T
  hat.l    <- traj$hat.l
  Ts <- T.values[T.values >= T.min.eval]
  
  score <- function(sel, variant, Tc) {
    ok <- !is.na(sel$T_entry) && !is.na(sel$T_exit)
    data.frame(
      rep = m, tau = tau, SNR = SNR, T = Tc, variant = variant,
      selected = ok,
      k = if (ok) sel$k else NA_real_,
      T_entry = if (ok) sel$T_entry else NA_real_,
      T_exit  = if (ok) sel$T_exit  else NA_real_,
      hat_l   = sel$hat_l_current,
      covered_strict = ok && sel$T_entry <  L && L <  sel$T_exit,
      covered_incl   = ok && sel$T_entry <= L && L <= sel$T_exit,
      stringsAsFactors = FALSE
    )
  }
  
  out <- vector("list", 2 * length(Ts))
  j <- 1L
  
  for (Tc in Ts) {
    
    S1 <- stage1_exit_candidates(
      T.values, hat.l, Tc, STAGE1_TOLERANCE
    )
    S2 <- stage2_all_candidates(
      T.values, hat.l, S1, SNR, L, STAGE2_TOLERANCE
    )
    
    sel.all <- select_estimate_candidate(S2, L)
    S2.k1   <- S2[S2$k >= 1, , drop = FALSE]
    sel.k1  <- select_estimate_candidate(S2.k1, L)
    # hat.l of the empty case when S2.k1 has no rows:
    if (nrow(S2.k1) == 0) sel.k1$hat_l_current <- sel.all$hat_l_current
    
    out[[j]]     <- score(sel.all, "all_k",  Tc)
    out[[j + 1]] <- score(sel.k1,  "k_ge_1", Tc)
    j <- j + 2L
  }
  
  do.call(rbind, out)
}

############################################################
# 3. RUN ALL SETTINGS (with checkpointing)
############################################################

run_setting <- function(tau, SNR) {
  file <- sprintf("mc_tau%.2f_SNR%g_M%d_%s.rds",
                  tau, SNR, M.rep, estimator)
  if (file.exists(file)) return(readRDS(file))
  
  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, ls(envir = .GlobalEnv),
                            envir = .GlobalEnv)
    reps <- parallel::parLapply(
      cl, seq_len(M.rep), one_replicate, tau = tau, SNR = SNR
    )
  } else {
    reps <- lapply(seq_len(M.rep), one_replicate,
                   tau = tau, SNR = SNR)
  }
  
  res <- do.call(rbind, reps)
  saveRDS(res, file)
  res
}

all.runs <- list()
for (tau in tau.values) {
  for (SNR in SNR.values) {
    cat("tau =", tau, " SNR =", SNR, "\n")
    all.runs[[length(all.runs) + 1L]] <- run_setting(tau, SNR)
  }
}
mc.raw <- do.call(rbind, all.runs)
rownames(mc.raw) <- NULL

############################################################
# 4. SUMMARY: q(T), p_sel(T), p_all(T)
############################################################

agg <- aggregate(
  cbind(selected, covered_strict, covered_incl) ~
    tau + SNR + variant + T,
  data = transform(mc.raw,
                   selected = as.numeric(selected),
                   covered_strict = as.numeric(covered_strict),
                   covered_incl = as.numeric(covered_incl)),
  FUN = sum
)

nrep <- tapply(mc.raw$rep, list(mc.raw$tau, mc.raw$SNR),
               function(x) length(unique(x)))
agg$M <- mapply(function(t, s) nrep[as.character(t), as.character(s)],
                agg$tau, agg$SNR)

se_p <- function(p, n) ifelse(n > 0, sqrt(p * (1 - p) / n), NA_real_)

agg$q            <- agg$selected / agg$M
agg$se_q         <- se_p(agg$q, agg$M)
agg$p_sel_strict <- ifelse(agg$selected > 0,
                           agg$covered_strict / agg$selected, NA_real_)
agg$se_sel_strict <- se_p(agg$p_sel_strict, agg$selected)
agg$p_sel_incl   <- ifelse(agg$selected > 0,
                           agg$covered_incl / agg$selected, NA_real_)
agg$p_all_strict <- agg$covered_strict / agg$M   # no interval = miss

agg <- agg[order(agg$variant, agg$tau, -agg$SNR, agg$T), ]
rownames(agg) <- NULL

write.csv(agg, sprintf("mc_coverage_summary_M%d_%s.csv",
                       M.rep, estimator), row.names = FALSE)
print(head(agg, 30))

############################################################
# 5. FIGURES
############################################################

snr.levels <- c(1000, 100, 1)
snr.cols   <- c("blue", "darkgreen", "red")
snr.ltys   <- c(1, 2, 4)
snr.pchs   <- c(16, 17, 15)
snr.labels <- c("SNR = 1000", "SNR = 100", "SNR = 1")

plot_mc_base <- function(variant.name, yvar, se.var = NULL,
                         ylab, file, show.band = TRUE,
                         T.range = NULL) {
  
  d <- agg[agg$variant == variant.name, ]
  if (!is.null(T.range)) {
    d <- d[d$T >= T.range[1] & d$T <= T.range[2], ]
  }
  taus <- sort(unique(d$tau))
  x.lim <- range(d$T)
  
  pdf(file, width = 10, height = 8)
  old.par <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1),
                 oma = c(0, 0, 2, 0))
  
  for (tau in taus) {
    
    plot(NA, xlim = x.lim, ylim = c(0, 1),
         xlab = expression(T), ylab = ylab,
         main = bquote(tau == .(tau)))
    abline(v = L, lty = 3)
    abline(h = c(0.5, 0.9, 0.95), col = "grey85", lty = 1)
    
    for (i in seq_along(snr.levels)) {
      
      dd <- d[d$tau == tau & d$SNR == snr.levels[i], ]
      dd <- dd[order(dd$T), ]
      if (nrow(dd) == 0) next
      
      y <- dd[[yvar]]
      
      if (show.band && !is.null(se.var)) {
        se <- dd[[se.var]]
        ok <- !is.na(y) & !is.na(se)
        if (any(ok)) {
          lo <- pmax(0, y - se)
          hi <- pmin(1, y + se)
          polygon(c(dd$T[ok], rev(dd$T[ok])),
                  c(lo[ok], rev(hi[ok])),
                  col = adjustcolor(snr.cols[i], alpha.f = 0.15),
                  border = NA)
        }
      }
      
      lines(dd$T, y, type = "o", col = snr.cols[i],
            lty = snr.ltys[i], pch = snr.pchs[i], cex = 0.7)
    }
    
    legend("bottomright", legend = snr.labels, col = snr.cols,
           lty = snr.ltys, pch = snr.pchs, bty = "n", cex = 0.8)
  }
  
  mtext(paste0("Selection variant: ", variant.name,
               "; M = ", M.rep, "; L = ", L),
        outer = TRUE, cex = 0.9)
  
  par(old.par)
  dev.off()
}

for (v in c("all_k", "k_ge_1")) {
  
  # Coverage among selected intervals
  plot_mc_base(
    v, "p_sel_strict", "se_sel_strict",
    ylab = expression(hat(p)[sel](T)),
    file = paste0("fig_mc_coverage_", v, ".pdf")
  )
  
  # Interval-selection rate
  plot_mc_base(
    v, "q", "se_q",
    ylab = expression(hat(q)(T)),
    file = paste0("fig_mc_selection_rate_", v, ".pdf")
  )
}

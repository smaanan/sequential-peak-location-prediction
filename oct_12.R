############################################################
# MONTE CARLO COVERAGE STUDY
# Rules:
#   narrowest     : argmin_{k in K(T)} w_k(T)
#   first_k       : min K(T)
#   k1_narrowest  : narrowest among k >= 1 in K(T);
#                   k = 0 if no k >= 1 is admissible
#   k1_first      : smallest k >= 1 in K(T);
#                   k = 0 if no k >= 1 is admissible
# Variants:
#   all_k  : all candidates k >= 0
#   k_ge_1 : only candidates k >= 1 (no fallback)
#
# Requires, already defined in the session:
#   stage1_exit_candidates(), stage2_all_candidates(),
#   the tolerances, and estimate_peak_curves()
#   (and estimate_peak_curves_h_equals_T() if used)
#
# All rules are applied to the SAME candidate tables.
# The true L is used only to score coverage.
############################################################

############################################################
# 0. SETTINGS
############################################################

M.rep      <- 1000          # test with 50; then use 1000
T.min.eval <- 4
estimator  <- "fixed_h"     # or "nested" (h = T version)
n.cores    <- 8             # set > 1 for parallel
seed.base  <- 280

rules          <- c("narrowest", "first_k",
                    "k1_narrowest", "k1_first")
fallback.rules <- c("k1_narrowest", "k1_first")  # all_k only

A <- 1; L <- 22; S <- 10; h <- 30
tau.values <- c(0.25, 0.50, 0.75, 0.95)
SNR.values <- c(1000, 100, 1)

############################################################
# 0b. SELECTION FUNCTIONS
############################################################

select_narrowest_candidate <- function(
    stage2.table,
    containment.column,
    L,
    rule = c("narrowest", "first_k",
             "k1_narrowest", "k1_first")
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
  
  # Fallback rules: restrict to k >= 1 when possible,
  # otherwise keep the admissible set (which is then {0}).
  if (rule %in% c("k1_narrowest", "k1_first")) {
    cand.k1 <- candidates[candidates$k >= 1, , drop = FALSE]
    if (nrow(cand.k1) > 0) candidates <- cand.k1
  }
  
  ord <- if (rule %in% c("narrowest", "k1_narrowest")) {
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
  
  score <- function(sel, variant, rule, Tc) {
    ok <- !is.na(sel$T_entry) && !is.na(sel$T_exit)
    data.frame(
      rep = m, tau = tau, SNR = SNR, T = Tc,
      variant = variant, rule = rule,
      selected = ok,
      k = if (ok) sel$k else NA_real_,
      T_entry = if (ok) sel$T_entry else NA_real_,
      T_exit  = if (ok) sel$T_exit  else NA_real_,
      width   = if (ok) sel$T_exit - sel$T_entry else NA_real_,
      hat_l   = sel$hat_l_current,
      covered_strict = ok && sel$T_entry <  L && L <  sel$T_exit,
      covered_incl   = ok && sel$T_entry <= L && L <= sel$T_exit,
      stringsAsFactors = FALSE
    )
  }
  
  out <- vector("list", 2L * length(rules) * length(Ts))
  j <- 1L
  
  for (Tc in Ts) {
    
    S1 <- stage1_exit_candidates(
      T.values, hat.l, Tc, STAGE1_TOLERANCE
    )
    S2 <- stage2_all_candidates(
      T.values, hat.l, S1, SNR, L, STAGE2_TOLERANCE
    )
    S2.k1 <- S2[S2$k >= 1, , drop = FALSE]
    
    for (rl in rules) {
      
      sel.all <- select_estimate_candidate(S2, L, rule = rl)
      out[[j]] <- score(sel.all, "all_k", rl, Tc)
      j <- j + 1L
      
      # The k >= 1 variant is meaningful only for the
      # non-fallback rules.
      if (!(rl %in% fallback.rules)) {
        sel.k1 <- select_estimate_candidate(S2.k1, L, rule = rl)
        if (nrow(S2.k1) == 0) {
          sel.k1$hat_l_current <- sel.all$hat_l_current
        }
        out[[j]] <- score(sel.k1, "k_ge_1", rl, Tc)
        j <- j + 1L
      }
    }
  }
  
  do.call(rbind, out[seq_len(j - 1L)])
}

############################################################
# 3. RUN ALL SETTINGS (with checkpointing)
############################################################

run_setting <- function(tau, SNR) {
  file <- sprintf("mc_rules4_tau%.2f_SNR%g_M%d_%s.rds",
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

write.csv(mc.raw,
          sprintf("mc_raw_rules4_M%d_%s.csv", M.rep, estimator),
          row.names = FALSE)

############################################################
# 4. SUMMARY: q(T), p_sel(T), p_all(T), width, k
############################################################

raw.num <- transform(mc.raw,
                     selected = as.numeric(selected),
                     covered_strict = as.numeric(covered_strict),
                     covered_incl = as.numeric(covered_incl))

agg <- aggregate(
  cbind(selected, covered_strict, covered_incl) ~
    tau + SNR + variant + rule + T,
  data = raw.num,
  FUN = sum
)

nrep <- tapply(mc.raw$rep, list(mc.raw$tau, mc.raw$SNR),
               function(x) length(unique(x)))
agg$M <- mapply(function(t, s) nrep[as.character(t), as.character(s)],
                agg$tau, agg$SNR)

se_p <- function(p, n) ifelse(n > 0, sqrt(p * (1 - p) / n), NA_real_)

agg$q             <- agg$selected / agg$M
agg$se_q          <- se_p(agg$q, agg$M)
agg$p_sel_strict  <- ifelse(agg$selected > 0,
                            agg$covered_strict / agg$selected, NA_real_)
agg$se_sel_strict <- se_p(agg$p_sel_strict, agg$selected)
agg$p_sel_incl    <- ifelse(agg$selected > 0,
                            agg$covered_incl / agg$selected, NA_real_)
agg$p_all_strict  <- agg$covered_strict / agg$M   # no interval = miss

# Mean width, mean k, proportion of k = 0 (selected only)
sel.rows <- mc.raw[mc.raw$selected, ]
if (nrow(sel.rows) > 0) {
  sel.rows$k0 <- as.numeric(sel.rows$k == 0)
  extra <- aggregate(
    cbind(width, k, k0) ~ tau + SNR + variant + rule + T,
    data = sel.rows, FUN = mean
  )
  names(extra)[names(extra) == "width"] <- "mean_width"
  names(extra)[names(extra) == "k"]     <- "mean_k"
  names(extra)[names(extra) == "k0"]    <- "prop_k0"
  agg <- merge(agg, extra,
               by = c("tau", "SNR", "variant", "rule", "T"),
               all.x = TRUE)
}

agg <- agg[order(agg$variant, agg$rule, agg$tau, -agg$SNR, agg$T), ]
rownames(agg) <- NULL

write.csv(agg,
          sprintf("mc_coverage_summary_rules4_M%d_%s.csv",
                  M.rep, estimator),
          row.names = FALSE)
print(head(agg, 30))

############################################################
# 4b. COMPARISON OF THE RULES (variant all_k)
############################################################

a.all <- agg[agg$variant == "all_k",
             c("tau", "SNR", "T", "rule", "q",
               "p_sel_strict", "mean_width", "prop_k0")]

cmp <- reshape(a.all,
               idvar = c("tau", "SNR", "T"),
               timevar = "rule", direction = "wide")

# CHECK: q(T) depends only on whether K(T) is empty, so it
# must be identical for all rules.
q.cols <- grep("^q\\.", names(cmp), value = TRUE)
q.ref  <- cmp[[q.cols[1]]]
same.q <- all(sapply(q.cols, function(cn)
  isTRUE(all.equal(cmp[[cn]], q.ref))))
if (!same.q) {
  warning("Selection rates differ between rules: check the code.")
} else {
  cat("\nCheck passed: selection rates are identical across rules.\n")
}

cmp <- cmp[order(cmp$tau, -cmp$SNR, cmp$T), ]
rownames(cmp) <- NULL

write.csv(cmp,
          sprintf("mc_comparison_rules4_M%d_%s.csv",
                  M.rep, estimator),
          row.names = FALSE)

cat("\nCOMPARISON AT T = 22 AND T = 30 (variant all_k)\n\n")
show.cols <- c("tau", "SNR", "T",
               paste0("p_sel_strict.", rules),
               paste0("mean_width.", rules))
print(cmp[cmp$T %in% c(22, 30), show.cols], digits = 3)

cat("\nProportion of k = 0 among selected intervals (all_k)\n\n")
print(cmp[cmp$T %in% c(22, 30),
          c("tau", "SNR", "T", paste0("prop_k0.", rules))],
      digits = 3)

print(
  knitr::kable(
    cmp[cmp$T == 30, c("tau", "SNR", "T",
                       paste0("p_sel_strict.", rules))],
    format = "latex", booktabs = TRUE, digits = 3
  )
)

############################################################
# 4c. T_est BY REPLICATE (the same for all rules)
############################################################

d.est <- mc.raw[mc.raw$rule == "narrowest" &
                  mc.raw$variant == "all_k" &
                  mc.raw$selected, ]
est <- aggregate(T ~ tau + SNR + rep, d.est, min)
names(est)[names(est) == "T"] <- "T_est"
M.all <- length(unique(mc.raw$rep))

T.est.summ <- do.call(rbind, lapply(
  split(est, list(est$tau, est$SNR), drop = TRUE),
  function(df) data.frame(
    tau = df$tau[1], SNR = df$SNR[1],
    never_selected  = 1 - nrow(df) / M.all,
    prop_T_est_lt_L = mean(df$T_est < L),
    median_T_est    = median(df$T_est),
    q10 = quantile(df$T_est, 0.1),
    q90 = quantile(df$T_est, 0.9)
  )))
rownames(T.est.summ) <- NULL
cat("\nDistribution of T_est (all_k, any rule):\n")
print(T.est.summ, digits = 3)

############################################################
# 5. FIGURES
############################################################

snr.levels <- c(1000, 100, 1)
snr.cols   <- c("blue", "darkgreen", "red")
snr.ltys   <- c(1, 2, 4)
snr.pchs   <- c(16, 17, 15)
snr.labels <- c("SNR = 1000", "SNR = 100", "SNR = 1")

plot_mc_base <- function(variant.name, rule.name, yvar,
                         se.var = NULL, ylab, file,
                         show.band = TRUE, T.range = NULL) {
  
  d <- agg[agg$variant == variant.name & agg$rule == rule.name, ]
  if (!is.null(T.range)) {
    d <- d[d$T >= T.range[1] & d$T <= T.range[2], ]
  }
  taus  <- sort(unique(d$tau))
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
           lty = snr.ltys, pch = snr.pchs, bty = "o", cex = 0.8)
  }
  
  # mtext(paste0("Variant: ", variant.name, "; rule: ", rule.name, "; M = ", M.rep, "; L = ", L), outer = TRUE, cex = 0.9)
  
  par(old.par)
  dev.off()
}

# Coverage figures: every (variant, rule) combination present
combos <- unique(agg[, c("variant", "rule")])

for (i in seq_len(nrow(combos))) {
  v  <- combos$variant[i]
  rl <- combos$rule[i]
  plot_mc_base(
    v, rl, "p_sel_strict", "se_sel_strict",
    ylab = expression(hat(p)[sel](T)),
    file = paste0("fig_mc_coverage_", v, "_", rl, ".pdf")
  )
}

# Selection-rate figure: identical for all rules, so one per
# variant, drawn with the first available rule.
for (v in unique(combos$variant)) {
  rl <- combos$rule[combos$variant == v][1]
  plot_mc_base(
    v, rl, "q", "se_q",
    ylab = expression(hat(q)(T)),
    file = paste0("fig_mc_selection_rate_", v, ".pdf")
  )
}

############################################################
# 6. DIRECT RULE COMPARISON (variant all_k), ONE FIGURE PER SNR
############################################################

rule.cols <- c(narrowest = "blue", first_k = "red",
               k1_narrowest = "darkgreen", k1_first = "orange")
rule.ltys <- c(narrowest = 1, first_k = 2,
               k1_narrowest = 4, k1_first = 3)
rule.pchs <- c(narrowest = 16, first_k = 17,
               k1_narrowest = 15, k1_first = 18)

for (snr in SNR.values) {
  
  dd <- a.all[a.all$SNR == snr, ]
  
  pdf(sprintf("fig_mc_rule_comparison_SNR%g.pdf", snr),
      width = 10, height = 8)
  old.par <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1),
                 oma = c(0, 0, 2, 0))
  
  for (tau in tau.values) {
    
    plot(NA, xlim = range(dd$T), ylim = c(0, 1),
         xlab = expression(T),
         ylab = expression(hat(p)[sel](T)),
         main = bquote(tau == .(tau)))
    abline(v = L, lty = 3)
    
    for (rl in rules) {
      x <- dd[dd$tau == tau & dd$rule == rl, ]
      x <- x[order(x$T), ]
      if (nrow(x) == 0) next
      lines(x$T, x$p_sel_strict, type = "o",
            col = rule.cols[rl], lty = rule.ltys[rl],
            pch = rule.pchs[rl], cex = 0.7)
    }
    
    legend("bottomright", legend = rules,
           col = rule.cols[rules], lty = rule.ltys[rules],
           pch = rule.pchs[rules], bty = "n", cex = 0.75)
  }
  
  mtext(sprintf("Coverage among selected, variant all_k, SNR = %g; M = %d",
                snr, M.rep), outer = TRUE, cex = 0.9)
  par(old.par)
  dev.off()
}

tstar <- function(tr, tol = 0) {
  ok <- abs(tr$hat.l - L) <= tol
  good <- sapply(seq_along(ok), function(i) all(ok[i:length(ok)]))
  if (any(good)) tr$T[min(which(good))] else Inf
}

res <- do.call(rbind, lapply(SNR.values, function(snr)
  do.call(rbind, lapply(seq_len(200), function(m) {
    tr <- get_trajectory(0.50, snr, seed.base + m)
    data.frame(SNR = snr, rep = m, Tstar = tstar(tr))
  }))))

aggregate(Tstar ~ SNR, res, function(x)
  c(median_finite = median(x[is.finite(x)]),
    prop_finite = mean(is.finite(x)),
    prop_gt_L = mean(x > L)))


res1 <- res[res$SNR == 1, ]
fin  <- res1[is.finite(res1$Tstar), ]; fin <- fin[order(fin$Tstar), ]
inf  <- res1[!is.finite(res1$Tstar), ]
head(fin, 10)   # finite T*: take the smallest one
head(inf, 10)   # infinite T*: any, ideally with large oscillation


tau.fig <- 0.50
rep.fin <- fin$rep[1]          # finite-T* example at SNR = 1
rep.inf <- inf$rep[1]          # infinite-T* example at SNR = 1
reps.hi <- 1:10

panel <- function(snr, reps, main, tol = 0) {
  plot(NA, xlim = c(4, h), ylim = c(1, h), xlab = "T",
       ylab = expression(hat(l)(T)), main = main)
  abline(h = L, lty = 2, lwd = 2, col = "red")
  for (m in reps) {
    tr <- get_trajectory(tau.fig, snr, seed.base + m)
    lines(tr$T, tr$hat.l, col = "steelblue", lwd = 0.9)
    ts <- tstar(tr, tol)
    if (is.finite(ts))
      points(ts, tr$hat.l[tr$T == ts], pch = 19, col = "black", cex = 0.9)
  }
}

pdf("fig_tstar_examples.pdf", width = 10, height = 8)
op <- par(mfrow = c(2, 2), mar = c(4.2, 5, 3, 1))
panel(1000, reps.hi, "SNR = 1000")
panel(100,  reps.hi, "SNR = 100")
panel(1, rep.fin,    expression("SNR = 1: finite " * T^"*"))
panel(1, rep.inf,    expression("SNR = 1: " * T^"*" == infinity))
par(op); dev.off()





tau.fig <- 0.50
osc <- sapply(inf$rep, function(m) {
  tr <- get_trajectory(tau.fig, 1, seed.base + m)
  diff(range(tr$hat.l[tr$T >= 20]))
})
rep.inf <- inf$rep[which.max(osc)]
rep.fin <- 138
c(rep.inf = rep.inf, osc = max(osc))

# check the finite example
tr138 <- get_trajectory(tau.fig, 1, seed.base + rep.fin)
print(data.frame(T = tr138$T, hat_l = tr138$hat.l)[tr138$T >= 18, ])
tstar(tr138)   # should return 24























M.str  <- 500
Ts.str <- c(22, 26, 30)
snr.s  <- 1

pairs_one <- function(tau, m) {
  tr <- get_trajectory(tau, snr.s, seed.base + m)
  out <- list()
  for (Tc in Ts.str) {
    S1 <- stage1_exit_candidates(tr$T, tr$hat.l, Tc, STAGE1_TOLERANCE)
    S2 <- stage2_all_candidates(tr$T, tr$hat.l, S1, snr.s, L,
                                STAGE2_TOLERANCE)
    S2 <- S2[!is.na(S2$T_entry), ]
    S2 <- S2[order(S2$k), ]
    n <- nrow(S2)
    if (n < 2) next
    a <- S2[-n, ]; b <- S2[-1, ]
    ok <- b$k == a$k + 1
    if (!any(ok)) next
    a <- a[ok, ]; b <- b[ok, ]
    out[[length(out) + 1L]] <- data.frame(
      tau = tau, rep = m, T = Tc, k = a$k,
      entry_k = a$T_entry, exit_k = a$T_exit,
      entry_k1 = b$T_entry, exit_k1 = b$T_exit)
  }
  if (length(out)) do.call(rbind, out) else NULL
}

pairs <- do.call(rbind, lapply(tau.values, function(tau)
  do.call(rbind, lapply(seq_len(M.str), function(m) pairs_one(tau, m)))))

pairs$exit_ok   <- pairs$exit_k1 < pairs$exit_k
pairs$entry_rev <- pairs$entry_k1 > pairs$entry_k
pairs$identity  <- pairs$entry_k == pairs$exit_k1
pairs$nested    <- pairs$entry_k1 >= pairs$entry_k &
  pairs$exit_k1  <= pairs$exit_k
pairs$disjoint  <- pairs$exit_k1 < pairs$entry_k

structure.summ <- do.call(rbind, lapply(
  split(pairs, list(pairs$tau, pairs$T), drop = TRUE),
  function(d) data.frame(
    tau = d$tau[1], T = d$T[1], n_pairs = nrow(d),
    exit_order  = mean(d$exit_ok),
    entry_rev   = mean(d$entry_rev),
    identity    = mean(d$identity),
    nested      = mean(d$nested),
    disjoint    = mean(d$disjoint))))
rownames(structure.summ) <- NULL
print(structure.summ, digits = 3)

# Number of candidates per (tau, T), averaged over replicates
ncand <- aggregate(k ~ tau + T + rep, pairs, function(x) length(x) + 1)
print(aggregate(k ~ tau + T, ncand, mean), digits = 3)





















############################################################
# EXIT / ENTRY TRAJECTORIES FOR A GIVEN SNR
# Required object: combined.diagnostics
# Required columns: tau_level, SNR_level, T_current, k,
#                   T_exit, T_entry
############################################################

snr.fig  <- 1                 # <- change here (1000, 100, 1)
k.show   <- 0:3               # candidates to display
out.file <- sprintf("fig_exit_entry_times_SNR%g.pdf", snr.fig)

req <- c("tau_level", "SNR_level", "T_current", "k",
         "T_exit", "T_entry")
miss <- setdiff(req, names(combined.diagnostics))
if (length(miss) > 0)
  stop("Missing columns: ", paste(miss, collapse = ", "))

diag.snr <- combined.diagnostics[
  combined.diagnostics$SNR_level == snr.fig &
    !is.na(combined.diagnostics$T_exit) &
    !is.na(combined.diagnostics$T_entry) &
    combined.diagnostics$k %in% k.show,
  , drop = FALSE]

if (nrow(diag.snr) == 0)
  stop("No valid observations for SNR = ", snr.fig)

diag.snr <- diag.snr[order(diag.snr$tau_level,
                           diag.snr$T_current, diag.snr$k), ]

tau.values <- sort(unique(diag.snr$tau_level))

L.plot <- if ("L" %in% names(combined.diagnostics) &&
              length(unique(na.omit(combined.diagnostics$L))) == 1)
  unique(na.omit(combined.diagnostics$L)) else 22

pch.symbols <- c(16, 17, 15, 18, 19, 8, 3, 4, 6, 7)

pdf(out.file, width = 8.5, height = 6.5)
old.par <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1))

for (tau in tau.values) {
  
  tmp <- diag.snr[diag.snr$tau_level == tau, , drop = FALSE]
  k.values <- sort(unique(tmp$k))
  
  plot(NA,
       xlim = range(tmp$T_current),
       ylim = range(c(tmp$T_exit, tmp$T_entry, L.plot), finite = TRUE),
       xlab = expression(T), ylab = "Time",
       main = bquote(tau == .(tau) * ", SNR = " * .(snr.fig)))
  
  abline(h = L.plot, lty = 2)
  
  for (i in seq_along(k.values)) {
    k  <- k.values[i]
    pc <- pch.symbols[(i - 1) %% length(pch.symbols) + 1]
    tk <- tmp[tmp$k == k, , drop = FALSE]
    lines(tk$T_current, tk$T_exit,  type = "o", pch = pc,
          lty = 1, col = "blue")
    lines(tk$T_current, tk$T_entry, type = "o", pch = pc,
          lty = 2, col = "red")
  }
  
  legend("topleft",
         legend = c("Exit", "Entry", expression(L),
                    paste0("k = ", k.values)),
         lty = c(1, 2, 2, rep(1, length(k.values))),
         pch = c(NA, NA, NA, pch.symbols[seq_along(k.values)]),
         col = c("blue", "red", "black",
                 rep("black", length(k.values))),
         bty = "o", bg = "white", cex = 0.8)
}

par(old.par)
dev.off()
cat("Saved:", out.file, "\n")
tstar_new <- function(tr, m = 3) {
  l <- tr$hat.l; n <- length(l)
  ok <- sapply(seq_len(n), function(i) all(l[i:n] == l[i]) && (n - i + 1) >= m)
  if (any(ok)) tr$T[min(which(ok))] else Inf
}

tau.fig <- 0.50
res2 <- do.call(rbind, lapply(c(1000, 100, 1), function(snr)
  do.call(rbind, lapply(seq_len(200), function(mm) {
    tr <- get_trajectory(tau.fig, snr, seed.base + mm)
    ts <- tstar_new(tr, m = 3)
    data.frame(SNR = snr, rep = mm, Tstar = ts,
               at_L = is.finite(ts) && tr$hat.l[tr$T == ts] == L)
  }))))

aggregate(cbind(finite = is.finite(Tstar), beyond_L = Tstar > L, at_L = at_L) ~ SNR, res2, mean)
aggregate(Tstar ~ SNR, res2[is.finite(res2$Tstar), ], median)

res1 <- res2[res2$SNR == 1, ]

# finite example: smallest finite T* among plateaus equal to L
fin <- res1[is.finite(res1$Tstar) & res1$at_L, ]
fin <- fin[order(fin$Tstar), ]
head(fin, 5)
rep.fin <- fin$rep[1]

# infinite example: typical oscillation over T >= 20
inf <- res1[!is.finite(res1$Tstar), ]
osc <- sapply(inf$rep, function(mm) {
  tr <- get_trajectory(tau.fig, 1, seed.base + mm)
  diff(range(tr$hat.l[tr$T >= 20]))
})
rep.inf <- inf$rep[which.min(abs(osc - median(osc)))]
c(rep.fin = rep.fin, rep.inf = rep.inf, median_osc = median(osc))

# check the finite example (if it is 138, T* should be 24)
trf <- get_trajectory(tau.fig, 1, seed.base + rep.fin)
print(data.frame(T = trf$T, hat_l = trf$hat.l)[trf$T >= 18, ])
tstar_new(trf)

reps.hi <- 1:10
panel <- function(snr, reps, main) {
  plot(NA, xlim = c(4, h), ylim = c(1, h), xlab = "T",
       ylab = expression(hat(l)(T)), main = main)
  abline(h = L, lty = 2, lwd = 2, col = "red")
  for (mm in reps) {
    tr <- get_trajectory(tau.fig, snr, seed.base + mm)
    lines(tr$T, tr$hat.l, col = "steelblue", lwd = 0.9)
    ts <- tstar_new(tr, m = 3)
    if (is.finite(ts))
      points(ts, tr$hat.l[tr$T == ts], pch = 19, col = "black", cex = 0.9)
  }
}

pdf("fig_tstar_examples.pdf", width = 10, height = 8)
op <- par(mfrow = c(2, 2), mar = c(4.2, 5, 3, 1))
panel(1000, reps.hi, "SNR = 1000")
panel(100,  reps.hi, "SNR = 100")
panel(1, rep.fin, expression("SNR = 1: finite " * T^"*"))
panel(1, rep.inf, expression("SNR = 1: " * T^"*" == infinity))
par(op); dev.off()

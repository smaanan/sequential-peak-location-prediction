M.chk <- 200

chk <- do.call(rbind, lapply(tau.values, function(tau) {
  do.call(rbind, lapply(SNR.values, function(snr) {
    do.call(rbind, lapply(seq_len(M.chk), function(m) {
      tr <- get_trajectory(tau, snr, seed.base + m)
      Tc <- max(tr$T)
      S1 <- stage1_exit_candidates(tr$T, tr$hat.l, Tc, STAGE1_TOLERANCE)
      S2 <- stage2_all_candidates(tr$T, tr$hat.l, S1, snr, L, STAGE2_TOLERANCE)
      e0 <- S2$T_entry[S2$k == 0]
      late <- tr$T >= 20
      data.frame(
        tau = tau, SNR = snr, rep = m,
        entry0 = e0,
        hat_l_T = tr$hat.l[tr$T == Tc],
        osc_last10 = diff(range(tr$hat.l[late])),
        err_after_entry = max(abs(tr$hat.l[tr$T >= e0] - L))
      )
    }))   # end inner lapply seq_len(M.chk)
  }))     # end lapply SNR.values
}))       # end lapply tau.values

summ <- aggregate(cbind(entry0, osc_last10, err_after_entry) ~ tau + SNR,
                  chk, median)
summ$prop_exact <- aggregate(I(hat_l_T == L) ~ tau + SNR, chk, mean)[[3]]
summ$prop_stable <- aggregate(I(osc_last10 == 0) ~ tau + SNR, chk, mean)[[3]]
print(summ, digits = 3)

# Step 2: Plot hat{l}(T) trajectories for three SNR values
# Base R only, save as PDF

tau.plot <- 0.50
SNR.plot <- c(1000, 100, 1)
M.plot <- 20
seed.offsets <- seq_len(M.plot)

# Generate trajectories
traj.list <- list()
idx <- 1L

for (snr in SNR.plot) {
  for (m in seed.offsets) {
    tr <- get_trajectory(
      tau.plot,
      snr,
      seed.base + m
    )
    
    traj.list[[idx]] <- data.frame(
      SNR = snr,
      rep = m,
      T = tr$T,
      hat_l = tr$hat.l
    )
    
    idx <- idx + 1L
  }
}

traj <- do.call(rbind, traj.list)

# Common plotting range for all panels
xlim <- range(traj$T, finite = TRUE)
ylim <- range(c(traj$hat_l, L), finite = TRUE)

# Save as PDF (inches). Matches roughly your previous PNG size:
pdf(file = "hatl_trajectories_tau050.pdf", width = 10, height = 1100/300)

# Increase left margin so the y-label (hat over l) is not clipped.
# mar = c(bottom, left, top, right)
op <- par(
  mfrow = c(1, 3),
  mar = c(4.2, 5.5, 3.2, 1),  # larger left margin
  oma = c(0, 0, 0, 0)
)

for (snr in SNR.plot) {
  
  d <- traj[traj$SNR == snr, ]
  
  plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    xlab = "T",
    ylab = expression(hat(l)(T)),
    main = paste("SNR =", snr),
    cex.lab = 1.1,
    cex.main = 1.2
  )
  
  # True peak location (red dashed line)
  abline(h = L, lty = 2, lwd = 2, col = "red")
  
  # One trajectory per replicate (blue)
  for (m in seed.offsets) {
    dm <- d[d$rep == m, ]
    dm <- dm[order(dm$T), ]
    
    lines(
      dm$T,
      dm$hat_l,
      lwd = 0.8,
      col = "steelblue"
    )
  }
}

# no mtext (removed as requested)

par(op)
dev.off()
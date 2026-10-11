############################################################
# 30. MONTE CARLO EMPIRICAL COVERAGE OF OPERATIONAL INTERVALS
#
# Uses the existing functions:
#   estimate_peak_curves()
#   run_entry_exit_at_T()
#   select_estimate_candidate()
#
# The true L is used only to assess coverage, never to select
# the operational candidate interval.
############################################################

run_monte_carlo_coverage <- function(
    M = 100,
    seed_start = 10000,
    tau.grid = c(0.25, 0.50, 0.75, 0.95),
    SNR.grid = c(1000, 100, 1),
    A.true = 1,
    L.true = 22,
    S.true = 10,
    h.max = 30
) {
  
  results <- list()
  counter <- 1L
  combo_id <- 0L
  
  for (tau in tau.grid) {
    for (SNR in SNR.grid) {
      
      combo_id <- combo_id + 1L
      
      cat(
        "\nStarting: tau =", tau,
        "| SNR =", SNR, "\n"
      )
      
      for (m in seq_len(M)) {
        
        # Each replicate gets a distinct seed.
        current_seed <- seed_start + combo_id * M + m
        
        # Generate the sequential estimation trajectory.
        res <- estimate_peak_curves(
          A = A.true,
          L = L.true,
          S = S.true,
          h = h.max,
          SNR = SNR,
          tau = tau,
          seed = current_seed
        )
        
        T.values <- res$T
        hat.l <- res$hat.l
        
        if (length(T.values) != length(hat.l)) {
          stop("T and hat.l have different lengths.")
        }
        
        ord <- order(T.values)
        T.values <- T.values[ord]
        hat.l <- hat.l[ord]
        
        if (anyDuplicated(T.values)) {
          stop("Duplicate observation horizons found.")
        }
        
        # Evaluate the operational procedure at every
        # available observation horizon.
        for (T.current in T.values) {
          
          result <- run_entry_exit_at_T(
            T.values = T.values,
            hat.l = hat.l,
            T.current = T.current,
            Data = SNR,
            L = L.true,
            tolerance.stage1 = STAGE1_TOLERANCE,
            tolerance.E = E_TOLERANCE,
            tolerance.stage2 = STAGE2_TOLERANCE
          )
          
          stage2.table <- result$stage2
          
          # IMPORTANT:
          # Selection is based on hat.l(T), not on true L.
          selected <- select_estimate_candidate(
            stage2.table = stage2.table,
            L = L.true
          )
          
          has_interval <-
            nrow(selected) == 1L &&
            is.finite(selected$T_entry[1]) &&
            is.finite(selected$T_exit[1]) &&
            isTRUE(selected$estimate.in.window[1])
          
          if (has_interval) {
            
            entry <- selected$T_entry[1]
            exit <- selected$T_exit[1]
            width <- exit - entry
            
            # Inclusive containment, consistent with
            # the interval [T_entry, T_exit].
            covers_L <-
              entry <= L.true && L.true <= exit
            
            k.selected <- selected$k[1]
            
          } else {
            
            entry <- NA_real_
            exit <- NA_real_
            width <- NA_real_
            covers_L <- FALSE
            k.selected <- NA_real_
          }
          
          results[[counter]] <- data.frame(
            replicate = m,
            seed = current_seed,
            tau = tau,
            SNR = SNR,
            T = T.current,
            hat_l = selected$hat_l_current[1],
            selected = has_interval,
            k = k.selected,
            T_entry = entry,
            T_exit = exit,
            width = width,
            covers_L = covers_L
          )
          
          counter <- counter + 1L
        }
        
        if (m %% 25 == 0L) {
          cat("  Completed", m, "of", M, "replicates\n")
        }
      }
    }
  }
  
  raw <- do.call(rbind, results)
  rownames(raw) <- NULL
  
  # Summarize each (tau, SNR, T) combination.
  groups <- split(
    raw,
    interaction(raw$tau, raw$SNR, raw$T, drop = TRUE)
  )
  
  summary.list <- lapply(groups, function(df) {
    
    n <- nrow(df)
    n.selected <- sum(df$selected)
    n.covered <- sum(df$covers_L)
    
    coverage.selected <- if (n.selected > 0) {
      mean(df$covers_L[df$selected])
    } else {
      NA_real_
    }
    
    mc.se <- if (n.selected > 0) {
      sqrt(
        coverage.selected * (1 - coverage.selected) /
          n.selected
      )
    } else {
      NA_real_
    }
    
    data.frame(
      tau = df$tau[1],
      SNR = df$SNR[1],
      T = df$T[1],
      M = n,
      n_selected = n.selected,
      n_covered = n.covered,
      selection_rate = mean(df$selected),
      coverage_selected = coverage.selected,
      coverage_all = mean(df$covers_L),
      MC_SE_selected = mc.se,
      mean_width = if (n.selected > 0) {
        mean(df$width[df$selected])
      } else {
        NA_real_
      }
    )
  })
  
  summary <- do.call(rbind, summary.list)
  rownames(summary) <- NULL
  
  summary <- summary[
    order(summary$tau, -summary$SNR, summary$T),
  ]
  
  list(
    raw = raw,
    summary = summary
  )
}


############################################################
# 31. RUN THE MONTE CARLO EXPERIMENT
############################################################

# First use M = 100 to test the code and estimate runtime.
# For the main experiment, increase M to 1000.

M.replicates <- 1000

mc.results <- run_monte_carlo_coverage(
  M = M.replicates,
  seed_start = 10000,
  tau.grid = tau.values,
  SNR.grid = SNR.values,
  A.true = A,
  L.true = L,
  S.true = S,
  h.max = h
)

mc.raw <- mc.results$raw
mc.summary <- mc.results$summary

# Save the full replicate-level results and the summary.
write.csv(
  mc.raw,
  "monte-carlo-interval-coverage-raw.csv",
  row.names = FALSE
)

write.csv(
  mc.summary,
  "monte-carlo-interval-coverage-summary.csv",
  row.names = FALSE
)

print(mc.summary)

############################################################
# HELPER FUNCTION: adjustAlpha
# Add this before the plotting code if you don't have it
############################################################

adjustAlpha <- function(col, alpha = 0.5) {
  rgb_vals <- col2rgb(col) / 255
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = alpha)
}




############################################################
# 32. FIGURE 1: COVERAGE AMONG SELECTED INTERVALS
############################################################

pdf(
  "monte-carlo-coverage-selected.pdf",
  width = 8.5,
  height = 6.5
)

old.par <- par(
  mfrow = c(2, 2),
  mar = c(4, 4.5, 3, 1)
)

snr.colors <- c("1000" = "red", "100" = "blue", "1" = "gray")
snr.pch <- c("1000" = 16, "100" = 17, "1" = 18)

for (tau in tau.values) {
  
  df.tau <- mc.summary[mc.summary$tau == tau, ]
  
  plot(
    NA,
    xlim = range(df.tau$T),
    ylim = c(0, 1),
    xlab = "Observation horizon T",
    ylab = "Empirical coverage",
    main = bquote(tau == .(tau))
  )
  
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  
  for (snr in SNR.values) {
    
    d <- df.tau[df.tau$SNR == snr, ]
    d <- d[order(d$T), ]
    
    ci_lower <- d$coverage_selected - 1.96 * d$MC_SE_selected
    ci_upper <- d$coverage_selected + 1.96 * d$MC_SE_selected
    
    ci_lower <- pmax(ci_lower, 0)
    ci_upper <- pmin(ci_upper, 1)
    
    polygon(
      c(d$T, rev(d$T)),
      c(ci_lower, rev(ci_upper)),
      col = adjustAlpha(snr.colors[as.character(snr)], 0.2),
      border = NA
    )
    
    lines(
      d$T,
      d$coverage_selected,
      type = "b",
      pch = snr.pch[as.character(snr)],
      col = snr.colors[as.character(snr)],
      lwd = 1.5
    )
  }
  
  legend(
    "topleft",
    legend = paste("SNR =", SNR.values),
    col = snr.colors[as.character(SNR.values)],
    lty = 1,
    pch = snr.pch[as.character(SNR.values)],
    cex = 0.8,
    bty = "o",
    bg = "white"
  )
}

par(old.par)
dev.off()


############################################################
# 33. FIGURE 2: INTERVAL-SELECTION RATE
############################################################

pdf(
  "monte-carlo-interval-selection-rate.pdf",
  width = 8.5,
  height = 6.5
)

old.par <- par(
  mfrow = c(2, 2),
  mar = c(4, 4.5, 3, 1)
)

snr.colors <- c("1000" = "red", "100" = "blue", "1" = "gray")
snr.pch <- c("1000" = 16, "100" = 17, "1" = 18)

for (tau in tau.values) {
  
  df.tau <- mc.summary[mc.summary$tau == tau, ]
  
  plot(
    NA,
    xlim = range(df.tau$T),
    ylim = c(0, 1),
    xlab = "Observation horizon T",
    ylab = "Interval-selection rate",
    main = bquote(tau == .(tau))
  )
  
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  
  for (snr in SNR.values) {
    
    d <- df.tau[df.tau$SNR == snr, ]
    d <- d[order(d$T), ]
    
    p <- d$selection_rate
    n <- d$M
    
    se <- sqrt(p * (1 - p) / n)
    ci_lower <- p - 1.96 * se
    ci_upper <- p + 1.96 * se
    
    ci_lower <- pmax(ci_lower, 0)
    ci_upper <- pmin(ci_upper, 1)
    
    polygon(
      c(d$T, rev(d$T)),
      c(ci_lower, rev(ci_upper)),
      col = adjustAlpha(snr.colors[as.character(snr)], 0.2),
      border = NA
    )
    
    lines(
      d$T,
      d$selection_rate,
      type = "b",
      pch = snr.pch[as.character(snr)],
      col = snr.colors[as.character(snr)],
      lwd = 1.5
    )
  }
  
  legend(
    "topleft",
    legend = paste("SNR =", SNR.values),
    col = snr.colors[as.character(SNR.values)],
    lty = 1,
    pch = snr.pch[as.character(SNR.values)],
    cex = 0.8,
    bty = "o",
    bg = "white"
  )
}

par(old.par)
dev.off()


############################################################
# 32. FIGURE 1: COVERAGE AMONG SELECTED INTERVALS
############################################################

pdf(
  "monte-carlo-coverage-selected.pdf",
  width = 8.5,
  height = 6.5
)

old.par <- par(
  mfrow = c(2, 2),
  mar = c(4, 4.5, 3, 1)
)

snr.colors <- c("1000" = "red", "100" = "blue", "1" = "gray")
snr.pch <- c("1000" = 16, "100" = 17, "1" = 18)

for (tau in tau.values) {
  
  df.tau <- mc.summary[mc.summary$tau == tau, ]
  
  plot(
    NA,
    xlim = range(df.tau$T),
    ylim = c(0, 1),
    xlab = expression(T),
    ylab = expression(hat(p)[sel](T)),
    main = bquote(tau == .(tau))
  )
  
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  
  for (snr in SNR.values) {
    
    d <- df.tau[df.tau$SNR == snr, ]
    d <- d[order(d$T), ]
    
    ci_lower <- d$coverage_selected - 1.96 * d$MC_SE_selected
    ci_upper <- d$coverage_selected + 1.96 * d$MC_SE_selected
    
    ci_lower <- pmax(ci_lower, 0)
    ci_upper <- pmin(ci_upper, 1)
    
    polygon(
      c(d$T, rev(d$T)),
      c(ci_lower, rev(ci_upper)),
      col = adjustAlpha(snr.colors[as.character(snr)], 0.2),
      border = NA
    )
    
    lines(
      d$T,
      d$coverage_selected,
      type = "b",
      pch = snr.pch[as.character(snr)],
      col = snr.colors[as.character(snr)],
      lwd = 1.5
    )
  }
  
  legend(
    "topleft",
    legend = paste("SNR =", SNR.values),
    col = snr.colors[as.character(SNR.values)],
    lty = 1,
    pch = snr.pch[as.character(SNR.values)],
    cex = 0.8,
    bty = "o",
    bg = "white"
  )
}

par(old.par)
dev.off()


############################################################
# 33. FIGURE 2: INTERVAL-SELECTION RATE
############################################################

pdf(
  "monte-carlo-interval-selection-rate.pdf",
  width = 8.5,
  height = 6.5
)

old.par <- par(
  mfrow = c(2, 2),
  mar = c(4, 4.5, 3, 1)
)

snr.colors <- c("1000" = "red", "100" = "blue", "1" = "gray")
snr.pch <- c("1000" = 16, "100" = 17, "1" = 18)

for (tau in tau.values) {
  
  df.tau <- mc.summary[mc.summary$tau == tau, ]
  
  plot(
    NA,
    xlim = range(df.tau$T),
    ylim = c(0, 1),
    xlab = expression(T),
    ylab = expression(hat(q)(T)),
    main = bquote(tau == .(tau))
  )
  
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  
  for (snr in SNR.values) {
    
    d <- df.tau[df.tau$SNR == snr, ]
    d <- d[order(d$T), ]
    
    p <- d$selection_rate
    n <- d$M
    
    se <- sqrt(p * (1 - p) / n)
    ci_lower <- p - 1.96 * se
    ci_upper <- p + 1.96 * se
    
    ci_lower <- pmax(ci_lower, 0)
    ci_upper <- pmin(ci_upper, 1)
    
    polygon(
      c(d$T, rev(d$T)),
      c(ci_lower, rev(ci_upper)),
      col = adjustAlpha(snr.colors[as.character(snr)], 0.2),
      border = NA
    )
    
    lines(
      d$T,
      d$selection_rate,
      type = "b",
      pch = snr.pch[as.character(snr)],
      col = snr.colors[as.character(snr)],
      lwd = 1.5
    )
  }
  
  legend(
    "topleft",
    legend = paste("SNR =", SNR.values),
    col = snr.colors[as.character(SNR.values)],
    lty = 1,
    pch = snr.pch[as.character(SNR.values)],
    cex = 0.8,
    bty = "o",
    bg = "white"
  )
}

par(old.par)
dev.off()
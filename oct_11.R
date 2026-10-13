############################################################
# 9'. SELECTION FUNCTIONS WITH A CHOICE OF RULE
#
#   rule = "narrowest":
#       k_hat(T) = argmin_{k in K(T)} w_k(T)
#       (ties: smaller k, smaller T_entry, smaller T_exit)
#
#   rule = "first_k":
#       k_hat(T) = min K(T)
#
# These definitions replace those of Sections 9-11.
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
  
  T.current <- unique(stage2.table$T_current)[1]
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
    order(
      candidates$w_k,
      candidates$k,
      candidates$T_entry,
      candidates$T_exit
    )
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

select_final_candidate <- function(
    stage2.table,
    L,
    rule = "narrowest"
) {
  select_narrowest_candidate(
    stage2.table = stage2.table,
    containment.column = "peak.in.strict",
    L = L,
    rule = rule
  )
}

select_estimate_candidate <- function(
    stage2.table,
    L,
    rule = "narrowest"
) {
  select_narrowest_candidate(
    stage2.table = stage2.table,
    containment.column = "estimate.in.window",
    L = L,
    rule = rule
  )
}


############################################################
# 30. MONTE CARLO EMPIRICAL COVERAGE OF OPERATIONAL INTERVALS
#
# For every replicate and every horizon T, the candidate
# table is computed ONCE, and each selection rule is applied
# to the same table. The true L is used only to assess
# coverage, never to select the interval.
############################################################

run_monte_carlo_coverage <- function(
    M = 100,
    seed_start = 10000,
    tau.grid = c(0.25, 0.50, 0.75, 0.95),
    SNR.grid = c(1000, 100, 1),
    A.true = 1,
    L.true = 22,
    S.true = 10,
    h.max = 30,
    rules = c("narrowest", "first_k")
) {
  
  results <- list()
  counter <- 1L
  combo_id <- 0L
  
  for (tau in tau.grid) {
    for (SNR in SNR.grid) {
      
      combo_id <- combo_id + 1L
      
      cat("\nStarting: tau =", tau, "| SNR =", SNR, "\n")
      
      for (m in seq_len(M)) {
        
        # Same seed scheme as before, so the trajectories
        # coincide with those of the previous simulation.
        current_seed <- seed_start + combo_id * M + m
        
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
          
          for (rule in rules) {
            
            # Selection uses hat.l(T), never the true L.
            selected <- select_estimate_candidate(
              stage2.table = stage2.table,
              L = L.true,
              rule = rule
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
              covers_L <- entry <= L.true && L.true <= exit
              k.selected <- selected$k[1]
            } else {
              entry <- NA_real_
              exit <- NA_real_
              width <- NA_real_
              covers_L <- FALSE
              k.selected <- NA_real_
            }
            
            results[[counter]] <- data.frame(
              rule = rule,
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
              covers_L = covers_L,
              stringsAsFactors = FALSE
            )
            
            counter <- counter + 1L
          }
        }
        
        if (m %% 25 == 0L) {
          cat("  Completed", m, "of", M, "replicates\n")
        }
      }
    }
  }
  
  raw <- do.call(rbind, results)
  rownames(raw) <- NULL
  
  # Summary for each (rule, tau, SNR, T).
  groups <- split(
    raw,
    interaction(raw$rule, raw$tau, raw$SNR, raw$T, drop = TRUE)
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
      sqrt(coverage.selected * (1 - coverage.selected) / n.selected)
    } else {
      NA_real_
    }
    
    data.frame(
      rule = df$rule[1],
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
      },
      mean_k = if (n.selected > 0) {
        mean(df$k[df$selected])
      } else {
        NA_real_
      },
      prop_k0 = if (n.selected > 0) {
        mean(df$k[df$selected] == 0)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  
  summary <- do.call(rbind, summary.list)
  rownames(summary) <- NULL
  
  summary <- summary[
    order(summary$rule, summary$tau, -summary$SNR, summary$T),
  ]
  rownames(summary) <- NULL
  
  list(raw = raw, summary = summary)
}


############################################################
# 31. RUN THE MONTE CARLO EXPERIMENT
############################################################

# Use M = 100 first to test the code and estimate runtime.

M.replicates <- 1000

mc.results <- run_monte_carlo_coverage(
  M = M.replicates,
  seed_start = 10000,
  tau.grid = tau.values,
  SNR.grid = SNR.values,
  A.true = A,
  L.true = L,
  S.true = S,
  h.max = h,
  rules = c("narrowest", "first_k")
)

mc.raw <- mc.results$raw
mc.summary <- mc.results$summary

write.csv(
  mc.raw,
  "monte-carlo-interval-coverage-raw-both-rules.csv",
  row.names = FALSE
)

write.csv(
  mc.summary,
  "monte-carlo-interval-coverage-summary-both-rules.csv",
  row.names = FALSE
)

# Summaries by rule (the "narrowest" file keeps the content
# of the previous summary file).
mc.summary.narrowest <- mc.summary[mc.summary$rule == "narrowest", ]
mc.summary.first <- mc.summary[mc.summary$rule == "first_k", ]

write.csv(
  mc.summary.narrowest,
  "monte-carlo-interval-coverage-summary.csv",
  row.names = FALSE
)

write.csv(
  mc.summary.first,
  "monte-carlo-interval-coverage-summary-first-k.csv",
  row.names = FALSE
)

print(head(mc.summary, 30))


############################################################
# 31b. COMPARISON OF THE TWO SELECTION RULES
############################################################

cmp <- merge(
  mc.summary.narrowest,
  mc.summary.first,
  by = c("tau", "SNR", "T"),
  suffixes = c("_narrow", "_first")
)

# CHECK: the selection rate depends only on whether K(T) is
# empty, so it must be identical under the two rules.
if (!all(cmp$n_selected_narrow == cmp$n_selected_first)) {
  warning(
    "Selection counts differ between the rules: check the code."
  )
} else {
  cat("\nCheck passed: selection rates are identical.\n")
}

comparison.table <- cmp[
  ,
  c(
    "tau", "SNR", "T",
    "selection_rate_narrow",
    "coverage_selected_narrow", "coverage_selected_first",
    "mean_width_narrow", "mean_width_first",
    "prop_k0_narrow", "prop_k0_first"
  )
]

comparison.table <- comparison.table[
  order(comparison.table$tau, -comparison.table$SNR,
        comparison.table$T),
]
rownames(comparison.table) <- NULL

write.csv(
  comparison.table,
  "monte-carlo-comparison-narrowest-vs-first-k.csv",
  row.names = FALSE
)

cat("\n========================================\n")
cat("COMPARISON AT T = 22 AND T = 30\n")
cat("========================================\n\n")

print(
  comparison.table[comparison.table$T %in% c(22, 30), ],
  digits = 3
)

# Compact LaTeX table at T = 30
print(
  knitr::kable(
    comparison.table[comparison.table$T == 30, ],
    format = "latex",
    booktabs = TRUE,
    digits = 3
  )
)


############################################################
# HELPER: adjustAlpha
############################################################

adjustAlpha <- function(col, alpha = 0.5) {
  rgb_vals <- col2rgb(col) / 255
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = alpha)
}


############################################################
# 32. GENERIC FIGURE FUNCTION (replaces the duplicated
#     sections 32 and 33)
#
# which = "coverage"  : coverage among selected intervals
# which = "selection" : interval-selection rate
############################################################

plot_mc_figure <- function(
    summary.table,
    which = c("coverage", "selection"),
    file
) {
  
  which <- match.arg(which)
  
  snr.colors <- c("1000" = "red", "100" = "blue", "1" = "gray")
  snr.pch <- c("1000" = 16, "100" = 17, "1" = 18)
  
  pdf(file, width = 8.5, height = 6.5)
  
  old.par <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1))
  
  for (tau in tau.values) {
    
    df.tau <- summary.table[summary.table$tau == tau, ]
    
    plot(
      NA,
      xlim = range(df.tau$T),
      ylim = c(0, 1),
      xlab = expression(T),
      ylab = if (which == "coverage") {
        expression(hat(p)[sel](T))
      } else {
        expression(hat(q)(T))
      },
      main = bquote(tau == .(tau))
    )
    
    grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
    
    for (snr in SNR.values) {
      
      d <- df.tau[df.tau$SNR == snr, ]
      d <- d[order(d$T), ]
      
      if (which == "coverage") {
        y <- d$coverage_selected
        se <- d$MC_SE_selected
      } else {
        y <- d$selection_rate
        se <- sqrt(y * (1 - y) / d$M)
      }
      
      ci_lower <- pmax(y - 1.96 * se, 0)
      ci_upper <- pmin(y + 1.96 * se, 1)
      
      # polygon() cannot handle NA values: draw only the
      # horizons where the band is defined.
      ok <- is.finite(ci_lower) & is.finite(ci_upper)
      
      if (any(ok)) {
        polygon(
          c(d$T[ok], rev(d$T[ok])),
          c(ci_lower[ok], rev(ci_upper[ok])),
          col = adjustAlpha(snr.colors[as.character(snr)], 0.2),
          border = NA
        )
      }
      
      lines(
        d$T,
        y,
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
}


############################################################
# 33. FIGURES FOR BOTH RULES
############################################################

# Narrowest-width rule (file names as in the manuscript)
plot_mc_figure(
  mc.summary.narrowest, "coverage",
  "monte-carlo-coverage-selected.pdf"
)

plot_mc_figure(
  mc.summary.narrowest, "selection",
  "monte-carlo-interval-selection-rate.pdf"
)

# First-admissible-k rule
plot_mc_figure(
  mc.summary.first, "coverage",
  "monte-carlo-coverage-selected-first-k.pdf"
)

# The selection-rate figure is identical for the two rules,
# so it is not drawn a second time.


############################################################
# 34. OPTIONAL: DIRECT COMPARISON FIGURE AT EACH T
#     (coverage and mean width, one panel per tau, SNR = 1)
############################################################

d1 <- comparison.table[comparison.table$SNR == 1, ]

pdf("monte-carlo-comparison-rules-SNR1.pdf",
    width = 8.5, height = 6.5)

old.par <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1))

for (tau in tau.values) {
  
  d <- d1[d1$tau == tau, ]
  d <- d[order(d$T), ]
  
  plot(
    NA,
    xlim = range(d$T),
    ylim = c(0, 1),
    xlab = expression(T),
    ylab = "Coverage among selected",
    main = bquote(tau == .(tau) * ", SNR = 1")
  )
  
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  
  lines(d$T, d$coverage_selected_narrow, type = "b",
        pch = 16, col = "blue", lwd = 1.5)
  lines(d$T, d$coverage_selected_first, type = "b",
        pch = 17, col = "red", lwd = 1.5, lty = 2)
  
  legend(
    "topleft",
    legend = c("narrowest", "first k"),
    col = c("blue", "red"),
    lty = c(1, 2),
    pch = c(16, 17),
    cex = 0.8, bty = "o", bg = "white"
  )
}

par(old.par)
dev.off()
############################################################
# STRUCTURAL ANALYSIS OF CONSECUTIVE CANDIDATE INTERVALS
# SNR = 1000
#
# Purpose:
#
#   1. Analyze the evolution of
#        T_exit^(k)(T)
#        T_entry^(k)(T)
#
#   2. Compare
#        T_entry^(k)(T)
#        with
#        T_exit^(k+1)(T)
#
#   3. Produce a compact summary by tau.
#
# Required object:
#
#   combined.diagnostics
#
# Required columns:
#
#   tau_level
#   SNR_level
#   T_current
#   k
#   T_exit
#   T_entry
#
# Main outputs:
#
#   1. fig_exit_entry_times_SNR1000.pdf
#      Combined exit/entry trajectories.
#
#   2. fig_entry_vs_next_exit_SNR1000.pdf
#      Difference
#      T_entry^(k)(T) - T_exit^(k+1)(T).
#
#   3. structural_summary_SNR1000.csv
#      Compact summary by tau.
#
# Additional working files:
#
#   4. relationship_table_SNR1000.csv
#   5. entry_exit_link_table_SNR1000.csv
#
############################################################


############################################################
# 0. CHECK INPUT
############################################################

required.columns <- c(
  "tau_level",
  "SNR_level",
  "T_current",
  "k",
  "T_exit",
  "T_entry"
)

missing.columns <- setdiff(
  required.columns,
  names(combined.diagnostics)
)

if (length(missing.columns) > 0) {
  stop(
    paste(
      "The following required columns are missing:",
      paste(missing.columns, collapse = ", ")
    )
  )
}


############################################################
# 1. KEEP SNR = 1000
############################################################

diagnostics1000 <- combined.diagnostics[
  combined.diagnostics$SNR_level == 1000 &
    !is.na(combined.diagnostics$T_exit) &
    !is.na(combined.diagnostics$T_entry),
  ,
  drop = FALSE
]

if (nrow(diagnostics1000) == 0) {
  stop("No valid observations found for SNR = 1000.")
}


############################################################
# 2. ORDER THE DATA
############################################################

diagnostics1000 <- diagnostics1000[
  order(
    diagnostics1000$tau_level,
    diagnostics1000$T_current,
    diagnostics1000$k
  ),
  ,
  drop = FALSE
]


############################################################
# 3. IDENTIFY CONSECUTIVE CANDIDATE PAIRS
############################################################
#
# For each fixed (tau,T), candidate k is matched with
# candidate k+1.
#
############################################################

consecutive1000 <- merge(
  diagnostics1000,
  diagnostics1000,
  by = c(
    "tau_level",
    "SNR_level",
    "T_current"
  ),
  suffixes = c("_k", "_kp1")
)

consecutive1000 <- consecutive1000[
  consecutive1000$k_kp1 ==
    consecutive1000$k_k + 1,
  ,
  drop = FALSE
]

consecutive1000 <- consecutive1000[
  order(
    consecutive1000$tau_level,
    consecutive1000$T_current,
    consecutive1000$k_k
  ),
  ,
  drop = FALSE
]

if (nrow(consecutive1000) == 0) {
  stop(
    "No consecutive candidate pairs were found for SNR = 1000."
  )
}


############################################################
# 4. COMPUTE THE STRUCTURAL DIFFERENCES
############################################################

# ----------------------------------------------------------
# 4.1 Difference between consecutive exit times
#
# D_exit^(k)(T)
#   = T_exit^(k)(T) - T_exit^(k+1)(T)
#
# By construction, this should be strictly positive whenever
# both candidates exist.
# ----------------------------------------------------------

consecutive1000$exit_difference <-
  consecutive1000$T_exit_k -
  consecutive1000$T_exit_kp1


# ----------------------------------------------------------
# 4.2 Entry time minus next exit time
#
# D_link^(k)(T)
#   = T_entry^(k)(T) - T_exit^(k+1)(T)
#
# The conjectured identity is equivalent to
#
# D_link^(k)(T) = 0.
# ----------------------------------------------------------

consecutive1000$entry_exit_difference <-
  consecutive1000$T_entry_k -
  consecutive1000$T_exit_kp1


# ----------------------------------------------------------
# 4.3 Difference between consecutive entry times
#
# D_entry^(k)(T)
#   = T_entry^(k)(T) - T_entry^(k+1)(T)
#
# Positive values correspond to strictly decreasing entry times.
# ----------------------------------------------------------

consecutive1000$entry_difference <-
  consecutive1000$T_entry_k -
  consecutive1000$T_entry_kp1


############################################################
# 5. BASIC INFORMATION
############################################################

cat("\n")
cat("============================================================\n")
cat("STRUCTURAL ANALYSIS: SNR = 1000\n")
cat("============================================================\n\n")

cat(
  "Number of valid observations:",
  nrow(diagnostics1000),
  "\n"
)

cat(
  "Number of consecutive candidate pairs:",
  nrow(consecutive1000),
  "\n\n"
)


############################################################
# 6. QUANTILE LEVELS
############################################################

tau.values <- sort(
  unique(diagnostics1000$tau_level)
)

cat(
  "Quantile levels:",
  paste(tau.values, collapse = ", "),
  "\n\n"
)


############################################################
# 7. CHECK EXIT ORDERING
############################################################
#
# The definition of the candidate exit times gives
#
# T_exit^(k+1)(T) < T_exit^(k)(T)
#
# whenever both candidates exist.
#
############################################################

exit.violations <- consecutive1000[
  consecutive1000$exit_difference <= 0,
  ,
  drop = FALSE
]

cat("------------------------------------------------------------\n")
cat("EXIT ORDERING\n")
cat("------------------------------------------------------------\n\n")

if (nrow(exit.violations) == 0) {
  
  cat(
    "No violations found.\n",
    "For every observed consecutive pair,\n",
    "T_exit^(k+1)(T) < T_exit^(k)(T).\n\n"
  )
  
} else {
  
  cat(
    "WARNING:",
    nrow(exit.violations),
    "violations found.\n\n"
  )
  
  print(exit.violations)
  
  cat("\n")
}


############################################################
# 8. CHECK ENTRY-NEXT-EXIT RELATIONSHIP
############################################################
#
# Main empirical question:
#
# T_entry^(k)(T) ?= T_exit^(k+1)(T)
#
############################################################

link.equal <- consecutive1000[
  consecutive1000$entry_exit_difference == 0,
  ,
  drop = FALSE
]

link.not.equal <- consecutive1000[
  consecutive1000$entry_exit_difference != 0,
  ,
  drop = FALSE
]

cat("------------------------------------------------------------\n")
cat("ENTRY - NEXT EXIT RELATIONSHIP\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Total consecutive pairs:",
  nrow(consecutive1000),
  "\n"
)

cat(
  "Exact equalities:",
  nrow(link.equal),
  "\n"
)

cat(
  "Non-equalities:",
  nrow(link.not.equal),
  "\n"
)

cat(
  "Proportion of exact equalities:",
  round(
    mean(
      consecutive1000$entry_exit_difference == 0
    ),
    4
  ),
  "\n\n"
)

if (nrow(link.not.equal) == 0) {
  
  cat(
    "The identity\n",
    "T_entry^(k)(T) = T_exit^(k+1)(T)\n",
    "holds for every observed consecutive pair.\n\n"
  )
  
} else {
  
  cat(
    "The identity is NOT universal.\n",
    "The range of the difference is:\n"
  )
  
  print(
    range(
      consecutive1000$entry_exit_difference
    )
  )
  
  cat("\n")
  
}


############################################################
# 9. CHECK ENTRY ORDERING
############################################################
#
# Test whether
#
# T_entry^(k+1)(T) < T_entry^(k)(T).
#
# Hence:
#
# entry_difference > 0
#
# means strict ordering,
#
# entry_difference = 0
#
# means equality.
#
############################################################

entry.strict.violations <- consecutive1000[
  consecutive1000$entry_difference <= 0,
  ,
  drop = FALSE
]

entry.ordering.nonpositive <- consecutive1000[
  consecutive1000$entry_difference < 0,
  ,
  drop = FALSE
]

entry.equal <- consecutive1000[
  consecutive1000$entry_difference == 0,
  ,
  drop = FALSE
]

cat("------------------------------------------------------------\n")
cat("ENTRY ORDERING\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Strict decreases:",
  sum(
    consecutive1000$entry_difference > 0
  ),
  "\n"
)

cat(
  "Equal entry times:",
  nrow(entry.equal),
  "\n"
)

cat(
  "Actual reversals:",
  nrow(entry.ordering.nonpositive),
  "\n\n"
)

if (nrow(entry.ordering.nonpositive) == 0) {
  
  cat(
    "No reversals found.\n",
    "Thus the entry times are empirically nonincreasing:\n",
    "T_entry^(k+1)(T) <= T_entry^(k)(T).\n\n"
  )
  
} else {
  
  cat(
    "WARNING: entry-time reversals were found.\n",
    nrow(entry.ordering.nonpositive),
    "pairs satisfy\n",
    "T_entry^(k+1)(T) > T_entry^(k)(T).\n\n"
  )
  
  print(
    entry.ordering.nonpositive[
      ,
      c(
        "tau_level",
        "T_current",
        "k_k",
        "T_entry_k",
        "T_entry_kp1",
        "entry_difference"
      ),
      drop = FALSE
    ]
  )
  
  cat("\n")
}


############################################################
# 10. COMPLETE WORKING RELATIONSHIP TABLE
############################################################
#
# This table is retained for detailed inspection, but is not
# intended for direct inclusion in the manuscript.
#
############################################################

relationship.table1000 <- consecutive1000[
  ,
  c(
    "tau_level",
    "T_current",
    "k_k",
    "T_exit_k",
    "T_entry_k",
    "T_exit_kp1",
    "T_entry_kp1",
    "exit_difference",
    "entry_exit_difference",
    "entry_difference"
  ),
  drop = FALSE
]

write.csv(
  relationship.table1000,
  "relationship_table_SNR1000.csv",
  row.names = FALSE
)


############################################################
# 11. COMPACT ENTRY-EXIT LINK TABLE
############################################################

link.table1000 <- consecutive1000[
  ,
  c(
    "tau_level",
    "T_current",
    "k_k",
    "T_entry_k",
    "T_exit_kp1",
    "entry_exit_difference"
  ),
  drop = FALSE
]

write.csv(
  link.table1000,
  "entry_exit_link_table_SNR1000.csv",
  row.names = FALSE
)


############################################################
# 12. COMPACT SUMMARY TABLE BY TAU
############################################################
#
# One row per tau.
#
# Columns:
#
#   tau
#   n_pairs
#   n_equal
#   proportion_equal
#   exit_ordering_holds
#   n_entry_reversals
#   entry_nonincreasing
#
############################################################

structural.summary1000 <- do.call(
  rbind,
  lapply(
    tau.values,
    function(tau) {
      
      tmp <- consecutive1000[
        consecutive1000$tau_level == tau,
        ,
        drop = FALSE
      ]
      
      exit.ordering.holds <- all(
        tmp$exit_difference > 0
      )
      
      n.entry.reversals <- sum(
        tmp$entry_difference < 0
      )
      
      entry.nonincreasing <- all(
        tmp$entry_difference >= 0
      )
      
      data.frame(
        tau = tau,
        n_pairs = nrow(tmp),
        n_equal =
          sum(
            tmp$entry_exit_difference == 0
          ),
        proportion_equal =
          mean(
            tmp$entry_exit_difference == 0
          ),
        exit_ordering_holds =
          exit.ordering.holds,
        n_entry_reversals =
          n.entry.reversals,
        entry_nonincreasing =
          entry.nonincreasing,
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(structural.summary1000) <- NULL


############################################################
# 13. DISPLAY AND SAVE SUMMARY TABLE
############################################################

cat("\n")
cat("============================================================\n")
cat("COMPACT STRUCTURAL SUMMARY: SNR = 1000\n")
cat("============================================================\n\n")

print(
  structural.summary1000,
  row.names = FALSE
)

write.csv(
  structural.summary1000,
  "structural_summary_SNR1000.csv",
  row.names = FALSE
)


############################################################
# 14. GLOBAL STRUCTURAL SUMMARY
############################################################

global.equality.proportion <- mean(
  consecutive1000$entry_exit_difference == 0
)

global.exit.ordering <- all(
  consecutive1000$exit_difference > 0
)

global.entry.nonincreasing <- all(
  consecutive1000$entry_difference >= 0
)

global.entry.strict <- all(
  consecutive1000$entry_difference > 0
)

cat("\n")
cat("------------------------------------------------------------\n")
cat("GLOBAL STRUCTURAL SUMMARY\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Exit ordering holds universally:",
  global.exit.ordering,
  "\n"
)

cat(
  "Proportion of exact entry-next-exit equalities:",
  round(
    global.equality.proportion,
    4
  ),
  "\n"
)

cat(
  "Entry times are nonincreasing:",
  global.entry.nonincreasing,
  "\n"
)

cat(
  "Entry times are strictly decreasing:",
  global.entry.strict,
  "\n\n"
)


############################################################
# 15. FIGURE 1:
#     COMBINED EXIT AND ENTRY TRAJECTORIES
#     SNR = 1000
############################################################
#
# One panel for each tau.
#
# Solid lines  : T_exit^(k)(T)  [Blue]
# Dashed lines : T_entry^(k)(T) [Red]
#
# Candidate index k is distinguished by plotting symbol.
#
# A horizontal line at L is added because the simulation uses
# L = 22.
#
############################################################

# Recover L from the diagnostics if available.
#
# If L is not present in combined.diagnostics, use 22, which
# corresponds to the simulation considered here.

if ("L" %in% names(combined.diagnostics)) {
  
  L.values <- unique(
    combined.diagnostics$L[
      !is.na(combined.diagnostics$L)
    ]
  )
  
  if (length(L.values) == 1) {
    L.plot <- L.values
  } else {
    L.plot <- 22
  }
  
} else {
  
  L.plot <- 22
  
}


# pdf(
#   "fig_exit_entry_times_SNR1000.pdf",
#   width = 10,
#   height = 8
# )
# 
# old.par <- par(
#   mfrow = c(2, 2),
#   mar = c(4, 4, 3, 1)
# )

pdf(
  "fig_exit_entry_times_SNR1000.pdf",
  width = 8.5,
  height = 6.5
)

old.par <- par(
  mfrow = c(2, 2),
  mar = c(4, 4.5, 3, 1)
)

# Define plotting symbols for each k value
pch.symbols <- c(16, 17, 15, 18, 19, 8, 3, 4, 6, 7)

for (tau in tau.values) {
  
  tmp <- diagnostics1000[
    diagnostics1000$tau_level == tau,
    ,
    drop = FALSE
  ]
  
  k.values <- sort(
    unique(tmp$k)
  )
  
  x.range <- range(
    tmp$T_current
  )
  
  y.range <- range(
    c(
      tmp$T_exit,
      tmp$T_entry,
      L.plot
    ),
    finite = TRUE
  )
  
  plot(
    NA,
    xlim = x.range,
    ylim = y.range,
    xlab = expression(T),
    ylab = "Time",
    main = bquote(tau == .(tau) * ", SNR = 1000")
  )
  
  # Peak location
  
  abline(
    h = L.plot,
    lty = 2
  )
  
  # Candidate curves
  
  for (i in seq_along(k.values)) {
    
    k <- k.values[i]
    pch.val <- pch.symbols[(i - 1) %% length(pch.symbols) + 1]
    
    tmp.k <- tmp[
      tmp$k == k,
      ,
      drop = FALSE
    ]
    
    # Exit time: solid line, blue
    
    lines(
      tmp.k$T_current,
      tmp.k$T_exit,
      type = "o",
      pch = pch.val,
      lty = 1,
      col = "blue"
    )
    
    # Entry time: dashed line, red
    
    lines(
      tmp.k$T_current,
      tmp.k$T_entry,
      type = "o",
      pch = pch.val,
      lty = 2,
      col = "red"
    )
  }
  
  # Create legend with framed background
  
  legend.labels <- c(
    "Exit",
    "Entry",
    expression(L),
    paste0("k = ", k.values)
  )
  
  legend.ltys <- c(
    1,
    2,
    3,
    rep(1, length(k.values))
  )
  
  legend.pchs <- c(
    NA,
    NA,
    NA,
    pch.symbols[1:length(k.values)]
  )
  
  legend.cols <- c(
    "blue",
    "red",
    "black",
    rep("black", length(k.values))
  )
  
  legend(
    "topleft",
    legend = legend.labels,
    lty = legend.ltys,
    pch = legend.pchs,
    col = legend.cols,
    bty = "o",
    bg = "white",
    cex = 0.8
  )
  
}

par(old.par)

dev.off()


############################################################
# 16. FIGURE 2:
#     ENTRY^(k)(T) - EXIT^(k+1)(T)
#     SNR = 1000
############################################################
#
# This is the principal diagnostic for the colleague's
# question.
#
# A value of zero means:
#
#   T_entry^(k)(T) = T_exit^(k+1)(T).
#
# Candidate index k is distinguished by line color:
#   k = 0: blue
#   k = 1: red
#   k = 2: black
#
############################################################

# pdf(
#   "fig_entry_vs_next_exit_SNR1000.pdf",
#   width = 10,
#   height = 8
# )
# 
# old.par <- par(
#   mfrow = c(2, 2),
#   mar = c(4, 5, 3, 1)
# )

pdf(
  "fig_entry_vs_next_exit_SNR1000.pdf",
  width = 8.5,
  height = 6.5
)

old.par <- par(
  mfrow = c(2, 2),
  mar = c(4, 4.5, 3, 1)
)

# Define colors for each k value
k.colors <- c("blue", "red", "black", "green", "purple", "orange", "brown", "pink", "gray", "cyan")

for (tau in tau.values) {
  
  tmp <- consecutive1000[
    consecutive1000$tau_level == tau,
    ,
    drop = FALSE
  ]
  
  if (nrow(tmp) == 0) {
    next
  }
  
  k.values <- sort(
    unique(tmp$k_k)
  )
  
  y.range <- range(
    c(
      tmp$entry_exit_difference,
      0
    ),
    finite = TRUE
  )
  
  # Avoid degenerate y-axis if all values are zero.
  
  if (diff(y.range) == 0) {
    
    y.range <- y.range + c(
      -1,
      1
    )
  }
  
  plot(
    NA,
    xlim = range(
      tmp$T_current
    ),
    ylim = y.range,
    xlab = expression(T),
    ylab = expression(
      T[entry]^{(k)} * (T) -
        T[exit]^{(k+1)} * (T)
    ),
    main = bquote(tau == .(tau) * ", SNR = 1000"),
    mgp = c(2.5, 1, 0)
  )
  
  # Equality line
  
  abline(
    h = 0,
    lty = 2
  )
  
  for (i in seq_along(k.values)) {
    
    k <- k.values[i]
    col.val <- k.colors[(i - 1) %% length(k.colors) + 1]
    
    tmp.k <- tmp[
      tmp$k_k == k,
      ,
      drop = FALSE
    ]
    
    lines(
      tmp.k$T_current,
      tmp.k$entry_exit_difference,
      type = "o",
      pch = 16,
      lty = 1,
      col = col.val
    )
  }
  
  legend(
    "bottomright",
    legend = paste0(
      "k = ",
      k.values
    ),
    lty = 1,
    pch = 16,
    col = k.colors[1:length(k.values)],
    bty = "o",
    bg = "white",
    cex = 0.8
  )
  
}

par(old.par)

dev.off()


############################################################
# 17. ADDITIONAL NUMERICAL DIAGNOSTICS
############################################################
#
# These are not manuscript figures.
# They are printed to the console to help interpret the
# structural pattern.
#
############################################################

cat("\n")
cat("============================================================\n")
cat("ADDITIONAL NUMERICAL DIAGNOSTICS\n")
cat("============================================================\n\n")


cat("Range of exit differences:\n")

print(
  range(
    consecutive1000$exit_difference
  )
)

cat("\n")


cat(
  "Unique entry - next exit differences:\n"
)

print(
  sort(
    unique(
      consecutive1000$entry_exit_difference
    )
  )
)

cat("\n")


cat(
  "Unique entry differences:\n"
)

print(
  sort(
    unique(
      consecutive1000$entry_difference
    )
  )
)

cat("\n")


############################################################
# 18. FREQUENCY TABLE OF ENTRY-NEXT-EXIT DIFFERENCE
############################################################

cat("------------------------------------------------------------\n")
cat("FREQUENCY OF ENTRY - NEXT EXIT DIFFERENCE\n")
cat("------------------------------------------------------------\n\n")

print(
  table(
    consecutive1000$entry_exit_difference
  )
)

cat("\n")


############################################################
# 19. FREQUENCY TABLE BY TAU
############################################################

cat("------------------------------------------------------------\n")
cat("ENTRY - NEXT EXIT DIFFERENCE BY TAU\n")
cat("------------------------------------------------------------\n\n")

print(
  with(
    consecutive1000,
    table(
      tau_level,
      entry_exit_difference
    )
  )
)

cat("\n")


############################################################
# 20. FREQUENCY TABLE OF ENTRY DIFFERENCES
############################################################

cat("------------------------------------------------------------\n")
cat("ENTRY DIFFERENCE BY TAU\n")
cat("------------------------------------------------------------\n\n")

print(
  with(
    consecutive1000,
    table(
      tau_level,
      entry_difference
    )
  )
)

cat("\n")


############################################################
# 21. FINAL INTERPRETATION
############################################################

cat("============================================================\n")
cat("AUTOMATIC INTERPRETATION\n")
cat("============================================================\n\n")


# ----------------------------------------------------------
# Exit ordering
# ----------------------------------------------------------

if (global.exit.ordering) {
  
  cat(
    "[1] EXIT ORDERING\n\n",
    "    Confirmed for every observed consecutive pair:\n\n",
    "    T_exit^(k+1)(T) < T_exit^(k)(T).\n\n",
    "    This ordering follows from the definition of the\n",
    "    exit candidates as ordered strict-decrease points.\n\n"
  )
  
} else {
  
  cat(
    "[1] EXIT ORDERING\n\n",
    "    Violations were detected in the numerical output.\n\n"
  )
}


# ----------------------------------------------------------
# Entry-next-exit identity
# ----------------------------------------------------------

if (
  global.equality.proportion == 1
) {
  
  cat(
    "[2] ENTRY-NEXT-EXIT RELATIONSHIP\n\n",
    "    Exact equality holds throughout the SNR = 1000\n",
    "    simulation:\n\n",
    "    T_entry^(k)(T) = T_exit^(k+1)(T).\n\n",
    "    This gives the simplified candidate structure:\n\n",
    "    R_k(T) = [T_exit^(k+1)(T), T_exit^(k)(T)].\n\n"
  )
  
} else {
  
  cat(
    "[2] ENTRY-NEXT-EXIT RELATIONSHIP\n\n",
    "    The identity is NOT universal at SNR = 1000.\n\n",
    "    Proportion of exact equalities = ",
    round(
      global.equality.proportion,
      4
    ),
    ".\n\n",
    "    Therefore the entry and next-exit times should not\n",
    "    be identified in general without an additional\n",
    "    theoretical argument.\n\n"
  )
}


# ----------------------------------------------------------
# Entry ordering
# ----------------------------------------------------------

if (global.entry.strict) {
  
  cat(
    "[3] ENTRY ORDERING\n\n",
    "    Entry times are strictly decreasing in k for every\n",
    "    observed consecutive pair:\n\n",
    "    T_entry^(k+1)(T) < T_entry^(k)(T).\n\n"
  )
  
} else if (global.entry.nonincreasing) {
  
  cat(
    "[3] ENTRY ORDERING\n\n",
    "    Entry times are empirically nonincreasing in k:\n\n",
    "    T_entry^(k+1)(T) <= T_entry^(k)(T).\n\n",
    "    Some consecutive candidates have equal entry times,\n",
    "    so strict decrease does not hold universally.\n\n"
  )
  
} else {
  
  cat(
    "[3] ENTRY ORDERING\n\n",
    "    Entry times are not universally ordered in the\n",
    "    observed SNR = 1000 simulation.\n\n",
    "    The number of actual reversals is ",
    nrow(entry.ordering.nonpositive),
    ".\n\n"
  )
}


############################################################
# 22. LIST OF GENERATED FILES
############################################################

cat("============================================================\n")
cat("GENERATED FILES\n")
cat("============================================================\n\n")

cat(
  "Main manuscript figures:\n",
  "  - fig_exit_entry_times_SNR1000.pdf\n",
  "  - fig_entry_vs_next_exit_SNR1000.pdf\n\n"
)

cat(
  "Compact manuscript table:\n",
  "  - structural_summary_SNR1000.csv\n\n"
)

cat(
  "Detailed working files:\n",
  "  - relationship_table_SNR1000.csv\n",
  "  - entry_exit_link_table_SNR1000.csv\n\n"
)

cat("Analysis completed.\n")

############################################################
# END
############################################################

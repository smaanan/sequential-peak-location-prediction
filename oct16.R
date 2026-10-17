main.rule <- "first_k"

plot_mc_base(
  "all_k", main.rule, "p_sel_incl", NULL,
  ylab = expression(hat(p)[sel](T)),
  file = "monte-carlo-coverage-selected.pdf",
  show.band = FALSE
)
plot_mc_base(
  "all_k", main.rule, "q", "se_q",
  ylab = expression(hat(q)(T)),
  file = "monte-carlo-interval-selection-rate.pdf"
)

agg$p_all_incl <- agg$covered_incl / agg$M

# Table tab:mc_coverage_T30 (rule first_k, all_k)
t30 <- agg[agg$variant == "all_k" & agg$rule == "first_k" &
             agg$T == 30,
           c("tau", "SNR", "q", "p_sel_incl", "p_all_incl",
             "mean_width")]
t30 <- t30[order(t30$tau, -t30$SNR), ]
print(t30, digits = 3, row.names = FALSE)
print(knitr::kable(t30, format = "latex", booktabs = TRUE, digits = 3))

# Table tab:rule_comparison: coverage, width, interval score
d <- mc.raw[mc.raw$variant == "all_k" & mc.raw$selected &
              mc.raw$rule %in% c("narrowest", "first_k") &
              mc.raw$SNR == 1 & mc.raw$T %in% c(22, 26, 30), ]
d$miss <- pmax(d$T_entry - L, L - d$T_exit, 0)

for (alpha in c(0.05, 0.10, 0.20)) {
  d$iscore <- d$width + (2 / alpha) * d$miss
  cat("\nalpha =", alpha, "\n")
  print(aggregate(cbind(covered_incl, width, iscore) ~ rule + tau + T,
                  d, mean), digits = 3)
}

agg[agg$rule=="first_k"&agg$"SNR"==1&agg$T<22,]

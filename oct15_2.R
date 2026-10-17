d <- mc.raw[mc.raw$variant == "all_k" & mc.raw$selected &
              mc.raw$rule %in% c("narrowest", "first_k") &
              mc.raw$SNR == 1 & mc.raw$T >= 22, ]

# Distance from L to the interval (0 if L is inside)
d$miss <- pmax(d$T_entry - L, L - d$T_exit, 0)

alpha <- 0.10   # for the interval score; report also 0.05
d$iscore <- d$width + (2 / alpha) * d$miss   # Winkler interval score

eff <- aggregate(cbind(covered_incl, width, miss, iscore) ~
                   rule + tau + T, d, mean)
eff$cov_per_width <- eff$covered_incl / eff$width
print(eff[eff$T %in% c(22, 26, 30), ], digits = 3)

# Paired comparison: same replicate, same T
w <- reshape(d[, c("rep","tau","T","rule","width","covered_incl","iscore")],
             idvar = c("rep","tau","T"), timevar = "rule",
             direction = "wide")
w$same <- w$width.narrowest == w$width.first_k
cat("\nShare of (rep,T) where both rules give the same interval:\n")
print(aggregate(same ~ tau + T, w[w$T %in% c(22,26,30), ], mean),
      digits = 3)
w2 <- w[!w$same, ]
cat("\nWhere they differ: mean coverage and width\n")
print(aggregate(cbind(covered_incl.narrowest, covered_incl.first_k,
                      width.narrowest, width.first_k,
                      iscore.narrowest, iscore.first_k) ~ tau,
                w2[w2$T == 30, ], mean), digits = 3)

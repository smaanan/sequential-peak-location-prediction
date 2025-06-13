
direct_quantile <- function(x, w, tau, grid_length = 500) {
  a_vals <- seq(min(x), max(x), length.out = grid_length)
  objective_vals <- sapply(a_vals, function(a) {
    sum(w * quantile_loss(x - a, tau))
  })
  a_vals[which.min(objective_vals)]
}

# Prepare sample and weights
set.seed(123)
x <- rnorm(10)
w <- runif(10)

# Methods to compare:
direct <- direct_quantile(x, w, tau = 0.5)
med_matrixStats <- matrixStats::weightedMedian(x, w)
med_isotone <- isotone::weighted.median(x, w)
med_laeken <- laeken::weightedMedian(x, w)
med_our <- med_weighted(x, l = 0, s = 10)  # special, not universal
med_new <- weighted_quantile(x, w, tau = 0.5)

cat("direct :", direct, "\n")
cat("matrixStats :", med_matrixStats, "\n")
cat("isotone :", med_isotone, "\n")
cat("laeken :", med_laeken, "\n")
cat("med_weighted :", med_our, "\n")
cat("weighted_quantile :", med_new, "\n")
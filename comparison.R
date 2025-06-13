library(matrixStats)
library(isotone)
library(laeken)


# ---------------------------------------------------------------------------
# 1. Prepare data
# ---------------------------------------------------------------------------
set.seed(123)  
x <- rnorm(10)  
w <- runif(10)  

# ---------------------------------------------------------------------------
# 2. Manual 50th percentile (`med_weighted`)
# ---------------------------------------------------------------------------
med_manual <- med_weighted(x, l = 0, s = 10) 
# Note: Here `l` and `s` are not used.

# ---------------------------------------------------------------------------
# 3. Generalized algorithm
# ---------------------------------------------------------------------------
med_general <- weighted_quantile(x, w, tau = 0.5)

# ---------------------------------------------------------------------------
# 4. Standard, from `matrixStats`
# ---------------------------------------------------------------------------
med_matrixStats <- matrixStats::weightedMedian(x, w)

# ---------------------------------------------------------------------------
# 5. Validation against `isotone`
# ---------------------------------------------------------------------------
med_isotone <- isotone::weighted.median(x, w)

# ---------------------------------------------------------------------------
# 6. Validation against `laeken`
# ---------------------------------------------------------------------------
med_laeken <- laeken::weightedMedian(x, w)

# ---------------------------------------------------------------------------
# 7. Summary of results
# ---------------------------------------------------------------------------
cat("Manual (med_weighted) :", med_manual, "\n")
cat("General (weighted_quantile) :", med_general, "\n")
cat("matrixStats::weightedMedian :", med_matrixStats, "\n")
cat("isotone::weighted.median :", med_isotone, "\n")
cat("laeken::weightedMedian :", med_laeken, "\n")
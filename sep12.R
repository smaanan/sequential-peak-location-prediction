###############################################################
## Sections 2.4 (large/medium SNR) and 2.5 (small SNR, h=30):
## estimator curves paired with their oscillation curves
##
## CURRENT MODEL:
##
##   I(t) = A exp(-(t-L)^2/S^2) + e(t)
##
## Quantile-loss criterion:
##
##   rho_tau(u) = u (tau - 1{u < 0})
##
## Noise:
##
##   tau-aligned asymmetric Laplace distribution
##
## SNR:
##
##   SNR = A / sqrt(Var(e))
##
## The amplitude for fixed (l,s) is obtained as a
## weighted tau-quantile of
##
##   x(t) = I(t) exp((t-l)^2/s^2)
##
## with weights
##
##   w(t) = exp(-(t-l)^2/s^2).
###############################################################


###############################################################
## Weighted tau-quantile
###############################################################

# weighted_quantile <- function(
#     x,
#     w,
#     tau
# ) {
#   
#   if (length(x) != length(w)) {
#     stop("x and w must have the same length.")
#   }
#   
#   if (tau <= 0 || tau >= 1) {
#     stop("tau must belong to (0,1).")
#   }
#   
#   if (any(!is.finite(x)) || any(!is.finite(w))) {
#     stop("x and w must be finite.")
#   }
#   
#   if (any(w < 0) || sum(w) <= 0) {
#     stop("Weights must be non-negative and have positive sum.")
#   }
#   
#   o <- order(x)
#   
#   x <- x[o]
#   w <- w[o]
#   
#   w <- w / sum(w)
#   
#   x[which(cumsum(w) >= tau)[1]]
# }

weighted_quantile <- function(x, w, tau) {
  
  if (length(x) != length(w)) {
    stop("x and w must have the same length.")
  }
  
  if (tau <= 0 || tau >= 1) {
    stop("tau must belong to (0,1).")
  }
  
  # Drop numerically-overflowed entries. These arise when l is far
  # outside 1:T; they carry weight exp(-large) ~ 0 and therefore do
  # not affect the weighted quantile.
  finite <- is.finite(x) & is.finite(w)
  
  if (!any(finite)) {
    return(NA_real_)
  }
  
  x <- x[finite]
  w <- w[finite]
  
  if (any(w < 0) || sum(w) <= 0) {
    stop("Weights must be non-negative and have positive sum.")
  }
  
  o <- order(x)
  x <- x[o]
  w <- w[o]
  w <- w / sum(w)
  
  x[which(cumsum(w) >= tau)[1]]
}

###############################################################
## Quantile check loss
###############################################################

rho_tau <- function(
    u,
    tau
) {
  
  u * (tau - (u < 0))
}



###############################################################
## tau-aligned asymmetric Laplace random generator
##
## Density:
##
##   f_{tau,b}(x)
##     = tau(1-tau)/b * exp(-rho_tau(x)/b)
##
## Hence:
##
##   P(e < 0) = tau.
###############################################################

rald_tau <- function(
    n,
    tau,
    b
) {
  
  if (tau <= 0 || tau >= 1) {
    stop("tau must belong to (0,1).")
  }
  
  if (b <= 0) {
    stop("b must be positive.")
  }
  
  u <- runif(n)
  
  e <- numeric(n)
  
  left <- u < tau
  right <- !left
  
  e[left] <-
    b / (1 - tau) *
    log(u[left] / tau)
  
  e[right] <-
    -b / tau *
    log(
      (1 - u[right]) / (1 - tau)
    )
  
  e
}



###############################################################
## Noise variance of the tau-aligned asymmetric Laplace law
##
## Var(e)
##   = b^2
##     [tau^2 + (1-tau)^2]
##     / [tau^2 (1-tau)^2]
###############################################################

ald_variance <- function(
    tau,
    b
) {
  
  b^2 *
    (
      tau^2 + (1 - tau)^2
    ) /
    (
      tau^2 * (1 - tau)^2
    )
}



###############################################################
## Calibrate b from the CURRENT SNR definition
##
##   SNR = A / sqrt(Var(e))
##
## Therefore
##
##   b =
##   (A/SNR)
##   tau(1-tau)
##   / sqrt(tau^2 + (1-tau)^2).
###############################################################

ald_scale_from_snr <- function(
    A,
    SNR,
    tau
) {
  
  if (A <= 0) {
    stop("A must be positive.")
  }
  
  if (SNR <= 0) {
    stop("SNR must be positive.")
  }
  
  if (tau <= 0 || tau >= 1) {
    stop("tau must belong to (0,1).")
  }
  
  (
    A / SNR
  ) *
    (
      tau * (1 - tau)
    ) /
    sqrt(
      tau^2 + (1 - tau)^2
    )
}



###############################################################
## Main sequential estimation function
###############################################################

estimate_peak_curves <- function(
    A = 1,
    L = 22,
    S = 10,
    h = 30,
    SNR = 1000,
    tau = 0.5,
    seed = 280,
    s.interval = c(1,20),
    tol = 1e-4
)
{
  
  
  #############################################################
  ## Basic checks
  #############################################################
  
  if (A <= 0) {
    stop("A must be positive.")
  }
  
  if (S <= 0) {
    stop("S must be positive.")
  }
  
  if (h < 1 || h != floor(h)) {
    stop("h must be a positive integer.")
  }
  
  if (L < 1 || L > h) {
    stop("L must satisfy 1 <= L <= h.")
  }
  
  if (tau <= 0 || tau >= 1) {
    stop("tau must belong to (0,1).")
  }
  
  
  #############################################################
  ## True Gaussian signal
  #############################################################
  
  signal <- function(
    A,
    L,
    S,
    h
  ) {
    
    t <- 1:h
    
    A * exp(
      -(t - L)^2 / S^2
    )
  }
  
  
  #############################################################
  ## Generate tau-aligned AL noise
  #############################################################
  
  set.seed(seed)
  
  
  b <- ald_scale_from_snr(
    A = A,
    SNR = SNR,
    tau = tau
  )
  
  
  noise <- rald_tau(
    n = h,
    tau = tau,
    b = b
  )
  
  
  y <-
    signal(
      A = A,
      L = L,
      S = S,
      h = h
    ) +
    noise
  
  
  #############################################################
  ## Check the achieved noise variance
  #############################################################
  
  theoretical.var <-
    ald_variance(
      tau = tau,
      b = b
    )
  
  
  #############################################################
  ## Quantile-loss objective for fixed (l,s,T)
  ##
  ## For fixed l and s, the minimizer with respect to a is the
  ## weighted tau-quantile of
  ##
  ##   y_t exp((t-l)^2/s^2)
  ##
  ## with weights
  ##
  ##   exp(-(t-l)^2/s^2).
  #############################################################
  
  objective_quantile <- function(
    l,
    s,
    T
  ) {
    
    tt <- 1:T
    
    
    weights <-
      exp(
        -(tt - l)^2 / s^2
      )
    
    
    z <-
      y[1:T] *
      exp(
        (tt - l)^2 / s^2
      )
    
    
    a <-
      weighted_quantile(
        x = z,
        w = weights,
        tau = tau
      )
    
    
    fit <-
      a *
      exp(
        -(tt - l)^2 / s^2
      )
    
    
    residuals <-
      y[1:T] - fit
    
    
    value <-
      mean(
        rho_tau(
          residuals,
          tau = tau
        )
      )
    
    
    list(
      value = value,
      a = a
    )
  }
  
  
  #############################################################
  ## Estimate parameters at one observation horizon T
  ##
  ## IMPORTANT:
  ## The current manuscript restricts l to 1,...,T.
  #############################################################
  
  estimate_T <- function(
    T
  ) {
    
    out <-
      matrix(
        NA_real_,
        nrow = h,
        ncol = 3
      )
    
    colnames(out) <-
      c(
        "a",
        "s",
        "obj"
      )
    
    
    for (l in 1:h) {
      
      f <- function(s) {
        
        objective_quantile(
          l = l,
          s = s,
          T = T
        )$value
      }
      
      
      opt <-
        optimize(
          f,
          interval = s.interval,
          tol = tol
        )
      
      
      tmp <-
        objective_quantile(
          l = l,
          s = opt$minimum,
          T = T
        )
      
      
      out[l,1] <- tmp$a
      out[l,2] <- opt$minimum
      out[l,3] <- tmp$value
    }
    
    
    best <-
      which.min(
        out[,3]
      )
    
    
    list(
      a = out[best,1],
      s = out[best,2],
      l = best,
      objective = out[best,3]
    )
  }
  
  
  #############################################################
  ## Sequential horizons
  #############################################################
  
  TT <- 4:h
  
  
  hat.a <-
    hat.l <-
    hat.s <-
    hat.obj <-
    numeric(
      length(TT)
    )
  
  
  #############################################################
  ## Sequential estimation
  #############################################################
  
  for (i in seq_along(TT)) {
    
    est <-
      estimate_T(
        TT[i]
      )
    
    
    hat.a[i] <- est$a
    hat.l[i] <- est$l
    hat.s[i] <- est$s
    hat.obj[i] <- est$objective
  }
  
  
  #############################################################
  ## Future-tail oscillation
  ##
  ## omega_theta(T,h)
  ##   = max_{T <= t1,t2 <= h}
  ##       |hat(theta)(t1)-hat(theta)(t2)|
  #############################################################
  
  tail_oscillation <- function(
    x
  ) {
    
    n <-
      length(x)
    
    omega <-
      numeric(n)
    
    
    for (i in seq_len(n)) {
      
      omega[i] <-
        max(
          x[i:n]
        ) -
        min(
          x[i:n]
        )
    }
    
    
    omega
  }
  
  
  omega.a <-
    tail_oscillation(
      hat.a
    )
  
  
  omega.l <-
    tail_oscillation(
      hat.l
    )
  
  
  omega.s <-
    tail_oscillation(
      hat.s
    )
  
  
  #############################################################
  ## Return all results
  #############################################################
  
  list(
    
    T = TT,
    
    hat.a = hat.a,
    hat.l = hat.l,
    hat.s = hat.s,
    hat.obj = hat.obj,
    
    omega.a = omega.a,
    omega.l = omega.l,
    omega.s = omega.s,
    
    A = A,
    L = L,
    S = S,
    h = h,
    
    SNR = SNR,
    tau = tau,
    
    b = b,
    theoretical.var = theoretical.var
  )
}



###############################################################
## Helper: one row = estimate curve | oscillation curve
###############################################################

plot_pair <- function(
    T,
    hat,
    omega,
    true.val,
    L,
    ylab.hat,
    ylab.omega,
    main
)
{
  
  plot(
    T,
    hat,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = ylab.hat,
    main = main
  )
  
  abline(
    h = true.val,
    lty = 2
  )
  
  abline(
    v = L,
    lty = 2,
    col = "red"
  )
  
  
  plot(
    T,
    omega,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = ylab.omega,
    main = main
  )
  
  abline(
    v = L,
    lty = 2,
    col = "red"
  )
}



###############################################################
## GLOBAL PARAMETERS
###############################################################

A <- 1
L <- 22
S <- 10

## The symmetric case tau = 0.5 reproduces the old
## symmetric-Laplace setting, but through the new
## tau-aligned AL parameterization.
tau <- 0.5



###############################################################
## Section 2.4: large SNR
##
## SNR = 1000, h = 30
###############################################################

# postscript(
#   "sec24_large_snr.eps",
#   onefile = FALSE,
#   paper = "special",
#   horizontal = FALSE,
#   width = 8,
#   height = 7,
#   family = "Helvetica",
#   pointsize = 11,
#   colormodel = "srgb"
# )
# 
# par(
#   mfrow = c(3,2),
#   mar = c(4,4,3,1),
#   mgp = c(2.2,0.6,0)
# )

png(
  "sec24_large_snr.png",
  width = 8,
  height = 7,
  units = "in",
  res = 300
)

par(
  mfrow = c(3,2),
  mar = c(4,4,3,1),
  mgp = c(2.2,0.6,0)
)

res.large <-
  estimate_peak_curves(
    A = A,
    L = L,
    S = S,
    h = 30,
    SNR = 1000,
    tau = tau,
    seed = 280
  )

plot_pair(
  res.large$T,
  res.large$hat.a,
  res.large$omega.a,
  A,
  L,
  expression(hat(a)(T)),
  expression(omega[a](T)),
  paste0(
    "SNR = 1000, ",
    expression(tau == 0.5)
  )
)

plot_pair(
  res.large$T,
  res.large$hat.l,
  res.large$omega.l,
  L,
  L,
  expression(hat(l)(T)),
  expression(omega[l](T)),
  paste0(
    "SNR = 1000, ",
    expression(tau == 0.5)
  )
)

plot_pair(
  res.large$T,
  res.large$hat.s,
  res.large$omega.s,
  S,
  L,
  expression(hat(s)(T)),
  expression(omega[s](T)),
  paste0(
    "SNR = 1000, ",
    expression(tau == 0.5)
  )
)

dev.off()



###############################################################
## Section 2.4: medium SNR
##
## SNR = 100, h = 30
###############################################################

# postscript(
#   "sec24_medium_snr.eps",
#   onefile = FALSE,
#   paper = "special",
#   horizontal = FALSE,
#   width = 8,
#   height = 7,
#   family = "Helvetica",
#   pointsize = 11,
#   colormodel = "srgb"
# )
# 
# par(
#   mfrow = c(3,2),
#   mar = c(4,4,3,1),
#   mgp = c(2.2,0.6,0)
# )

png(
  "sec24_medium_snr.png",
  width = 8,
  height = 7,
  units = "in",
  res = 300
)

par(
  mfrow = c(3,2),
  mar = c(4,4,3,1),
  mgp = c(2.2,0.6,0)
)

res.medium <-
  estimate_peak_curves(
    A = A,
    L = L,
    S = S,
    h = 30,
    SNR = 100,
    tau = tau,
    seed = 280
  )

plot_pair(
  res.medium$T,
  res.medium$hat.a,
  res.medium$omega.a,
  A,
  L,
  expression(hat(a)(T)),
  expression(omega[a](T)),
  paste0(
    "SNR = 100, ",
    expression(tau == 0.5)
  )
)

plot_pair(
  res.medium$T,
  res.medium$hat.l,
  res.medium$omega.l,
  L,
  L,
  expression(hat(l)(T)),
  expression(omega[l](T)),
  paste0(
    "SNR = 100, ",
    expression(tau == 0.5)
  )
)

plot_pair(
  res.medium$T,
  res.medium$hat.s,
  res.medium$omega.s,
  S,
  L,
  expression(hat(s)(T)),
  expression(omega[s](T)),
  paste0(
    "SNR = 100, ",
    expression(tau == 0.5)
  )
)

dev.off()



###############################################################
## Section 2.5: small SNR
##
## SNR = 40, h = 30
###############################################################

# postscript(
#   "sec25_small_snr_h30.eps",
#   onefile = FALSE,
#   paper = "special",
#   horizontal = FALSE,
#   width = 8,
#   height = 7,
#   family = "Helvetica",
#   pointsize = 11,
#   colormodel = "srgb"
# )
# 
# par(
#   mfrow = c(3,2),
#   mar = c(4,4,3,1),
#   mgp = c(2.2,0.6,0)
# )

png(
  "sec25_small_snr_h30.png",
  width = 8,
  height = 7,
  units = "in",
  res = 300
)

par(
  mfrow = c(3,2),
  mar = c(4,4,3,1),
  mgp = c(2.2,0.6,0)
)

res.small <-
  estimate_peak_curves(
    A = A,
    L = L,
    S = S,
    h = 30,
    SNR = 40,
    tau = tau,
    seed = 280
  )

plot_pair(
  res.small$T,
  res.small$hat.a,
  res.small$omega.a,
  A,
  L,
  expression(hat(a)(T)),
  expression(omega[a](T)),
  paste0(
    "SNR = 40, ",
    expression(tau == 0.5)
  )
)

plot_pair(
  res.small$T,
  res.small$hat.l,
  res.small$omega.l,
  L,
  L,
  expression(hat(l)(T)),
  expression(omega[l](T)),
  paste0(
    "SNR = 40, ",
    expression(tau == 0.5)
  )
)

plot_pair(
  res.small$T,
  res.small$hat.s,
  res.small$omega.s,
  S,
  L,
  expression(hat(s)(T)),
  expression(omega[s](T)),
  paste0(
    "SNR = 40, ",
    expression(tau == 0.5)
  )
)

dev.off()



###############################################################
## Compute all results once
###############################################################

results <- list(
  
  large =
    res.large,
  
  medium =
    res.medium,
  
  small =
    res.small
)



###############################################################
## Display results
###############################################################

results
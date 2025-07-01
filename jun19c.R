# # ---------------------------------------------------------------------------
# # Figures 3 & 4 — Estimation Paths with Asymmetric Laplace Noise (z overlay)
# # Author: Said Maanan | Date: 2025-06-12
# # ---------------------------------------------------------------------------
# 
# # --- Utility Functions ------------------------------------------------------
# 
# asy_laplace_variance <- function(tau, z) {
#   stopifnot(tau > 0, tau < 1, abs(z) < 1)
#   mu <- tau / (1 + z) - (1 - tau) / (1 - z)
#   ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
#   ex2 - mu^2
# }
# 
# compute_noise_scale <- function(A, snr, tau, z) {
#   sqrt((A / snr)^2 / asy_laplace_variance(tau, z))
# }
# 
# r_asy_laplace <- function(n, tau, scale = 1) {
#   u <- runif(n)
#   ifelse(u < tau, -rexp(n, rate = 1 / scale), rexp(n, rate = 1 / scale))
# }
# 
# quantile_loss <- function(u, tau) {
#   u * (tau - (u < 0))
# }
# 
# med_weighted_tau <- function(I, l, s, tau) {
#   T <- length(I)
#   x <- I * exp(((1:T) - l)^2 / s^2)
#   w <- exp(-((1:T) - l)^2 / s^2)
#   o <- order(x)
#   cum_w <- cumsum(w[o])
#   cutoff <- tau * sum(w)
#   x[o][min(which(cum_w >= cutoff))]
# }
# 
# compute_a <- function(I, l, s, tau) {
#   med_weighted_tau(I, l, s, tau)
# }
# 
# objective_function <- function(I, a, l, s, tau) {
#   t <- seq_along(I)
#   fitted <- a * exp(-((t - l)^2) / s^2)
#   mean(quantile_loss(I - fitted, tau))
# }
# 
# search_best_s <- function(I, l, tau) {
#   f <- function(s) objective_function(I, compute_a(I, l, s, tau), l, s, tau)
#   optimize(f, c(1, 20), tol = 1e-4)$minimum
# }
# 
# estimate_params <- function(I, tau, h) {
#   losses <- sapply(1:h, function(l) {
#     s <- search_best_s(I, l, tau)
#     a <- compute_a(I, l, s, tau)
#     objective_function(I, a, l, s, tau)
#   })
#   l_hat <- which.min(losses)
#   s_hat <- search_best_s(I, l_hat, tau)
#   a_hat <- compute_a(I, l_hat, s_hat, tau)
#   list(l_hat = l_hat, s_hat = s_hat, a_hat = a_hat)
# }
# 
# # --- Simulation Settings ----------------------------------------------------
# 
# A_true <- 1
# L_true <- 22
# S_true <- 10
# snr <- 1000
# seed <- 280
# h <- 30
# T_vals <- 4:h
# taus <- c(0.25, 0.5, 0.75, 0.95)
# z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
# colors <- c("orange", "red", "black", "blue", "darkgreen")
# 
# est_paths <- list()
# 
# for (tau in taus) {
#   est_paths[[as.character(tau)]] <- list(
#     l_hat = list(),
#     a_hat = list(),
#     s_hat = list()
#   )
#   
#   for (z in z_vals) {
#     set.seed(seed)
#     signal <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
#     noise <- r_asy_laplace(h, tau, scale = compute_noise_scale(A_true, snr, tau, z))
#     I <- signal + noise
#     
#     path_l <- numeric(length(T_vals))
#     path_a <- numeric(length(T_vals))
#     path_s <- numeric(length(T_vals))
#     
#     for (i in seq_along(T_vals)) {
#       T <- T_vals[i]
#       p <- estimate_params(I[1:T], tau, h)
#       path_l[i] <- p$l_hat
#       path_a[i] <- p$a_hat
#       path_s[i] <- p$s_hat
#     }
#     
#     est_paths[[as.character(tau)]]$l_hat[[as.character(z)]] <- path_l
#     est_paths[[as.character(tau)]]$a_hat[[as.character(z)]] <- path_a
#     est_paths[[as.character(tau)]]$s_hat[[as.character(z)]] <- path_s
#   }
# }
# 
# get_ylim <- function(est_name) {
#   all_vals <- unlist(lapply(taus, function(tau)
#     unlist(est_paths[[as.character(tau)]][[est_name]])
#   ))
#   range(all_vals)
# }
# 
# expand_range <- function(x, margin = 0.05) {
#   r <- range(x)
#   delta <- diff(r) * margin
#   r + c(-delta, delta)
# }
# 
# ylim_l <- expand_range(get_ylim("l_hat"))
# ylim_a <- expand_range(get_ylim("a_hat"))
# ylim_s <- expand_range(get_ylim("s_hat"))
# 
# 
# pdf(file = "fig3_4_superposed_3x3.pdf", width = 8, height = 9)
# par(mfrow = c(4, 3), mar = c(2, 2, 1, 1), oma = c(3, 3, 2.5, 1), bg = "white")
# 
# for (tau in taus) {
#   for (j in seq_along(c("l_hat", "a_hat", "s_hat"))) {
#     est_name <- c("l_hat", "a_hat", "s_hat")[j]
#     
#     # Compute the appropriate ylim for this estimator
#     ylim <- switch(est_name,
#                    l_hat = ylim_l,
#                    a_hat = ylim_a,
#                    s_hat = ylim_s)
#     
#     plot(NA, xlim = range(T_vals), ylim = ylim,
#          xlab = "", ylab = "", axes = FALSE)
#     
#     abline(h = switch(est_name, l_hat = L_true, a_hat = A_true, s_hat = S_true),
#            col = "gray40", lty = 2)
#     grid()
#     box()
#     axis(1); axis(2)
#     
#     for (i in seq_along(z_vals)) {
#       z <- z_vals[i]
#       set.seed(seed)
#       signal <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
#       noise <- r_asy_laplace(h, tau, scale = compute_noise_scale(A_true, snr, tau, z))
#       I <- signal + noise
#       
#       est_path <- sapply(T_vals, function(T) {
#         p <- estimate_params(I[1:T], tau, h)
#         p[[est_name]]
#       })
#       
#       lines(T_vals, est_path, type = "b", col = colors[i], pch = 19)
#     }
#     
#     legend("topright", legend = paste("z =", z_vals),
#            col = colors, lwd = 2, pch = 19, cex = 0.75, bty = "o")
#     
#     if (tau == taus[1]) {
#       mtext(switch(est_name,
#                    l_hat = expression(hat(l)),
#                    a_hat = expression(hat(a)),
#                    s_hat = expression(hat(s))),
#             side = 3, line = 1.2)
#     }
#     
#     if (j == 1) {
#       mtext(bquote(tau == .(tau)), side = 2, line = 3)
#     }
#   }
# }
# 
# dev.off()
# cat("Cleaned superposed 3x3 plot saved as 'fig3_4_superposed_3x3.pdf'\n")


# ---------------------------------------------------------------------------
# Figures 3 & 4 — Estimation Paths with Asymmetric Laplace Noise (z overlay) - UPDATED
# ---------------------------------------------------------------------------

# --- Utility Functions ------------------------------------------------------

asy_laplace_variance <- function(tau, z) {
  stopifnot(tau > 0, tau < 1, abs(z) < 1)
  mu <- tau / (1 + z) - (1 - tau) / (1 - z)
  ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
  ex2 - mu^2
}

compute_noise_scale <- function(A, snr, tau, z) {
  sqrt((A / snr)^2 / asy_laplace_variance(tau, z))
}

# UPDATED noise generator: includes asymmetry parameter z and scale b
r_asy_laplace <- function(n, tau, z, b = 1) {
  u <- runif(n)
  ifelse(u < tau,
         -rexp(n, rate = (1 + z) / b),
         rexp(n, rate = (1 - z) / b))
}

quantile_loss <- function(u, tau) {
  u * (tau - (u < 0))
}

med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  o <- order(x)
  cum_w <- cumsum(w[o])
  cutoff <- tau * sum(w)
  x[o][min(which(cum_w >= cutoff))]
}

compute_a <- function(I, l, s, tau) {
  med_weighted_tau(I, l, s, tau)
}

objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  fitted <- a * exp(-((t - l)^2) / s^2)
  mean(quantile_loss(I - fitted, tau))
}

search_best_s <- function(I, l, tau) {
  f <- function(s) objective_function(I, compute_a(I, l, s, tau), l, s, tau)
  optimize(f, c(1, 20), tol = 1e-4)$minimum
}

estimate_params <- function(I, tau, h) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau)
    a <- compute_a(I, l, s, tau)
    objective_function(I, a, l, s, tau)
  })
  l_hat <- which.min(losses)
  s_hat <- search_best_s(I, l_hat, tau)
  a_hat <- compute_a(I, l_hat, s_hat, tau)
  list(l_hat = l_hat, s_hat = s_hat, a_hat = a_hat)
}

# --- Simulation Settings ----------------------------------------------------

A_true <- 1
L_true <- 22
S_true <- 10
snr <- 1000
seed <- 280
h <- 30
T_vals <- 4:h
taus <- c(0.25, 0.5, 0.75, 0.95)
z_vals <- c(-0.9, -0.5, 0, 0.5, 0.9)
colors <- c("orange", "red", "black", "blue", "darkgreen")

est_paths <- list()

for (tau in taus) {
  est_paths[[as.character(tau)]] <- list(
    l_hat = list(),
    a_hat = list(),
    s_hat = list()
  )
  
  for (z in z_vals) {
    set.seed(seed)
    signal <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
    b <- compute_noise_scale(A_true, snr, tau, z)  # compute scale b for each (tau, z)
    noise <- r_asy_laplace(h, tau, z, b)          # pass z and b here
    I <- signal + noise
    
    path_l <- numeric(length(T_vals))
    path_a <- numeric(length(T_vals))
    path_s <- numeric(length(T_vals))
    
    for (i in seq_along(T_vals)) {
      T <- T_vals[i]
      p <- estimate_params(I[1:T], tau, h)
      path_l[i] <- p$l_hat
      path_a[i] <- p$a_hat
      path_s[i] <- p$s_hat
    }
    
    est_paths[[as.character(tau)]]$l_hat[[as.character(z)]] <- path_l
    est_paths[[as.character(tau)]]$a_hat[[as.character(z)]] <- path_a
    est_paths[[as.character(tau)]]$s_hat[[as.character(z)]] <- path_s
  }
}

get_ylim <- function(est_name) {
  all_vals <- unlist(lapply(taus, function(tau)
    unlist(est_paths[[as.character(tau)]][[est_name]])
  ))
  range(all_vals)
}

expand_range <- function(x, margin = 0.05) {
  r <- range(x)
  delta <- diff(r) * margin
  r + c(-delta, delta)
}

ylim_l <- expand_range(get_ylim("l_hat"))
ylim_a <- expand_range(get_ylim("a_hat"))
ylim_s <- expand_range(get_ylim("s_hat"))


pdf(file = "fig3_4_superposed_3x3.pdf", width = 8, height = 9)
par(mfrow = c(4, 3), mar = c(2, 2, 1, 1), oma = c(3, 3, 2.5, 1), bg = "white")

for (tau in taus) {
  for (j in seq_along(c("l_hat", "a_hat", "s_hat"))) {
    est_name <- c("l_hat", "a_hat", "s_hat")[j]
    
    # Compute the appropriate ylim for this estimator
    ylim <- switch(est_name,
                   l_hat = ylim_l,
                   a_hat = ylim_a,
                   s_hat = ylim_s)
    
    plot(NA, xlim = range(T_vals), ylim = ylim,
         xlab = "", ylab = "", axes = FALSE)
    
    abline(h = switch(est_name, l_hat = L_true, a_hat = A_true, s_hat = S_true),
           col = "gray40", lty = 2)
    grid()
    box()
    axis(1); axis(2)
    
    for (i in seq_along(z_vals)) {
      z <- z_vals[i]
      set.seed(seed)
      signal <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
      b <- compute_noise_scale(A_true, snr, tau, z)
      noise <- r_asy_laplace(h, tau, z, b)
      I <- signal + noise
      
      est_path <- sapply(T_vals, function(T) {
        p <- estimate_params(I[1:T], tau, h)
        p[[est_name]]
      })
      
      lines(T_vals, est_path, type = "b", col = colors[i], pch = 19)
    }
    
    legend("topright", legend = paste("z =", z_vals),
           col = colors, lwd = 2, pch = 19, cex = 0.75, bty = "o")
    
    if (tau == taus[1]) {
      mtext(switch(est_name,
                   l_hat = expression(hat(l)),
                   a_hat = expression(hat(a)),
                   s_hat = expression(hat(s))),
            side = 3, line = 1.2)
    }
    
    if (j == 1) {
      mtext(bquote(tau == .(tau)), side = 2, line = 3)
    }
  }
}

dev.off()
cat("Cleaned superposed 3x3 plot saved as 'fig3_4_superposed_3x3.pdf'\n")

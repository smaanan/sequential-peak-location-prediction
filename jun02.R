# -------------------------------------------------------
# Algorithm 2 + Confidence Intervals for log(I(T))
# COVID-19 Data - China, T = 14, 18, 21
# Author: [Your Name], Date: [Current Date]
# -------------------------------------------------------

# --- Libraries ---
library(dplyr)
library(quantreg)
library(COVID19)
library(pracma)
library(gridExtra)
library(grid)
library(ggplot2)

# --- Load and preprocess COVID-19 data for China ---
china <- covid19("China", start = "2020-01-22", end = "2020-04-01", level = 1, verbose = FALSE) %>%
  arrange(date) %>%
  mutate(new_confirmed = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_confirmed)) %>%
  mutate(log_cases = log(pmax(new_confirmed, 1)),
         t = row_number())

# --- Weighted median function ---
wmedian <- function(x, w) {
  valid <- !is.na(x) & !is.na(w)
  x <- x[valid]
  w <- w[valid]
  o <- order(x)
  x <- x[o]; w <- w[o] / sum(w)
  cw <- cumsum(w)
  x[which(cw >= 0.5)[1]]
}

# --- Algorithm 2: Estimate (l_hat, s_hat, a_hat) ---
estimate_peak_params <- function(log_I, T, h = 40) {
  objective <- function(s, l) {
    weights <- exp(-((1:T) - l)^2 / s^2)
    z <- log_I[1:T] + ((1:T) - l)^2 / s^2
    log_a <- wmedian(z, weights)
    sum(abs(log_I[1:T] - log_a + ((1:T) - l)^2 / s^2))
  }
  
  best <- list(val = Inf)
  for (l in 1:h) {
    opt <- optimize(function(s) objective(s, l), interval = c(1, 20))
    s_opt <- opt$minimum
    val <- opt$objective
    weights <- exp(-((1:T) - l)^2 / s_opt^2)
    z <- log_I[1:T] + ((1:T) - l)^2 / s_opt^2
    log_a <- wmedian(z, weights)
    if (val < best$val) best <- list(val = val, l = l, s = s_opt, log_a = log_a)
  }
  
  list(l_hat = best$l, s_hat = best$s, a_hat = exp(best$log_a))
}

# --- Delta-method-based CI for log(I(T)) ---
ci_log_I_T <- function(log_I, T) {
  df <- data.frame(log_cases = log_I[1:T], t = 1:T, t_sq = (1:T)^2)
  model <- rq(log_cases ~ t + t_sq, tau = 0.5, data = df)
  beta <- coef(model); names(beta) <- c("b0", "b1", "b2")
  residuals <- model$residuals
  f0 <- density(residuals, from = 0, to = 0, n = 1, kernel = "gaussian")$y
  omega2 <- (0.5 * 0.5) / (f0^2)
  
  t_vec <- 1:T
  scaling <- c(sqrt(T), sqrt(sum(t_vec^2)), sqrt(sum(t_vec^4)))
  H_inv <- solve(hilbert(3))
  Sigma_beta <- omega2 * diag(1 / scaling) %*% H_inv %*% diag(1 / scaling)
  
  grad <- c(1, T, T^2)
  var_log_I_T <- as.numeric(t(grad) %*% Sigma_beta %*% grad)
  se <- sqrt(var_log_I_T)
  pred <- sum(beta * grad)
  
  list(pred = pred, CI = pred + c(-1.96, 1.96) * se, se = se)
}

# --- Run Algorithm 2 and CI estimation for T = 14, 18, 21 ---
log_I <- china$log_cases
T_vals <- c(14, 18, 21)
results <- lapply(T_vals, function(T_val) {
  params <- estimate_peak_params(log_I, T_val)
  ci <- ci_log_I_T(log_I, T_val)
  list(T = T_val, l_hat = params$l_hat, s_hat = params$s_hat,
       a_hat = params$a_hat, log_I_pred = ci$pred,
       CI_log_I = ci$CI, se_log_I = ci$se)
})

# --- Tabulate and export results ---
results_df <- do.call(rbind, lapply(results, function(res) {
  data.frame(
    T = res$T,
    l_hat = sprintf("%.2f", res$l_hat),
    s_hat = sprintf("%.2f", res$s_hat),
    a_hat = sprintf("%.2f", res$a_hat),
    log_I_actual = sprintf("%.3f", log_I[res$T]),
    log_I_pred = sprintf("%.3f", res$log_I_pred),
    CI_lower = sprintf("%.3f", res$CI_log_I[1]),
    CI_upper = sprintf("%.3f", res$CI_log_I[2])
  )
}))
table_grob <- tableGrob(results_df, rows = NULL)
png("algorithm2_china_confidence_table.png", width = 800, height = 150, res = 120)
grid.draw(table_grob)
dev.off()

# --- Predict log-median trajectory ---
predict_log_median_curve <- function(T_val, l_hat, s_hat, a_hat) {
  t_vals <- T_val:35
  log_median <- log(a_hat) - ((t_vals - l_hat)^2) / s_hat^2
  data.frame(T = T_val, t = t_vals, log_median = log_median)
}
curves <- do.call(rbind, lapply(results, function(res) {
  predict_log_median_curve(res$T, res$l_hat, res$s_hat, res$a_hat)
}))

# --- Function to generate LaTeX table for predicted curves ---
latex_rowwise_curve_table_split <- function(df) {
  t_vals <- df$t
  log_vals <- sprintf("%.2f", df$log_median)
  split <- ceiling(length(t_vals) / 2)
  make_tab <- function(ts, vs) {
    paste0(
      "\\begin{tabular}{|", paste(rep("r", length(ts) + 1), collapse = ""), "|}\n\\hline\n",
      "$t$ & ", paste(ts, collapse = " & "), " \\\\\n\\hline\n",
      "$\\log(\\hat{a}) - \\frac{(t - \\hat{\\ell})^2}{\\hat{s}^2}$ & ", paste(vs, collapse = " & "), " \\\\\n\\hline\n\\end{tabular}\n"
    )
  }
  paste0(
    "\\begin{table}[H]\n\\centering\n",
    make_tab(t_vals[1:split], log_vals[1:split]), "\n\\vspace{0.4cm}\n",
    make_tab(t_vals[(split + 1):length(t_vals)], log_vals[(split + 1):length(t_vals)]), "\n",
    sprintf("\\caption{Predicted log-median values from $t = %d$ to $35$ using Algorithm~\\ref{algo2}.}\n", df$T[1]),
    sprintf("\\label{tab:log_median_T%d_rowwise}\n", df$T[1]),
    "\\end{table}"
  )
}

# --- Export LaTeX tables ---
cat(latex_rowwise_curve_table_split(subset(curves, T == 14)), "\n\n")
cat(latex_rowwise_curve_table_split(subset(curves, T == 18)), "\n\n")
cat(latex_rowwise_curve_table_split(subset(curves, T == 21)), "\n\n")

# --- Delta-method-based CI at arbitrary t with optional adjustment ---
ci_log_I_at <- function(log_I, T_fit, t_eval, adjust = c("none", "max", "range")) {
  adjust <- match.arg(adjust)
  
  df <- data.frame(log_cases = log_I[1:T_fit], t = 1:T_fit, t_sq = (1:T_fit)^2)
  model <- rq(log_cases ~ t + t_sq, tau = 0.5, data = df)
  beta <- coef(model); names(beta) <- c("b0", "b1", "b2")
  residuals <- model$residuals
  
  f0 <- density(residuals, from = 0, to = 0, n = 1)$y
  omega2 <- 0.25 / f0^2
  
  t_vec <- 1:T_fit
  scaling <- c(sqrt(T_fit), sqrt(sum(t_vec^2)), sqrt(sum(t_vec^4)))
  H_inv <- solve(hilbert(3))
  Sigma <- omega2 * diag(1 / scaling) %*% H_inv %*% diag(1 / scaling)
  
  grad <- c(1, t_eval, t_eval^2)
  pred <- sum(beta * grad)
  se <- sqrt(as.numeric(t(grad) %*% Sigma %*% grad))
  CI <- pred + c(-1.96, 1.96) * se
  
  if (adjust == "max") {
    CI[2] <- CI[2] + max(residuals)
  } else if (adjust == "range") {
    CI[2] <- CI[2] + (max(residuals) - min(residuals))
  }
  
  list(
    T_fit = T_fit,
    t_eval = t_eval,
    pred = pred,
    CI = CI,
    inside = log_I[t_eval] >= CI[1] && log_I[t_eval] <= CI[2],
    se = se
  )
}

# Evaluate prediction at t = 22 using chosen adjustment ("none", "max", or "range")
adjust_method <- "none"  # or "max" or "range"
log_I_obs_22 <- log_I[22]

preds_at_22 <- lapply(T_vals, function(T_fit) {
  ci_log_I_at(log_I, T_fit, 22, adjust = adjust_method)
})

# Assemble result table
preds_df <- do.call(rbind, lapply(preds_at_22, function(res) {
  data.frame(
    T = res$T_fit,
    t = res$t_eval,
    pred_log_I = sprintf("%.3f", res$pred),
    CI_lower = sprintf("%.3f", res$CI[1]),
    CI_upper = sprintf("%.3f", res$CI[2]),
    log_I_obs = sprintf("%.3f", log_I_obs_22),
    inside_CI = ifelse(res$inside, "Yes", "No")
  )
}))

print(preds_df)

# Optional: export as PNG
table_grob <- tableGrob(preds_df, rows = NULL)
png("table_logI22_prediction.png", width = 1000, height = 200, res = 120)
grid.draw(table_grob)
dev.off()


# Function to compute log differences
compute_log_diff <- function(df, Tmin) {
  df <- subset(df, t >= Tmin & t <= 35)
  diffs <- diff(df$log_median)
  t_vals <- df$t[-1]
  return(data.frame(t = t_vals, diff_log = diffs))
}

# Compute all three differences
diff_14 <- compute_log_diff(subset(curves, T == 14), 14)
diff_18 <- compute_log_diff(subset(curves, T == 18), 18)
diff_21 <- compute_log_diff(subset(curves, T == 21), 21)

# Determine the full t range for all curves
t_min <- min(diff_14$t)
t_max <- max(diff_14$t)
y_min <- min(c(diff_14$diff_log, diff_18$diff_log, diff_21$diff_log))
y_max <- max(c(diff_14$diff_log, diff_18$diff_log, diff_21$diff_log))

png("log_diff_unified_plot.png", width = 900, height = 600, res = 120)

# Base plot with updated expression
plot(diff_14$t, diff_14$diff_log, type = "b", pch = 16, col = "blue",
     xlim = c(t_min, t_max), ylim = c(y_min, y_max),
     xlab = "t",
     ylab = expression(log(hat(I)[T~","~t+1]) - log(hat(I)[T~","~t])),
     main = expression(paste("Predicted differences ", log(hat(I)[T~","~t+1]) - log(hat(I)[T~","~t]))))

# Add lines for other T values
lines(diff_18$t, diff_18$diff_log, col = "red", pch = 16, type = "b")
lines(diff_21$t, diff_21$diff_log, col = "grey60", pch = 16, type = "b")

# Reference line and grid
abline(h = 0, lty = 2, col = "gray")
grid()

# Legend
legend("topright", legend = c("14", "18", "21"),
       col = c("blue", "red", "grey60"),
       pch = 16, lty = 1,
       title = expression("Fitted at T ="),
       bty = "o")

dev.off()

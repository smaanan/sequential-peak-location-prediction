# -------------------------------------------------------
# Algorithm 2 + Confidence Intervals for log(I(T))
# COVID-19 Data - China, T = 14, 18, 21
# -------------------------------------------------------

# --- Libraries ---
library(dplyr)
library(quantreg)
library(COVID19)
library(pracma)

# --- Step 1: Load and prepare COVID-19 data for China ---
china <- covid19("China", start = "2020-01-22", end = "2020-04-01", level = 1, verbose = FALSE) %>%
  arrange(date) %>%
  mutate(new_confirmed = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_confirmed)) %>%
  mutate(log_cases = log(pmax(new_confirmed, 1)),
         t = row_number())

# --- Define weighted median function ---
wmedian <- function(x, w) {
  valid <- !is.na(x) & !is.na(w)
  x <- x[valid]
  w <- w[valid]
  o <- order(x)
  x <- x[o]
  w <- w[o]
  w <- w / sum(w)
  cw <- cumsum(w)
  idx <- which(cw >= 0.5)[1]
  return(x[idx])
}

# --- Step 2: Algorithm 2 to estimate l_hat, s_hat, a_hat ---
estimate_peak_params <- function(log_I, T, h = 40) {
  objective <- function(s, l) {
    weights <- exp(-((1:T) - l)^2 / s^2)
    z <- log_I[1:T] + ((1:T) - l)^2 / s^2
    log_a <- wmedian(z, weights)
    sum(abs(log_I[1:T] - log_a + ((1:T) - l)^2 / s^2))
  }

  # Grid over l and optimize s
  best_val <- Inf
  best_l <- NA
  best_s <- NA
  best_log_a <- NA

  for (l in 1:h) {
    opt <- optimize(f = function(s) objective(s, l), interval = c(1, 20))
    s_opt <- opt$minimum
    val <- opt$objective
    weights <- exp(-((1:T) - l)^2 / s_opt^2)
    z <- log_I[1:T] + ((1:T) - l)^2 / s_opt^2
    log_a <- wmedian(z, weights)

    if (val < best_val) {
      best_val <- val
      best_l <- l
      best_s <- s_opt
      best_log_a <- log_a
    }
  }

  return(list(l_hat = best_l, s_hat = best_s, a_hat = exp(best_log_a)))
}

# --- Step 3: Compute 95% CI for log(I(T)) using quantile regression ---
ci_log_I_T <- function(log_I, T) {
  df <- data.frame(log_cases = log_I[1:T], t = 1:T, t_sq = (1:T)^2)
  model <- rq(log_cases ~ t + t_sq, tau = 0.5, data = df)
  beta <- coef(model)
  names(beta) <- c("b0", "b1", "b2")

  # Estimate f(0) from residuals
  residuals <- model$residuals
  f0 <- density(residuals, from = 0, to = 0, n = 1, kernel = "gaussian")$y
  omega2 <- (0.5 * 0.5) / (f0^2)

  # Covariance estimation (Hilbert-based)
  t_vec <- 1:T
  scaling <- c(sqrt(T), sqrt(sum(t_vec^2)), sqrt(sum(t_vec^4)))
  H_inv <- solve(hilbert(3))
  Sigma_beta <- omega2 * diag(1 / scaling) %*% H_inv %*% diag(1 / scaling)

  # Delta method for Var(log(I(T)))
  grad <- c(1, T, T^2)
  var_log_I_T <- as.numeric(t(grad) %*% Sigma_beta %*% grad)
  se_log_I_T <- sqrt(var_log_I_T)
  pred_log_I_T <- as.numeric(beta["b0"] + beta["b1"] * T + beta["b2"] * T^2)

  CI <- pred_log_I_T + c(-1.96, 1.96) * se_log_I_T
  return(list(pred = pred_log_I_T, CI = CI, se = se_log_I_T))
}

# --- Step 4: Run for T = 14, 18, 21 ---
results <- lapply(c(14, 18, 21), function(T_val) {
  log_I <- china$log_cases
  params <- estimate_peak_params(log_I, T = T_val)
  ci <- ci_log_I_T(log_I, T = T_val)

  list(
    T = T_val,
    l_hat = params$l_hat,
    s_hat = params$s_hat,
    a_hat = params$a_hat,
    log_I_pred = ci$pred,
    CI_log_I = ci$CI,
    se_log_I = ci$se
  )
})

# --- Step 5: Print results ---
cat("=== Results for Algorithm 2 and Confidence Intervals ===\\n")
for (res in results) {
  cat(sprintf("T = %d\\n", res$T))
  cat(sprintf("  l_hat = %.2f, s_hat = %.2f, a_hat = %.2f\\n", res$l_hat, res$s_hat, res$a_hat))
  cat(sprintf("  log(I(T)) = %.3f, 95%% CI = [%.3f, %.3f], SE = %.3f\\n\\n",
              res$log_I_pred, res$CI_log_I[1], res$CI_log_I[2], res$se_log_I))
}


library(gridExtra)
library(grid)
library(ggplot2)


# Build results data frame from the 'results' list
results_df <- do.call(rbind, lapply(results, function(res) {
  data.frame(
    T = res$T,
    l_hat = res$l_hat,
    s_hat = res$s_hat,
    a_hat = res$a_hat,
    log_I_actual = china$log_cases[res$T],  # Actual log value at T
    log_I_pred = res$log_I_pred,
    CI_lower = res$CI_log_I[1],
    CI_upper = res$CI_log_I[2]
  )
}))

# Format numbers (2 decimals for parameters, 3 for log values)
formatted_df <- transform(
  results_df,
  l_hat = sprintf("%.2f", l_hat),
  s_hat = sprintf("%.2f", s_hat),
  a_hat = sprintf("%.2f", a_hat),
  log_I_actual = sprintf("%.3f", log_I_actual),
  log_I_pred = sprintf("%.3f", log_I_pred),
  CI_lower = sprintf("%.3f", CI_lower),
  CI_upper = sprintf("%.3f", CI_upper)
)

# Generate and save table
table_grob <- tableGrob(formatted_df, rows = NULL)
png("algorithm2_china_confidence_table.png", width = 800, height = 150, res = 120)
grid.draw(table_grob)
dev.off()


# Function to generate predicted log-median curves from T to 35
predict_log_median_curve <- function(T_val, l_hat, s_hat, a_hat) {
  t_vals <- T_val:35
  log_a <- log(a_hat)
  pred_log_median <- log_a - ((t_vals - l_hat)^2) / s_hat^2
  data.frame(T = T_val, t = t_vals, log_median = pred_log_median)
}

# Use the fitted results
curves <- do.call(rbind, lapply(results, function(res) {
  predict_log_median_curve(res$T, res$l_hat, res$s_hat, res$a_hat)
}))


# Split curves by T value
curve_14 <- subset(curves, T == 14)
curve_18 <- subset(curves, T == 18)
curve_21 <- subset(curves, T == 21)


# Split LaTeX row-wise table into two halves stacked vertically
latex_rowwise_curve_table_split <- function(df) {
  t_vals <- df$t
  log_vals <- sprintf("%.2f", df$log_median)
  
  n <- length(t_vals)
  split_idx <- ceiling(n / 2)
  
  # First block
  t_row_1 <- paste0("$t$ & ", paste(t_vals[1:split_idx], collapse = " & "), " \\\\")
  log_row_1 <- paste0("$\\log(\\hat{a}) - \\frac{(t - \\hat{\\ell})^2}{\\hat{s}^2}$ & ",
                      paste(log_vals[1:split_idx], collapse = " & "), " \\\\")
  
  # Second block
  t_row_2 <- paste0("$t$ & ", paste(t_vals[(split_idx + 1):n], collapse = " & "), " \\\\")
  log_row_2 <- paste0("$\\log(\\hat{a}) - \\frac{(t - \\hat{\\ell})^2}{\\hat{s}^2}$ & ",
                      paste(log_vals[(split_idx + 1):n], collapse = " & "), " \\\\")
  
  n_cols_1 <- split_idx + 1
  n_cols_2 <- n - split_idx + 1
  
  paste0(
    "\\begin{table}[H]\n\\centering\n",
    "\\begin{tabular}{|", paste(rep("r", n_cols_1), collapse = ""), "|}\n\\hline\n",
    t_row_1, "\n\\hline\n",
    log_row_1, "\n\\hline\n",
    "\\end{tabular}\n\n\\vspace{0.4cm}\n\n",
    "\\begin{tabular}{|", paste(rep("r", n_cols_2), collapse = ""), "|}\n\\hline\n",
    t_row_2, "\n\\hline\n",
    log_row_2, "\n\\hline\n",
    "\\end{tabular}\n\n",
    sprintf("\\caption{Predicted log-median values from $t = %d$ to $35$ using Algorithm~\\ref{algo2}.}\n", df$T[1]),
    sprintf("\\label{tab:log_median_T%d_rowwise}\n", df$T[1]),
    "\\end{table}"
  )
}


cat(latex_rowwise_curve_table_split(curve_14), "\n\n")
cat(latex_rowwise_curve_table_split(curve_18), "\n\n")
cat(latex_rowwise_curve_table_split(curve_21), "\n")

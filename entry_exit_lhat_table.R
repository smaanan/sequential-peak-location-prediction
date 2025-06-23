# -----------------------------------------------------------------------------
# Entry and Exit Time into Peak Region — Table-Style Output (Like Table 2)
# -----------------------------------------------------------------------------

# --- Noise Generator and Helpers --------------------------------------------
rgen_asy_laplace <- function(n, tau, z, b = 1) {
  u <- runif(n)
  ifelse(u < tau,
         -rexp(n, rate = (1 + z) / b),
         rexp(n, rate = (1 - z) / b))
}

asy_laplace_variance <- function(tau, z) {
  mu <- tau / (1 + z) - (1 - tau) / (1 - z)
  ex2 <- 2 * tau / (1 + z)^2 + 2 * (1 - tau) / (1 - z)^2
  ex2 - mu^2
}

compute_noise_scale <- function(A, snr, tau, z) {
  sqrt((A / snr)^2 / asy_laplace_variance(tau, z))
}

quantile_loss <- function(u, tau) u * (tau - (u < 0))

med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  o <- order(x)
  cum_w <- cumsum(w[o])
  cutoff <- tau * sum(w)
  x[o][min(which(cum_w >= cutoff))]
}

compute_a <- function(I, l, s, tau) med_weighted_tau(I, l, s, tau)

objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  mean(quantile_loss(I - a * exp(-((t - l)^2) / s^2), tau))
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

# --- Settings ---------------------------------------------------------------
A_true <- 1
L_true <- 22
S_true <- 10
snr <- 40
seed <- 280
h <- 35
T_vals <- 4:h
tau <- 0.5
z <- 0.0

threshold <- 4
threshold_prime <- 0

# --- Simulate l_hat(T) path -------------------------------------------------
set.seed(seed)
signal <- A_true * exp(-((1:h) - L_true)^2 / S_true^2)
b <- compute_noise_scale(A_true, snr, tau, z)
noise <- rgen_asy_laplace(h, tau, z, b)
I <- signal + noise
l_path <- sapply(T_vals, function(T) estimate_params(I[1:T], tau, h)$l_hat)

# --- Print table-like layout ------------------------------------------------
cat("T:
", paste(T_vals[1:16], collapse = " "), "
")
cat("l_hat:
", paste(l_path[1:16], collapse = " "), "

")
cat("T:
", paste(T_vals[17:32], collapse = " "), "
")
cat("l_hat:
", paste(l_path[17:32], collapse = " "), "

")

# --- Compute T_entry and T_exit --------------------------------------------
find_stable_time <- function(l_path, T_vals, tol) {
  for (i in seq_along(T_vals)) {
    future_vals <- l_path[i:length(l_path)]
    if (max(future_vals) - min(future_vals) <= tol) {
      return(T_vals[i])
    }
  }
  return(NA)
}

T_entry <- find_stable_time(l_path, T_vals, threshold)
T_exit  <- find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)

cat(sprintf("T_entry (|Δl̂| ≤ %d): %d
", threshold, T_entry))
cat(sprintf("T_exit  (|Δl̂| ≤ %d): %d
", threshold_prime, T_exit))




# --- Required Libraries ---
library(dplyr)
library(gt)
library(tibble)

# --- Prepare Table Data ---
table_data <- tibble(
  T_1 = T_vals[1:16],
  lhat_1 = l_path[1:16],
  T_2 = T_vals[17:32],
  lhat_2 = l_path[17:32]
)

# --- Create gt Table ---
gt_table <- table_data %>%
  gt() %>%
  tab_header(
    title = md("**Table: Sequence T →  l̂(T)**"),
    subtitle = md("Simulated under SNR = 40, with τ = 0.5, z = 0")
  ) %>%
  cols_label(
    T_1 = md("**T**"), lhat_1 = md("**l̂(T)**"),
    T_2 = md("**T**"), lhat_2 = md("**l̂(T)**")
  ) %>%
  fmt_number(columns = everything(), decimals = 0) %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )

# --- Export Table as PNG ---
gtsave(gt_table, "table8_like_lhat_path.png")
cat("Saved as table8_like_lhat_path.png\n")

# --- Required Libraries ---
library(dplyr)
library(gt)
library(tibble)
library(glue)  # for cleaner string interpolation

# --- Transpose and Prepare Table ---
transposed_matrix <- t(table_data)

# Dummy column names just to avoid tibble error
colnames(transposed_matrix) <- paste0("V", seq_len(ncol(transposed_matrix)))

# Convert to tibble with row names
transposed_data <- as.data.frame(transposed_matrix) %>%
  rownames_to_column(var = "Label") %>%
  as_tibble()

# --- Rename Labels T_1/T_2 → "T", lhat_1/lhat_2 → "l̂(T)" ---
transposed_data$Label <- recode(transposed_data$Label,
                                "T_1" = "T", "T_2" = "T",
                                "lhat_1" = "l̂(T)", "lhat_2" = "l̂(T)")

# --- Format subtitle with dynamic entry/exit info ---
subtitle_text <- glue(
  "Simulated under SNR = {snr}, with τ = {tau}, z = {z}. ",
  "Entry at **T_entry = {T_entry}** (ε = {threshold}), ",
  "exit at **T_exit = {T_exit}** (ε' = {threshold_prime})."
)

# --- Create gt Table ---
gt_table <- transposed_data %>%
  gt() %>%
  tab_header(
    title = md("**Table: Sequence T →  l̂(T)**"),
    subtitle = md(subtitle_text)
  ) %>%
  cols_label(.list = setNames(rep("", ncol(transposed_data)), names(transposed_data))) %>%
  tab_options(
    column_labels.hidden = TRUE,
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )

# --- Save Table as PNG ---
gtsave(gt_table, "table8_transposed_lhat_path_with_entry_exit.png")
cat("Saved as table8_transposed_lhat_path_with_entry_exit.png\n")

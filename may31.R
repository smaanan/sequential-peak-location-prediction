# -------------------------------------------
# Forecasting Peak Location L(T) from COVID-19 Data (China)
# Reproduce Table 8 with h = 60
# -------------------------------------------

# Load libraries
library(dplyr)
library(zoo)
library(quantreg)
library(pracma)
library(COVID19)

# -------------------------------------------
# Step 1: Prepare COVID-19 data for China
# -------------------------------------------
china <- covid19("China", start = "2020-01-22", end = "2020-04-01", level = 1, verbose = FALSE) %>%
  arrange(date) %>%
  mutate(new_confirmed = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_confirmed)) %>%
  mutate(
    log_cases = log(pmax(new_confirmed, 1)),
    t = row_number(),
    t_sq = t^2
  )

# -------------------------------------------
# Step 2: Estimate L̂(T) for T from 4 to 60
# -------------------------------------------
results <- data.frame(T = integer(), L_hat = numeric())

for (T_val in 4:60) {
  df_T <- filter(china, t <= T_val)
  model <- rq(log_cases ~ t + t_sq, tau = 0.95, data = df_T, method = "fn")
  beta <- coef(model)
  b1 <- beta["t"]
  b2 <- beta["t_sq"]
  # Only compute L̂ if the parabola is concave
  L_hat <- if (b2 >= 0) NA else -b1 / (2 * b2)
  results <- rbind(results, data.frame(T = T_val, L_hat = round(L_hat, 2)))
}

# -------------------------------------------
# Step 3: Define LaTeX Table Builder
# -------------------------------------------
build_latex_table_block <- function(T_vals, L_vals) {
  T_line <- paste0("$T$ & ", paste(T_vals, collapse = " & "), " \\\\")
  L_line <- paste0("$\\hat{L}$ & ", paste(ifelse(is.na(L_vals), "--", L_vals), collapse = " & "), " \\\\")
  return(paste(
    "\\begin{tabular}{|", paste(rep("r", length(T_vals) + 1), collapse = ""), "|}",
    "\\hline", T_line, "\\hline", L_line, "\\\\ \\hline",
    "\\end{tabular}",
    sep = "\n"
  ))
}

# -------------------------------------------
# Step 4: Format into 3 LaTeX Table Blocks
# -------------------------------------------
block1 <- results[1:19, ]
block2 <- results[20:38, ]
block3 <- results[39:57, ]

# -------------------------------------------
# Step 5: Output Final LaTeX Table
# -------------------------------------------
cat("%% LaTeX Table for Estimated L(T), h = 60\n\n")
cat("\\begin{table}[H]\n\\centering\n")

cat(build_latex_table_block(block1$T, block1$L_hat), "\n\n\\vspace{0.4cm}\n\n")
cat(build_latex_table_block(block2$T, block2$L_hat), "\n\n\\vspace{0.4cm}\n\n")
cat(build_latex_table_block(block3$T, block3$L_hat), "\n")

cat("\\caption{Estimated peak location $\\hat{L}(T)$ for China COVID-19 data with $h=60$ using quantile regression ($\\tau = 0.95$).}\n")
cat("\\label{tab:china_Lhat_h60}\n\\end{table}\n")

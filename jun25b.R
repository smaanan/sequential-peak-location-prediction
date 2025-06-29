# --- Libraries ---
library(dplyr)
library(tidyr)
library(gt)
library(tibble)
library(glue)
library(readr)
library(zoo)

# --- Parameters ---
taus <- c(0.25, 0.5, 0.75, 0.95)
threshold <- 15
threshold_prime <- 2
h <- 160
T_vals <- 40:75  # From Table 11 in the paper

# --- Load ECDC Germany COVID-19 Data (to get t in [1,160]) ---
get_de_cases_ecdc <- function(start_date = "2020-01-28", end_date = "2020-07-05") {
  url <- "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
  
  read_csv(url, col_types = cols()) %>%
    filter(countriesAndTerritories == "Germany") %>%
    mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
    filter(date >= as.Date(start_date) & date <= as.Date(end_date)) %>%
    arrange(date) %>%
    select(date, daily_cases = cases) %>%
    mutate(daily_cases = ifelse(daily_cases < 0, NA, daily_cases),
           daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE))
}

# --- Load Data ---
de <- get_de_cases_ecdc()
I_full <- de$daily_cases

# --- Estimation Functions ---
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
  optimize(f, c(1, 40), tol = 1e-4)$minimum
}

estimate_params <- function(I, tau, h) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau)
    a <- compute_a(I, l, s, tau)
    objective_function(I, a, l, s, tau)
  })
  l_hat <- which.min(losses)
  list(l_hat = l_hat)
}

find_stable_time <- function(l_path, T_vals, tol) {
  for (i in seq_along(T_vals)) {
    future_vals <- l_path[i:length(l_path)]
    if (max(future_vals) - min(future_vals) <= tol) {
      return(T_vals[i])
    }
  }
  return(NA)
}

# --- Store objective function values ---
# min_objectives <- list()
# 
# # --- Main Loop Over τ Values ---
# for (tau in taus) {
#   l_path <- sapply(T_vals, function(T) {
#     I <- I_full[1:T]
#     estimate_params(I, tau, h)$l_hat
#   })
#   
#   T_max <- max(T_vals)
#   I_Tmax <- I_full[1:T_max]
#   res <- estimate_params(I_Tmax, tau, h)
#   s_opt <- search_best_s(I_Tmax, res$l_hat, tau)
#   a_opt <- compute_a(I_Tmax, res$l_hat, s_opt, tau)
#   obj_val <- objective_function(I_Tmax, a_opt, res$l_hat, s_opt, tau)
#   min_objectives[[as.character(tau)]] <- obj_val
#   
#   T_entry <- find_stable_time(l_path, T_vals, threshold)
#   T_exit  <- find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)
#   
#   table_data <- tibble(
#     T_1 = T_vals[1:18],
#     lhat_1 = l_path[1:18],
#     T_2 = T_vals[19:36],
#     lhat_2 = l_path[19:36]
#   )
#   
#   transposed_matrix <- t(table_data)
#   colnames(transposed_matrix) <- paste0("V", seq_len(ncol(transposed_matrix)))
#   transposed_data <- as.data.frame(transposed_matrix) %>%
#     rownames_to_column(var = "Label") %>%
#     as_tibble()
#   
#   transposed_data$Label <- recode(transposed_data$Label,
#                                   "T_1" = "T", "T_2" = "T",
#                                   "lhat_1" = "l̂(T)", "lhat_2" = "l̂(T)")
#   
#   subtitle_text <- glue(
#     "COVID-19 in Germany, τ = {tau}. ",
#     "Entry at **T_entry = {T_entry}** (ε = {threshold}), ",
#     "exit at **T_exit = {T_exit}** (ε' = {threshold_prime})."
#   )
#   
#   gt_table <- transposed_data %>%
#     gt() %>%
#     tab_header(
#       title = md("**Table: Sequence T → l̂(T)**"),
#       subtitle = md(subtitle_text)
#     ) %>%
#     cols_label(.list = setNames(rep("", ncol(transposed_data)), names(transposed_data))) %>%
#     tab_options(
#       column_labels.hidden = TRUE,
#       table.font.size = "small",
#       heading.title.font.size = 14,
#       heading.subtitle.font.size = 11,
#       data_row.padding = px(3)
#     )
#   
#   # Save outputs
#   png_file <- glue("germany_table8_lhat_path_tau{tau}.png")
#   tex_file <- glue("germany_table8_lhat_path_tau{tau}.tex")
#   gtsave(gt_table, png_file)
#   cat(glue("Saved: {png_file}\n"))
#   cat(as.character(as_latex(gt_table)), file = tex_file)
#   cat(glue("Saved: {tex_file}\n"))
# }
# 
# # --- Report the tau that minimizes the objective function ---
# best_tau <- names(which.min(unlist(min_objectives)))
# cat(glue("The τ that minimizes the objective function at T = {max(T_vals)} is τ = {best_tau}\n"))




# --- Create Combined Table Output --------------------------------------------------

all_blocks <- list()

for (tau in taus) {
  l_path <- sapply(T_vals, function(T) {
    I <- I_full[1:T]
    estimate_params(I, tau, h)$l_hat
  })
  
  T_max <- max(T_vals)
  I_Tmax <- I_full[1:T_max]
  res <- estimate_params(I_Tmax, tau, h)
  l_hat_Tmax <- res$l_hat
  s_hat <- search_best_s(I_Tmax, l_hat_Tmax, tau)
  a_hat <- compute_a(I_Tmax, l_hat_Tmax, s_hat, tau)
  obj_val <- objective_function(I_Tmax, a_hat, l_hat_Tmax, s_hat, tau)
  min_objectives[[as.character(tau)]] <- obj_val
  
  T_entry <- find_stable_time(l_path, T_vals, threshold)
  T_exit <- find_stable_time(l_path[T_vals >= T_entry], T_vals[T_vals >= T_entry], threshold_prime)
  
  # Split into chunks of 18 (Germany: 36 values = 2 chunks)
  T_chunks <- split(T_vals, ceiling(seq_along(T_vals) / 18))
  l_chunks <- split(l_path, ceiling(seq_along(l_path) / 18))
  
  for (i in seq_along(T_chunks)) {
    if (i == 1) {
      title <- glue("$\\tau = {tau}$ | T_{{entry}} = {T_entry}, T_{{exit}} = {T_exit}, \\varepsilon = {threshold}, \\varepsilon' = {threshold_prime}")
      all_blocks[[length(all_blocks) + 1]] <- list(Label = title, is_title = TRUE)
    }
    all_blocks[[length(all_blocks) + 1]] <- list(Label = "T", is_title = FALSE, values = T_chunks[[i]])
    all_blocks[[length(all_blocks) + 1]] <- list(Label = "\\hat{l}(T)", is_title = FALSE, values = l_chunks[[i]])
  }
}

# Format as tibble for gt preview (PNG only)
max_cols <- max(sapply(all_blocks, function(row) if (!row$is_title) length(row$values) else 0))
table_df <- tibble(Label = sapply(all_blocks, function(x) x$Label))
for (i in 1:max_cols) {
  table_df[[paste0("V", i)]] <- sapply(all_blocks, function(row) {
    if (!row$is_title && length(row$values) >= i) row$values[i] else ""
  })
}

# Save PNG using gt (visual preview only)
gt_table <- table_df %>%
  gt() %>%
  tab_header(
    title = md("**Combined Table: Sequence T → \\u0302l(T)**"),
    subtitle = md("Germany, with \\( \\tau \\in \\{0.25, 0.5, 0.75, 0.95\\} \\)")
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold", align = "center")),
    locations = cells_body(rows = grepl("^\\$", Label))
  ) %>%
  cols_label(.list = setNames(rep("", ncol(table_df)), names(table_df))) %>%
  tab_options(
    column_labels.hidden = TRUE,
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )

gtsave(gt_table, "germany_table8_lhat_paths_combined.png")

# --- Custom LaTeX Export with Centered Subheaders -------------------------

tabular_line <- glue("\\begin{{tabular*}}{{\\linewidth}}{{@{{\\extracolsep{{\\fill}}}}l{paste(rep('r', max_cols), collapse = '')}}}")

latex_lines <- c(
  "\\begin{table}[t]",
  "\\caption*{",
  "{\\large \\textbf{Combined Table: Sequence T → \\( \\hat{l}(T) \\)}} \\\\",
  "{\\small Germany, with \\( \\tau \\in \\{0.25, 0.5, 0.75, 0.95\\} \\)}",
  "}",
  tabular_line,
  "\\toprule"
)

for (row in all_blocks) {
  if (row$is_title) {
    latex_lines <- c(
      latex_lines,
      "\\midrule\\addlinespace[2.5pt]",
      glue("\\multicolumn{{{1 + max_cols}}}{{c}}{{\\textbf{{{row$Label}}}}} \\\\"),
      "\\addlinespace[2.5pt]",
      "\\midrule"
    )
  } else {
    latex_row <- c(row$Label, row$values)
    latex_row_fmt <- paste(latex_row, collapse = " & ")
    latex_lines <- c(latex_lines, paste0(latex_row_fmt, " \\\\"))
  }
}

latex_lines <- c(latex_lines, "\\bottomrule", "\\end{tabular*}", "\\end{table}")

writeLines(latex_lines, "germany_table8_lhat_paths_combined.tex")
cat("✓ Saved PNG and custom LaTeX table with centered subheaders.\n")

# --- Report τ with minimum objective ---
best_tau <- names(which.min(unlist(min_objectives)))
cat(glue("The τ that minimizes the objective at T = {T_max} is τ = {best_tau}\n"))

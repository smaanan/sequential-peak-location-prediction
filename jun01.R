# --- Libraries ---
library(dplyr)
library(zoo)
library(COVID19)

# --- Load COVID-19 data ---
get_country_cases <- function(country, start_date, end_date) {
  if (country == "France") {
    data <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv",
                     na.strings = "", fileEncoding = "UTF-8-BOM") %>%
      filter(countriesAndTerritories == "France") %>%
      mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
      filter(date >= as.Date(start_date) & date <= as.Date(end_date)) %>%
      arrange(date) %>%
      mutate(daily_cases = ifelse(cases < 0, NA, cases),
             daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE)) %>%
      pull(daily_cases)
  } else {
    data <- covid19(country, start = start_date, end = end_date, level = 1, verbose = FALSE) %>%
      arrange(date) %>%
      mutate(new_confirmed = c(NA, diff(confirmed))) %>%
      filter(!is.na(new_confirmed)) %>%
      pull(new_confirmed)
  }
  return(data)
}

# --- Algorithm 1: Tmin function ---
Tmin <- function(data, T_initial, h, tau_1, tau_2, tau_3) {
  objf <- function(x, l, T) sum(abs(x[1] * exp(-((1:T) - l)^2 / (x[2])^2) - data[1:T])) / T
  
  med_weighted <- function(l, s, T) {
    x <- data[1:T] * exp(((1:T) - l)^2 / s^2)
    o <- order(x)
    w <- exp(-((1:T) - l)^2 / s^2)
    wo <- w[o]
    W <- sum(w)
    y <- which(cumsum(wo) <= (W / 2))
    return(x[o[max(y) + 1]])
  }
  
  objfw <- function(l, s, T) objf(c(med_weighted(l, s, T), s), l, T)
  
  hat_l <- function(T) {
    out <- matrix(0, h, 3)
    for (l in 1:h) {
      smin <- optimize(function(s) objfw(l, s, T), c(1, 20), tol = 1e-4)
      out[l, ] <- c(med_weighted(l, smin$minimum, T), smin$minimum, smin$objective)
    }
    list(out, which.min(out[, 3]))
  }
  
  premier_indice <- function(x) {
    for (i in 1:(length(x) - tau_1)) {
      if (all(abs(x[i:(i + tau_1)] - x[i]) <= tau_2)) return(i)
    }
    return(NA)
  }
  
  deuxieme_indice <- function(x) {
    pi <- premier_indice(x)
    if (is.na(pi)) return(NA)
    for (i in (pi + 1):(length(x) - tau_1)) {
      if (all(abs(x[i:(length(x) - tau_1)] - x[i]) <= tau_3)) return(i)
    }
    return(NA)
  }
  
  x <- sapply(T_initial:h, function(T) hat_l(T)[[2]])
  list(x, T_initial - 1 + premier_indice(x), T_initial - 1 + deuxieme_indice(x))
}

# --- LaTeX Table Builder ---
build_latex_table_block <- function(T_vals, L_vals) {
  T_line <- paste0("$T$ & ", paste(T_vals, collapse = " & "), " \\")
  L_line <- paste0("$\\hat{l}(T)$ & ", paste(ifelse(is.na(L_vals), "--", L_vals), collapse = " & "), " \\")
  paste("\\begin{tabular}{|", paste(rep("r", length(T_vals) + 1), collapse = ""), "|}",
        "\\hline", T_line, "\\hline", L_line, "\\\\ \\hline",
        "\\end{tabular}", sep = "\n")
}

# --- Run for a given country and parameters ---
run_estimation_and_print_table <- function(country, start_date, end_date, h, caption_label) {
  data <- get_country_cases(country, start_date, end_date)
  result <- Tmin(data, T_initial = 4, h = h, tau_1 = 3, tau_2 = 1, tau_3 = 1)
  
  T_vals <- 4:h
  L_vals <- result[[1]]
  T_vals_trimmed <- T_vals[2:length(T_vals)]
  L_vals_trimmed <- L_vals[2:length(L_vals)]
  
  block1 <- T_vals_trimmed[1:32];  L1 <- L_vals_trimmed[1:32]
  block2 <- T_vals_trimmed[33:64]; L2 <- L_vals_trimmed[33:64]
  block3 <- T_vals_trimmed[65:96]; L3 <- L_vals_trimmed[65:96]
  
  cat(paste0("%% LaTeX Table for Estimated $\\hat{l}(T)$ from Algorithm 1, ", country, " COVID-19, h = ", h, "\n\n"))
  cat("\\begin{table}[H]\n\\centering\n")
  cat(build_latex_table_block(block1, L1), "\n\n\\vspace{0.4cm}\n\n")
  cat(build_latex_table_block(block2, L2), "\n\n\\vspace{0.4cm}\n\n")
  cat(build_latex_table_block(block3, L3), "\n")
  cat(paste0("\\caption{Estimated peak location $\\hat{l}(T)$ for ", country, " COVID-19 data with $h=", h, "$ using Algorithm 1.}\n"))
  cat(paste0("\\label{tab:lhat_", tolower(country), "_h", h, "}\n\\end{table}\n"))
}

# --- Run for China and France ---
run_estimation_and_print_table("China", "2020-01-22", "2020-04-01", 60, "lhat_china_h60")
run_estimation_and_print_table("France", "2020-01-15", "2020-05-31", 100, "lhat_france_h100")

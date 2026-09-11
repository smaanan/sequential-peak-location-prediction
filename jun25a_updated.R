############################################################
# COVID-19 FRANCE — QUANTILE PEAK-LOCATION ESTIMATION
#
# Updated version:
# - Candidate peak locations are restricted to l = 1,...,T
#   at each observation horizon T.
# - Only the stabilization/entry time T_entry is computed.
# - T_exit and threshold_prime have been removed.
#
# Data:
#   ECDC daily COVID-19 cases
#   January 15 to May 31, 2020
#
# Observation horizons:
#   T = 50,...,89
#
# Quantile levels:
#   tau = 0.25, 0.50, 0.75, 0.95
#
# Stabilization tolerance:
#   epsilon = 5
############################################################


# --- Libraries -------------------------------------------------------------

library(dplyr)
library(tidyr)
library(gt)
library(tibble)
library(glue)
library(readr)
library(zoo)


# --- Parameters ------------------------------------------------------------

taus <- c(0.25, 0.50, 0.75, 0.95)

# Stabilization tolerance
threshold <- 5

# Observation horizons
T_vals <- 50:89

# Store objective values at T_max
min_objectives <- list()


# --- Load ECDC France COVID-19 Data ----------------------------------------

get_france_cases_ecdc <- function(
    start_date = "2020-01-15",
    end_date   = "2020-05-31"
) {
  
  url <- paste0(
    "https://opendata.ecdc.europa.eu/",
    "covid19/casedistribution/csv"
  )
  
  read_csv(
    url,
    col_types = cols()
  ) %>%
    filter(
      countriesAndTerritories == "France"
    ) %>%
    mutate(
      date = as.Date(
        sprintf(
          "%04d-%02d-%02d",
          year,
          month,
          day
        )
      )
    ) %>%
    filter(
      date >= as.Date(start_date),
      date <= as.Date(end_date)
    ) %>%
    arrange(date) %>%
    select(
      date,
      daily_cases = cases
    ) %>%
    mutate(
      daily_cases = ifelse(
        daily_cases < 0,
        NA,
        daily_cases
      ),
      daily_cases = zoo::na.locf(
        daily_cases,
        na.rm = FALSE
      )
    )
}


france <- get_france_cases_ecdc()

# Raw daily cases
I_full <- france$daily_cases


# --- Quantile Loss ----------------------------------------------------------

quantile_loss <- function(u, tau) {
  u * (tau - (u < 0))
}


# --- Weighted Quantile for the Amplitude -----------------------------------

med_weighted_tau <- function(I, l, s, tau) {
  
  T <- length(I)
  
  x <- I * exp(
    ((1:T) - l)^2 / s^2
  )
  
  w <- exp(
    -((1:T) - l)^2 / s^2
  )
  
  o <- order(x)
  
  cum_w <- cumsum(
    w[o]
  )
  
  cutoff <- tau * sum(w)
  
  x[o][
    min(
      which(
        cum_w >= cutoff
      )
    )
  ]
}


compute_a <- function(I, l, s, tau) {
  
  med_weighted_tau(
    I,
    l,
    s,
    tau
  )
}


# --- Quantile-Loss Objective ------------------------------------------------

objective_function <- function(
    I,
    a,
    l,
    s,
    tau
) {
  
  t <- seq_along(I)
  
  fitted <- a * exp(
    -((t - l)^2) / s^2
  )
  
  mean(
    quantile_loss(
      I - fitted,
      tau
    )
  )
}


# --- Optimization over the Width -------------------------------------------

search_best_s <- function(
    I,
    l,
    tau
) {
  
  f <- function(s) {
    
    a <- compute_a(
      I,
      l,
      s,
      tau
    )
    
    objective_function(
      I,
      a,
      l,
      s,
      tau
    )
  }
  
  optimize(
    f,
    c(1, 40),
    tol = 1e-4
  )$minimum
}


# --- Estimate Peak Location -------------------------------------------------
#
# At each observation horizon T, only the currently observed locations
# l = 1,...,T are considered.
# ---------------------------------------------------------------------------

estimate_params <- function(
    I,
    tau
) {
  
  T <- length(I)
  
  losses <- sapply(
    1:T,
    function(l) {
      
      s <- search_best_s(
        I,
        l,
        tau
      )
      
      a <- compute_a(
        I,
        l,
        s,
        tau
      )
      
      objective_function(
        I,
        a,
        l,
        s,
        tau
      )
    }
  )
  
  l_hat <- which.min(
    losses
  )
  
  list(
    l_hat = l_hat
  )
}


# --- Stabilization Time -----------------------------------------------------

find_stable_time <- function(
    l_path,
    T_vals,
    tol
) {
  
  for (i in seq_along(T_vals)) {
    
    future_vals <- l_path[
      i:length(l_path)
    ]
    
    if (
      max(future_vals) -
      min(future_vals) <= tol
    ) {
      
      return(
        T_vals[i]
      )
    }
  }
  
  return(NA)
}


# --- Create Combined Table Output ------------------------------------------

all_blocks <- list()


for (tau in taus) {
  
  # ---------------------------------------------------------
  # Estimation path
  # ---------------------------------------------------------
  
  l_path <- sapply(
    T_vals,
    function(T) {
      
      I <- I_full[
        1:T
      ]
      
      estimate_params(
        I,
        tau
      )$l_hat
    }
  )
  
  
  # ---------------------------------------------------------
  # Objective at T_max
  # ---------------------------------------------------------
  
  T_max <- max(T_vals)
  
  I_Tmax <- I_full[
    1:T_max
  ]
  
  res <- estimate_params(
    I_Tmax,
    tau
  )
  
  l_hat_Tmax <- res$l_hat
  
  s_hat <- search_best_s(
    I_Tmax,
    l_hat_Tmax,
    tau
  )
  
  a_hat <- compute_a(
    I_Tmax,
    l_hat_Tmax,
    s_hat,
    tau
  )
  
  obj_val <- objective_function(
    I_Tmax,
    a_hat,
    l_hat_Tmax,
    s_hat,
    tau
  )
  
  min_objectives[[as.character(tau)]] <- obj_val
  
  
  # ---------------------------------------------------------
  # Stabilization time
  # ---------------------------------------------------------
  
  T_entry <- find_stable_time(
    l_path,
    T_vals,
    threshold
  )
  
  
  # ---------------------------------------------------------
  # Split rows into chunks of 20
  # ---------------------------------------------------------
  
  T_chunks <- split(
    T_vals,
    ceiling(
      seq_along(T_vals) / 20
    )
  )
  
  l_chunks <- split(
    l_path,
    ceiling(
      seq_along(l_path) / 20
    )
  )
  
  
  # ---------------------------------------------------------
  # Build table blocks
  # ---------------------------------------------------------
  
  for (i in seq_along(T_chunks)) {
    
    if (i == 1) {
      
      title <- glue(
        "$\\tau = {tau}$ | ",
        "T_{{entry}} = {T_entry}, ",
        "\\varepsilon = {threshold}"
      )
      
      all_blocks[[length(all_blocks)+1]] <- list(Label = title,is_title = TRUE)
    }
    
    
    all_blocks[[length(all_blocks)+1]] <- list(Label = "T",is_title = FALSE,values = T_chunks[[i]])
    
    all_blocks[[length(all_blocks)+1]] <- list(Label = "\\hat{l}(T)",is_title = FALSE,values = l_chunks[[i]])
  }
}


# --- Format as Tibble -------------------------------------------------------

max_cols <- max(
  sapply(
    all_blocks,
    function(row)
      if (!row$is_title)
        length(row$values)
    else
      0
  )
)


table_df <- tibble(
  Label = sapply(
    all_blocks,
    function(x) x$Label
  )
)


for (i in 1:max_cols) {
  
  table_df[[paste0("V", i)]] <- sapply(
    all_blocks,
    function(row) {
      
      if (
        !row$is_title &&
        length(row$values) >= i
      ) {
        
        row$values[i]
        
      } else {
        
        ""
      }
    }
  )
}


# --- Save PNG ---------------------------------------------------------------

gt_table <- table_df %>%
  gt() %>%
  tab_header(
    title = md(
      "**Combined Table: Sequence T → \\hat{l}(T)**"
    ),
    subtitle = md(
      "France, with \\( \\tau \\in \\{0.25, 0.5, 0.75, 0.95\\} \\)"
    )
  ) %>%
  tab_style(
    style = list(
      cell_text(
        weight = "bold",
        align = "center"
      )
    ),
    locations = cells_body(
      rows = grepl(
        "^\\$",
        Label
      )
    )
  ) %>%
  cols_label(
    .list = setNames(
      rep("", ncol(table_df)),
      names(table_df)
    )
  ) %>%
  tab_options(
    column_labels.hidden = TRUE,
    table.font.size = "small",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 11,
    data_row.padding = px(3)
  )


gtsave(
  gt_table,
  "france_table8_lhat_paths_combined.png"
)


# --- Custom LaTeX Export ----------------------------------------------------

tabular_line <- glue(
  "\\begin{{tabular*}}{{\\linewidth}}",
  "{{@{{\\extracolsep{{\\fill}}}}l",
  "{paste(rep('r', max_cols), collapse = '')}}}"
)


latex_lines <- c(
  
  "\\begin{table}[t]",
  
  "\\caption*{",
  
  "{\\large \\textbf{Combined Table: ",
  "Sequence T → \\( \\hat{l}(T) \\)}} \\\\",
  
  "{\\small France, with ",
  "\\( \\tau \\in \\{0.25, 0.5, 0.75, 0.95\\} \\)}",
  
  "}",
  
  tabular_line,
  
  "\\toprule"
)


for (row in all_blocks) {
  
  if (row$is_title) {
    
    latex_lines <- c(
      latex_lines,
      
      "\\midrule\\addlinespace[2.5pt]",
      
      glue(
        "\\multicolumn{{{1 + max_cols}}}",
        "{{c}}{{\\textbf{{{row$Label}}}}} \\\\"
      ),
      
      "\\addlinespace[2.5pt]",
      
      "\\midrule"
    )
    
  } else {
    
    latex_row <- c(
      row$Label,
      row$values
    )
    
    latex_row_fmt <- paste(
      latex_row,
      collapse = " & "
    )
    
    latex_lines <- c(
      latex_lines,
      paste0(
        latex_row_fmt,
        " \\\\"
      )
    )
  }
}


latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular*}",
  "\\end{table}"
)


writeLines(
  latex_lines,
  "france_table8_lhat_paths_combined.tex"
)


cat(
  "Saved France PNG and LaTeX table.\n"
)


# --- Best tau ---------------------------------------------------------------

best_tau <- names(
  which.min(
    unlist(min_objectives)
  )
)


cat(
  glue(
    "The tau that minimizes the objective ",
    "at T = {T_max} is tau = {best_tau}\n"
  )
)
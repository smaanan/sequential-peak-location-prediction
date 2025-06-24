# ---------------------------------------------------------------------------
# Entry/Exit Time Estimation under GAL Log-Likelihood Objective
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(gt)
library(ggplot2)
library(purrr)
library(tibble)
library(glue)

# --- GAL density (used in log-likelihood) ---
dgal <- function(x, tau, z, b) {
  ifelse(x < 0,
         tau * (1 + z) / b * exp((1 + z) * x / b),
         (1 - tau) * (1 - z) / b * exp(-(1 - z) * x / b))
}

log_likelihood <- function(I, a, l, s, tau, z, b) {
  t <- seq_along(I)
  u <- I - a * exp(-((t - l)^2) / s^2)
  ll <- log(dgal(u, tau, z, b))
  -sum(ll) / length(I)  # normalized negative log-likelihood
}

# --- Estimation Functions ---
search_best_s <- function(I, l, tau, z, b) {
  f <- function(s) {
    a <- compute_a(I, l, s, tau)
    log_likelihood(I, a, l, s, tau, z, b)
  }
  optimize(f, c(1, 20), tol = 1e-4)$minimum
}

compute_a <- function(I, l, s, tau) {
  T <- length(I)
  x <- I * exp(((1:T) - l)^2 / s^2)
  w <- exp(-((1:T) - l)^2 / s^2)
  o <- order(x)
  cum_w <- cumsum(w[o])
  cutoff <- tau * sum(w)
  x[o][min(which(cum_w >= cutoff))]
}

estimate_params <- function(I, tau, z, b, h) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau, z, b)
    a <- compute_a(I, l, s, tau)
    log_likelihood(I, a, l, s, tau, z, b)
  })
  l_hat <- which.min(losses)
  s_hat <- search_best_s(I, l_hat, tau, z, b)
  a_hat <- compute_a(I, l_hat, s_hat, tau)
  list(l_hat = l_hat, s_hat = s_hat, a_hat = a_hat)
}

# --- Noise generator ---
rgen_gal <- function(n, tau, z, b) {
  u <- runif(n)
  ifelse(u < tau,
         -rexp(n, rate = (1 + z) / b),
         rexp(n, rate = (1 - z) / b))
}

# --- Stability Detection ---
find_stable_time <- function(l_path, T_vals, tol) {
  for (i in seq_along(T_vals)) {
    future_vals <- l_path[i:length(l_path)]
    if (max(future_vals) - min(future_vals) <= tol) {
      return(T_vals[i])
    }
  }
  return(NA)
}

# --- Main Grid Simulation ---
run_sim_grid <- function(A = 1, L = 22, S = 10, snr = 40,
                         taus = c(0.5, 0.75, 0.95),
                         z_vals = c(-0.9, -0.5, 0, 0.5, 0.9),
                         b_vals = c(0.25, 0.5, 1, 2),
                         T_range = 4:35,
                         threshold = 4, threshold_prime = 0,
                         h = 35, seed = 42) {
  expand.grid(tau = taus, z = z_vals, b = b_vals) %>%
    mutate(data = pmap(list(tau, z, b), function(tau, z, b) {
      set.seed(seed)
      t_full <- 1:h
      signal <- A * exp(-((t_full - L)^2) / S^2)
      noise <- rgen_gal(h, tau, z, b)
      I <- signal + noise
      l_path <- sapply(T_range, function(T) estimate_params(I[1:T], tau, z, b, h)$l_hat)
      T_entry <- find_stable_time(l_path, T_range, threshold)
      T_exit  <- find_stable_time(l_path[T_range >= T_entry], T_range[T_range >= T_entry], threshold_prime)
      list(l_path = l_path, T_vals = T_range, T_entry = T_entry, T_exit = T_exit)
    }))
}


# --- Heatmap (Saved) ---
entry_exit_long <- run_sim_grid() %>%
  unnest_wider(data) %>%
  mutate(b = as.factor(b)) %>%
  pivot_longer(cols = c(T_entry, T_exit),
               names_to = "type", values_to = "T_val")

heatmap_plot <- ggplot(entry_exit_long, aes(x = factor(z), y = b, fill = T_val)) +
  geom_tile(color = "white") +
  facet_grid(type ~ tau, labeller = label_both) +
  scale_fill_viridis_c(option = "B", direction = -1, na.value = "grey90") +
  labs(
    title = "Entry and Exit Times under GAL Log-Likelihood Objective",
    x = "Asymmetry (z)", y = "Scale (b)", fill = "T",
    caption = "Each facet shows Entry (T_entry) or Exit (T_exit) for a given τ"
  ) +
  theme_minimal(base_size = 12)

ggsave("heatmap_entry_exit_loglik_tau_z_b.png", plot = heatmap_plot, width = 10, height = 6, dpi = 300)
cat("Saved: heatmap_entry_exit_loglik_tau_z_b.png\n")


# --- Horizontal gt tables by b ---
summary_table <- run_sim_grid() %>%
  unnest_wider(data) %>%
  select(tau, z, b, T_entry, T_exit)

b_vals <- unique(summary_table$b)

for (b_val in b_vals) {
  df_b <- summary_table %>%
    filter(b == b_val) %>%
    pivot_longer(cols = c(T_entry, T_exit), names_to = "type", values_to = "Tval") %>%
    mutate(
      tau_lab = paste0("τ = ", tau),
      colname = paste0(tau_lab, "__T_", type)
    ) %>%
    select(z, colname, Tval) %>%
    pivot_wider(names_from = colname, values_from = Tval) %>%
    arrange(z)
  
  col_info <- tibble(
    name = names(df_b)[-1]
  ) %>%
    separate(name, into = c("tau", "metric"), sep = "__T_", remove = FALSE)
  
  gt_table <- df_b %>%
    gt(rowname_col = "z") %>%
    {
      gt_obj <- .
      for (tau_val in unique(col_info$tau)) {
        columns_in_tau <- col_info %>% filter(tau == tau_val) %>% pull(name)
        gt_obj <- gt_obj %>%
          tab_spanner(label = tau_val, columns = all_of(columns_in_tau))
      }
      gt_obj
    } %>%
    cols_label(.list = setNames(
      rep(c("T_entry", "T_exit"), times = length(unique(col_info$tau))),
      col_info$name
    )) %>%
    tab_header(
      title = md(glue("**Entry and Exit Times for b = {b_val}**")),
      subtitle = md("Rows = z values; Columns = (τ, T_entry / T_exit)")
    ) %>%
    fmt_number(everything(), decimals = 1) %>%
    tab_options(
      table.font.size = "small",
      heading.title.font.size = 14,
      heading.subtitle.font.size = 11,
      data_row.padding = px(3)
    )
  
  gtsave(gt_table, glue("gt_table_loglik_entry_exit_b{b_val}.png"))
  cat(glue("Saved: gt_table_loglik_entry_exit_b{b_val}.png\n"))
}

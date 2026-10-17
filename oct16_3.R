library(dplyr); library(readr); library(zoo)

############################################################
# 1. DATA
############################################################
get_ecdc <- function(country, start_date, end_date,
                     file = "6YbQxT6Z", fix_negative = FALSE) {
  dat <- readr::read_csv(file, col_types = readr::cols(),
                         progress = FALSE) %>%
    filter(countriesAndTerritories == country) %>%
    mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
    filter(date >= as.Date(start_date), date <= as.Date(end_date)) %>%
    arrange(date) %>% select(date, daily_cases = cases)
  if (fix_negative)
    dat <- dat %>% mutate(
      daily_cases = ifelse(daily_cases < 0, NA, daily_cases),
      daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE))
  dat$daily_cases
}

I_france  <- get_ecdc("France",  "2020-01-15", "2020-05-31")
I_germany <- get_ecdc("Germany", "2020-01-28", "2020-07-05",
                      fix_negative = TRUE)
# China: I_full as in your China script (COVID19 package)
I_china <- I_full

############################################################
# 2. ESTIMATOR (s.max set per country)
############################################################
quantile_loss <- function(u, tau) u * (tau - (u < 0))

med_weighted_tau <- function(I, l, s, tau) {
  T <- length(I); tt <- 1:T
  x <- I * exp((tt - l)^2 / s^2); w <- exp(-(tt - l)^2 / s^2)
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  if (any(w < 0) || sum(w) <= 0) return(NA_real_)
  o <- order(x)
  x[o][min(which(cumsum(w[o]) >= tau * sum(w)))]
}
compute_a <- function(I, l, s, tau) med_weighted_tau(I, l, s, tau)

objective_function <- function(I, a, l, s, tau) {
  t <- seq_along(I)
  mean(quantile_loss(I - a * exp(-((t - l)^2) / s^2), tau))
}

search_best_s <- function(I, l, tau, s.max) {
  f <- function(s) objective_function(I, compute_a(I, l, s, tau),
                                      l, s, tau)
  optimize(f, c(1, s.max), tol = 1e-4)$minimum
}

estimate_l <- function(I, tau, h, s.max) {
  losses <- sapply(1:h, function(l) {
    s <- search_best_s(I, l, tau, s.max)
    objective_function(I, compute_a(I, l, s, tau), l, s, tau)
  })
  which.min(losses)
}

estimate_path <- function(I_full, T_vals, tau, s.max) {
  h <- max(T_vals)
  data.frame(T = T_vals,
             hat.l = sapply(T_vals, function(T)
               estimate_l(I_full[1:T], tau, h, s.max)))
}

############################################################
# 3. RUN ONE COUNTRY (paths cached)
############################################################
taus  <- c(0.25, 0.50, 0.75, 0.95)
rules <- c("narrowest", "first_k")

run_country <- function(name, I_full, T_vals, s.max, L_obs = NA) {
  cands <- list(); sels <- list(); paths <- list()
  for (tau in taus) {
    cat(name, "tau =", tau, "\n")
    f <- sprintf("%s_path_tau%.2f.rds", tolower(name), tau)
    path <- if (file.exists(f)) readRDS(f) else {
      p <- estimate_path(I_full, T_vals, tau, s.max); saveRDS(p, f); p }
    paths[[length(paths) + 1L]] <- transform(path, tau = tau)
    for (Tc in path$T) {
      s2 <- run_entry_exit_at_T(path$T, path$hat.l, Tc, Data = name,
                                L = L_obs)$stage2
      s2$tau <- tau
      cands[[length(cands) + 1L]] <- s2
      for (rl in rules) {
        sel <- select_estimate_candidate(s2, L_obs, rule = rl)
        ok <- !is.na(sel$T_entry) && !is.na(sel$T_exit)
        sels[[length(sels) + 1L]] <- data.frame(
          country = name, tau = tau, rule = rl, T = Tc,
          hat_l = sel$hat_l_current, selected = ok,
          k = if (ok) sel$k else NA_real_,
          T_entry = if (ok) sel$T_entry else NA_real_,
          T_exit  = if (ok) sel$T_exit  else NA_real_)
      }
    }
  }
  list(cands = do.call(rbind, cands), sel = do.call(rbind, sels),
       paths = do.call(rbind, paths))
}

############################################################
# 4. SIZE OF K(T)
############################################################
analyze_K <- function(res, name) {
  cc <- res$cands
  cc$adm <- !is.na(cc$estimate.in.window) & cc$estimate.in.window
  sz <- do.call(rbind, lapply(
    split(cc, list(cc$tau, cc$T_current), drop = TRUE),
    function(d) data.frame(
      tau = d$tau[1], T = d$T_current[1], n_cand = nrow(d),
      n_K = sum(d$adm), K = paste(d$k[d$adm], collapse = ","),
      min_K = if (any(d$adm)) min(d$k[d$adm]) else NA,
      hat_l = d$hat_l_current[1])))
  sz <- sz[order(sz$tau, sz$T), ]; rownames(sz) <- NULL
  
  cat("\n====", name, ": K(T) at the terminal horizon ====\n")
  print(sz[sz$T == max(sz$T), ], row.names = FALSE)
  
  cat("\n====", name, ": summary by tau (all horizons) ====\n")
  print(do.call(rbind, lapply(split(sz, sz$tau), function(d)
    data.frame(tau = d$tau[1],
               prop_empty = mean(d$n_K == 0),
               prop_single = mean(d$n_K == 1),
               prop_multi = mean(d$n_K >= 2),
               max_K = max(d$n_K),
               mean_K_nonempty = mean(d$n_K[d$n_K > 0]),
               mean_n_cand = mean(d$n_cand)))), row.names = FALSE)
  
  cat("\n====", name, ": |K| counts (0,1,2,3+) ====\n")
  print(table(sz$tau, pmin(sz$n_K, 3)))
  
  d2 <- res$sel[res$sel$selected, ]
  w <- reshape(d2[, c("tau", "T", "rule", "k", "T_entry", "T_exit")],
               idvar = c("tau", "T"), timevar = "rule",
               direction = "wide")
  w$differ <- w$k.narrowest != w$k.first_k
  cat("\n====", name, ": share of horizons where rules differ ====\n")
  print(aggregate(differ ~ tau, w, mean))
  print(w[w$differ, ], row.names = FALSE)
  
  write.csv(sz, sprintf("%s_K_sizes.csv", tolower(name)),
            row.names = FALSE)
  invisible(sz)
}

############################################################
# 5. RUN
############################################################
res_france  <- run_country("France",  I_france,  50:89, s.max = 40)
res_germany <- run_country("Germany", I_germany, 40:75, s.max = 40)
sz_france   <- analyze_K(res_france,  "France")
sz_germany  <- analyze_K(res_germany, "Germany")

# China, to compare on the same footing (s in [1,20], T = 4:35)
res_china <- run_country("China", I_china, 4:35, s.max = 20)
sz_china  <- analyze_K(res_china, "China")

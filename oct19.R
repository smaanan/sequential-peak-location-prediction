############################################################
# REAL-DATA ANALYSIS OF THE ADMISSIBLE SETS K(T)
# (Section "Admissible candidates K(T) on the real data")
#
# Reproduces, for China, France and Germany:
#   Table  real_K_terminal      K(h) at the terminal horizon
#   Table  real_K_props         size of K(T) over all horizons
#   Text   rule disagreement    min K(T) vs narrowest candidate
#   Text   position vs L        pooled covering counts at L-3, L+3, L+10
#   Table  real_data_intervals  illustrative (non-operational) candidates
#   Table  real_data_minK       operational rule  k_hat(T) = min K(T)
#   Text   coverage by k        share of admissible candidates covering L
#   (extra) union / envelope diagnostics, not reported in the manuscript
#
# REQUIREMENTS (run the three country scripts first, in the same
# R session; they source fixed_h_core.R and create these objects):
#   I_china, I_france, I_germany            daily case series
#   res_china_fh, res_france_fh, res_germany_fh
#                                           output of run_country_fixed_h()
#
# NOTE: this script replaces the earlier ad hoc code. The per-horizon
# candidate table is rebuilt here (build_candidates) from the
# sequential paths stored in res_*_fh$paths; the column that used to
# be called `estimate.in.window` is now called `admissible`.
############################################################

suppressPackageStartupMessages(library(dplyr))


############################################################
# 0. CONFIGURATION
############################################################

OUT_DIR <- "real_K_output"
dir.create(OUT_DIR, showWarnings = FALSE)

TAUS <- c(0.25, 0.50, 0.75, 0.95)

# Observed peak days used in the manuscript (checked against the data below)
L_MANUSCRIPT <- c(China = 22, France = 78, Germany = 61)

# Horizons examined relative to the observed peak L
OFFSETS <- c(before = -3, after_near = 3, after_far = 10)

# Width cap used for the "illustrative candidates" table
MAX_WIDTH_ILLUSTRATIVE <- 10

# Extra diagnostics (union / envelope intervals); not in the manuscript
RUN_EXTRA_DIAGNOSTICS <- TRUE

# Rounding that matches the manuscript tables (round half UP;
# base R rounds exact halves to even, e.g. round(0.625, 2) = 0.62)
round_half_up <- function(x, digits = 2)
  floor(x * 10^digits + 0.5 + 1e-9) / 10^digits

round_df <- function(df, digits = 2, skip = c("tau", "max_size")) {
  num <- vapply(df, is.numeric, logical(1)) & !(names(df) %in% skip)
  df[num] <- lapply(df[num], round_half_up, digits = digits)
  df
}

show <- function(title, x) {
  cat("\n==== ", title, " ====\n", sep = "")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}

save_csv <- function(x, name)
  write.csv(x, file.path(OUT_DIR, name), row.names = FALSE)


############################################################
# 1. PREREQUISITES
############################################################

if (!exists("run_entry_exit")) source("fixed_h_core.R")

required <- c("run_entry_exit", "I_china", "I_france", "I_germany",
              "res_china_fh", "res_france_fh", "res_germany_fh")
missing_objs <- required[!vapply(required, exists, logical(1))]
if (length(missing_objs) > 0)
  stop("Run the country scripts first. Missing objects: ",
       paste(missing_objs, collapse = ", "))

country_fh <- list(China   = res_china_fh,
                   France  = res_france_fh,
                   Germany = res_germany_fh)

# Observed peak day = day of the maximum daily count (analysis window only)
L_obs <- c(China   = as.integer(which.max(I_china[1:35])),
           France  = as.integer(which.max(I_france[1:89])),
           Germany = as.integer(which.max(I_germany[1:75])))
print(L_obs)
if (!isTRUE(all.equal(unname(L_obs), unname(L_MANUSCRIPT))))
  warning("Observed peak days differ from the manuscript values (",
          paste(L_MANUSCRIPT, collapse = ", "), ").")


############################################################
# 2. CANDIDATES AT EVERY CURRENT HORIZON
#
# For each country, tau and current horizon T, the entry-exit
# construction is applied to the sequential path restricted to
# horizons <= T. Candidate k is ADMISSIBLE at T when it contains the
# current estimate:  T_entry <= l_hat(T) <= T_exit.
# At T = h this coincides with the terminal construction.
############################################################

build_candidates <- function(res_fh, country) {
  bind_rows(lapply(names(res_fh$paths), function(nm) {
    path <- res_fh$paths[[nm]]
    bind_rows(lapply(seq_len(nrow(path)), function(i) {
      st2 <- run_entry_exit(path$T[1:i], path$hat.l[1:i])$stage2
      data.frame(country = country, tau = as.numeric(nm),
                 T_current = path$T[i], hat_l_current = path$hat.l[i],
                 k = st2$step, T_entry = st2$T_entry, T_exit = st2$T_exit,
                 stringsAsFactors = FALSE)
    }))
  })) %>%
    mutate(width = T_exit - T_entry,
           admissible = !is.na(T_entry) &
             T_entry <= hat_l_current & hat_l_current <= T_exit)
}

# Sanity check: at T = h the rebuilt candidates must equal the
# terminal construction stored in res_*_fh$stage2
check_terminal_consistency <- function(cands, res_fh) {
  h <- max(cands$T_current)
  a <- cands %>% filter(T_current == h) %>%
    select(tau, k, T_entry, T_exit) %>% arrange(tau, k)
  b <- res_fh$stage2 %>%
    transmute(tau, k = step, T_entry, T_exit) %>% arrange(tau, k)
  stopifnot(isTRUE(all.equal(as.data.frame(a), as.data.frame(b),
                             check.attributes = FALSE)))
}

cands_list <- Map(build_candidates, country_fh, names(country_fh))
invisible(Map(check_terminal_consistency, cands_list, country_fh))

cands <- bind_rows(cands_list) %>%
  group_by(country) %>% mutate(h = max(T_current)) %>% ungroup()

adm <- cands %>% filter(admissible)          # admissible candidates only
save_csv(cands, "candidates_all_horizons.csv")


############################################################
# 3. SIZE OF K(T)
############################################################

# 3a. K(h) at the terminal horizon                   [real_K_terminal]
K_terminal <- adm %>%
  filter(T_current == h) %>%
  group_by(country, tau) %>%
  summarise(K = paste0("{", paste(sort(k), collapse = ","), "}"),
            .groups = "drop")
show("K(h) at the terminal horizon (long)", K_terminal)
cat("\nWide layout (rows: country, columns: tau)\n")
print(with(K_terminal, tapply(K, list(country, tau), identity)))
save_csv(K_terminal, "K_terminal.csv")

# 3b. Proportion of horizons with K(T) empty / singleton / >= 2
#     (every horizon of the grid is counted, including empty ones)   [real_K_props]
K_sizes <- cands %>%
  group_by(country, tau, T_current) %>%
  summarise(n = sum(admissible), .groups = "drop")

K_props <- K_sizes %>%
  group_by(country, tau) %>%
  summarise(n_horizons = n(),
            empty     = mean(n == 0),
            singleton = mean(n == 1),
            ge2       = mean(n >= 2),
            max_size  = max(n),
            mean_size_nonempty = mean(n[n > 0]),
            .groups = "drop")
show("Size of K(T) over horizons (rounded half up)",
     round_df(K_props, 2, skip = c("tau", "max_size", "n_horizons")))
save_csv(K_props, "K_props.csv")

# 3c. Disagreement between  min K(T)  and the narrowest candidate
#     (ties in width broken by the smaller k), among horizons with an
#     interval selected.
rule_disagreement <- adm %>%
  group_by(country, tau, T_current) %>%
  summarise(k_min    = min(k),
            k_narrow = k[order(width, k)][1],
            .groups = "drop") %>%
  group_by(country, tau) %>%
  summarise(pct_differ = 100 * mean(k_min != k_narrow), .groups = "drop")
show("Rule disagreement (% of selected horizons)", rule_disagreement)
show("Range over tau (text: 0-38%, 0-38%, 6-50%)",
     rule_disagreement %>% group_by(country) %>%
       summarise(min_pct = min(pct_differ), max_pct = max(pct_differ),
                 .groups = "drop"))
save_csv(rule_disagreement, "rule_disagreement.csv")


############################################################
# 4. POSITION OF THE CANDIDATES RELATIVE TO L
############################################################

# Nearest available horizon to L + offset (ties -> smaller horizon)
pick_horizon <- function(grid, target) grid[which.min(abs(grid - target))]

horizons <- bind_rows(lapply(names(L_obs), function(ct) {
  grid <- sort(unique(cands$T_current[cands$country == ct]))
  data.frame(country = ct, case = names(OFFSETS),
             T_current = vapply(OFFSETS, function(o)
               pick_horizon(grid, L_obs[[ct]] + o), numeric(1)),
             stringsAsFactors = FALSE)
}))
show("Horizons examined (L-3, L+3, L+10)", horizons)

# Every (country, case, tau) combination, so that empty K(T) shows up as NA
design <- expand.grid(country = names(L_obs), case = names(OFFSETS),
                      tau = TAUS, stringsAsFactors = FALSE) %>%
  inner_join(horizons, by = c("country", "case"))

# Admissible candidates at those horizons, with coverage of L
#   inclusive: T_entry <= L <= T_exit   (manuscript convention)
#   strict   : T_entry <  L <  T_exit   (old convention; check only)
tab_K <- adm %>%
  inner_join(horizons, by = c("country", "T_current")) %>%
  mutate(L = unname(L_obs[country]),
         covers_incl   = T_entry <= L & L <= T_exit,
         covers_strict = T_entry <  L & L <  T_exit) %>%
  arrange(country, case, tau, k)
save_csv(tab_K, "real_data_K_intervals.csv")

# 4a. Pooled counts: admissible candidates covering L           [text + last col. of
#     (all four tau pooled)                                      real_data_intervals]
pooled_counts <- horizons %>%
  left_join(tab_K %>% group_by(country, case) %>%
              summarise(n_admissible = n(),
                        n_cover_incl = sum(covers_incl),
                        n_cover_strict = sum(covers_strict),
                        .groups = "drop"),
            by = c("country", "case")) %>%
  mutate(across(c(n_admissible, n_cover_incl, n_cover_strict),
                ~ ifelse(is.na(.x), 0L, .x)))
show("Pooled counts of admissible candidates covering L", pooled_counts)
cat("\nBefore L: ",
    sum(pooled_counts$n_cover_incl[pooled_counts$case == "before"]), " of ",
    sum(pooled_counts$n_admissible[pooled_counts$case == "before"]),
    " candidates contain L (structurally impossible for T < L)\n", sep = "")
save_csv(pooled_counts, "pooled_counts.csv")

# 4b. Illustrative candidates: narrowest k >= 1, width <= cap, covering L;
#     ties: smaller k, then smaller tau. Uses L: NOT operational.
#                                                    [real_data_intervals]
illustrative <- tab_K %>%
  filter(case != "before", k >= 1, covers_incl,
         width <= MAX_WIDTH_ILLUSTRATIVE) %>%
  arrange(country, case, width, k, tau) %>%
  group_by(country, case) %>% slice(1) %>% ungroup() %>%
  left_join(pooled_counts %>% select(country, case, n_admissible, n_cover_incl),
            by = c("country", "case")) %>%
  select(country, case, T_current, tau, k, hat_l_current,
         T_entry, T_exit, width, n_cover_incl, n_admissible)
show("Illustrative candidates (uses L; not operational)", illustrative)
save_csv(illustrative, "illustrative_candidates.csv")

# 4c. Operational rule  k_hat(T) = min K(T)                  [real_data_minK]
minK <- tab_K %>%
  arrange(country, case, tau, k) %>%
  group_by(country, case, tau) %>% slice(1) %>% ungroup() %>%
  select(country, case, T_current, tau, k, hat_l_current,
         T_entry, T_exit, width, covers_incl, covers_strict)

minK_full <- design %>%
  left_join(minK, by = c("country", "case", "T_current", "tau")) %>%
  arrange(country, case, tau)           # NA row = K(T) empty (dash in table)
show("Operational rule min K(T) (NA = K(T) empty)",
     minK_full %>% filter(case != "before"))
save_csv(minK_full, "operational_minK.csv")

operational_summary <- minK_full %>%
  filter(case != "before") %>%
  group_by(case) %>%
  summarise(n_selected = sum(!is.na(k)),
            n_covering = sum(covers_incl, na.rm = TRUE),
            n_empty    = sum(is.na(k)), .groups = "drop")
show("Operational rule: coverage of L (text: 4 of 9 and 8 of 12)",
     operational_summary)


############################################################
# 5. COVERAGE OF L BY CANDIDATE INDEX k  (all horizons T > L)
#    Manuscript reports k = 1, 2, 3, >=4 (k = 0 is shown here too).
############################################################

coverage_by_k <- adm %>%
  mutate(L = unname(L_obs[country])) %>%
  filter(T_current > L) %>%
  mutate(covers = T_entry <= L & L <= T_exit,
         k_group = pmin(k, 4)) %>%
  group_by(country, k_group) %>%
  summarise(n = n(), prop_covering = mean(covers), .groups = "drop")
show("Share of admissible candidates covering L, by k (k_group 4 = k>=4)",
     round_df(coverage_by_k, 2, skip = c("k_group", "n")))
save_csv(coverage_by_k, "coverage_by_k.csv")


############################################################
# 6. SELF-CHECK AGAINST THE MANUSCRIPT VALUES
############################################################

expected_pooled <- data.frame(
  country = rep(c("China", "France", "Germany"), each = 2),
  case    = rep(c("after_near", "after_far"), times = 3),
  n_admissible = c(5, 7, 3, 9, 9, 9),
  n_cover_incl = c(3, 5, 1, 6, 6, 6),
  stringsAsFactors = FALSE)

chk <- expected_pooled %>%
  inner_join(pooled_counts, by = c("country", "case"),
             suffix = c("_ms", "_new"))
ok_pooled <- all(chk$n_admissible_ms == chk$n_admissible_new &
                   chk$n_cover_incl_ms == chk$n_cover_incl_new)

ok_oper <- with(operational_summary,
                all(n_selected[case == "after_near"] == 9,
                    n_covering[case == "after_near"] == 4,
                    n_selected[case == "after_far"]  == 12,
                    n_covering[case == "after_far"]  == 8))

cat("\n==== Self-check against manuscript ====\n")
cat("Pooled counts (Table real_data_intervals, last column): ",
    ifelse(ok_pooled, "OK", "MISMATCH"), "\n", sep = "")
cat("Operational rule counts (4/9 and 8/12):                 ",
    ifelse(ok_oper, "OK", "MISMATCH"), "\n", sep = "")


############################################################
# 7. EXTRA DIAGNOSTICS (not reported in the manuscript)
############################################################

if (RUN_EXTRA_DIAGNOSTICS) {
  
  # 7a. Union interval over ALL candidates with a valid entry time,
  #     for horizons T >= L:  [min_k T_entry, max_k T_exit].
  #     first_T_covering = earliest T >= L from which the union contains
  #     L at every later horizon (corrected: the earlier version simply
  #     returned the first horizon T >= L).
  union_check <- cands %>%
    filter(!is.na(T_entry)) %>%
    mutate(L = unname(L_obs[country])) %>%
    filter(T_current >= L) %>%
    group_by(country, tau, T_current, L) %>%
    summarise(lo = min(T_entry), hi = max(T_exit), .groups = "drop") %>%
    mutate(width = hi - lo, covers = lo <= L & L <= hi) %>%
    arrange(country, tau, T_current) %>%
    group_by(country, tau) %>%
    summarise(prop_covered = mean(covers),
              mean_width   = mean(width),
              first_T_covering = {
                ok <- rev(cumprod(rev(covers))) == 1
                if (any(ok)) min(T_current[ok]) else NA_real_
              },
              .groups = "drop")
  show("Union interval over all candidates (T >= L)", union_check)
  save_csv(union_check, "extra_union_check.csv")
  
  # 7b. Envelope of the ADMISSIBLE candidates:
  #     [min T_entry, max T_exit] over k in K(T), at horizons where K(T)
  #     is non-empty.
  envelope <- adm %>%
    mutate(L = unname(L_obs[country])) %>%
    group_by(country, tau, T_current, L) %>%
    summarise(entry_min = min(T_entry), exit_max = max(T_exit),
              .groups = "drop") %>%
    mutate(covers_L = entry_min <= L & L <= exit_max)
  
  show("Envelope covers L: share of horizons with K(T) non-empty",
       envelope %>% group_by(country, tau) %>%
         summarise(share = mean(covers_L), .groups = "drop"))
  show("Envelope covers L: share over horizons T >= L",
       envelope %>% filter(T_current >= L) %>%
         group_by(country, tau) %>%
         summarise(share = mean(covers_L), .groups = "drop"))
  show("Envelope at the terminal horizon",
       envelope %>% group_by(country, tau) %>%
         filter(T_current == max(T_current)) %>% ungroup())
  save_csv(envelope, "extra_envelope.csv")
}


############################################################
# 8. RECORD THE SESSION
############################################################

writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nAll tables written to ", normalizePath(OUT_DIR), "\n", sep = "")
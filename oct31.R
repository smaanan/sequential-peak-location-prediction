############################################################
# REAL DATA: LINK BETWEEN T_stab AND T_entry(eps), T_exit(eps')
#
# Tail oscillation of the estimated peak-location path, at the
# terminal horizon h:
#     omega_l(t, h) = max_{t <= u, u' <= h} | l_hat(u) - l_hat(u') |
#
# Entry and exit times (first paper, window extended to the horizon):
#     T_entry(eps)  = min{ t : omega_l(t, h) <= eps }
#     T_exit(eps')  = min{ t > T_entry(eps) : omega_l(t, h) <= eps' },
#     eps > eps'.
# T_stab(eps) of the manuscript is exactly T_entry(eps).
#
# Data-driven thresholds: for a candidate k of the two-stage
# construction at T = h, take
#     eps = eps_k   and   eps' = eps'_k.
#
# For each country and tau, and for two choices of k
#   - "retrospective": narrowest candidate containing L
#                      (ties: smaller k), as in the manuscript tables;
#   - "operational"  : k_hat(h) = min K(h);
# the script reports
#   eps_k, eps'_k, the two-stage window [T_entry^(k), T_exit^(k)],
#   T_entry(eps_k), T_exit(eps'_k), the former T_stab (eps = 3, 5, 15),
#   T_entry(0) (= T*, start of the final plateau), and whether the
#   windows contain L.
#
# It also checks the exact "ladder" identity, for every candidate k:
#     T_entry(eps'_k) = T_exit^(k+1) + 1.
#
# Requires:
#   - oct_09.R (only its function definitions, sections 0-11,
#     are loaded: the simulation part is NOT executed)
#   - {china,france,germany}_path_tau{0.25,0.50,0.75,0.95}.rds
#     (data frames with columns T and hat.l)
# Run in a fresh R session (do not load fixed_h_core.R).
#
# Outputs:
#   real_data_tentry_eps_table.csv
#   real_data_tentry_eps_table.tex   (rows for the manuscript)
#   real_data_ladder.csv             (eps'_k, T_exit^(k), identity check)
############################################################


############################################################
# 0. LOAD ONLY THE FUNCTION DEFINITIONS OF oct_09.R
############################################################

src <- readLines("oct_09.R")
cut <- grep("^# 12\\. SIMULATION PARAMETERS", src)

if (length(cut) != 1L) {
  stop("Could not find the start of Section 12 in oct_09.R.")
}

eval(parse(text = src[seq_len(cut - 1L)]), envir = globalenv())

stopifnot(exists("run_entry_exit_at_T"))


############################################################
# 1. SETTINGS
############################################################

countries <- c("China", "France", "Germany")
taus      <- c(0.25, 0.50, 0.75, 0.95)

# Observed peak locations (manuscript)
L.values <- c(China = 22, France = 78, Germany = 61)

# Thresholds previously used for T_stab in the manuscript
EPS_OLD  <- c(China = 3, France = 5, Germany = 15)

TOL <- 1e-12


############################################################
# 2. TAIL OSCILLATION, T_entry(eps), T_exit(eps')
############################################################

# omega_l(t, h) for every observation horizon t of the path
tail_oscillation <- function(T.values, hat.l) {
  o <- order(T.values)
  T.values <- T.values[o]
  hat.l    <- hat.l[o]
  rng <- rev(cummax(rev(hat.l))) - rev(cummin(rev(hat.l)))
  data.frame(t = T.values, omega = rng)
}

# T_entry(eps) = min{t : omega(t, h) <= eps}
t_entry_eps <- function(om, eps) {
  ok <- om$t[om$omega <= eps + TOL]
  if (length(ok) > 0L) min(ok) else NA_real_
}

# T_exit(eps') = min{t > T_entry(eps) : omega(t, h) <= eps'}
t_exit_epsp <- function(om, eps, epsp) {
  te <- t_entry_eps(om, eps)
  if (is.na(te)) return(NA_real_)
  ok <- om$t[om$t >= te & om$omega <= epsp + TOL]
  if (length(ok) > 0L) min(ok) else NA_real_
}


############################################################
# 3. ANALYSIS FOR ONE COUNTRY AND TAU (terminal horizon T = h)
############################################################

analyse <- function(country, tau) {
  
  filename <- sprintf("%s_path_tau%.2f.rds", tolower(country), tau)
  if (!file.exists(filename)) stop("Missing cached path file: ", filename)
  
  path <- readRDS(filename)
  if (!all(c("T", "hat.l") %in% names(path))) {
    stop("File ", filename, " must contain columns T and hat.l.")
  }
  path <- path[order(path$T), , drop = FALSE]
  
  h <- max(path$T)
  L <- L.values[[country]]
  
  # Stage-2 candidates at T = h (L is not used by the construction)
  s2 <- run_entry_exit_at_T(
    T.values  = path$T,
    hat.l     = path$hat.l,
    T.current = h,
    Data      = country,
    L         = NA_real_
  )$stage2
  s2 <- s2[order(s2$k), , drop = FALSE]
  
  om <- tail_oscillation(path$T, path$hat.l)
  
  # ---- ladder identity: T_entry(eps'_k) = T_exit^(k+1) + 1 ----
  n <- nrow(s2)
  ladder <- data.frame(
    country = country, tau = tau, k = s2$k,
    epsp_k = s2$epsilon_prime, T_exit_k = s2$T_exit,
    T_entry_of_epsp = vapply(s2$epsilon_prime,
                             function(e) t_entry_eps(om, e), numeric(1)),
    stringsAsFactors = FALSE
  )
  ladder$T_exit_next_plus1 <- c(s2$T_exit[-1L] + 1, NA_real_)
  ladder$identity_ok <- ifelse(is.na(ladder$T_exit_next_plus1), NA,
                               ladder$T_entry_of_epsp ==
                                 ladder$T_exit_next_plus1)
  
  # ---- candidate selection ----
  has.entry <- !is.na(s2$T_entry)
  covers    <- has.entry & s2$T_entry <= L & L <= s2$T_exit
  width     <- s2$T_exit - s2$T_entry
  
  i.cov   <- which(covers)
  idx.retro <- if (length(i.cov) > 0L) {
    i.cov[order(width[i.cov], s2$k[i.cov])][1L]
  } else NA_integer_
  
  i.adm   <- which(!is.na(s2$estimate.in.window) & s2$estimate.in.window)
  idx.oper <- if (length(i.adm) > 0L) {
    i.adm[which.min(s2$k[i.adm])]
  } else NA_integer_
  
  make_row <- function(i, rule) {
    if (is.na(i)) return(NULL)
    
    eps  <- s2$epsilon[i]
    epsp <- s2$epsilon_prime[i]
    te   <- t_entry_eps(om, eps)
    tx   <- t_exit_epsp(om, eps, epsp)
    
    data.frame(
      country = country, tau = tau, rule = rule, L = L,
      k       = s2$k[i],
      eps_k   = eps,
      epsp_k  = epsp,
      eps_gt_epsp = eps > epsp,               # strict: eps > eps'
      T_entry_k = s2$T_entry[i],              # two-stage window
      T_exit_k  = s2$T_exit[i],
      T_entry_eps = te,                       # T_entry(eps_k)
      T_exit_epsp = tx,                       # T_exit(eps'_k)
      T_stab_old  = t_entry_eps(om, EPS_OLD[[country]]),
      T_star      = t_entry_eps(om, 0),       # start of the final plateau
      twostage_has_L = s2$T_entry[i] <= L & L <= s2$T_exit[i],
      hybrid_has_L   = !is.na(te) & !is.na(tx) & te <= L & L <= tx,
      same_entry = !is.na(te) & te == s2$T_entry[i],
      same_exit  = !is.na(tx) & tx == s2$T_exit[i],
      stringsAsFactors = FALSE
    )
  }
  
  list(rows = rbind(make_row(idx.retro, "retrospective"),
                    make_row(idx.oper,  "operational")),
       ladder = ladder)
}


############################################################
# 4. RUN FOR ALL COUNTRIES AND TAU
############################################################

results <- lapply(countries, function(cn)
  lapply(taus, function(tt) analyse(cn, tt)))

tab    <- do.call(rbind, lapply(results, function(r)
  do.call(rbind, lapply(r, `[[`, "rows"))))
ladder <- do.call(rbind, lapply(results, function(r)
  do.call(rbind, lapply(r, `[[`, "ladder"))))
rownames(tab) <- NULL
rownames(ladder) <- NULL

write.csv(tab,    "real_data_tentry_eps_table.csv", row.names = FALSE)
write.csv(ladder, "real_data_ladder.csv",           row.names = FALSE)


############################################################
# 5. REPORT
############################################################

options(width = 200)

for (rule in c("retrospective", "operational")) {
  cat("\n==== Selection rule:", rule, "====\n")
  d <- tab[tab$rule == rule,
           c("country", "tau", "k", "eps_k", "epsp_k", "T_entry_k",
             "T_exit_k", "T_entry_eps", "T_exit_epsp", "T_stab_old",
             "T_star", "L", "twostage_has_L", "hybrid_has_L")]
  print(d, row.names = FALSE)
}

cat("\n==== Ladder identity T_entry(eps'_k) = T_exit^(k+1) + 1 ====\n")
cat("checked:", sum(!is.na(ladder$identity_ok)),
    " violations:", sum(ladder$identity_ok == FALSE, na.rm = TRUE), "\n")

cat("\n==== Summary (retrospective rule) ====\n")
r <- tab[tab$rule == "retrospective", ]
cat("eps_k > eps'_k (strict)      :", sum(r$eps_gt_epsp), "of", nrow(r), "\n")
cat("T_entry(eps_k) = T_entry^(k) :", sum(r$same_entry), "of", nrow(r), "\n")
cat("T_exit(eps'_k) = T_exit^(k)  :", sum(r$same_exit),  "of", nrow(r), "\n")
cat("two-stage window contains L  :", sum(r$twostage_has_L), "of", nrow(r), "\n")
cat("[T_entry(eps_k), T_exit(eps'_k)] contains L:",
    sum(r$hybrid_has_L), "of", nrow(r), "\n")


############################################################
# 6. LaTeX ROWS FOR THE MANUSCRIPT
#    columns: country & tau & k & eps_k & eps'_k &
#             [T_entry^(k), T_exit^(k)] & [T_entry(eps_k), T_exit(eps'_k)]
#             & T_stab(old) & L
############################################################

fmt <- function(x) ifelse(is.na(x), "---", as.character(x))

latex_rows <- function(d) {
  sprintf("%s & %.2f & %s & %s & %s & [%s, %s] & [%s, %s] & %s & %s \\\\",
          d$country, d$tau, fmt(d$k), fmt(d$eps_k), fmt(d$epsp_k),
          fmt(d$T_entry_k), fmt(d$T_exit_k),
          fmt(d$T_entry_eps), fmt(d$T_exit_epsp),
          fmt(d$T_stab_old), fmt(d$L))
}

tex <- c("% retrospective selection",
         latex_rows(tab[tab$rule == "retrospective", ]),
         "",
         "% operational selection k_hat(h) = min K(h)",
         latex_rows(tab[tab$rule == "operational", ]))

writeLines(tex, "real_data_tentry_eps_table.tex")

cat("\nSaved: real_data_tentry_eps_table.csv, real_data_tentry_eps_table.tex,",
    "real_data_ladder.csv\n")

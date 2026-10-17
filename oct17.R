L_obs <- c(China   = which.max(I_china[1:35]),
           France  = which.max(I_france[1:89]),
           Germany = which.max(I_germany[1:75]))
print(L_obs)          # replace by the manuscript values if they differ

show_K <- function(res, name, L, offsets = c(near = 3, far = 10)) {
  cc <- res$cands
  Tgrid <- sort(unique(cc$T_current))
  pick <- function(target) Tgrid[which.min(abs(Tgrid - target))]
  Ts <- c(before     = pick(L - offsets["near"]),
          after_near = pick(L + offsets["near"]),
          after_far  = pick(L + offsets["far"]))
  Ts <- pmin(Ts, max(Tgrid))
  out <- do.call(rbind, lapply(names(Ts), function(lbl) {
    d <- cc[cc$T_current == Ts[lbl] & !is.na(cc$estimate.in.window) &
              cc$estimate.in.window, ]
    if (nrow(d) == 0) return(data.frame(country = name, case = lbl,
                                        T = Ts[lbl], tau = NA, k = NA, hat_l = NA, entry = NA, exit = NA,
                                        width = NA, L = L, covers_incl = NA, covers_strict = NA))
    data.frame(country = name, case = lbl, T = Ts[lbl], tau = d$tau,
               k = d$k, hat_l = d$hat_l_current, entry = d$T_entry,
               exit = d$T_exit, width = d$T_exit - d$T_entry, L = L,
               covers_incl = d$T_entry <= L & L <= d$T_exit,
               covers_strict = d$T_entry < L & L < d$T_exit)
  }))
  out[order(out$case, out$tau, out$k), ]
}

tabK <- rbind(
  show_K(res_china,   "China",   L_obs["China"]),
  show_K(res_france,  "France",  L_obs["France"]),
  show_K(res_germany, "Germany", L_obs["Germany"]))
print(tabK, row.names = FALSE)
write.csv(tabK, "real_data_K_intervals.csv", row.names = FALSE)

# Share of k in K(T) whose interval covers L, over all T > L
cc_all <- do.call(rbind, lapply(
  list(China = res_china, France = res_france, Germany = res_germany),
  function(r) r$cands[!is.na(r$cands$estimate.in.window) &
                        r$cands$estimate.in.window, ]))


cc_all$country <- rep(c("China", "France", "Germany"),
                      times = sapply(list(res_china, res_france, res_germany),
                                     function(r) sum(!is.na(r$cands$estimate.in.window) &
                                                       r$cands$estimate.in.window)))
cc_all$L <- L_obs[cc_all$country]
a <- cc_all[cc_all$T_current > cc_all$L & !is.na(cc_all$T_entry), ]
a$cov <- a$T_entry <= a$L & a$L <= a$T_exit
a$kk <- pmin(a$k, 4)
print(aggregate(cov ~ country + kk, a, function(x) c(n = length(x), p = mean(x))))

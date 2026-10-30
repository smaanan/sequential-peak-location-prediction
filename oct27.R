############################################################
# REAL DATA: ENTRY / EXIT CURVES, THREE COUNTRIES IN ONE FIGURE
#
# Layout: 3 rows (China, France, Germany) x 4 columns
#         (tau = 0.25, 0.50, 0.75, 0.95).
#
#   T <  L :  T -> T_entry^(0)(T)
#   T >= L :  T -> T_entry^(1)(T)  and  T -> T_exit^(1)(T)
#
# The curves before L (k = 0) and from L on (k = 1) refer to
# different candidates; they are plotted whether or not the
# candidates belong to K(T). Each row has its own y-axis range
# (the three countries have very different horizon ranges).
#
# Input (written by real_data_epsilon_windows_v2.R):
#   real_data_epsilon_windows.csv
#
# Output:
#   fig_real_entry_exit_3x4.pdf
############################################################


############################################################
# 1. SETTINGS AND DATA
############################################################

csv <- "real_data_epsilon_windows.csv"
if (!file.exists(csv)) {
  stop("File ", csv, " not found. Run real_data_epsilon_windows_v2.R first.")
}

curves <- read.csv(csv)

needed <- c("country", "tau", "L", "T", "entry0", "entry1", "exit1")
stopifnot(all(needed %in% names(curves)))

countries <- c("China", "France", "Germany")
taus      <- c(0.25, 0.50, 0.75, 0.95)

# Dashed reference lines at T = L (vertical) and at the value L
# (horizontal). Set to FALSE to draw strictly the three curves.
SHOW_L_REF <- TRUE

col.entry0 <- "black"
col.window <- "darkgreen"

# Series actually plotted
curves$entry0_plot <- ifelse(curves$T <  curves$L, curves$entry0, NA_real_)
curves$entry1_plot <- ifelse(curves$T >= curves$L, curves$entry1, NA_real_)
curves$exit1_plot  <- ifelse(curves$T >= curves$L, curves$exit1,  NA_real_)

y_range <- function(...) {
  x <- c(...)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(0, 1))
  r <- range(x)
  if (r[1] == r[2]) r <- r + c(-1, 1)
  r + c(-1, 1) * 0.04 * diff(r)
}


############################################################
# 2. FIGURE
############################################################

graphics.off()

pdf("fig_real_entry_exit_3x4.pdf", width = 12, height = 12)

par(mfrow = c(3, 4),
    oma = c(3.5, 0.5, 1.5, 0.5),
    mar = c(3.6, 4.6, 2.2, 0.8))

# legend content
leg <- expression(T[entry]^(0) * (T) ~~ "(T < L)",
                  T[entry]^(1) * (T) ~~ "(T >= L)",
                  T[exit]^(1)  * (T) ~~ "(T >= L)")
col    <- c(col.entry0, col.window, col.window)
pch    <- c(24, 24, 25)
pt.bg  <- c(col.entry0, col.window, col.window)
lty    <- c(1, 1, 2)

if (SHOW_L_REF) {
  leg   <- c(leg, expression("T = L", "value L"))
  col   <- c(col, "grey40", "grey40")
  pch   <- c(pch, NA, NA)
  pt.bg <- c(pt.bg, NA, NA)
  lty   <- c(lty, 2, 3)
}

for (ci in seq_along(countries)) {
  
  country <- countries[ci]
  dc  <- curves[curves$country == country, ]
  L.c <- unique(dc$L)
  stopifnot(length(L.c) == 1L)
  
  xr <- range(dc$T)
  yr <- y_range(dc$entry0_plot, dc$entry1_plot, dc$exit1_plot,
                if (SHOW_L_REF) L.c)
  
  for (ti in seq_along(taus)) {
    
    tau <- taus[ti]
    d <- dc[abs(dc$tau - tau) < 1e-8, ]
    d <- d[order(d$T), ]
    
    plot(NA, xlim = xr, ylim = yr,
         xlab = if (ci == length(countries)) expression(T) else "",
         ylab = if (ti == 1) paste0(country, ": time index") else "",
         main = if (ci == 1) bquote(tau == .(tau)) else "",
         font.lab = 2)
    
    if (SHOW_L_REF) {
      abline(v = L.c, lty = 2, lwd = 1.3, col = "grey40")
      abline(h = L.c, lty = 3, lwd = 1.3, col = "grey40")
    }
    
    lines(d$T, d$entry0_plot, type = "o", col = col.entry0,
          pch = 24, bg = col.entry0, cex = 0.6)
    lines(d$T, d$entry1_plot, type = "o", col = col.window,
          pch = 24, bg = col.window, cex = 0.6)
    lines(d$T, d$exit1_plot,  type = "o", col = col.window,
          pch = 25, bg = col.window, cex = 0.6, lty = 2)
    
    # legend inside the uppermost subplot only
    if (ci == 1) {
      legend("topleft",
             legend = leg,
             col = col, pch = pch, pt.bg = pt.bg, lty = lty,
             bty = "o", cex = 0.75)
    }
  }
}

dev.off()

cat("Saved: fig_real_entry_exit_3x4.pdf\n")

############################################################
# REAL DATA: ENTRY / EXIT CURVES, THREE COUNTRIES IN ONE FIGURE
#
# Layout: 3 rows (China, France, Germany) x 4 columns
#         (tau = 0.25, 0.50, 0.75, 0.95).
#
#   all T :  T -> T_entry^(0)(T)
#   T >= L:  T -> T_entry^(1)(T)  and  T -> T_exit^(1)(T)
#
# T_entry^(0)(T) is drawn over the whole horizon range (below and
# above L). The k = 1 curves are drawn only for T >= L. All curves
# are plotted whether or not the candidates belong to K(T).
# Each row has its own y-axis range.
#
# The PDF is created AT PRINT SIZE (FIG_W x FIG_H inches, small
# pointsize), so that it is included at scale 1 and the fonts stay
# readable. Set FIG_W to the text width of the paper:
#   FIG_W = \the\textwidth / 72.27   (inches)
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

# ---- figure geometry (inches) and base font size (points) ----
# 7 x 5.5 gives panels with a plotting area of about 1.2 x 1.3 in
# (close to square). Use FIG_H = 6.5 for taller panels, or
# FIG_W = 9, FIG_H = 6.2, PS = 9 for a landscape page (sidewaysfigure).
FIG_W <- 7
FIG_H <- 5.5
PS    <- 8

# Panel (row, column) that carries the legend: (1, 1) = China, tau = 0.25.
# Use NULL for no legend.
LEGEND_PANEL <- c(1, 1)
LEGEND_CEX   <- 0.7

# Dashed reference lines at T = L (vertical) and at the value L
# (horizontal). Set to FALSE to draw strictly the three curves.
SHOW_L_REF <- TRUE

col.entry0 <- "black"
col.window <- "darkgreen"

# Series actually plotted
curves$entry0_plot <- curves$entry0                       # all T
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
# 2. LEGEND CONTENT (short labels; "k=1 only for T >= L"
#    is explained in the caption)
############################################################

leg   <- expression(T[entry]^(0), T[entry]^(1), T[exit]^(1))
col   <- c(col.entry0, col.window, col.window)
pch   <- c(24, 24, 25)
pt.bg <- c(col.entry0, col.window, col.window)
lty   <- c(1, 1, 2)

if (SHOW_L_REF) {
  leg   <- c(leg, expression("T = L", "value L"))
  col   <- c(col, "grey40", "grey40")
  pch   <- c(pch, NA, NA)
  pt.bg <- c(pt.bg, NA, NA)
  lty   <- c(lty, 2, 3)
}


############################################################
# 3. FIGURE
############################################################

graphics.off()

pdf("fig_real_entry_exit_3x4.pdf",
    width = FIG_W, height = FIG_H, pointsize = PS)

# margins in text lines (1 line = 1.2 * PS points)
par(mfrow = c(3, 4),
    oma = c(0, 0, 0, 0),
    mar = c(2.6, 3.4, 1.4, 0.4),
    mgp = c(1.9, 0.6, 0),
    tcl = -0.3,
    cex.main = 1)

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
      abline(v = L.c, lty = 2, lwd = 1, col = "grey40")
      abline(h = L.c, lty = 3, lwd = 1, col = "grey40")
    }
    
    lines(d$T, d$entry0_plot, type = "o", col = col.entry0,
          pch = 24, bg = col.entry0, cex = 0.5, lwd = 0.8)
    lines(d$T, d$entry1_plot, type = "o", col = col.window,
          pch = 24, bg = col.window, cex = 0.5, lwd = 0.8)
    lines(d$T, d$exit1_plot,  type = "o", col = col.window,
          pch = 25, bg = col.window, cex = 0.5, lwd = 0.8, lty = 2)
    
    # legend inside ONE panel only
    if (!is.null(LEGEND_PANEL) && ci == LEGEND_PANEL[1] &&
        ti == LEGEND_PANEL[2]) {
      legend("topleft", legend = leg,
             col = col, pch = pch, pt.bg = pt.bg, lty = lty,
             bty = "o", bg = "white", box.lwd = 0.5,
             cex = LEGEND_CEX, pt.cex = 0.7,
             seg.len = 1.2, x.intersp = 0.5, y.intersp = 0.9,
             inset = 0.02)
    }
  }
}

dev.off()

cat("Saved: fig_real_entry_exit_3x4.pdf (", FIG_W, "x", FIG_H, "in )\n")
############################################################
# REAL DATA: ENTRY / EXIT CURVES
# Generate either:
#   - one 3x4 combined figure across all countries and taus
#   - one 2x2 figure per country
# or both, depending on the generate vector below.
############################################################

# ---- choose output(s) ----
generate <- c("combined", "by_country")
# valid values: "combined", "by_country"

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

# ---- options ----
SHOW_L_REF <- TRUE
col.entry0 <- "black"
col.window <- "darkgreen"

# figure options for combined figure
FIG_W <- 7
FIG_H <- 5.5
PS    <- 8
LEGEND_PANEL <- c(1, 1)
LEGEND_CEX   <- 0.7

# Set to TRUE to draw the reference lines at T = L and value L
curves$entry0_plot <- curves$entry0
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

# legend for combined figure
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
# 2. PLOTTING FUNCTIONS
############################################################

draw_combined_figure <- function(curves, countries, taus) {
  graphics.off()
  
  pdf("fig_real_entry_exit_3x4.pdf",
      width = FIG_W, height = FIG_H, pointsize = PS)
  
  par(mfrow = c(3, 4),
      oma = c(0, 0, 0, 0),
      mar = c(2.6, 3.4, 1.4, 0.4),
      mgp = c(1.9, 0.6, 0),
      tcl = -0.3,
      cex.main = 1)
  
  for (ci in seq_along(countries)) {
    country <- countries[ci]
    dc <- curves[curves$country == country, ]
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
      lines(d$T, d$exit1_plot, type = "o", col = col.window,
            pch = 25, bg = col.window, cex = 0.5, lwd = 0.8, lty = 2)
      
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
}

draw_country_figures <- function(curves, countries, taus) {
  graphics.off()
  files <- character(0)
  
  for (country in countries) {
    dc <- curves[curves$country == country, ]
    L.c <- unique(dc$L)
    stopifnot(length(L.c) == 1L)
    
    xr <- range(dc$T)
    yr <- y_range(dc$entry0_plot, dc$entry1_plot, dc$exit1_plot,
                  if (SHOW_L_REF) L.c)
    
    f <- sprintf("fig_real_entry_exit_curves_%s.pdf", tolower(country))
    files <- c(files, f)
    
    pdf(f, width = 10, height = 8)
    par(mfrow = c(2, 2), oma = c(0, 0, 2.5, 0), mar = c(4.2, 4.8, 3, 1))
    
    for (tau in taus) {
      d <- dc[abs(dc$tau - tau) < 1e-8, ]
      d <- d[order(d$T), ]
      
      plot(NA, xlim = xr, ylim = yr,
           xlab = expression(T), ylab = "time index",
           main = bquote(tau == .(tau)))
      
      if (SHOW_L_REF) {
        abline(v = L.c, lty = 2, lwd = 1.3, col = "grey40")
        abline(h = L.c, lty = 3, lwd = 1.3, col = "grey40")
      }
      
      lines(d$T, d$entry0_plot, type = "o", col = col.entry0,
            pch = 24, bg = col.entry0, cex = 0.6)
      lines(d$T, d$entry1_plot, type = "o", col = col.window,
            pch = 24, bg = col.window, cex = 0.6)
      lines(d$T, d$exit1_plot, type = "o", col = col.window,
            pch = 25, bg = col.window, cex = 0.6, lty = 2)
      
      legend("topleft",
             legend = expression(T[entry]^(0) * (T) ~~ "(all T)",
                                 T[entry]^(1) * (T) ~~ "(T >= L)",
                                 T[exit]^(1)  * (T) ~~ "(T >= L)"),
             col = c(col.entry0, col.window, col.window),
             pch = c(24, 24, 25),
             pt.bg = c(col.entry0, col.window, col.window),
             lty = c(1, 1, 2),
             bty = "o", bg = "white", cex = 0.7)
    }
    
    mtext(sprintf("%s (L = %s)", country, L.c),
          outer = TRUE, line = 0.5, font = 2)
    
    dev.off()
    cat("Figure drawn for", country, "\n")
  }
  
  cat("\nFinished. Files written:\n")
  cat(paste0("  ", files, "\n"), sep = "")
}

############################################################
# 3. RUN
############################################################

for (mode in generate) {
  if (mode == "combined") {
    draw_combined_figure(curves, countries, taus)
  } else if (mode == "by_country") {
    draw_country_figures(curves, countries, taus)
  } else {
    stop("Unknown mode: ", mode, ". Use 'combined' and/or 'by_country'.")
  }
}
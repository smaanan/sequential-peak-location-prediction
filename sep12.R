###############################################################
## GLOBAL PARAMETERS
###############################################################

A <- 1
L <- 22
S <- 10

## Quantile levels considered in the new manuscript
tau.values <- c(
  0.25,
  0.50,
  0.75,
  0.95
)

## SNR values corresponding to the accepted-article
## simulation design
SNR.values <- c(
  1000,
  100,
  40
)

## Common simulation parameters
h <- 30
seed <- 280


###############################################################
## Function to generate the complete 3 x 2 figure
## for one value of tau
###############################################################

run_tau_experiment <- function(
    tau.value,
    A = 1,
    L = 22,
    S = 10,
    h = 30,
    seed = 280,
    SNR.values = c(1000, 100, 40),
    filename
) {
  
  
  #############################################################
  ## Open figure
  ##
  ## Layout:
  ##
  ##   row 1: amplitude estimate | amplitude oscillation
  ##   row 2: location estimate  | location oscillation
  ##   row 3: width estimate     | width oscillation
  #############################################################
  
  png(
    filename,
    width = 8,
    height = 7,
    units = "in",
    res = 300
  )
  
  par(
    mfrow = c(3,2),
    mar = c(4,4,3,1),
    mgp = c(2.2,0.6,0)
  )
  
  
  #############################################################
  ## Run the three SNR experiments
  #############################################################
  
  results <- vector(
    "list",
    length(SNR.values)
  )
  
  names(results) <- paste0(
    "SNR_",
    SNR.values
  )
  
  
  for (j in seq_along(SNR.values)) {
    
    SNR.value <- SNR.values[j]
    
    results[[j]] <-
      estimate_peak_curves(
        A = A,
        L = L,
        S = S,
        h = h,
        SNR = SNR.value,
        tau = tau.value,
        seed = seed
      )
  }
  
  
  #############################################################
  ## Plot amplitude curves
  ##
  ## Large SNR
  #############################################################
  
  plot(
    results[[1]]$T,
    results[[1]]$hat.a,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(hat(a)(T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 1000"
    )
  )
  
  abline(
    h = A,
    lty = 2
  )
  
  abline(
    v = L,
    lty = 2,
    col = "red"
  )
  
  
  plot(
    results[[1]]$T,
    results[[1]]$omega.a,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(omega[a](T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 1000"
    )
  )
  
  abline(
    v = L,
    lty = 2,
    col = "red"
  )
  
  
  #############################################################
  ## Plot location curves
  #############################################################
  
  plot(
    results[[1]]$T,
    results[[1]]$hat.l,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(hat(l)(T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 1000"
    )
  )
  
  abline(
    h = L,
    lty = 2
  )
  
  abline(
    v = L,
    lty = 2,
    col = "red"
  )
  
  
  plot(
    results[[1]]$T,
    results[[1]]$omega.l,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(omega[l](T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 1000"
    )
  )
  
  abline(
    v = L,
    lty = 2,
    col = "red"
  )
  
  
  #############################################################
  ## Plot width curves
  #############################################################
  
  plot(
    results[[1]]$T,
    results[[1]]$hat.s,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(hat(s)(T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 1000"
    )
  )
  
  abline(
    h = S,
    lty = 2
  )
  
  abline(
    v = L,
    lty = 2,
    col = "red"
  )
  
  
  plot(
    results[[1]]$T,
    results[[1]]$omega.s,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(omega[s](T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 1000"
    )
  )
  
  abline(
    v = L,
    lty = 2,
    col = "red"
  )
  
  
  #############################################################
  ## The above six panels are for SNR = 1000.
  ##
  ## For the medium and small SNR cases, create separate
  ## figures below.
  #############################################################
  
  dev.off()
  
  
  #############################################################
  ## Medium SNR figure
  #############################################################
  
  filename.medium <- sub(
    "\\.png$",
    "_SNR100.png",
    filename
  )
  
  png(
    filename.medium,
    width = 8,
    height = 7,
    units = "in",
    res = 300
  )
  
  par(
    mfrow = c(3,2),
    mar = c(4,4,3,1),
    mgp = c(2.2,0.6,0)
  )
  
  
  res.medium <- results[[2]]
  
  
  plot(
    res.medium$T,
    res.medium$hat.a,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(hat(a)(T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 100"
    )
  )
  
  abline(h = A, lty = 2)
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.medium$T,
    res.medium$omega.a,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(omega[a](T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 100"
    )
  )
  
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.medium$T,
    res.medium$hat.l,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(hat(l)(T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 100"
    )
  )
  
  abline(h = L, lty = 2)
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.medium$T,
    res.medium$omega.l,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(omega[l](T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 100"
    )
  )
  
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.medium$T,
    res.medium$hat.s,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(hat(s)(T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 100"
    )
  )
  
  abline(h = S, lty = 2)
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.medium$T,
    res.medium$omega.s,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(omega[s](T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 100"
    )
  )
  
  abline(v = L, lty = 2, col = "red")
  
  
  dev.off()
  
  
  #############################################################
  ## Small SNR figure
  #############################################################
  
  filename.small <- sub(
    "\\.png$",
    "_SNR40.png",
    filename
  )
  
  png(
    filename.small,
    width = 8,
    height = 7,
    units = "in",
    res = 300
  )
  
  par(
    mfrow = c(3,2),
    mar = c(4,4,3,1),
    mgp = c(2.2,0.6,0)
  )
  
  
  res.small <- results[[3]]
  
  
  plot(
    res.small$T,
    res.small$hat.a,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(hat(a)(T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 40"
    )
  )
  
  abline(h = A, lty = 2)
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.small$T,
    res.small$omega.a,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(omega[a](T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 40"
    )
  )
  
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.small$T,
    res.small$hat.l,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(hat(l)(T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 40"
    )
  )
  
  abline(h = L, lty = 2)
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.small$T,
    res.small$omega.l,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(omega[l](T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 40"
    )
  )
  
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.small$T,
    res.small$hat.s,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(hat(s)(T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 40"
    )
  )
  
  abline(h = S, lty = 2)
  abline(v = L, lty = 2, col = "red")
  
  
  plot(
    res.small$T,
    res.small$omega.s,
    type = "b",
    pch = 19,
    lwd = 2,
    xlab = "T",
    ylab = expression(omega[s](T)),
    main = bquote(
      tau == .(tau.value) ~ ", SNR = 40"
    )
  )
  
  abline(v = L, lty = 2, col = "red")
  
  
  dev.off()
  
  
  #############################################################
  ## Return all results
  #############################################################
  
  list(
    tau = tau.value,
    results = results
  )
}



###############################################################
## Run the complete experiment for all four tau values
###############################################################

all.results <- vector(
  "list",
  length(tau.values)
)

names(all.results) <- paste0(
  "tau_",
  tau.values
)


for (i in seq_along(tau.values)) {
  
  tau.value <- tau.values[i]
  
  filename <- paste0(
    "tau_",
    sprintf("%.2f", tau.value),
    ".png"
  )
  
  all.results[[i]] <-
    run_tau_experiment(
      tau.value = tau.value,
      A = A,
      L = L,
      S = S,
      h = h,
      seed = seed,
      SNR.values = SNR.values,
      filename = filename
    )
}



###############################################################
## Display all results
###############################################################

all.results
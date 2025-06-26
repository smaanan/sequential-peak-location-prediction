# --- Libraries ---
library(COVID19)
library(readr)
library(dplyr)
library(ggplot2)
library(zoo)
library(patchwork)  # for combining plots

# --- 1. Load China Data ---
china <- covid19("China", level = 1, start = "2020-01-22", end = "2020-03-22", verbose = FALSE) %>%
  arrange(date) %>%
  mutate(new_cases = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_cases)) %>%
  mutate(country = "China")

# --- 2. Load France Data ---
france <- read_csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", col_types = cols()) %>%
  filter(countriesAndTerritories == "France") %>%
  mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
  filter(date >= as.Date("2020-01-15") & date <= as.Date("2020-05-31")) %>%
  arrange(date) %>%
  mutate(daily_cases = ifelse(cases < 0, NA, cases),
         daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE),
         country = "France") %>%
  select(date, daily_cases, country) %>%
  rename(new_cases = daily_cases)

# --- 3. Load Germany Data ---
germany <- read_csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", col_types = cols()) %>%
  filter(countriesAndTerritories == "Germany") %>%
  mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))) %>%
  filter(date >= as.Date("2020-01-28") & date <= as.Date("2020-07-05")) %>%
  arrange(date) %>%
  mutate(daily_cases = ifelse(cases < 0, NA, cases),
         daily_cases = zoo::na.locf(daily_cases, na.rm = FALSE),
         country = "Germany") %>%
  select(date, daily_cases, country) %>%
  rename(new_cases = daily_cases)


# Save plot as high-quality PDF
pdf("covid_series_1x3.pdf", width = 12, height = 4)

# Set up 1x3 layout
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

# Plot for China
plot(china$date, china$new_cases,
     type = "l", col = "steelblue", lwd = 2,
     main = "China", xlab = "Date", ylab = "Daily Cases")
grid()

# Plot for France
plot(france$date, france$new_cases,
     type = "l", col = "darkred", lwd = 2,
     main = "France", xlab = "Date", ylab = "")
grid()

# Plot for Germany
plot(de$date, de$daily_cases,
     type = "l", col = "darkgreen", lwd = 2,
     main = "Germany", xlab = "Date", ylab = "")
grid()

# Common title
mtext("Daily Confirmed COVID-19 Cases", outer = TRUE, cex = 1.5)

# Close PDF device
dev.off()

cat("Saved: covid_series_1x3.pdf\n")

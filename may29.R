# Load required libraries
library(tidyr)
library(dplyr)
library(COVID19)

x <- covid19(country = "China", level = 1, start = "2020-01-22", end = "2020-04-01", verbose = FALSE)

# Compute daily new confirmed cases
x <- x %>%
  arrange(date) %>%
  mutate(new_confirmed = c(NA, diff(confirmed))) %>%
  filter(!is.na(new_confirmed)) %>%
  select(date, new_confirmed)

# View the result
head(x)
tail(x)

png("china.png", width = 600, height = 400, res = 120)
plot(x$date, x$new_confirmed, type = "l", col = "blue",
     xlab = "Date", ylab = "New Confirmed Cases",
     main = "Daily New COVID-19 Cases in China")
dev.off()
x <- x %>%
  mutate(
    log_cases = log(pmax(new_confirmed, 1)),
    t = row_number()
  )
data_T18 <- x %>% filter(t <= 18)
data_T21 <- x %>% filter(t <= 21)

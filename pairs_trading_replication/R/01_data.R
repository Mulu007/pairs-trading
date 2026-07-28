library(quantmod)
library(dplyr)

tickers <- c("COP", "CVX", "XOM", "XLE", "SPY")

prices_list <- lapply(tickers, function(tk) {
  message("Pulling ", tk)
  x <- getSymbols(tk,
                  src = "yahoo",
                  from = "1999-12-31",
                  to   = "2026-07-26"(),
                  auto.assign = FALSE)
  Ad(x)                      # adjusted close only
})

prices <- do.call(merge, prices_list)
colnames(prices) <- tickers
prices <- na.omit(prices)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
saveRDS(prices, "data/raw/prices_raw.rds")
write.csv(as.data.frame(prices), "data/raw/prices_raw.csv")

dim(prices)
head(prices, 3)
tail(prices, 3)
range(index(prices))
colSums(is.na(prices))

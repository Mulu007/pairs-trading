library(quantmod)
library(ggplot2)
library(tidyr)
library(dplyr)

prices <- readRDS("data/raw/prices_raw.rds")

returns <- diff(log(prices))
returns <- na.omit(returns)
colnames(returns) <- paste0(colnames(prices), "_r")

# dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(returns, "data/processed/returns.rds")

dim(returns) # -> differencing costs one obs therefore 6679


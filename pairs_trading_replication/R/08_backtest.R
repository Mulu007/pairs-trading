# =============================================================
# 10_backtest_fixed.R
# Vectorised backtest of the COP-XOM mean-reversion strategy with
# fixed in-sample OU parameters, replicating the original study.
#
# Strategy:
#   z_t = (spread_t - mu_e) / sigma_eq
#   enter SHORT spread when z >  sd_mult   (spread too high)
#   enter LONG  spread when z < -sd_mult   (spread too low)
#   exit to flat when z crosses back through 0
#   execution lagged by `slip` days to approximate slippage
#
# Reads:  data/processed/spread_fixed.rds
#         data/processed/ou_params.rds
#         data/raw/prices_raw.rds
# Writes: output/tables/backtest_performance.tex
#         output/tables/backtest_optimization.tex
#         output/figures/backtest_equity_curves.png
#         output/figures/backtest_position_spread.png
#         data/processed/backtest_results.rds
# =============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)
library(zoo)

dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

spread_df <- readRDS("data/processed/spread_fixed.rds")
ou        <- readRDS("data/processed/ou_params.rds")
prices    <- readRDS("data/raw/prices_raw.rds")
price_df  <- as.data.frame(prices)

dates  <- spread_df$date
spread <- spread_df$spread
n      <- length(spread)

# Log returns of the two legs, for computing spread-position P&L
r_cop <- c(0, diff(log(price_df$COP)))
r_xom <- c(0, diff(log(price_df$XOM)))

# Benchmark: XLE total return, and buy-and-hold of the pair
r_xle <- c(0, diff(log(price_df$XLE)))

# Period masks
IN_END   <- ou$in_end
OOS_END  <- as.Date("2018-12-31")   # original's out-of-sample end
in_mask  <- dates <= IN_END
oos_mask <- dates > IN_END & dates <= OOS_END
ext_mask <- dates > OOS_END          # the novel extension, 2019-2026
tot_mask <- rep(TRUE, n)

# -------------------------------------------------------------
# Core backtester.
#
# Returns a daily strategy return series for given threshold multiple
# and slippage. The spread position is +1 (long spread: long COP,
# short beta*XOM) or -1 (short spread) or 0 (flat).
#
# Long spread earns   r_cop - beta * r_xom
# Short spread earns  -(r_cop - beta * r_xom)
# -------------------------------------------------------------
run_backtest <- function(sd_mult, slip, beta,
                         mu_e, sigma_eq, max_hold = Inf) {
  z <- (spread - mu_e) / sigma_eq
  
  pos <- numeric(n)          # desired position each day
  cur <- 0
  hold <- 0                  # days in current position
  for (t in seq_len(n)) {
    if (cur == 0) {
      if (z[t] >  sd_mult) { cur <- -1; hold <- 0 }   # short the spread
      else if (z[t] < -sd_mult) { cur <-  1; hold <- 0 }  # long the spread
    } else {
      hold <- hold + 1
      crossed <- (cur == 1 && z[t] >= 0) || (cur == -1 && z[t] <= 0)
      timed_out <- hold >= max_hold
      if (crossed || timed_out) cur <- 0
    }
    pos[t] <- cur
  }
  
  # Execute with slippage: position effective `slip` days after signal
  pos_exec <- c(rep(0, slip), head(pos, n - slip))
  
  spread_ret <- r_cop - beta * r_xom
  strat_ret  <- pos_exec * spread_ret
  
  list(pos = pos_exec, strat_ret = strat_ret, z = z)
}

# -------------------------------------------------------------
# Performance statistics on a return series over a mask.
# -------------------------------------------------------------
perf_stats <- function(ret, mask, ann = 252) {
  r <- ret[mask]
  r <- r[is.finite(r)]
  if (length(r) < 2 || sd(r) == 0) {
    return(c(CAGR = NA, Vol = NA, Sharpe = NA, Sortino = NA,
             MaxDD = NA, Pct_pos = NA))
  }
  cum      <- cumprod(1 + r)
  years    <- length(r) / ann
  cagr     <- cum[length(cum)]^(1 / years) - 1
  vol      <- sd(r) * sqrt(ann)
  sharpe   <- (mean(r) * ann) / vol
  downside <- sqrt(mean(pmin(r, 0)^2)) * sqrt(ann)
  sortino  <- (mean(r) * ann) / downside
  peak     <- cummax(cum)
  maxdd    <- max(1 - cum / peak)
  pct_pos  <- mean(r > 0)
  c(CAGR = cagr, Vol = vol, Sharpe = sharpe, Sortino = sortino,
    MaxDD = maxdd, Pct_pos = pct_pos)
}

# -------------------------------------------------------------
# Grid search over sd multiple and slippage, selecting on in-sample
# Sharpe, exactly as the original does. Reported as its own table.
# -------------------------------------------------------------
sd_grid   <- seq(1, 5, by = 0.5)
slip_grid <- 1:3

grid <- expand.grid(sd_mult = sd_grid, slip = slip_grid)
grid_res <- lapply(seq_len(nrow(grid)), function(i) {
  bt <- run_backtest(grid$sd_mult[i], grid$slip[i],
                     ou$beta_e, ou$mu_e, ou$sigma_eq)
  s <- perf_stats(bt$strat_ret, in_mask)
  data.frame(sd_mult = grid$sd_mult[i], slip = grid$slip[i],
             in_sharpe = s["Sharpe"], in_cagr = s["CAGR"],
             in_vol = s["Vol"])
})
grid_res <- do.call(rbind, grid_res)

best <- grid_res[which.max(grid_res$in_sharpe), ]
cat("=== Best in-sample parameters ===\n")
print(best, row.names = FALSE, digits = 4)

# -------------------------------------------------------------
# Run the selected specification and compute stats on every period.
# -------------------------------------------------------------
bt <- run_backtest(best$sd_mult, best$slip, ou$beta_e, ou$mu_e, ou$sigma_eq)

# Benchmarks
bh_ret <- r_cop - ou$beta_e * r_xom   # always-long-spread buy & hold

periods <- list("In-sample" = in_mask, "Out-of-sample" = oos_mask,
                "Extension (2019-2026)" = ext_mask, "Total" = tot_mask)

perf_tbl <- lapply(names(periods), function(pn) {
  m  <- periods[[pn]]
  st <- perf_stats(bt$strat_ret, m)
  bh <- perf_stats(bh_ret, m)
  bm <- perf_stats(r_xle, m)
  data.frame(
    Period = pn,
    Metric = c("Strategy", "Buy \\& hold", "Benchmark (XLE)"),
    CAGR    = c(st["CAGR"],  bh["CAGR"],  bm["CAGR"]),
    Vol     = c(st["Vol"],   bh["Vol"],   bm["Vol"]),
    Sharpe  = c(st["Sharpe"],bh["Sharpe"],bm["Sharpe"]),
    Sortino = c(st["Sortino"],bh["Sortino"],bm["Sortino"]),
    MaxDD   = c(st["MaxDD"], bh["MaxDD"], bm["MaxDD"])
  )
})
perf_tbl <- do.call(rbind, perf_tbl)

cat("\n=== Performance by period ===\n")
print(perf_tbl, row.names = FALSE, digits = 3)

saveRDS(list(bt = bt, best = best, perf = perf_tbl,
             bh_ret = bh_ret, r_xle = r_xle, dates = dates,
             periods = periods), "data/processed/backtest_results.rds")

# -------------------------------------------------------------
# Time-stop diagnostic: force-close after 2 half-lives.
# Isolates whether losses come from bad entries or slow exits.
# -------------------------------------------------------------
max_hold <- round(2 * ou$half_life)
bt_stop  <- run_backtest(best$sd_mult, best$slip, ou$beta_e,
                         ou$mu_e, ou$sigma_eq, max_hold = max_hold)

cat(sprintf("\n=== With time stop at %d days (2 half-lives) ===\n", max_hold))
stop_tbl <- do.call(rbind, lapply(names(periods), function(pn) {
  s <- perf_stats(bt_stop$strat_ret, periods[[pn]])
  data.frame(Period = pn, Sharpe = s["Sharpe"], CAGR = s["CAGR"],
             MaxDD = s["MaxDD"])
}))
print(stop_tbl, row.names = FALSE, digits = 3)

# -------------------------------------------------------------
# Export performance table
# -------------------------------------------------------------
perf_out <- perf_tbl
perf_out$CAGR    <- sprintf("%.2f\\%%", 100 * perf_out$CAGR)
perf_out$Vol     <- sprintf("%.2f\\%%", 100 * perf_out$Vol)
perf_out$Sharpe  <- sprintf("%.3f", perf_out$Sharpe)
perf_out$Sortino <- sprintf("%.3f", perf_out$Sortino)
perf_out$MaxDD   <- sprintf("%.2f\\%%", 100 * perf_out$MaxDD)

print(
  xtable(perf_out,
         caption = paste("Backtest performance of the fixed-parameter COP-XOM",
                         "strategy against buy-and-hold of the spread and the XLE",
                         "benchmark, across the in-sample, out-of-sample, extension",
                         "and total periods. Parameters selected on in-sample",
                         "Sharpe ratio."),
         label = "tab:backtest_perf"),
  file = "output/tables/backtest_performance.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "H",
  sanitize.text.function = identity,
  sanitize.colnames.function = identity
)

# Optimization grid table
grid_out <- grid_res
grid_out$in_sharpe <- sprintf("%.3f", grid_out$in_sharpe)
grid_out$in_cagr   <- sprintf("%.2f\\%%", 100 * grid_out$in_cagr)
grid_out$in_vol    <- sprintf("%.2f\\%%", 100 * grid_out$in_vol)
names(grid_out) <- c("SD multiple", "Slippage (days)",
                     "In-sample Sharpe", "In-sample CAGR", "In-sample vol")

print(
  xtable(grid_out,
         caption = paste("In-sample grid search over the standard deviation",
                         "threshold multiple and slippage. The combination",
                         "maximising the in-sample Sharpe ratio is carried to the",
                         "out-of-sample and extension periods."),
         label = "tab:backtest_opt"),
  file = "output/tables/backtest_optimization.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "H",
  sanitize.text.function = identity,
  sanitize.colnames.function = identity
)

# -------------------------------------------------------------
# Equity curves
# -------------------------------------------------------------
eq_df <- data.frame(
  date     = dates,
  Strategy = cumprod(1 + ifelse(is.finite(bt$strat_ret), bt$strat_ret, 0)),
  Benchmark = cumprod(1 + ifelse(is.finite(r_xle), r_xle, 0)),
  BuyHold  = cumprod(1 + ifelse(is.finite(bh_ret), bh_ret, 0))
) %>%
  pivot_longer(-date, names_to = "series", values_to = "value")

p_eq <- ggplot(eq_df, aes(date, value, colour = series)) +
  geom_vline(xintercept = c(IN_END, OOS_END),
             linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.55) +
  scale_y_log10() +
  scale_colour_manual(values = c("Strategy" = "#C0392B",
                                 "Benchmark" = "#3B4CC0",
                                 "BuyHold" = "#1BA37A")) +
  labs(x = NULL, y = "Cumulative growth (log scale)", colour = NULL,
       title = "Strategy equity curve against benchmark and buy-and-hold",
       subtitle = paste0("Fixed parameters, SD multiple ", best$sd_mult,
                         ", slippage ", best$slip,
                         " day. Dotted lines mark period boundaries.")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("output/figures/backtest_equity_curves.png", p_eq,
       width = 10, height = 5.5, dpi = 300, bg = "white")

# -------------------------------------------------------------
# Spread with position overlay
# -------------------------------------------------------------
pos_df <- data.frame(date = dates, z = bt$z, pos = bt$pos)

p_pos <- ggplot(pos_df, aes(date, z)) +
  geom_hline(yintercept = c(-best$sd_mult, 0, best$sd_mult),
             linetype = c("dashed","solid","dashed"),
             colour = c("#1BA37A","grey40","#C0392B")) +
  geom_line(colour = "#3B4CC0", linewidth = 0.35) +
  labs(x = NULL, y = "Spread z-score",
       title = "Spread z-score with entry thresholds",
       subtitle = paste0("Entry at \u00b1", best$sd_mult,
                         " standard deviations, exit at the mean")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

ggsave("output/figures/backtest_position_spread.png", p_pos,
       width = 10, height = 4.5, dpi = 300, bg = "white")

message("Backtest complete.")

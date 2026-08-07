# =============================================================
# 07_rolling_cointegration.R
# Rolling window Engle and Granger cointegration tests.
#
# A single full sample test returns one verdict averaged over 26 years,
# which conceals any change in the relationship. Re-estimating on
# overlapping subsamples returns a path, which can date a breakdown.
#
# Reads:  data/raw/prices_raw.rds
# Writes: output/figures/rolling_cointegration_pairs.png
#         output/figures/rolling_cointegration_window_robustness.png
#         output/figures/rolling_hedge_ratio.png
#         output/tables/rolling_cointegration_summary.tex
#         output/tables/rolling_cointegration_regime.tex
#         data/processed/rolling_coint.rds
# =============================================================

library(urca)
library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)
library(zoo)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

prices   <- readRDS("data/raw/prices_raw.rds")
price_df <- as.data.frame(prices)

# -------------------------------------------------------------
# Date handling.
#
# The saved xts object stores its index as numeric seconds since epoch
# while carrying a tclass attribute of "Date", so zoo::index returns raw
# numbers rather than Dates. Converting through POSIXct first and then
# dropping the time component recovers the correct dates. Applying
# as.Date with origin = "1970-01-01" directly would treat those seconds
# as day counts and produce dates several million years in the future.
# -------------------------------------------------------------
dates <- as.Date(as.POSIXct(zoo::index(prices),
                            origin = "1970-01-01", tz = "UTC"))

stopifnot(inherits(dates, "Date"), length(dates) == nrow(price_df))
cat("Date range:", format(min(dates)), "to", format(max(dates)),
    "|", length(dates), "observations\n")

# Log prices. Cointegration is tested on logs rather than raw levels so
# that the cointegrating coefficient is interpretable as an elasticity
# and the spread is invariant to the scale of either security.
lp <- log(price_df)

# -------------------------------------------------------------
# Engle and Granger step one on a single window.
#
# Regress log y on log x, extract the residual, and test it for a unit
# root. Rejection of the unit root null indicates a stationary residual,
# which is to say the two series are cointegrated.
#
# The critical values are NOT the standard Dickey and Fuller values.
# Because the residual is estimated rather than observed, the statistic
# has a different limiting distribution and the correct critical values
# are those of Engle and Granger, which are more negative. Using
# standard ADF values would over-reject substantially.
# -------------------------------------------------------------
EG_CV <- c("1pct" = -3.90, "5pct" = -3.34, "10pct" = -3.04)  # n=2, no trend

eg_test <- function(y, x) {
  fit <- lm(y ~ x)
  r   <- residuals(fit)
  
  # No drift or trend, since OLS residuals are mean zero by construction
  adf <- ur.df(r, type = "none", lags = 4, selectlags = "AIC")
  
  list(stat = as.numeric(adf@teststat[1]),
       beta = unname(coef(fit)[2]))
}

# -------------------------------------------------------------
# Roll the test across the sample.
#
# Stepping monthly rather than daily cuts compute by a factor of twenty
# at no cost to interpretation, since adjacent daily windows share all
# but one observation and are almost perfectly correlated.
# -------------------------------------------------------------
rolling_eg <- function(y_name, x_name, data, dates, width = 504, by = 21) {
  
  starts <- seq(1, nrow(data) - width + 1, by = by)
  
  res <- lapply(starts, function(s) {
    idx <- s:(s + width - 1)
    r   <- eg_test(data[[y_name]][idx], data[[x_name]][idx])
    data.frame(
      y          = y_name,
      x          = x_name,
      window_end = dates[s + width - 1],
      adf_stat   = r$stat,
      beta       = r$beta,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, res)
  # rbind on data frames can strip the Date class; restore it explicitly
  out$window_end <- as.Date(out$window_end, origin = "1970-01-01")
  out
}

PAIRS <- list(
  c("COP", "XOM"), c("XOM", "COP"),
  c("CVX", "XOM"), c("XOM", "CVX"),
  c("COP", "CVX"), c("XOM", "XLE"),
  c("CVX", "XLE"), c("COP", "XLE")
)

WIDTH <- 504   # roughly two trading years
STEP  <- 21    # roughly one trading month

message("Rolling ", length(PAIRS), " pairs at width ", WIDTH, " ...")

roll <- do.call(rbind, lapply(PAIRS, function(p) {
  rolling_eg(p[1], p[2], lp, dates, width = WIDTH, by = STEP)
}))

roll <- roll %>%
  mutate(
    pair_label     = paste0(y, " on ", x),
    cointegrated_5 = adf_stat < EG_CV["5pct"],
    cointegrated_1 = adf_stat < EG_CV["1pct"]
  )

# Verify the date column survived before anything depends on it
stopifnot(inherits(roll$window_end, "Date"))
cat("Window end range:", format(min(roll$window_end)), "to",
    format(max(roll$window_end)), "\n")
cat("Windows ending before March 2020:",
    sum(roll$window_end < as.Date("2020-03-01")), "of", nrow(roll), "\n\n")

saveRDS(roll, "data/processed/rolling_coint.rds")

# -------------------------------------------------------------
# Summary
# -------------------------------------------------------------
summ <- roll %>%
  group_by(pair_label) %>%
  summarise(
    n_windows     = n(),
    pct_coint_5   = 100 * mean(cointegrated_5),
    pct_coint_1   = 100 * mean(cointegrated_1),
    last_reject_5 = if (any(cointegrated_5))
      format(max(window_end[cointegrated_5])) else "never",
    mean_beta     = mean(beta),
    sd_beta       = sd(beta),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_coint_5))

cat("=== Rolling cointegration summary, width", WIDTH, "days ===\n")
print(as.data.frame(summ), digits = 4)

print(
  xtable(as.data.frame(
    summ %>% transmute(
      Pair                   = pair_label,
      Windows                = n_windows,
      `Reject at 5\\%`       = pct_coint_5,
      `Reject at 1\\%`       = pct_coint_1,
      `Last rejection`       = last_reject_5,
      `Mean $\\hat{\\beta}$` = mean_beta,
      `SD $\\hat{\\beta}$`   = sd_beta)),
    digits = c(0, 0, 0, 1, 1, 0, 4, 4),
    caption = paste("Rolling Engle and Granger cointegration tests on a", WIDTH,
                    "day window stepped monthly. Percentages give the share of",
                    "windows in which the unit root null is rejected for the",
                    "cointegrating residual, evaluated against Engle and",
                    "Granger critical values. Under the null a correctly sized",
                    "test rejects in 5 percent of windows."),
    label = "tab:rolling_coint"),
  file = "output/tables/rolling_cointegration_summary.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "H",
  sanitize.text.function = identity,
  sanitize.colnames.function = identity
)

# -------------------------------------------------------------
# Regime comparison
# -------------------------------------------------------------
regime <- roll %>%
  mutate(regime = ifelse(window_end < as.Date("2020-03-01"), "pre", "post")) %>%
  group_by(pair_label, regime) %>%
  summarise(pct = 100 * mean(cointegrated_5),
            stat = mean(adf_stat), .groups = "drop") %>%
  pivot_wider(names_from = regime, values_from = c(pct, stat)) %>%
  mutate(pct_change = pct_post - pct_pre) %>%
  arrange(pct_change)

cat("\n=== Rejection rate before and after March 2020 ===\n")
print(as.data.frame(regime), digits = 4)

print(
  xtable(as.data.frame(
    regime %>% transmute(
      Pair                = pair_label,
      `Reject pre \\%`    = pct_pre,
      `Reject post \\%`   = pct_post,
      `Change`            = pct_change,
      `Mean stat pre`     = stat_pre,
      `Mean stat post`    = stat_post)),
    digits = c(0, 0, 1, 1, 1, 3, 3),
    caption = paste("Share of rolling windows rejecting the unit root null",
                    "before and after March 2020, with the mean test",
                    "statistic in each regime. The Engle and Granger critical",
                    "value at 5 percent is $-3.34$."),
    label = "tab:rolling_regime"),
  file = "output/tables/rolling_cointegration_regime.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "H",
  sanitize.text.function = identity,
  sanitize.colnames.function = identity
)

# -------------------------------------------------------------
# Figures
# -------------------------------------------------------------
BREAKS <- data.frame(date = as.Date(c("2020-03-01", "2024-05-01")))

focus_pairs <- c("COP on XOM", "XOM on COP", "CVX on XOM", "XOM on CVX")

p_main <- roll %>%
  filter(pair_label %in% focus_pairs) %>%
  ggplot(aes(x = window_end, y = adf_stat)) +
  geom_vline(data = BREAKS, aes(xintercept = date),
             linetype = "dotted", colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = EG_CV["5pct"], colour = "#C0392B",
             linetype = "dashed", linewidth = 0.45) +
  geom_hline(yintercept = EG_CV["1pct"], colour = "#7B241C",
             linetype = "dotdash", linewidth = 0.4) +
  geom_line(colour = "#3B4CC0", linewidth = 0.6) +
  facet_wrap(~ pair_label, ncol = 2) +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
  labs(
    x = "Window end date",
    y = "ADF statistic on cointegrating residual",
    title = "Rolling Engle and Granger cointegration tests",
    subtitle = paste0(WIDTH, " day window stepped monthly. Values below the ",
                      "dashed line reject the unit root null at 5 percent."),
    caption = paste(
      "Dotted vertical lines mark March 2020 and May 2024. Critical values are",
      "those of Engle and Granger rather than\nstandard Dickey and Fuller,",
      "since the tested residual is estimated rather than observed."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, colour = "grey40", size = 8))

ggsave("output/figures/rolling_cointegration_pairs.png", p_main,
       width = 10, height = 6.5, dpi = 300, bg = "white")

# Window length robustness
message("Running window robustness ...")

widths <- c("1 year" = 252, "2 years" = 504, "3 years" = 756)

robust <- do.call(rbind, lapply(names(widths), function(w) {
  r <- rolling_eg("COP", "XOM", lp, dates, width = widths[[w]], by = STEP)
  r$window <- w
  r
}))
robust$window_end <- as.Date(robust$window_end, origin = "1970-01-01")
robust$window     <- factor(robust$window, levels = names(widths))

p_robust <- ggplot(robust, aes(x = window_end, y = adf_stat, colour = window)) +
  geom_vline(data = BREAKS, aes(xintercept = date),
             linetype = "dotted", colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = EG_CV["5pct"], colour = "#C0392B",
             linetype = "dashed", linewidth = 0.45) +
  geom_line(linewidth = 0.55) +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
  scale_colour_manual(values = c("1 year"  = "#8E9BE8",
                                 "2 years" = "#3B4CC0",
                                 "3 years" = "#1A2270")) +
  labs(x = "Window end date",
       y = "ADF statistic on cointegrating residual",
       colour = "Window length",
       title = "ConocoPhillips on Exxon Mobil: sensitivity to window length",
       subtitle = paste("Shorter windows detect a change sooner at the cost",
                        "of greater sampling noise")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("output/figures/rolling_cointegration_window_robustness.png", p_robust,
       width = 9.5, height = 5.5, dpi = 300, bg = "white")

# Hedge ratio stability
p_beta <- roll %>%
  filter(pair_label %in% focus_pairs) %>%
  ggplot(aes(x = window_end, y = beta)) +
  geom_vline(data = BREAKS, aes(xintercept = date),
             linetype = "dotted", colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_line(colour = "#1BA37A", linewidth = 0.6) +
  facet_wrap(~ pair_label, ncol = 2, scales = "free_y") +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
  labs(x = "Window end date",
       y = expression(hat(beta)),
       title = "Rolling cointegrating coefficient",
       subtitle = paste0("Estimated on a ", WIDTH, " day window. A coefficient ",
                         "fixed once in sample is held constant across all of ",
                         "this variation."),
       caption = paste("The grey line marks zero. A negative coefficient",
                       "implies the regression assigns the two securities",
                       "opposing exposures.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, colour = "grey40", size = 8))

ggsave("output/figures/rolling_hedge_ratio.png", p_beta,
       width = 10, height = 6, dpi = 300, bg = "white")

message("Rolling cointegration complete.")

# =============================================================
# 08_full_sample_eg.R
# Full sample Engle and Granger cointegration tests, for the
# candidate pairs carried from the Granger screening.
#
# Provides the single sample verdict that the rolling analysis of
# 07_rolling_cointegration.R then decomposes through time.
#
# Reads:  data/raw/prices_raw.rds
# Writes: output/tables/eg_full_sample.tex
# =============================================================

library(urca)
library(xtable)
library(zoo)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

prices   <- readRDS("data/raw/prices_raw.rds")
price_df <- as.data.frame(prices)
lp       <- log(price_df)

EG_CV <- c("1pct" = -3.90, "5pct" = -3.34, "10pct" = -3.04)

# Engle and Granger step one, full sample
eg_full <- function(y_name, x_name, data) {
  fit <- lm(data[[y_name]] ~ data[[x_name]])
  r   <- residuals(fit)
  adf <- ur.df(r, type = "none", lags = 4, selectlags = "AIC")
  data.frame(
    caused      = y_name,
    causing     = x_name,
    beta        = unname(coef(fit)[2]),
    adf_stat    = as.numeric(adf@teststat[1]),
    coint_5     = as.numeric(adf@teststat[1]) < EG_CV["5pct"],
    coint_1     = as.numeric(adf@teststat[1]) < EG_CV["1pct"],
    stringsAsFactors = FALSE
  )
}

# The pairs carried from the Granger shortlist, both directions where
# both survived. COP-XOM is retained regardless as the original's choice.
PAIRS <- list(
  c("COP", "XOM"), c("XOM", "COP"),
  c("CVX", "XOM"), c("XOM", "CVX"),
  c("COP", "CVX"), c("CVX", "COP"),
  c("XOM", "XLE"), c("CVX", "XLE"),
  c("COP", "XLE")
)

res <- do.call(rbind, lapply(PAIRS, function(p) eg_full(p[1], p[2], lp)))
res <- res[order(res$adf_stat), ]

cat("=== Full sample Engle and Granger, 1999 to 2026 ===\n")
cat("EG 5% critical value:", EG_CV["5pct"], "\n\n")
print(res, digits = 4, row.names = FALSE)

# Compare against the original's training window, 1999-12-31 to 2013-04-30,
# to show the pair tested as cointegrated then and does not now.
train_end <- which(as.Date(as.POSIXct(zoo::index(prices),
                                      origin = "1970-01-01", tz = "UTC")) <= as.Date("2013-04-30"))
train_end <- max(train_end)

lp_train <- lp[1:train_end, ]

cat("\n=== COP-XOM on the original training window only (to 2013-04-30) ===\n")
print(eg_full("COP", "XOM", lp_train), digits = 4, row.names = FALSE)
print(eg_full("XOM", "COP", lp_train), digits = 4, row.names = FALSE)

# Export the full sample table
out <- data.frame(
  Pair       = paste0(res$caused, " on ", res$causing),
  `$\\hat{\\beta}$` = res$beta,
  `ADF stat` = res$adf_stat,
  `Cointegrated at 5\\%` = ifelse(res$coint_5, "Yes", "No"),
  check.names = FALSE
)

print(
  xtable(out,
         digits = c(0, 0, 4, 3, 0),
         caption = paste("Full sample Engle and Granger cointegration tests,",
                         "31 December 1999 to 24 July 2026, on log price levels.",
                         "The residual ADF statistic is evaluated against the",
                         "Engle and Granger 5 percent critical value of",
                         "$-3.34$. No pair rejects the unit root null over the",
                         "full sample."),
         label = "tab:eg_full"),
  file = "output/tables/eg_full_sample.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "H",
  sanitize.text.function = identity,
  sanitize.colnames.function = identity
)

message("Full sample EG complete. Report the console COP-XOM training-window ",
        "statistics in the chapter where marked.")


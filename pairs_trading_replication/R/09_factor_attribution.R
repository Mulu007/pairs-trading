# =============================================================
# 11_factor_attribution.R
# Factor attribution of the strategy returns against the
# Fama-French five-factor model plus momentum.
#
# Tests whether the strategy's risk-adjusted return is genuine alpha
# or compensation for exposure to common factors, in particular the
# low-volatility and value tilts a market-neutral energy spread would
# plausibly carry.
#
# Reads:  data/processed/backtest_results.rds
# Writes: output/tables/factor_attribution.tex
#         output/figures/factor_loadings.png
#         data/processed/factor_data.rds
#
# NOTE ON DATA SOURCE
# The Fama-French factors are obtained with the frenchdata package,
# which downloads directly from Kenneth French's data library. If that
# package or network access is unavailable, download the daily files
#   F-F_Research_Data_5_Factors_2x3_daily.CSV
#   F-F_Momentum_Factor_daily.CSV
# manually from the library and point read_ff_local() at them.
# =============================================================
install.packages("frenchdata")
library(dplyr)
library(tidyr)
library(ggplot2)
library(xtable)
library(sandwich)   # Newey-West standard errors
library(lmtest)     # coeftest

dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

bt <- readRDS("data/processed/backtest_results.rds")

strat <- data.frame(date = bt$dates, ret = bt$bt$strat_ret) %>%
  filter(is.finite(ret))

# -------------------------------------------------------------
# Load Fama-French daily factors.
# -------------------------------------------------------------
get_ff <- function() {
  if (requireNamespace("frenchdata", quietly = TRUE)) {
    ff5 <- frenchdata::download_french_data(
      "Fama/French 5 Factors (2x3) [Daily]")$subsets$data[[1]]
    mom <- frenchdata::download_french_data(
      "Momentum Factor (Mom) [Daily]")$subsets$data[[1]]
    
    ff5 <- ff5 %>%
      mutate(date = as.Date(as.character(date), "%Y%m%d")) %>%
      rename(MktRF = `Mkt-RF`)
    mom <- mom %>%
      mutate(date = as.Date(as.character(date), "%Y%m%d")) %>%
      rename(MOM = Mom)
    
    out <- inner_join(ff5, mom, by = "date")
    # French factors are in percent; convert to decimal
    num <- c("MktRF","SMB","HML","RMW","CMA","RF","MOM")
    out[num] <- out[num] / 100
    return(out)
  }
  stop("frenchdata not available; use read_ff_local() with manual files.")
}

read_ff_local <- function(ff5_path, mom_path) {
  # Fallback loader for manually downloaded CSVs. The French CSVs carry
  # header and footer text that must be stripped; adjust skip/nrows to
  # the actual file. Left as a template rather than hard-coded.
  stop("Point this function at your local FF CSVs and set skip/nrows.")
}

ff <- get_ff()
saveRDS(ff, "data/processed/factor_data.rds")

# -------------------------------------------------------------
# Merge and form excess strategy return.
# The strategy is self-financing and market-neutral, but the excess
# return convention is retained for comparability with the factor model.
# -------------------------------------------------------------
dat <- strat %>%
  inner_join(ff, by = "date") %>%
  mutate(exc = ret - RF)

cat("Merged observations:", nrow(dat), "\n")
cat("Date range:", format(min(dat$date)), "to", format(max(dat$date)), "\n\n")

# -------------------------------------------------------------
# Regressions.
#
# Three nested specifications:
#   (1) CAPM: excess return on market only
#   (2) FF5:  five factors
#   (3) FF5 + momentum
#
# Standard errors are Newey-West to account for the residual
# autocorrelation that the long holding periods induce.
# -------------------------------------------------------------
run_reg <- function(formula, data, lags = 10) {
  m  <- lm(formula, data = data)
  nw <- coeftest(m, vcov = NeweyWest(m, lag = lags, prewhite = FALSE))
  list(model = m, coef = nw)
}

m1 <- run_reg(exc ~ MktRF, dat)
m2 <- run_reg(exc ~ MktRF + SMB + HML + RMW + CMA, dat)
m3 <- run_reg(exc ~ MktRF + SMB + HML + RMW + CMA + MOM, dat)

# annualise alpha (daily -> annual)
ann_alpha <- function(nw) {
  a  <- nw["(Intercept)", "Estimate"]
  se <- nw["(Intercept)", "Std. Error"]
  p  <- nw["(Intercept)", "Pr(>|t|)"]
  c(alpha_ann = a * 252, alpha_daily = a, se = se, p = p)
}

cat("=== CAPM ===\n");    print(m1$coef)
cat("\n=== FF5 ===\n");   print(m2$coef)
cat("\n=== FF5 + MOM ===\n"); print(m3$coef)

cat("\n=== Annualised alpha ===\n")
cat("CAPM:      ", sprintf("%.2f%% (p=%.3f)\n",
                           100*ann_alpha(m1$coef)["alpha_ann"], ann_alpha(m1$coef)["p"]))
cat("FF5:       ", sprintf("%.2f%% (p=%.3f)\n",
                           100*ann_alpha(m2$coef)["alpha_ann"], ann_alpha(m2$coef)["p"]))
cat("FF5 + MOM: ", sprintf("%.2f%% (p=%.3f)\n",
                           100*ann_alpha(m3$coef)["alpha_ann"], ann_alpha(m3$coef)["p"]))

# -------------------------------------------------------------
# Build the results table (FF5 + MOM as primary specification)
# -------------------------------------------------------------
build_row <- function(nw, name) {
  if (!name %in% rownames(nw)) return(c("", ""))
  est <- nw[name, "Estimate"]
  p   <- nw[name, "Pr(>|t|)"]
  star <- ifelse(p<0.01,"***",ifelse(p<0.05,"**",ifelse(p<0.10,"*","")))
  c(sprintf("%.4f%s", est, star), sprintf("(%.3f)", p))
}

factors <- c("(Intercept)","MktRF","SMB","HML","RMW","CMA","MOM")
labels  <- c("Alpha (daily)","Market","SMB (size)","HML (value)",
             "RMW (profitability)","CMA (investment)","MOM (momentum)")

tbl <- data.frame(
  Factor = labels,
  CAPM   = sapply(factors, function(f) build_row(m1$coef, f)[1]),
  FF5    = sapply(factors, function(f) build_row(m2$coef, f)[1]),
  FF5_MOM= sapply(factors, function(f) build_row(m3$coef, f)[1]),
  stringsAsFactors = FALSE
)
names(tbl) <- c("Factor","CAPM","FF5","FF5 + MOM")

# add R-squared row
r2row <- data.frame(Factor = "$R^2$",
                    CAPM = sprintf("%.3f", summary(m1$model)$r.squared),
                    FF5  = sprintf("%.3f", summary(m2$model)$r.squared),
                    `FF5 + MOM` = sprintf("%.3f", summary(m3$model)$r.squared),
                    check.names = FALSE)
tbl <- rbind(tbl, r2row)

print(
  xtable(tbl,
         caption = paste("Factor attribution of daily strategy excess returns.",
                         "Each column is a regression on the indicated factor set. Newey-West",
                         "standard errors with 10 lags; $p$ values in the discussion. Stars:",
                         "*** $p<0.01$, ** $p<0.05$, * $p<0.10$. A daily alpha of 0.0001 is",
                         "approximately 2.5 percent annualised."),
         label = "tab:factor_attr"),
  file = "output/tables/factor_attribution.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "H",
  sanitize.text.function = identity,
  sanitize.colnames.function = identity
)

# -------------------------------------------------------------
# Coefficient plot with confidence intervals (FF5 + MOM)
# -------------------------------------------------------------
ci_df <- as.data.frame(m3$coef[,c("Estimate","Std. Error")])
ci_df$factor <- rownames(ci_df)
ci_df <- ci_df %>%
  filter(factor != "(Intercept)") %>%
  mutate(lo = Estimate - 1.96*`Std. Error`,
         hi = Estimate + 1.96*`Std. Error`,
         factor = factor(factor,
                         levels = c("MktRF","SMB","HML","RMW","CMA","MOM")))

p_load <- ggplot(ci_df, aes(x = factor, y = Estimate)) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2,
                colour = "#3B4CC0") +
  geom_point(size = 2.5, colour = "#C0392B") +
  labs(x = NULL, y = "Factor loading",
       title = "Strategy factor loadings",
       subtitle = "Fama-French five factors plus momentum, 95 percent intervals") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("output/figures/factor_loadings.png", p_load,
       width = 8, height = 5, dpi = 300, bg = "white")

message("Factor attribution complete.")

# =============================================================
# 11b_factor_periods.R
# Per-period factor attribution: does the alpha decay across the
# extension period the way the raw returns did?
#
# Run AFTER 11_factor_attribution.R (reuses factor_data.rds).
#
# Reads:  data/processed/backtest_results.rds
#         data/processed/factor_data.rds
# Writes: output/tables/factor_attribution_periods.tex
#         output/figures/factor_alpha_by_period.png
# =============================================================

library(dplyr)
library(ggplot2)
library(xtable)
library(sandwich)
library(lmtest)

bt <- readRDS("data/processed/backtest_results.rds")
ff <- readRDS("data/processed/factor_data.rds")

strat <- data.frame(date = bt$dates, ret = bt$bt$strat_ret) %>%
  filter(is.finite(ret))

dat <- strat %>%
  inner_join(ff, by = "date") %>%
  mutate(exc = ret - RF)

IN_END  <- as.Date("2013-04-30")
OOS_END <- as.Date("2018-12-31")

dat <- dat %>%
  mutate(period = case_when(
    date <= IN_END  ~ "In-sample",
    date <= OOS_END ~ "Out-of-sample",
    TRUE            ~ "Extension"
  ))
dat$period <- factor(dat$period,
                     levels = c("In-sample","Out-of-sample","Extension"))

# -------------------------------------------------------------
# FF5 + MOM within each period, plus the full sample for reference.
# Newey-West standard errors throughout.
# -------------------------------------------------------------
fit_period <- function(d, lab) {
  m  <- lm(exc ~ MktRF + SMB + HML + RMW + CMA + MOM, data = d)
  nw <- coeftest(m, vcov = NeweyWest(m, lag = 10, prewhite = FALSE))
  a   <- nw["(Intercept)","Estimate"]
  se  <- nw["(Intercept)","Std. Error"]
  p   <- nw["(Intercept)","Pr(>|t|)"]
  mom <- nw["MOM","Estimate"]
  momp<- nw["MOM","Pr(>|t|)"]
  data.frame(
    Period      = lab,
    n           = nrow(d),
    alpha_ann   = 252 * a,
    alpha_lo    = 252 * (a - 1.96*se),
    alpha_hi    = 252 * (a + 1.96*se),
    alpha_p     = p,
    mom_load    = mom,
    mom_p       = momp,
    r2          = summary(m)$r.squared
  )
}

res <- bind_rows(
  fit_period(dat %>% filter(period=="In-sample"),     "In-sample"),
  fit_period(dat %>% filter(period=="Out-of-sample"), "Out-of-sample"),
  fit_period(dat %>% filter(period=="Extension"),     "Extension"),
  fit_period(dat,                                     "Full sample")
)

cat("=== Alpha and momentum loading by period (FF5 + MOM) ===\n")
print(as.data.frame(res), digits = 4, row.names = FALSE)

# -------------------------------------------------------------
# Table
# -------------------------------------------------------------
star <- function(p) ifelse(p<0.01,"***",ifelse(p<0.05,"**",ifelse(p<0.10,"*","")))

out <- res %>% transmute(
  Period                       = Period,
  `Obs.`                       = n,
  `Alpha (ann.)`               = sprintf("%.2f\\%%%s", 100*alpha_ann, star(alpha_p)),
  `95\\% CI`                   = sprintf("[%.1f, %.1f]", 100*alpha_lo, 100*alpha_hi),
  `MOM loading`                = sprintf("%.3f%s", mom_load, star(mom_p)),
  `$R^2$`                      = sprintf("%.3f", r2)
)

print(
  xtable(as.data.frame(out),
         caption = paste("Factor-adjusted alpha and momentum loading by period,",
                         "from the Fama-French five-factor plus momentum model estimated",
                         "separately within each period. Newey-West standard errors, 10 lags.",
                         "Confidence intervals are for the annualised alpha. Stars: *** $p<0.01$,",
                         "** $p<0.05$, * $p<0.10$."),
         label = "tab:factor_periods"),
  file = "output/tables/factor_attribution_periods.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "H",
  sanitize.text.function = identity,
  sanitize.colnames.function = identity
)

# -------------------------------------------------------------
# Alpha-by-period plot with confidence intervals
# -------------------------------------------------------------
plot_df <- res %>%
  filter(Period != "Full sample") %>%
  mutate(Period = factor(Period,
                         levels = c("In-sample","Out-of-sample","Extension")))

p_alpha <- ggplot(plot_df, aes(x = Period, y = 100*alpha_ann)) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
  geom_errorbar(aes(ymin = 100*alpha_lo, ymax = 100*alpha_hi),
                width = 0.15, colour = "#3B4CC0") +
  geom_point(size = 3, colour = "#C0392B") +
  labs(x = NULL, y = "Annualised alpha (%)",
       title = "Factor-adjusted alpha by period",
       subtitle = "FF5 + momentum, 95 percent confidence intervals") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("output/figures/factor_alpha_by_period.png", p_alpha,
       width = 8, height = 5, dpi = 300, bg = "white")

message("Per-period factor attribution complete.")

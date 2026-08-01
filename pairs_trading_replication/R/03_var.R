#install.packages("vars")
library(vars)
library (urca)
library(tidyverse)

prices  <- readRDS("data/raw/prices_raw.rds")
returns <- readRDS("data/processed/returns.rds")
ret_df  <- as.data.frame(returns)
price_df <- as.data.frame(prices)

adf_table <- function(x, series_name) {
  # Schwert rule of thumb for max lag
  n <- length(x)
  pmax <- floor(12 * (n/100)^(1/4))
  
  test <- ur.df(x, type = "drift", lags = pmax, selectlags = "AIC")
  
  # count z.diff.lag terms actually retained
  coef_names <- rownames(test@testreg$coefficients)
  n_lags <- sum(grepl("z.diff.lag", coef_names))
  
  data.frame(
    series   = series_name,
    adf_stat = test@teststat[1],
    cv_1pct  = test@cval[1, 1],
    cv_5pct  = test@cval[1, 2],
    cv_10pct = test@cval[1, 3],
    lags     = n_lags
  )
}

# returns
adf_returns <- do.call(rbind, lapply(colnames(ret_df), function(nm)
  adf_table(ret_df[[nm]], nm)))
print(adf_returns)

# adf on price levels
price_df <- as.data.frame(prices)
adf_levels <- do.call(rbind, lapply(colnames(price_df), function(nm)
  adf_table(price_df[[nm]], nm)))
print(adf_levels)

# kpss on price levels
kpss_results <- do.call(rbind, lapply(colnames(price_df), function(nm) {
  k <- ur.kpss(price_df[[nm]], type = "mu", lags = "short")
  data.frame(
    series   = nm,
    kpss_stat = k@teststat,
    cv_10pct = k@cval[1], cv_5pct = k@cval[2],
    cv_2.5pct = k@cval[3], cv_1pct = k@cval[4]
  )
}))
kpss_results

# Lag selection
ret_mat <- as.matrix(ret_df)
lag_sel <- VARselect(ret_mat, lag.max = 20, type = "const")
lag_sel$selection
lag_sel$criteria

# Visualisation of all lag selection criterias
library(ggplot2)
library(tidyr)

crit_df <- as.data.frame(t(lag_sel$criteria))
crit_df$lag <- 1:nrow(crit_df)

crit_long <- pivot_longer(crit_df, cols = c("AIC(n)", "HQ(n)", "SC(n)", "FPE(n)"),
                          names_to = "criterion", values_to = "value")

# FPE is on a different scale (raw, not log), so plot it separately
crit_long_log <- crit_long[crit_long$criterion != "FPE(n)", ]

selected <- data.frame(
  criterion = c("AIC(n)", "HQ(n)", "SC(n)"),
  lag = c(which.min(crit_df[["AIC(n)"]]),
          which.min(crit_df[["HQ(n)"]]),
          which.min(crit_df[["SC(n)"]]))
)
selected$value <- mapply(function(c, l) crit_df[l, c], selected$criterion, selected$lag)

p <- ggplot(crit_long_log, aes(x = lag, y = value, colour = criterion)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  geom_point(data = selected, aes(x = lag, y = value, colour = criterion),
             size = 3.5, shape = 21, fill = "white", stroke = 1.2) +
  scale_x_continuous(breaks = 1:20) +
  labs(x = "Lag order p", y = "Criterion value",
       colour = NULL,
       title = "Information criteria as a function of lag order",
       subtitle = "Open circles mark each criterion's selected minimum") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/lag_criteria.png", p, width = 8, height = 5, dpi = 300)

# Visualisation of fpe
p_fpe <- ggplot(crit_df, aes(x = lag, y = `FPE(n)`)) +
  geom_line(colour = "#2a78d6", linewidth = 0.6) +
  geom_point(size = 1.5, colour = "#2a78d6") +
  scale_x_continuous(breaks = 1:20) +
  labs(x = "Lag order p", y = "FPE",
       title = "Final prediction error as a function of lag order") +
  theme_minimal(base_size = 11)

ggsave("output/figures/fpe_criteria.png", p_fpe, width = 7, height = 4.5, dpi = 300)

# Table Exportaion
library(xtable)

crit_export <- as.data.frame(round(t(lag_sel$criteria), 6))
crit_export <- cbind(Lag = 1:nrow(crit_export), crit_export)
colnames(crit_export) <- c("Lag", "AIC", "HQ", "SC", "FPE")

print(
  xtable(crit_export, digits = c(0, 0, 6, 6, 6, 10),
         caption = "Information criteria by lag order, $p = 1$ to $20$",
         label = "tab:lag_criteria_full"),
  file = "output/tables/lag_criteria_full.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "ht"
)

# Estimate
var_model <- VAR(ret_mat, p = 2, type = "const")
summary(var_model)

# stability
roots(var_model)
all(roots(var_model) < 1)

serial.test(var_model, lags.pt = 16, type = "PT.asymptotic")  # residual autocorrelation
arch.test(var_model, lags.multi = 5)                           # conditional heteroskedasticity
normality.test(var_model)                                      # residual normality

library(xtable)
library(vars)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# ---- helper: write tabular-only fragment ----
export_tbl <- function(obj, file, caption, label, digits = 4, rownames = TRUE) {
  print(
    xtable(obj, digits = digits, caption = caption, label = label),
    file = file, booktabs = TRUE, include.rownames = rownames,
    caption.placement = "top", table.placement = "ht"
  )
}

# ---- 1. ADF levels ----
adf_lev <- adf_levels
adf_lev$series <- gsub("_r$", "", adf_lev$series)
export_tbl(adf_lev,
           "output/tables/adf_levels.tex",
           "Augmented Dickey and Fuller tests on price levels. Lag length selected by AIC subject to the Schwert bound of 34.",
           "tab:adf_levels", digits = c(0,0,4,2,2,2,0), rownames = FALSE)

# ---- 2. ADF returns ----
adf_ret <- adf_returns
adf_ret$series <- gsub("_r$", "", adf_ret$series)
export_tbl(adf_ret,
           "output/tables/adf_returns.tex",
           "Augmented Dickey and Fuller tests on daily log returns.",
           "tab:adf_returns", digits = c(0,0,4,2,2,2,0), rownames = FALSE)

# ---- 3. KPSS levels ----
kpss_out <- kpss_results
kpss_out$series <- gsub("_r$", "", kpss_out$series)
export_tbl(kpss_out,
           "output/tables/kpss_levels.tex",
           "Kwiatkowski, Phillips, Schmidt and Shin tests on price levels. The null hypothesis is stationarity.",
           "tab:kpss_levels", digits = c(0,0,3,3,3,3,3), rownames = FALSE)

# ---- 4. VAR coefficient table ----
# Build a compact coefficient matrix: rows = regressors, cols = equations
# with significance stars, plus a fit-statistics block at the bottom.

star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01,  "**",
                ifelse(p < 0.05,  "*",
                       ifelse(p < 0.10,  ".", ""))))
}

eqs <- names(var_model$varresult)
regressors <- rownames(summary(var_model)$varresult[[1]]$coefficients)

coef_mat <- sapply(eqs, function(eq) {
  cf <- summary(var_model)$varresult[[eq]]$coefficients
  paste0(formatC(cf[, "Estimate"], format = "f", digits = 4),
         star(cf[, "Pr(>|t|)"]))
})
rownames(coef_mat) <- regressors
colnames(coef_mat) <- gsub("_r$", "", eqs)

# fit stats block
fit_block <- sapply(eqs, function(eq) {
  s <- summary(var_model)$varresult[[eq]]
  c(formatC(s$r.squared, format = "f", digits = 4),
    formatC(s$adj.r.squared, format = "f", digits = 4),
    formatC(s$fstatistic[1], format = "f", digits = 2),
    formatC(s$sigma, format = "f", digits = 5))
})
rownames(fit_block) <- c("$R^2$", "Adj. $R^2$", "$F$ statistic", "Resid. s.e.")
colnames(fit_block) <- colnames(coef_mat)

var_tbl <- rbind(coef_mat, fit_block)

# clean regressor labels for LaTeX (underscores)
rownames(var_tbl) <- gsub("_r\\.l", ".L", rownames(var_tbl))

print(
  xtable(var_tbl,
         caption = "VAR(2) coefficient estimates. Significance: *** $p<0.001$, ** $p<0.01$, * $p<0.05$, . $p<0.10$.",
         label = "tab:var_coefs"),
  file = "output/tables/var_coefficients.tex",
  booktabs = TRUE, include.rownames = TRUE,
  caption.placement = "top", table.placement = "ht",
  sanitize.text.function = identity,
  sanitize.rownames.function = identity
)

# ---- 5. Residual correlation matrix ----
resid_cor <- cor(residuals(var_model))
rownames(resid_cor) <- colnames(resid_cor) <- gsub("_r$", "", colnames(resid_cor))
export_tbl(resid_cor,
           "output/tables/var_resid_correlation.tex",
           "Contemporaneous correlation matrix of VAR(2) residuals.",
           "tab:var_resid_cor", digits = 4)

message("Exported 5 tables to output/tables/")

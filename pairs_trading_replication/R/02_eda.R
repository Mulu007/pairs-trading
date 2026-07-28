library(quantmod)
library(ggplot2)
library(tidyr)
library(dplyr)

prices <- readRDS("data/raw/prices_raw.rds")

# log return transformations
returns <- diff(log(prices))
returns <- na.omit(returns)
colnames(returns) <- paste0(colnames(prices), "_r")

saveRDS(returns, "data/processed/returns.rds")

dim(returns) # -> differencing costs one obs therefore 6679

# --------- Summary Statistics ---------
library(moments)

ret_df <- as.data.frame(returns)

summary_stats <- data.frame(
  count = sapply(ret_df, length),
  mean  = sapply(ret_df, mean),
  sd    = sapply(ret_df, sd),
  min   = sapply(ret_df, min),
  q25   = sapply(ret_df, quantile, 0.25),
  median= sapply(ret_df, median),
  q75   = sapply(ret_df, quantile, 0.75),
  max   = sapply(ret_df, max),
  skew  = sapply(ret_df, skewness),
  kurt  = sapply(ret_df, function(x) kurtosis(x) - 3)
)

round(summary_stats, 5)

library(tseries)

# Formal Normality tests Jarque Bera similar to D'Agostino-Pearson
normality_tests <- data.frame(
  jb_stat = sapply(ret_df, function(x) jarque.bera.test(x)$statistic),
  jb_pval = sapply(ret_df, function(x) jarque.bera.test(x)$p.value),
  ad_skew_pval = sapply(ret_df, function(x) agostino.test(x)$p.value),
  ansc_kurt_pval = sapply(ret_df, function(x) anscombe.test(x)$p.value)
)

normality_tests <- round(normality_tests, 6)
normality_tests$normal <- ifelse(normality_tests$jb_pval > 0.05, "Yes", "No")
normality_tests
rownames(normality_tests) <- gsub("_r\\.X-squared$", "", rownames(normality_tests))
rownames(summary_stats)   <- gsub("_r$", "", rownames(summary_stats))

# Standardization 
ret_std <- scale(ret_df)
ret_std_df <- as.data.frame(ret_std)
ret_std_df$std_normal <- rnorm(nrow(ret_std_df))

# Plots
library(ggplot2)

qq_data <- ret_std_df %>%
  pivot_longer(everything(), names_to = "series", values_to = "value")

p <- ggplot(qq_data, aes(sample = value)) +
  stat_qq(size = 0.5, alpha = 0.4) +
  stat_qq_line(colour = "red") +
  facet_wrap(~ series, ncol = 3) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "Normal QQ plots, standardized daily log returns") +
  theme_minimal()

ggsave("output/figures/qq_plots.png", p, width = 9, height = 6, dpi = 300)

# Density and correlation
p2 <- ggplot(qq_data, aes(x = value, colour = series)) +
  geom_density() +
  xlim(-6, 6) +
  labs(x = "Standardized return", y = "Density") +
  theme_minimal()

ggsave("output/figures/density.png", p2, width = 8, height = 5, dpi = 300)

cor_mat <- cor(ret_df)
round(cor_mat, 4)


library(xtable)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# Summary statistics table
summary_t <- t(summary_stats)

print(
  xtable(summary_t, digits = 6,
         caption = "Descriptive statistics of daily log returns, 2000 to 2026",
         label = "tab:summary_stats"),
  file = "output/tables/summary_stats.tex",
  booktabs = TRUE, include.rownames = TRUE, caption.placement = "top"
)
summary_t <- t(summary_stats[, !colnames(summary_stats) %in% "count"])

# Normality test table
print(
  xtable(normality_tests,
         digits = 6,
         caption = "Normality tests on daily log returns",
         label = "tab:normality"),
  file = "output/tables/normality_tests.tex",
  booktabs = TRUE,
  include.rownames = TRUE,
  caption.placement = "top"
)

# Correlation matrix
cor_mat <- cor(ret_df)
rownames(cor_mat) <- colnames(cor_mat) <- gsub("_r$", "", colnames(cor_mat))
print(xtable(cor_mat, digits = 4, caption = "Correlation matrix of daily log returns",
             label = "tab:correlation"),
      file = "output/tables/correlation.tex", booktabs = TRUE,
      caption.placement = "top")

write_bordermatrix <- function(m, file, digits = 3) {
  lbl <- colnames(m)
  hdr <- paste("", paste(lbl, collapse = " & "), sep = " & ")
  rows <- sapply(seq_len(nrow(m)), function(i) {
    paste(lbl[i], paste(round(m[i, ], digits), collapse = " & "), sep = " & ")
  })
  writeLines(
    c("\\[", "\\bordermatrix{",
      paste0(hdr, " \\cr"),
      paste0(rows, " \\cr"),
      "}", "\\]"),
    file
  )
}

write_bordermatrix(cor_mat, "output/tables/correlation_matrix.tex")

# Tail figure
# Tail comparison: Gaussian vs Student-t on log scale
# Produces output/figures/tail_comparison.png

library(ggplot2)

x <- seq(0, 16, length.out = 2000)

tail_df <- data.frame(
  x = rep(x, 3),
  density = c(dnorm(x),
              dt(x, df = 3),
              dt(x, df = 5)),
  dist = rep(c("Normal(0,1)", "Student t, 3 df", "Student t, 5 df"),
             each = length(x))
)

# observed extremes, in absolute standardised units
extremes <- data.frame(
  k = c(9.1, 9.6, 12.4, 13.6, 14.5),
  label = c("XOM", "SPY", "XLE", "COP", "CVX")
)

p <- ggplot(tail_df, aes(x = x, y = density, colour = dist)) +
  geom_line(linewidth = 0.7) +
  geom_vline(data = extremes, aes(xintercept = k),
             linetype = "dotted", colour = "grey40", linewidth = 0.3) +
  geom_text(data = extremes,
            aes(x = k, y = 1e-2, label = label),
            inherit.aes = FALSE, angle = 90, vjust = -0.4,
            size = 2.6, colour = "grey30") +
  scale_y_log10(
    limits = c(1e-50, 1),
    breaks = 10^seq(-50, 0, by = 10),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(
    x = "Standard deviations from mean",
    y = "Density (log scale)",
    colour = NULL,
    title = "Tail decay under Gaussian and power law specifications",
    subtitle = "Dotted lines mark the worst observed day in each series"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/tail_comparison.png", p,
       width = 8, height = 5.5, dpi = 300)

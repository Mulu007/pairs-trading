# =============================================================
# 06_granger.R
# Pairwise Granger causality across all ordered security pairs,
# with a multiple comparison correction the original study omits.
#
# Reads:  data/processed/returns.rds
# Writes: output/tables/granger_pairwise.tex
#         output/tables/granger_shortlist.tex
#         output/figures/granger_heatmap.png
# =============================================================

library(vars)
library(ggplot2)
library(dplyr)
library(xtable)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

returns <- readRDS("data/processed/returns.rds")
ret_df  <- as.data.frame(returns)

P <- 2  # matching the baseline VAR specification

# -------------------------------------------------------------
# Why bivariate VARs rather than vars::causality on the full system.
#
# causality() tests one variable against ALL others jointly, which
# answers "does X predict the rest of the system." The pairs analysis
# needs the narrower question "does X predict Y specifically," so each
# pair is fitted as its own two-variable VAR. This also matches the
# structure of the table in the original study.
#
# Note the cost: a bivariate VAR omits the other three securities, so
# a rejection may reflect a common driver rather than a direct link.
# That limitation is inherent to the pairwise framing and is recorded
# in the write-up rather than corrected here.
# -------------------------------------------------------------
pairwise_gc <- function(y_name, x_name, data, p = P) {
  sub <- as.matrix(data[, c(y_name, x_name)])
  m   <- VAR(sub, p = p, type = "const")
  tst <- causality(m, cause = x_name)
  data.frame(
    caused  = gsub("_r$", "", y_name),
    causing = gsub("_r$", "", x_name),
    f_stat  = unname(tst$Granger$statistic),
    df1     = unname(tst$Granger$parameter[1]),
    df2     = unname(tst$Granger$parameter[2]),
    p_value = unname(tst$Granger$p.value),
    stringsAsFactors = FALSE
  )
}

series <- colnames(ret_df)
grid <- expand.grid(y = series, x = series, stringsAsFactors = FALSE)
grid <- grid[grid$y != grid$x, ]   # 20 ordered pairs from 5 securities

gc_results <- do.call(
  rbind,
  Map(function(y, x) pairwise_gc(y, x, ret_df), grid$y, grid$x)
)
rownames(gc_results) <- NULL

# -------------------------------------------------------------
# Multiple comparison correction.
#
# Twenty simultaneous tests at alpha = 0.10 yields two expected
# false rejections by chance alone. The original study reports ten
# surviving pairs at the 10 percent level with no adjustment, which
# materially overstates the evidence.
#
# Benjamini-Hochberg controls the false discovery rate, the expected
# proportion of rejections that are false. Bonferroni is reported
# alongside as the conservative bound: it controls the probability of
# ANY false rejection, which is a stricter and usually too-blunt
# criterion for a screening exercise like this one.
# -------------------------------------------------------------
gc_results$p_bh   <- p.adjust(gc_results$p_value, method = "BH")
gc_results$p_bonf <- p.adjust(gc_results$p_value, method = "bonferroni")

gc_results <- gc_results %>%
  mutate(
    sig_raw  = p_value < 0.10,
    sig_bh   = p_bh    < 0.10,
    sig_bonf = p_bonf  < 0.10,
    pair     = paste0(causing, " \u2192 ", caused)
  ) %>%
  arrange(p_value)

cat("\n=== All 20 ordered pairs ===\n")
print(gc_results %>%
        select(pair, f_stat, p_value, p_bh, p_bonf, sig_raw, sig_bh, sig_bonf),
      digits = 4)

cat("\n=== Rejection counts at alpha = 0.10 ===\n")
cat("Unadjusted:      ", sum(gc_results$sig_raw),  "of 20\n")
cat("Benjamini-Hochberg:", sum(gc_results$sig_bh),   "of 20\n")
cat("Bonferroni:      ", sum(gc_results$sig_bonf), "of 20\n")

# how many survive raw screening but fail after correction
lost <- sum(gc_results$sig_raw & !gc_results$sig_bh)
cat("Pairs lost to the FDR correction:", lost, "\n\n")

# -------------------------------------------------------------
# Full results table
# -------------------------------------------------------------
star <- function(p) {
  ifelse(p < 0.01, "***",
         ifelse(p < 0.05, "**",
                ifelse(p < 0.10, "*", "")))
}

full_tbl <- gc_results %>%
  transmute(
    Pair       = paste0(causing, " $\\rightarrow$ ", caused),
    `$F$`      = f_stat,
    `$p$`      = p_value,
    `$p$ (BH)` = p_bh,
    Sig        = star(p_bh)
  )

print(
  xtable(as.data.frame(full_tbl),
         digits = c(0, 0, 3, 4, 4, 0),
         caption = "Pairwise Granger causality tests across all twenty ordered pairs, bivariate VAR(2). The null is that lagged values of the causing series carry no predictive information for the caused series. Stars reflect Benjamini-Hochberg adjusted $p$ values: *** $p<0.01$, ** $p<0.05$, * $p<0.10$.",
         label = "tab:granger_full"),
  file = "output/tables/granger_pairwise.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "htbp",
  sanitize.text.function = identity,
  sanitize.colnames.function = identity
)

# -------------------------------------------------------------
# Shortlist: pairs surviving the FDR correction, which become the
# candidate set carried into the cointegration analysis.
# -------------------------------------------------------------
shortlist <- gc_results %>%
  filter(sig_bh) %>%
  transmute(
    Pair       = paste0(causing, " $\\rightarrow$ ", caused),
    `$F$`      = f_stat,
    `$p$`      = p_value,
    `$p$ (BH)` = p_bh
  )

print(
  xtable(as.data.frame(shortlist),
         digits = c(0, 0, 3, 4, 4),
         caption = "Ordered pairs surviving Granger causality screening after Benjamini-Hochberg correction at the 10 percent level.",
         label = "tab:granger_shortlist"),
  file = "output/tables/granger_shortlist.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "htbp",
  sanitize.text.function = identity,
  sanitize.colnames.function = identity
)

# -------------------------------------------------------------
# Heatmap of adjusted p values.
# Rows are the caused series, columns the causing series.
# -------------------------------------------------------------
heat_df <- gc_results %>%
  select(causing, caused, p_bh)

lvl <- c("COP", "CVX", "XOM", "XLE", "SPY")
heat_df$causing <- factor(heat_df$causing, levels = lvl)
heat_df$caused  <- factor(heat_df$caused,  levels = rev(lvl))

p_heat <- ggplot(heat_df, aes(x = causing, y = caused, fill = p_bh)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", p_bh)), size = 3, colour = "grey15") +
  scale_fill_gradientn(
    colours = c("#c0392b", "#e8b84b", "#f2f2f2"),
    values  = scales::rescale(c(0, 0.10, 1)),
    limits  = c(0, 1),
    name    = "BH-adjusted\n p value"
  ) +
  labs(
    x = "Causing series", y = "Caused series",
    title = "Pairwise Granger causality",
    subtitle = "Benjamini-Hochberg adjusted p values; red indicates rejection of no causality"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 9))

ggsave("output/figures/granger_heatmap.png", p_heat,
       width = 7.5, height = 6, dpi = 300)

saveRDS(gc_results, "data/processed/granger_results.rds")
message("Granger causality analysis complete.")


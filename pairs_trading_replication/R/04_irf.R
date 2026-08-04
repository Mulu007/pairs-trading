# =============================================================
# 04b_irf_plots.R
# Impulse response figures.
#
# Reads:  data/processed/returns.rds
#         data/processed/irf_tidy.rds
# Writes: output/figures/irf_marginal_responses.png
#         output/figures/irf_cumulative_responses.png
#         output/figures/irf_decay_dynamics_excluding_impact.png
#         output/figures/irf_ordering_comparison_energy_pairs.png
#         output/figures/irf_impact_response_ordering_contrast.png
#         output/tables/irf_longrun_cumulative.tex
# =============================================================

library(vars)
library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# Horizon shown in figures. The original study plots five days; six is
# used here so the reader sees one period of flat line confirming the
# responses have settled. The underlying data extend to ten and the decay
# table in the text still draws on the full range.
MAX_H <- 6

irf_all <- readRDS("data/processed/irf_tidy.rds")

ORD_A_PATH <- "SPY \u2192 XLE \u2192 XOM \u2192 CVX \u2192 COP"
ORD_B_PATH <- "COP \u2192 CVX \u2192 XOM \u2192 XLE \u2192 SPY"

BLUE   <- "#3B4CC0"
BLUEFI <- "#8E9BE8"
RED    <- "#C0392B"

# -------------------------------------------------------------
# Panel labelling.
#
# Each panel carries the full transmission it depicts, in the convention
# of the original study: CVX_r -> XOM_r reads as the response of Exxon
# Mobil returns to a Chevron shock. This is more legible than separate
# row and column strips, which require cross referencing two margins to
# identify any given panel.
# -------------------------------------------------------------
SERIES <- c("COP", "CVX", "XOM", "XLE", "SPY")

add_pair_label <- function(df) {
  lv <- as.vector(t(outer(SERIES, SERIES,
                          function(i, j) paste0(i, "_r \u2192 ", j, "_r"))))
  df %>%
    mutate(pair_label = factor(paste0(impulse, "_r \u2192 ", response, "_r"),
                               levels = lv))
}

theme_irf <- function(base_size = 9) {
  theme_grey(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "#EDEDED", colour = NA),
      panel.grid.major = element_line(colour = "white", linewidth = 0.35),
      panel.grid.minor = element_line(colour = "white", linewidth = 0.18),
      panel.spacing    = unit(0.3, "lines"),
      strip.background = element_rect(fill = "#D6D6D6", colour = NA),
      strip.text       = element_text(size = base_size - 1, face = "bold",
                                      colour = "grey15",
                                      margin = margin(2.5, 2, 2.5, 2)),
      plot.title       = element_text(face = "bold", size = base_size + 4),
      plot.subtitle    = element_text(size = base_size + 1, colour = "grey30"),
      plot.caption     = element_text(size = base_size - 1, colour = "grey40",
                                      hjust = 0),
      axis.text        = element_text(size = base_size - 2.5),
      legend.position  = "bottom",
      legend.title     = element_blank()
    )
}

# -------------------------------------------------------------
# Generic panel plotter
# -------------------------------------------------------------
plot_panels <- function(df, title, subtitle, caption,
                        drop_impact = FALSE, ncol = 5) {
  d <- df %>% filter(horizon <= MAX_H)
  if (drop_impact) d <- d %>% filter(horizon >= 1)
  
  ggplot(d, aes(x = horizon, y = estimate)) +
    geom_hline(yintercept = 0, colour = "grey35",
               linetype = "dashed", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = BLUEFI, alpha = 0.55) +
    geom_line(aes(y = lower), colour = BLUE,
              linetype = "dotted", linewidth = 0.25) +
    geom_line(aes(y = upper), colour = BLUE,
              linetype = "dotted", linewidth = 0.25) +
    geom_line(colour = BLUE, linewidth = 0.65) +
    facet_wrap(~ pair_label, ncol = ncol, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, MAX_H, by = 2)) +
    labs(x = "Horizon (trading days after shock)",
         y = "Response",
         title = title, subtitle = subtitle, caption = caption) +
    theme_irf()
}

# =============================================================
# 1. Marginal responses, Ordering A
# =============================================================
mar_a <- irf_all %>% filter(ordering == "Ordering A") %>% add_pair_label()

p_marginal <- plot_panels(
  mar_a,
  title    = "Orthogonalised impulse responses",
  subtitle = paste0("Ordering A: ", ORD_A_PATH),
  caption  = paste(
    "Each panel is labelled impulse followed by response, so that",
    "CVX_r \u2192 XOM_r depicts the response of Exxon Mobil returns to a",
    "Chevron shock.\nShaded regions are bootstrapped 95 percent confidence",
    "bands from 500 replications. Vertical scales are free."
  )
)

ggsave("output/figures/irf_marginal_responses.png", p_marginal,
       width = 10.5, height = 9, dpi = 300, bg = "white")

# =============================================================
# 2. Cumulative responses
#
# The cumulative response at horizon h is the running sum of marginal
# responses from 0 to h. Because returns are log price differences, this
# sum is the response of the log price LEVEL, and its limit is the
# permanent effect of the shock on price.
#
# Bands come from re-running the bootstrap with cumulative = TRUE rather
# than from summing the marginal bands, since the confidence interval of
# a sum is not the sum of the confidence intervals.
# =============================================================
returns <- readRDS("data/processed/returns.rds")
ret_mat <- as.matrix(as.data.frame(returns))
ord_a   <- c("SPY_r", "XLE_r", "XOM_r", "CVX_r", "COP_r")

var_a <- VAR(ret_mat[, ord_a], p = 2, type = "const")

set.seed(42)
irf_cum_raw <- irf(var_a, n.ahead = 10, ortho = TRUE, cumulative = TRUE,
                   boot = TRUE, runs = 500, ci = 0.95)

tidy_irf <- function(obj, ordering_label) {
  do.call(rbind, lapply(names(obj$irf), function(imp) {
    pt <- as.data.frame(obj$irf[[imp]])
    lo <- as.data.frame(obj$Lower[[imp]])
    up <- as.data.frame(obj$Upper[[imp]])
    do.call(rbind, lapply(colnames(pt), function(resp) {
      data.frame(impulse  = gsub("_r$", "", imp),
                 response = gsub("_r$", "", resp),
                 horizon  = seq_len(nrow(pt)) - 1,
                 estimate = pt[[resp]], lower = lo[[resp]], upper = up[[resp]],
                 ordering = ordering_label, stringsAsFactors = FALSE)
    }))
  }))
}

cum_a <- tidy_irf(irf_cum_raw, "Ordering A")
saveRDS(cum_a, "data/processed/irf_cumulative_tidy.rds")

p_cumulative <- plot_panels(
  add_pair_label(cum_a),
  title    = "Cumulative orthogonalised impulse responses",
  subtitle = paste0("Ordering A: ", ORD_A_PATH),
  caption  = paste(
    "The cumulative response is the running sum of the marginal responses.",
    "Because returns are log price differences,\nthis sum traces the response",
    "of the log price level, and the level at which each panel settles is the",
    "permanent effect of the shock on price."
  )
)

ggsave("output/figures/irf_cumulative_responses.png", p_cumulative,
       width = 10.5, height = 9, dpi = 300, bg = "white")

# -------------------------------------------------------------
# Long run cumulative effects, own shocks
# -------------------------------------------------------------
longrun <- cum_a %>%
  filter(horizon == 10, impulse == response) %>%
  transmute(Series = impulse,
            `Permanent effect` = estimate,
            `Lower 95\\%` = lower,
            `Upper 95\\%` = upper,
            `Percent` = 100 * estimate) %>%
  arrange(desc(`Permanent effect`))

cat("\n=== Permanent effect of own shock on log price level ===\n")
print(as.data.frame(longrun), digits = 4)

print(
  xtable(as.data.frame(longrun),
         digits = c(0, 0, 5, 5, 5, 2),
         caption = "Permanent effect of a one standard deviation own shock on the log price level, given by the cumulative orthogonalised response at a ten day horizon under Ordering A. Convergence to a nonzero constant rather than to zero is the impulse response counterpart of the unit root established in Section~\\ref{sec:var}.",
         label = "tab:irf_longrun"),
  file = "output/tables/irf_longrun_cumulative.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "htbp",
  sanitize.colnames.function = identity
)

# =============================================================
# 3. Decay dynamics, impact horizon removed
# =============================================================
p_decay <- plot_panels(
  mar_a, drop_impact = TRUE,
  title    = "Impulse responses excluding the impact horizon",
  subtitle = paste0("Ordering A: ", ORD_A_PATH),
  caption  = paste(
    "With horizon zero omitted the vertical scale is set by the decay path",
    "rather than by the contemporaneous spike,\nand the confidence bands",
    "become legible. Band width is roughly constant across horizons while",
    "point estimates fall by an order of magnitude."
  )
)

ggsave("output/figures/irf_decay_dynamics_excluding_impact.png", p_decay,
       width = 10.5, height = 9, dpi = 300, bg = "white")

# =============================================================
# 4. Ordering comparison. Requires both orderings: this figure is the
#    evidence that impact responses are assumption driven.
# =============================================================
focus <- c("COP", "CVX", "XOM")

comp <- irf_all %>%
  filter(impulse %in% focus, response %in% focus,
         impulse != response, horizon <= MAX_H) %>%
  add_pair_label()

p_comp <- ggplot(comp, aes(x = horizon, y = estimate,
                           colour = ordering, fill = ordering)) +
  geom_hline(yintercept = 0, colour = "grey35",
             linetype = "dashed", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.1) +
  facet_wrap(~ pair_label, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = c("Ordering A" = BLUE, "Ordering B" = RED),
                      labels = c(paste0("Ordering A  (", ORD_A_PATH, ")"),
                                 paste0("Ordering B  (", ORD_B_PATH, ")"))) +
  scale_fill_manual(values = c("Ordering A" = BLUE, "Ordering B" = RED),
                    labels = c(paste0("Ordering A  (", ORD_A_PATH, ")"),
                               paste0("Ordering B  (", ORD_B_PATH, ")"))) +
  scale_x_continuous(breaks = seq(0, MAX_H, by = 2)) +
  labs(x = "Horizon (trading days after shock)", y = "Response",
       title = "Cholesky ordering sensitivity",
       subtitle = "Cross responses among the individual securities under both orderings",
       caption = paste(
         "Separation at horizon zero is attributable to the ordering",
         "assumption rather than to the data.\nFrom the first day onward the",
         "two orderings coincide."
       )) +
  theme_irf() +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

ggsave("output/figures/irf_ordering_comparison_energy_pairs.png", p_comp,
       width = 10, height = 7, dpi = 300, bg = "white")

# =============================================================
# 5. Impact contrast
# =============================================================
impact <- irf_all %>%
  filter(horizon == 0, impulse %in% focus, response %in% focus,
         impulse != response) %>%
  mutate(pair = paste0(impulse, "_r \u2192 ", response, "_r"))

p_impact <- ggplot(impact, aes(x = reorder(pair, estimate),
                               y = estimate, fill = ordering)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65,
           colour = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.75),
                width = 0.2, linewidth = 0.4, colour = "grey25") +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = c("Ordering A" = BLUE, "Ordering B" = RED)) +
  labs(x = NULL, y = "Impact response (horizon 0)",
       title = "Impact responses reverse entirely when the ordering is reversed",
       subtitle = "Every pair is exactly zero under one ordering and nonzero under the other",
       caption = paste(
         "Zeros are structural rather than estimated: a lower triangular",
         "Cholesky factor defines the contemporaneous\ncoefficient running",
         "from a later ordered variable to an earlier ordered one to be zero."
       )) +
  theme_irf(base_size = 10) +
  theme(panel.grid.major.y = element_blank())

ggsave("output/figures/irf_impact_response_ordering_contrast.png", p_impact,
       width = 9, height = 5.5, dpi = 300, bg = "white")

message("Five figures and one table written.")


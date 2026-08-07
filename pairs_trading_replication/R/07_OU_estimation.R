# =============================================================
# 09_ou_estimation.R
# Ornstein-Uhlenbeck parameter estimation on the COP-XOM spread,
# replicating the original study's approach with a fixed
# cointegrating coefficient estimated on the in-sample period.
#
# Reads:  data/raw/prices_raw.rds
# Writes: output/tables/ou_parameters.tex
#         data/processed/spread_fixed.rds
# =============================================================

library(urca)
library(xtable)
library(zoo)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

prices <- readRDS("data/raw/prices_raw.rds")
dates  <- as.Date(as.POSIXct(zoo::index(prices),
                             origin = "1970-01-01", tz = "UTC"))
lp     <- log(as.data.frame(prices))

# -------------------------------------------------------------
# In-sample period matches the original: 1999-12-31 to 2013-04-30.
# The cointegrating coefficient and OU parameters are estimated here
# and then held fixed, exactly as the original does.
# -------------------------------------------------------------
IN_END <- as.Date("2013-04-30")
in_idx <- dates <= IN_END

# -------------------------------------------------------------
# Step 1: cointegrating coefficient on the in-sample window.
#
# COP on XOM, matching the pair and direction the original selected.
# Recall from Chapter 5 that this residual does NOT reject the null of
# no cointegration at the Engle-Granger 5% level even on this window;
# the strategy is built regardless, to demonstrate the consequence.
# -------------------------------------------------------------
fit    <- lm(lp$COP[in_idx] ~ lp$XOM[in_idx])
tau_e  <- unname(coef(fit)[1])
beta_e <- unname(coef(fit)[2])

# Full-sample spread using the fixed in-sample coefficient
spread <- lp$COP - (tau_e + beta_e * lp$XOM)

spread_df <- data.frame(date = dates, spread = spread,
                        in_sample = in_idx)
saveRDS(spread_df, "data/processed/spread_fixed.rds")

# -------------------------------------------------------------
# Step 2: OU parameters via AR(1) on the in-sample spread.
#
# The OU process de_t = -theta (e_t - mu) dt + sigma dW_t has exact
# discrete solution  e_t = C + B e_{t-1} + eps,  from which:
#   B     = e^{-theta * tau}          (AR1 coefficient)
#   theta = -ln(B) / tau              (speed of reversion)
#   mu    = C / (1 - B)               (equilibrium level)
#   sigma_OU^2 = 2 theta Var(eps) / (1 - B^2)
#   sigma_eq   = sigma_OU / sqrt(2 theta)   (stationary sd)
#   half-life  = ln(2) / theta         (time to close half a deviation)
#
# tau = 1 trading day.
# -------------------------------------------------------------
tau <- 1
e   <- spread[in_idx]
e_lag <- e[-length(e)]
e_now <- e[-1]

ar1 <- lm(e_now ~ e_lag)
C   <- unname(coef(ar1)[1])
B   <- unname(coef(ar1)[2])
resid_var <- summary(ar1)$sigma^2

theta    <- -log(B) / tau
mu_e     <- C / (1 - B)
sigma_OU <- sqrt(2 * theta * resid_var / (1 - B^2))
sigma_eq <- sigma_OU / sqrt(2 * theta)
half_life <- log(2) / theta

# -------------------------------------------------------------
# Report
# -------------------------------------------------------------
cat("=== Fixed cointegrating coefficient, in-sample ===\n")
cat(sprintf("  tau_e (intercept): %.4f\n", tau_e))
cat(sprintf("  beta_e (COP on XOM): %.4f\n\n", beta_e))

cat("=== OU parameters, in-sample ===\n")
cat(sprintf("  AR(1) coefficient B: %.5f\n", B))
cat(sprintf("  theta (reversion speed): %.5f per day\n", theta))
cat(sprintf("  mu_e (equilibrium): %.5f\n", mu_e))
cat(sprintf("  sigma_eq (stationary sd): %.5f\n", sigma_eq))
cat(sprintf("  HALF-LIFE: %.1f days  (%.2f trading years)\n",
            half_life, half_life / 252))

ou_tbl <- data.frame(
  Parameter = c("Cointegrating coefficient $\\hat{\\beta}$",
                "AR(1) coefficient $B$",
                "Reversion speed $\\theta$ (per day)",
                "Equilibrium level $\\mu_e$",
                "Stationary sd $\\sigma_{eq}$",
                "Half-life (days)",
                "Half-life (trading years)"),
  Value = c(sprintf("%.4f", beta_e),
            sprintf("%.5f", B),
            sprintf("%.5f", theta),
            sprintf("%.5f", mu_e),
            sprintf("%.5f", sigma_eq),
            sprintf("%.1f", half_life),
            sprintf("%.2f", half_life / 252)),
  check.names = FALSE
)

print(
  xtable(ou_tbl,
         caption = paste("Ornstein-Uhlenbeck parameters for the COP-XOM spread,",
                         "estimated on the in-sample period 1999 to April 2013 with",
                         "the cointegrating coefficient held fixed thereafter. The",
                         "half-life is the time for the process to close half of a",
                         "deviation from equilibrium."),
         label = "tab:ou_params"),
  file = "output/tables/ou_parameters.tex",
  booktabs = TRUE, include.rownames = FALSE,
  caption.placement = "top", table.placement = "H",
  sanitize.text.function = identity
)

# Store parameters for the backtester
saveRDS(list(tau_e = tau_e, beta_e = beta_e,
             theta = theta, mu_e = mu_e, sigma_eq = sigma_eq,
             half_life = half_life, in_end = IN_END),
        "data/processed/ou_params.rds")

message("OU estimation complete. Note the half-life before proceeding.")

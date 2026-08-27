# eff_moments.R
# Compute every empirical moment used by the paper from the clean EFF 2022 data,
# with full multiply-imputed inference (bootstrap replicate-weight variance
# combined across implicates by Rubin's rules), and write a single target table.
#
# Each row of the output states the estimate, its standard error, the sample it
# is computed on, the role it plays (calibration target, model input, stylized
# fact, or context), and the model object it disciplines. Reading the table top
# to bottom is the audit that the model is calibrated against the right objects.
#
# Input : data/clean/hogar_imp{1..5}.rds        (from build_eff.R)
#         data/raw/replicate_weights_2022.csv    (optional; enables SEs)
# Output: output/eff_targets.csv
#
# Run from the pipeline root:  Rscript eff_moments.R

suppressPackageStartupMessages(library(data.table))
source("eff_functions.R")

# --- Configuration -----------------------------------------------------------
IMP_DIR  <- "data/clean"
WEIGHTS  <- "data/raw/replicate_weights_2022.csv"
OUT      <- "output/eff_targets.csv"
dir.create("output", showWarnings = FALSE, recursive = TRUE)

HTM_DIV             <- 24     # hand-to-mouth threshold: liquid <= income / 24
WEALTH_DIV          <- 24     # wealthy-HtM threshold:   illiquid > income / 24
SHOCK               <- 0.04   # mechanical stress: +400 bp, full pass-through
REF_EURIBOR_1Y_2022 <- 1.09   # ECB Data Portal, 1-year Euribor, 2022 average (%)

# --- Load and derive ---------------------------------------------------------
annuity <- function(M, r_m, n) fifelse(r_m == 0, M / n, M * r_m / (1 - (1 + r_m)^(-n)))

prep <- function(d) {
  d <- as.data.table(d)
  d[, illiquid   := riquezanet - liquid]
  d[, iliq_income := fifelse(income > 0, illiquid / income, NA_real_)]
  d[, nw_income   := fifelse(income > 0, riquezanet / income, NA_real_)]
  d[, htm         := income_ok & (liquid <= income / HTM_DIV)]
  d[, htm_cat     := fifelse(!income_ok, NA_character_,
                      fifelse(!htm, "no",
                       fifelse(illiquid > income / WEALTH_DIV, "wealthy", "poor")))]
  # Mechanical +400 bp stress on the variable-rate balance.
  d[, r0    := p2_13_1 / 100]
  d[, n_m   := p2_17_1 * 12]
  d[, dpay  := fifelse(M_var > 0 & !is.na(n_m) & n_m > 0 & !is.na(r0),
                       annuity(M_var, (r0 + SHOCK) / 12, n_m) -
                       annuity(M_var,  r0          / 12, n_m), NA_real_)]
  d[, dsr_pre  := mortgage_pay * 12 / income]
  d[, dsr_post := (mortgage_pay * 12 + dpay * 12) / income]
  d
}

imp <- lapply(1:5, function(m) prep(readRDS(file.path(IMP_DIR, sprintf("hogar_imp%d.rds", m)))))
N   <- nrow(imp[[1]])

HAVE_SE <- file.exists(WEIGHTS)
if (HAVE_SE) {
  message("Loading replicate weights for standard errors (this reads a large file)...")
  repw <- load_replicate_weights(WEIGHTS)
} else {
  message("Replicate weights not found; reporting point estimates without SEs.")
}

# --- Inference wrappers ------------------------------------------------------
# A moment is an estimator est(d, w) -> scalar. The point estimate averages it
# over implicates on the population weight; the SE adds bootstrap sampling
# variance (replicate weights) and Rubin between-implicate variance.
estimate <- function(est) {
  if (HAVE_SE) {
    r <- rubin_apply(imp, function(d) replicate_var(d, repw, est))
    c(value = r$estimate, se = r$se)
  } else {
    c(value = mean(vapply(imp, function(d) est(d, d$facine3), numeric(1))), se = NA_real_)
  }
}

ratio_est <- function(num, den) function(d, w) {
  ok <- is.finite(d[[num]]) & is.finite(d[[den]]) & is.finite(w)
  sum(w[ok] * d[[num]][ok]) / sum(w[ok] * d[[den]][ok])
}

# --- Output table ------------------------------------------------------------
ROWS <- list()
emit <- function(key, value, se = NA, n = NA, unit = "", role = "", object = "", def = "") {
  ROWS[[length(ROWS) + 1L]] <<- data.table(key, value = round(value, 4),
    se = round(se, 4), n = n, unit = unit, role = role, model_object = object, definition = def)
  cat(sprintf("  %-34s %9.4f %s%-8s  [%s] %s\n", key, value,
              ifelse(is.na(se), "", sprintf("(%.4f) ", se)), "", role, object))
}
emit_scalar <- function(key, est, n, unit, role, object, def) {
  r <- estimate(est); emit(key, r[["value"]], r[["se"]], n, unit, role, object, def)
}
banner <- function(x) cat("\n===== ", x, " =====\n", sep = "")

# =============================================================================
banner("Calibration targets (household model)")
# Three of the four SMM targets are EFF moments; the fourth (annual MPC = 0.25)
# is external (empirical MPC literature) and is recorded here for completeness.
emit_scalar("liquid_income_p50",
  function(d, w) { ok <- d$income_ok; w_quantile(d$liq_income[ok], w[ok], .5) },
  n = sum(imp[[1]]$income_ok), unit = "ratio", role = "target",
  object = "liquid-asset target (chi1 / beta split)",
  def = "median liquid assets / annual income")
emit_scalar("illiquid_income_p50",
  function(d, w) w_quantile(d$iliq_income, w, .5),
  n = sum(imp[[1]]$income_ok), unit = "ratio", role = "target",
  object = "illiquid-asset target (ra)",
  def = "median illiquid net wealth / annual income")
emit_scalar("htm_share",
  function(d, w) { ok <- d$income_ok; w_share((d$liquid <= d$income / HTM_DIV)[ok], w[ok]) },
  n = sum(imp[[1]]$income_ok), unit = "share", role = "target",
  object = "hand-to-mouth share (beta heterogeneity)",
  def = "liquid <= income/24")
emit("mpc_annual", 0.25, NA, NA, "ratio", "target (external)",
     "sigma_eps in SMM",
     "NOT from the EFF: HANK calibration norm (Kaplan-Violante, 0.15-0.25), consistent with the 22% HtM; HFCS direct survey is higher (~0.5). See docs/sources.md")

# =============================================================================
banner("Mortgage structure (model inputs)")
emit_scalar("mortgage_share",
  function(d, w) w_share(d$M_tot > 0, w),
  n = N, unit = "share", role = "input", object = "mortgage prevalence",
  def = "households with dwelling debt")
emit_scalar("variable_share_among_mortgaged",
  function(d, w) { s <- d$M_tot > 0; w_share((d$M_var > 0)[s], w[s]) },
  n = sum(imp[[1]]$M_tot > 0), unit = "share", role = "input", object = "PVAR (variable-rate share)",
  def = "variable-rate share of the mortgaged stock")
emit_scalar("Mvar_income_p50_exposed",
  function(d, w) { s <- d$M_var > 0; w_quantile(d$exp_income[s], w[s], .5) },
  n = sum(imp[[1]]$M_var > 0), unit = "ratio", role = "input", object = "MV (exposure)",
  def = "median variable balance / income among exposed")
emit_scalar("mortgage_rate_level_p50",
  function(d, w) { s <- d$mortgage_type %in% c("variable", "mixed"); w_quantile(d$p2_13_1[s], w[s], .5) },
  n = sum(imp[[1]]$mortgage_type %in% c("variable", "mixed")), unit = "% p.a.",
  role = "input", object = "mortgage rate level (NOT a spread)",
  def = "median current interest rate, variable/mixed loans (p2_13)")
emit_scalar("mortgage_term_p50",
  function(d, w) { s <- d$mortgage_type %in% c("variable", "mixed"); w_quantile(d$p2_17_1[s], w[s], .5) },
  n = sum(imp[[1]]$mortgage_type %in% c("variable", "mixed")), unit = "years",
  role = "input", object = "residual maturity / amortization",
  def = "median remaining term, variable/mixed loans (p2_17)")

# =============================================================================
banner("Portfolio composition (two-asset justification)")
emit_scalar("net_wealth_income_p50",
  function(d, w) w_quantile(d$nw_income, w, .5),
  n = sum(imp[[1]]$income_ok), unit = "ratio", role = "context", object = "total wealth target",
  def = "median net wealth / income")
emit_scalar("share_real_over_gross_wealth",
  ratio_est("actreales", "riquezabr"), n = N, unit = "share", role = "context",
  object = "illiquid = housing-heavy", def = "real assets / gross wealth")
emit_scalar("share_financial_over_gross_wealth",
  ratio_est("actfinanc", "riquezabr"), n = N, unit = "share", role = "context",
  object = "financial wealth", def = "financial assets / gross wealth")

# =============================================================================
banner("Hand-to-mouth split (Kaplan-Violante)")
emit_scalar("htm_poor_share",
  function(d, w) { ok <- !is.na(d$htm_cat); w_share((d$htm_cat == "poor")[ok], w[ok]) },
  n = sum(imp[[1]]$income_ok), unit = "share", role = "fact", object = "poor-HtM",
  def = "HtM with low illiquid wealth")
emit_scalar("htm_wealthy_share",
  function(d, w) { ok <- !is.na(d$htm_cat); w_share((d$htm_cat == "wealthy")[ok], w[ok]) },
  n = sum(imp[[1]]$income_ok), unit = "share", role = "fact", object = "wealthy-HtM (dominant -> two assets)",
  def = "HtM with sizeable illiquid wealth")

# =============================================================================
banner("Amplification object: extensive margin, not fragility")
# HtM by households (H0) vs weighted by euros of variable balance (H1). Ratio
# near 1 means exposed euros are not concentrated in fragile hands.
emit_scalar("htm_share_exposed_H0",
  function(d, w) { s <- d$M_var > 0 & d$income_ok; w_share(d$htm[s], w[s]) },
  n = sum(imp[[1]]$M_var > 0 & imp[[1]]$income_ok), unit = "share", role = "fact",
  object = "HtM among exposed (by households)", def = "H0")
emit_scalar("htm_share_exposed_H1",
  function(d, w) { s <- d$M_var > 0 & d$income_ok
                   sum(w[s] * d$htm[s] * d$M_var[s]) / sum(w[s] * d$M_var[s]) },
  n = sum(imp[[1]]$M_var > 0 & imp[[1]]$income_ok), unit = "share", role = "fact",
  object = "HtM among exposed (by euros of Mvar)", def = "H1; ratio H1/H0 ~ 1 -> extensive margin")
emit_scalar("neg_liquid_share",
  function(d, w) { ok <- d$income_ok; w_share((d$liquid < 0)[ok], w[ok]) },
  n = sum(imp[[1]]$income_ok), unit = "share", role = "fact", object = "borrowing marginal (b_low)",
  def = "households with negative net liquid (revolving debt)")

# =============================================================================
banner("Mechanical stress test (+400 bp, first round)")
stress_ok <- function(d) d$M_var > 0 & d$income_ok & is.finite(d$dsr_post)
emit_scalar("dsr_over40_pre",
  function(d, w) { s <- stress_ok(d); w_share((d$dsr_pre > 0.40)[s], w[s]) },
  n = sum(stress_ok(imp[[1]])), unit = "share", role = "fact", object = "debt-service tail (pre)",
  def = "share of exposed with DSR > 40% before shock")
emit_scalar("dsr_over40_post",
  function(d, w) { s <- stress_ok(d); w_share((d$dsr_post > 0.40)[s], w[s]) },
  n = sum(stress_ok(imp[[1]])), unit = "share", role = "fact", object = "debt-service tail (post)",
  def = "share of exposed with DSR > 40% after +400 bp")

# =============================================================================
banner("Point-only diagnostics (no SE)")
pt <- function(est) mean(vapply(imp, function(d) est(d, d$facine3), numeric(1)))
# Spread over the reference rate: resolves the p2_13 mislabel.
rate50 <- ROWS[[which(vapply(ROWS, function(r) r$key, "") == "mortgage_rate_level_p50")]]$value
emit("spread_over_euribor1Y_2022", rate50 - REF_EURIBOR_1Y_2022, NA, NA, "pp",
     "diagnostic", "mortgage spread kappa",
     sprintf("rate_p50 - 2022 avg 1Y Euribor (%.2f)", REF_EURIBOR_1Y_2022))
# Liquid/income dispersion (targets the spread of beta).
for (p in c(.10, .25, .50, .75, .90)) {
  v <- pt(function(d, w) { ok <- d$income_ok; w_quantile(d$liq_income[ok], w[ok], p) })
  emit(sprintf("liq_income_pctile_p%02d", 100 * p), v, NA, NA, "ratio", "check",
       "beta dispersion", "liquid/income percentile")
}
# Partial correlation of exposure and buffer, controlling income (imp1, unweighted).
d1 <- imp[[1]][M_var > 0 & income_ok == TRUE]
rx <- resid(lm(log(M_var) ~ log(income), d1))
ry <- resid(lm(liq_income ~ log(income), d1))
emit("partial_corr_Mvar_liqincome", cor(rx, ry, method = "spearman"), NA, nrow(d1),
     "corr", "fact", "Cov(MPC,exposure) ~ 0", "partial rank corr | income")

# --- Write -------------------------------------------------------------------
out <- rbindlist(ROWS)
fwrite(out, OUT)
cat(sprintf("\nWrote %d moments -> %s\n", nrow(out), OUT))

# =============================================================================
banner("Cross-checks (hard assertions)")
c1 <- imp[[1]]
# (1) Constructed variable balance vs the BdE mortgage aggregate deuhipv.
# Informative, not asserted: under M_DEF = "dwelling" the exposed balance also
# includes non-mortgage dwelling loans, so it need not equal strictly-mortgage
# debt. Under M_DEF = "mortgage" the two coincide.
vv     <- c1[mortgage_type == "variable" & is.finite(deuhipv)]
rho_mv <- cor(vv$M_var, vv$deuhipv)
cat(sprintf("  M_var vs deuhipv (variable): corr = %.4f  (informative; depends on M_DEF)\n", rho_mv))
# (2) Mortgage type only assigned where the household reports dwelling debt.
stopifnot(c1[!is.na(mortgage_type), all(deuhipv > 0 | M_tot > 0, na.rm = TRUE)])
# (3) Liquid assets reconcile with the BdE gross series (add card debt back).
rho <- cor(c1$liquid + fifelse(is.na(c1$p8_5a), 0, c1$p8_5a),
           c1$salcuentas + c1$p4_7_3, use = "complete.obs")
cat(sprintf("  corr(liquid_gross, BdE salcuentas + p4_7_3) = %.4f\n", rho))
stopifnot(rho > 0.95)
# (4) Implicate stability: variable-loan counts vary little across implicates.
vc <- vapply(imp, function(d) sum(d$mortgage_type == "variable", na.rm = TRUE), integer(1))
cat(sprintf("  variable count by implicate: %s\n", paste(vc, collapse = ", ")))
stopifnot((max(vc) - min(vc)) / mean(vc) < 0.05)
cat("\nAll cross-checks passed.\n")

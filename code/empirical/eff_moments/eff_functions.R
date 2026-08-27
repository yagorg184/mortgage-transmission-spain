# eff_functions.R
# Shared helpers for the EFF 2022 pipeline: recoding, weighted statistics, and
# multiply-imputed inference (bootstrap replicate-weight variance combined
# across implicates by Rubin's rules). Sourced by build_eff.R and eff_moments.R.

suppressPackageStartupMessages(library(data.table))

# --- Recoding ----------------------------------------------------------------
# The EFF encodes non-response as negative values; these are not economic
# quantities and are set to NA before any computation.
neg_to_na <- function(x, codes = NULL) {
  if (is.null(codes)) x[x < 0] <- NA else x[x %in% codes] <- NA
  x
}

clean_neg <- function(df, cols) {
  df[cols] <- lapply(df[cols], neg_to_na)
  df
}

# --- Weighted statistics -----------------------------------------------------
# All population statistics use the main household weight (facine3) unless a
# replicate weight is passed in its place.
w_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  sum(x[ok] * w[ok]) / sum(w[ok])
}

w_share <- function(cond, w) {
  ok <- !is.na(cond) & is.finite(w) & w > 0
  sum(w[ok] * cond[ok]) / sum(w[ok])
}

w_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  approx(cumsum(w) / sum(w), x, xout = probs, rule = 2, ties = "ordered")$y
}

w_median <- function(x, w) w_quantile(x, w, 0.5)

# --- Multiply-imputed inference ----------------------------------------------
# Point estimates use the population weight; the within-implicate sampling
# variance comes from the 1,000 BdE bootstrap replicate weights; the five
# implicates are combined by Rubin's rules.
REP_R <- 1000L

load_replicate_weights <- function(path) {
  cols <- c("h_2022", paste0("wt3r_", seq_len(REP_R)))
  fread(path, sep = ";", select = cols, showProgress = FALSE)
}

# Bootstrap sampling variance of a weighted estimator within one implicate.
# est_fun(d, w) returns the statistic given the implicate and a weight vector.
replicate_var <- function(d, repw, est_fun, w_main = "facine3") {
  d <- merge(as.data.table(d), repw, by = "h_2022", all.x = TRUE)
  theta <- est_fun(d, d[[w_main]])
  reps  <- vapply(seq_len(REP_R),
                  function(r) est_fun(d, d[[paste0("wt3r_", r)]]),
                  numeric(1))
  list(estimate = theta, var = mean((reps - theta)^2))
}

# Rubin (1987): combine M point estimates and within-imputation variances.
rubin_combine <- function(estimates, variances) {
  M    <- length(estimates)
  qbar <- mean(estimates)
  Ubar <- mean(variances)
  B    <- if (M > 1) stats::var(estimates) else 0
  Tvar <- Ubar + (1 + 1 / M) * B
  list(estimate = qbar, se = sqrt(Tvar), within = Ubar, between = B, M = M)
}

# Apply an estimator returning list(estimate, var) to each implicate, combine.
rubin_apply <- function(imp_list, fun) {
  res <- lapply(imp_list, fun)
  rubin_combine(vapply(res, `[[`, numeric(1), "estimate"),
                vapply(res, `[[`, numeric(1), "var"))
}

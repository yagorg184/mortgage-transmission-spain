# build_eff.R
# Construct analysis-ready household data from the raw EFF 2022 files.
#
# Input : data/raw/otras_secciones_2022_imp{1..5}.csv
#         data/raw/databol{1..5}.csv
# Output: data/clean/hogar_imp{1..5}.rds   (one per implicate)
#         data/clean/codebook.csv
#
# Run from the pipeline root:  Rscript build_eff.R

suppressPackageStartupMessages(library(data.table))
source("eff_functions.R")

RAW   <- "data/raw"
CLEAN <- "data/clean"
dir.create(CLEAN, showWarnings = FALSE, recursive = TRUE)

# --- Design switch: definition of the exposed balance M ----------------------
# "dwelling": every loan secured on the main dwelling (all reprice on reset).
# "mortgage": only loans with a formal mortgage guarantee (p2_9 == 1).
M_DEF <- "dwelling"

LOANS <- 1:4
lc <- function(stub) paste0(stub, "_", LOANS)   # e.g. p2_12_1 .. p2_12_4

# Per-loan blocks come from the "other sections" file; household aggregates and
# portfolio totals come from the derived "databol" file.
vars_otras <- c("h_2022", "facine3", "hogarpanel", "p2_8a",
                lc("p2_9"), lc("p2_12"), lc("p2_14"), lc("p2_18"),
                "p2_13_1", "p2_17_1",
                "p4_7_1", "p4_7_2", "p4_7_3", "p8_5a",
                "p9_1", "p9_2")

vars_databol <- c("h_2022", "renthog21_eur22", "bage", "percrent", "percriq",
                  "deuhipv", "dpdtehipo", "salcuentas", "ptmos_tarj", "pagodeuda",
                  "nodur", "alim", "riquezanet",
                  "actreales", "actfinanc", "riquezabr")

# Columns whose negative codes are non-response and must become NA.
neg_cols <- c(lc("p2_9"), lc("p2_12"), lc("p2_14"), lc("p2_18"),
              "p2_13_1", "p2_17_1", "p4_7_1", "p4_7_2", "p4_7_3",
              "p8_5a", "p9_1", "p9_2")

message(sprintf("Building five implicates [M_DEF = %s]", M_DEF))

for (m in 1:5) {
  os <- fread(file.path(RAW, sprintf("otras_secciones_2022_imp%d.csv", m)),
              sep = ";", select = vars_otras, showProgress = FALSE)
  db <- fread(file.path(RAW, sprintf("databol%d.csv", m)),
              sep = ";", select = vars_databol, showProgress = FALSE)
  d <- merge(os, db, by = "h_2022", all.x = TRUE)
  d <- clean_neg(as.data.frame(d), intersect(neg_cols, names(d)))
  setDT(d)

  # --- Exposure from every dwelling loan -------------------------------------
  # p2_12 outstanding balance, p2_14 rate-reset flag (11 variable / 22 fixed),
  # p2_9 loan type (mortgage guarantee), p2_18 monthly payment.
  bal <- lc("p2_12"); rev <- lc("p2_14"); gar <- lc("p2_9"); pay <- lc("p2_18")
  B  <- as.matrix(d[, ..bal]); B[is.na(B)]  <- 0
  PY <- as.matrix(d[, ..pay]); PY[is.na(PY)] <- 0
  RV <- as.matrix(d[, ..rev])
  G  <- as.matrix(d[, ..gar])

  keep <- matrix(TRUE, nrow(d), length(LOANS))
  if (M_DEF == "mortgage") keep <- (G == 1)
  keep[is.na(keep)] <- FALSE
  is_var <- (RV == 11) & keep; is_var[is.na(is_var)] <- FALSE
  is_fix <- (RV == 22) & keep; is_fix[is.na(is_fix)] <- FALSE

  d[, M_var          := rowSums(B * is_var)]   # variable-rate dwelling balance
  d[, M_fix          := rowSums(B * is_fix)]   # fixed-rate dwelling balance
  d[, M_tot          := rowSums(B * keep)]     # total dwelling balance
  d[, mortgage_pay   := rowSums(PY * keep)]    # total monthly dwelling payment
  d[, mortgage_type  := fifelse(M_var > 0 & M_fix > 0, "mixed",
                         fifelse(M_var > 0, "variable",
                          fifelse(M_fix > 0, "fixed", NA_character_)))]

  # --- Liquidity, income, ratios ---------------------------------------------
  # Liquid assets are deposits and cash net of revolving card debt (p8_5a).
  d[, liquid     := rowSums(cbind(p4_7_1, p4_7_2, p4_7_3), na.rm = TRUE) -
                    fifelse(is.na(p8_5a), 0, p8_5a)]
  d[, income     := renthog21_eur22]
  d[, income_ok  := income > 0]
  d[, liq_income := fifelse(income > 0, liquid / income, NA_real_)]
  d[, exp_income := fifelse(income > 0, M_var  / income, NA_real_)]

  saveRDS(d, file.path(CLEAN, sprintf("hogar_imp%d.rds", m)))
  message(sprintf("  implicate %d: %d households", m, nrow(d)))
}

# --- Living codebook from the Stata label files ------------------------------
build_codebook <- function() {
  dofiles <- list.files(RAW, pattern = "^etiquetas_.*\\.do$", full.names = TRUE)
  rows <- lapply(dofiles, function(f) {
    ln <- grep("^label var", readLines(f, warn = FALSE), value = TRUE)
    mm <- regmatches(ln, regexec('^label var\\s+(\\S+)\\s+"(.*)"', ln))
    do.call(rbind, lapply(mm, function(x)
      if (length(x) == 3) data.frame(var = x[2], label = x[3]) else NULL))
  })
  cb <- do.call(rbind, rows)
  fwrite(cb, file.path(CLEAN, "codebook.csv"))
  message(sprintf("Codebook: %d labelled variables", nrow(cb)))
}
tryCatch(build_codebook(),
         error = function(e) message("Codebook not written: ", conditionMessage(e)))

message("Build complete -> ", CLEAN)

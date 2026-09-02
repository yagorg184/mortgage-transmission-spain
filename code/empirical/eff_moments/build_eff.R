# build_eff.R
# Construct analysis-ready household data from the raw EFF 2022 files.
#
# Provides build_implicates(): raw CSVs -> a list of five cleaned implicate
# tables, in memory. eff_moments.R calls it directly, so there is no manual
# "clean" step: drop the raw files in data/raw/ and run eff_moments.R. Running
# this file on its own additionally caches the tables as data/clean/*.rds.
#
# Input : data/raw/otras_secciones_2022_imp{1..5}.csv
#         data/raw/databol{1..5}.csv
# Output (only when run standalone): data/clean/hogar_imp{1..5}.rds, codebook.csv
#
#   Rscript build_eff.R                      # build and cache to data/clean/
#   source("build_eff.R"); imp <- build_implicates()   # integrated, in memory

suppressPackageStartupMessages(library(data.table))
if (!exists("clean_neg")) source("eff_functions.R")

# Raw location. Override with the EFF_RAW_DIR environment variable; otherwise
# data/raw relative to the pipeline root (just unzip the raw folder there).
RAW_DIR_DEFAULT <- Sys.getenv("EFF_RAW_DIR", "data/raw")

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


# Build one implicate (m = 1..5) from the raw CSVs into a cleaned data.table.
build_one <- function(m, raw) {
  os <- fread(file.path(raw, sprintf("otras_secciones_2022_imp%d.csv", m)),
              sep = ";", select = vars_otras, showProgress = FALSE)
  db <- fread(file.path(raw, sprintf("databol%d.csv", m)),
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
  d
}

# --- Living codebook from the Stata label files ------------------------------
build_codebook <- function(raw, clean) {
  dofiles <- list.files(raw, pattern = "^etiquetas_.*\\.do$", full.names = TRUE)
  rows <- lapply(dofiles, function(f) {
    ln <- grep("^label var", readLines(f, warn = FALSE), value = TRUE)
    mm <- regmatches(ln, regexec('^label var\\s+(\\S+)\\s+"(.*)"', ln))
    do.call(rbind, lapply(mm, function(x)
      if (length(x) == 3) data.frame(var = x[2], label = x[3]) else NULL))
  })
  cb <- do.call(rbind, rows)
  fwrite(cb, file.path(clean, "codebook.csv"))
  message(sprintf("Codebook: %d labelled variables", nrow(cb)))
}

# Build the five implicates in memory. If clean_dir is given, also cache them as
# rds (the old data/clean/ artefacts) plus the codebook.
build_implicates <- function(raw_dir = RAW_DIR_DEFAULT, clean_dir = NULL) {
  message(sprintf("Building five implicates from %s [M_DEF = %s]", raw_dir, M_DEF))
  imp <- lapply(1:5, function(m) {
    d <- build_one(m, raw_dir)
    message(sprintf("  implicate %d: %d households", m, nrow(d)))
    d
  })
  if (!is.null(clean_dir)) {
    dir.create(clean_dir, showWarnings = FALSE, recursive = TRUE)
    for (m in 1:5) saveRDS(imp[[m]], file.path(clean_dir, sprintf("hogar_imp%d.rds", m)))
    tryCatch(build_codebook(raw_dir, clean_dir),
             error = function(e) message("Codebook not written: ", conditionMessage(e)))
    message("Cached -> ", clean_dir)
  }
  imp
}

# Run standalone -> build and cache to data/clean/. When sourced by eff_moments.R
# (which sets EFF_LIB_ONLY) this guard skips the auto-run so only the function is
# defined.
if (!exists("EFF_LIB_ONLY")) invisible(build_implicates(clean_dir = "data/clean"))

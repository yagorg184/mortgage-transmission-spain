#!/usr/bin/env bash
# Master reproduction script.
#
# Steps 2-4 run out of the box from the committed calibration and public data.
# Step 1 (the EFF moments) requires the restricted EFF 2022 microdata; see
# DATA_AVAILABILITY.md. Place the raw files under data/eff_raw/ and the cleaned
# implicates under data/eff_clean/, then uncomment the block below.
set -e
cd "$(dirname "$0")"

# ---- Step 1: EFF moments (R) — requires the restricted EFF microdata ----------
# (cd code/empirical/eff_moments \
#    && Rscript build_eff.R \       # data/eff_raw  -> data/eff_clean
#    && Rscript eff_moments.R)      # data/eff_clean -> calibration/eff_targets.csv
echo "Step 1 (EFF moments): skipped — needs the EFF microdata. Using committed calibration/eff_targets.csv."

# ---- Step 2: household calibration (Python) ----------------------------------
# The full SMM is ~1.5 h ('python calibrate.py'); --refine reproduces the frozen
# calibration from eff_targets.csv in about a minute.
echo "Step 2: calibrating the household..."
( cd code/model && python calibrate.py --refine )

# ---- Step 3: structural headline results (Python) ----------------------------
echo "Step 3: structural results..."
( cd code/model && python results.py )

# ---- Step 4: macro validation SVAR (Python) ----------------------------------
echo "Step 4: SVAR macro validation..."
( cd code/empirical/macro_svar && python svar.py )

echo "Done."

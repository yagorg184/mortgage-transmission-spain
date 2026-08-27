#!/usr/bin/env python3
# run_all.py — cross-platform reproduction driver (Windows, macOS, Linux).
#
# Runs the household calibration (fast refine), the structural results, and the
# SVAR from the committed data and public inputs. Regenerating the EFF moments
# needs the restricted microdata (see DATA_AVAILABILITY.md) and is left to the R
# scripts in code/empirical/eff_moments/.
#
# Usage:  python run_all.py     (on Windows, if `python` is not found, try `py run_all.py`)

import os
import sys
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable  # reuse the interpreter that is running this script


def run(rel_dir, script, *args):
    print(f"\n>>> {script}  ({rel_dir})")
    subprocess.run([PY, script, *args], cwd=os.path.join(HERE, rel_dir), check=True)


def main():
    print("Step 1 (EFF moments): skipped — needs the EFF microdata.")
    print("                      Using the committed calibration/eff_targets.csv.")
    print("Step 2: calibrating the household (fast refine)...")
    run(os.path.join("code", "model"), "calibrate.py", "--refine")
    print("Step 3: structural headline results...")
    run(os.path.join("code", "model"), "results.py")
    print("Step 4: SVAR macro validation...")
    run(os.path.join("code", "empirical", "macro_svar"), "svar.py")
    print("\nDone.")


if __name__ == "__main__":
    main()

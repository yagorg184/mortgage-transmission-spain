#!/usr/bin/env python3
# run_all.py
# Reproduces every number, table and figure of the paper from the files in the
# repository. The EFF moments themselves need the restricted microdata and are
# produced separately (see README.md); their output is committed under
# calibration/.
#
# Usage:
#   python run_all.py               from the committed data and calibrations
#   python run_all.py --calibrate   also re-solve the baseline and alternative calibrations
#   python run_all.py --download    first download the public data again
#   python run_all.py --tests       finish with the test suite
#
# On Windows, `py run_all.py` if `python` is not on the path.

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join("code", "model")
EMPIRICAL = os.path.join("code", "empirical")


def run(folder, script, *args):
    t0 = time.time()
    print(f"\n>>> {os.path.join(folder, script)} {' '.join(args)}", flush=True)
    subprocess.run([sys.executable, script, *args], cwd=os.path.join(HERE, folder), check=True)
    print(f"    done in {time.time() - t0:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate", action="store_true", help="re-solve the calibrations (about 30 minutes)")
    ap.add_argument("--download", action="store_true", help="download the public data again")
    ap.add_argument("--tests", action="store_true", help="run the tests at the end")
    args = ap.parse_args()

    if args.download:
        run(os.path.join(EMPIRICAL, "public_data"), "download.py")
    run(os.path.join(EMPIRICAL, "io"), "import_content.py")
    run(os.path.join(EMPIRICAL, "macro_svar"), "svar.py")
    run(os.path.join(EMPIRICAL, "macro_svar"), "svar_annual.py")
    run(os.path.join(EMPIRICAL, "lp"), "lp.py")

    if args.calibrate:
        run(MODEL, "calibrate.py")
        run(MODEL, "calibrate.py", "--beta-patient", "0.962")
        run(MODEL, "calibrate.py", "--beta-patient", "0.965")
        run(MODEL, "calibrate.py", "--terciles")
    else:
        run(MODEL, "calibrate.py", "--check")

    run(MODEL, "results.py")
    run(MODEL, "robustness.py")
    run(MODEL, "numerics.py")
    run(os.path.join("code", "figures"), "make_figures.py")

    if args.tests:
        for t in ("test_mortgages.py", "test_equilibrium.py", "test_nesting.py"):
            run(os.path.join(MODEL, "tests"), t)
    print("\nAll steps finished. Model numbers: code/model/output/{results,robustness,numerics}.json; "
          "figures: paper/figures/.")


if __name__ == "__main__":
    main()

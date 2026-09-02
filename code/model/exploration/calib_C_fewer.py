
# calib_C_fewer.py
# Experiment C: fewer free parameters. Fix the two parameters that external
# evidence can pin -- the illiquid return r_a and the patient discount factor --
# and re-estimate only {beta_modal, pi_modal, chi1, sig_eps} against the four
# usual moments. Four free parameters vs four targets: exactly identified, the
# properly-disciplined version of the baseline (which had six free parameters).
#
# The two fixed values are set here and are meant to be changed and re-run:
#   RA_FIX      external real return on illiquid wealth
#   BETA_PAT_FIX  a "very patient" saver's discount factor (~1/(1+r))
#
#   python calib_C_fewer.py --quick
#   python calib_C_fewer.py
#   python calib_C_fewer.py --workers -1

import argparse
from scipy.optimize import differential_evolution
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # find _calib_common in code/model/
import _calib_common as CC

TAG = "fewer"

# ---- externally fixed (edit these and re-run to test sensitivity) --------------
RA_FIX = 0.035
BETA_PAT_FIX = 0.98

# free: {beta_modal, pi_modal, chi1, sig_eps}
BOUNDS = [(0.93, 0.975), (0.60, 1.00), (2.0, 8.0), (0.20, 0.45)]
PARAM_NAMES = ["beta_modal", "pi_modal", "chi1", "sig_eps"]

# four targets -> exactly identified with four free parameters. MPC kept here so
# the contrast with A is clean (A drops it; C keeps it under a leaner parameter set).
TARGETS = [("liquid_income_p50", 1.0), ("illiquid_income_p50", 0.7),
           ("htm_share", 1.0), ("mpc_annual", 0.8)]

REPORT = ["liquid_income_p50", "illiquid_income_p50", "htm_share", "mpc_annual",
          "htm_poor_share", "htm_wealthy_share", "net_wealth_income_p50",
          "liq_income_pctile_p10", "liq_income_pctile_p90"]


def unpack(x):
    bm, pim, chi1, se = x
    if BETA_PAT_FIX < bm:            # modal must stay more impatient than the fixed patient type
        return None
    return ([bm, BETA_PAT_FIX], [pim, 1 - pim], chi1, se, RA_FIX)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxiter", type=int, default=120)
    ap.add_argument("--popsize", type=int, default=18)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--every", type=int, default=15)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    if a.quick:
        a.maxiter, a.popsize, a.every = 5, 4, 1

    loss, _ = CC.make_loss(TARGETS, unpack)
    cb = CC.make_callback(TAG, loss, unpack, TARGETS, REPORT, every=a.every)
    print(f"[{TAG}] fixed: ra={RA_FIX}  beta_patient={BETA_PAT_FIX}")
    print(f"[{TAG}] free ({len(PARAM_NAMES)}): {PARAM_NAMES}  vs {len(TARGETS)} targets -> exactly identified")
    print(f"[{TAG}] maxiter={a.maxiter} popsize={a.popsize} tol={a.tol} workers={a.workers}\n")

    kw = dict(maxiter=a.maxiter, popsize=a.popsize, tol=a.tol, seed=a.seed,
              polish=True, disp=False, callback=cb, mutation=(0.5, 1.0), recombination=0.7)
    if a.workers != 1:
        kw.update(workers=a.workers, updating="deferred")
    res = differential_evolution(loss, BOUNDS, **kw)
    CC.final_report(TAG, res.x, loss, unpack, TARGETS, REPORT, PARAM_NAMES)


if __name__ == "__main__":
    main()


# calib_A_mpcfree.py
# Experiment A: drop the MPC from the loss. Recalibrate the six parameters to the
# three EFF *wealth* moments only (median liquid/income, median illiquid/income,
# hand-to-mouth share) and report where the annual MPC lands as an untargeted
# outcome. If it lands near 0.25, the MPC can leave the loss and become a
# validation number; if it drifts, it stays as a (soft) target.
#
#   python calib_A_mpcfree.py --quick          fast smoke test (~minutes)
#   python calib_A_mpcfree.py                  full run (tune --maxiter/--popsize)
#   python calib_A_mpcfree.py --workers -1     use all cores
#
# Progress is printed and written to output/calib_mpcfree_* every --every gens.

import argparse
from scipy.optimize import differential_evolution
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # find _calib_common in code/model/
import _calib_common as CC

TAG = "mpcfree"

# free: {beta_modal, beta_patient, pi_modal, chi1, sig_eps, ra}
BOUNDS = [(0.93, 0.975), (0.96, 0.995), (0.60, 1.00), (2.0, 8.0), (0.20, 0.45), (0.020, 0.045)]
PARAM_NAMES = ["beta_modal", "beta_patient", "pi_modal", "chi1", "sig_eps", "ra"]

# targeted moments (csv_key, weight) — no MPC
TARGETS = [("liquid_income_p50", 1.0), ("illiquid_income_p50", 0.7), ("htm_share", 1.0)]

# everything to print/save (targeted + untargeted checks, MPC first among checks)
REPORT = ["liquid_income_p50", "illiquid_income_p50", "htm_share",
          "mpc_annual", "htm_poor_share", "htm_wealthy_share",
          "net_wealth_income_p50", "liq_income_pctile_p10", "liq_income_pctile_p90"]


def unpack(x):
    bm, bp, pim, chi1, se, ra = x
    if bp < bm:                      # keep the patient type more patient
        return None
    return ([bm, bp], [pim, 1 - pim], chi1, se, ra)


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
    print(f"[{TAG}] targets: {[k for k, _ in TARGETS]}  (MPC untargeted)")
    print(f"[{TAG}] maxiter={a.maxiter} popsize={a.popsize} tol={a.tol} workers={a.workers}\n")

    kw = dict(maxiter=a.maxiter, popsize=a.popsize, tol=a.tol, seed=a.seed,
              polish=True, disp=False, callback=cb, mutation=(0.5, 1.0), recombination=0.7)
    if a.workers != 1:
        kw.update(workers=a.workers, updating="deferred")
    res = differential_evolution(loss, BOUNDS, **kw)
    CC.final_report(TAG, res.x, loss, unpack, TARGETS, REPORT, PARAM_NAMES)


if __name__ == "__main__":
    main()

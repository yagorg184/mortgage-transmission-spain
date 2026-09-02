
# calib_E_5v5.py
# Final two-type design, exactly identified (5 free vs 5 targets), ra external.
# The redundant net-wealth *level* is dropped; two wealth-*concentration* moments
# (top-10% and top-1% net wealth share, measured on the EFF) identify the patient
# type -- how many are rich (pi_modal -> top-10%) vs how rich (beta_patient -> top-1%).
#
#   sigma_eps    liquidity demand      -> median liquid/income
#   beta_modal   middle wealth level   -> median illiquid/income
#   chi1         liquid/illiquid wedge -> hand-to-mouth share
#   pi_modal     size of the rich tail -> top-10% net-wealth share
#   beta_patient intensity of the tail -> top-1%  net-wealth share
# External: ra (real illiquid return). MPC untargeted (reported).
#
#   python calib_E_5v5.py --quick
#   python calib_E_5v5.py --workers -1

import argparse
from scipy.optimize import differential_evolution
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # find _calib_common in code/model/
import _calib_common as CC

TAG = "5v5"

RA_FIX = 0.035

# Wealth-concentration targets measured on the EFF 2022 microdata (mean over the
# five implicates; net wealth = riquezanet, weight = facine3). Re-run eff_moments.R
# for the Rubin SEs; the point estimates are what the calibration needs.
CC._T_ALL.loc["net_wealth_top10_share"] = 0.5362
CC._T_ALL.loc["net_wealth_top1_share"] = 0.2114

# free: {beta_modal, beta_patient, pi_modal, chi1, sig_eps}
BOUNDS = [(0.93, 0.975), (0.96, 0.995), (0.40, 0.98), (2.0, 8.0), (0.20, 0.45)]
PARAM_NAMES = ["beta_modal", "beta_patient", "pi_modal", "chi1", "sig_eps"]

TARGETS = [("liquid_income_p50", 1.0),
           ("illiquid_income_p50", 0.8),
           ("htm_share", 1.0),
           ("net_wealth_top10_share", 1.0),
           ("net_wealth_top1_share", 0.7)]

REPORT = ["liquid_income_p50", "illiquid_income_p50", "htm_share",
          "htm_poor_share", "htm_wealthy_share",
          "net_wealth_income_p50", "net_wealth_top10_share", "net_wealth_top1_share",
          "mpc_annual"]


def unpack(x):
    bm, bp, pim, chi1, se = x
    if bp < bm:
        return None
    return ([bm, bp], [pim, 1 - pim], chi1, se, RA_FIX)


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
    print(f"[{TAG}] fixed: ra={RA_FIX}")
    print(f"[{TAG}] free ({len(PARAM_NAMES)}): {PARAM_NAMES}")
    print(f"[{TAG}] targets ({len(TARGETS)}): {[k for k, _ in TARGETS]}  (MPC untargeted)")
    print(f"[{TAG}] maxiter={a.maxiter} popsize={a.popsize} tol={a.tol} workers={a.workers}\n")

    kw = dict(maxiter=a.maxiter, popsize=a.popsize, tol=a.tol, seed=a.seed,
              polish=True, disp=False, callback=cb, mutation=(0.5, 1.0), recombination=0.7)
    if a.workers != 1:
        kw.update(workers=a.workers, updating="deferred")
    res = differential_evolution(loss, BOUNDS, **kw)
    CC.final_report(TAG, res.x, loss, unpack, TARGETS, REPORT, PARAM_NAMES)


if __name__ == "__main__":
    main()

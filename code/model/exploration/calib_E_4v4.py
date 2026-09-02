
# calib_E_4v4.py
# Final two-type design, exactly identified (4 free vs 4 targets). ra and the
# patient discount factor are external; the redundant net-wealth level is dropped
# in favour of one reachable wealth-concentration moment, the top-10% net-wealth
# share (measured on the EFF).
#
#   sigma_eps    liquidity demand      -> median liquid/income
#   beta_modal   middle wealth level   -> median illiquid/income
#   chi1         liquid/illiquid wedge -> hand-to-mouth share
#   pi_modal     size of the rich tail -> top-10% net-wealth share
# External: ra (real illiquid return); beta_patient (the accumulating type of the
#   beta-heterogeneity literature, ~0.98, verified to reproduce the top-10% share).
# MPC untargeted. Two extreme-tail moments are reported but NOT targeted, because
# a two-type model are not targeted: poor-HtM (bottom) and top-1% (very top).
#
#   python calib_E_4v4.py --quick
#   python calib_E_4v4.py --workers -1

import argparse
from scipy.optimize import differential_evolution
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # find _calib_common in code/model/
import _calib_common as CC

TAG = "4v4"

RA_FIX = 0.035
BETA_PAT_FIX = 0.981      # accumulating patient type; beta*(1+ra) ~ 1.016

# top-10% net wealth share, EFF 2022 microdata (mean over implicates).
CC._T_ALL.loc["net_wealth_top10_share"] = 0.5362
CC._T_ALL.loc["net_wealth_top1_share"] = 0.2114     # reported only (reported only)

# free: {beta_modal, pi_modal, chi1, sig_eps}
# Bounds narrowed around the known good region: illiquid wealth is steep in
# beta_modal, so a wide box makes differential evolution fall off the knife-edge.
BOUNDS = [(0.950, 0.958), (0.92, 0.985), (5.0, 8.0), (0.28, 0.42)]
PARAM_NAMES = ["beta_modal", "pi_modal", "chi1", "sig_eps"]

TARGETS = [("liquid_income_p50", 1.0),
           ("illiquid_income_p50", 0.8),
           ("htm_share", 1.0),
           ("net_wealth_top10_share", 1.0)]

REPORT = ["liquid_income_p50", "illiquid_income_p50", "htm_share",
          "htm_poor_share", "htm_wealthy_share",
          "net_wealth_income_p50", "net_wealth_top10_share", "net_wealth_top1_share",
          "mpc_annual"]


def unpack(x):
    bm, pim, chi1, se = x
    if BETA_PAT_FIX < bm:
        return None
    return ([bm, BETA_PAT_FIX], [pim, 1 - pim], chi1, se, RA_FIX)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxiter", type=int, default=120)
    ap.add_argument("--popsize", type=int, default=16)
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
    print(f"[{TAG}] targets: {[k for k, _ in TARGETS]}  (MPC, poor-HtM, top-1% reported only)")
    print(f"[{TAG}] maxiter={a.maxiter} popsize={a.popsize} tol={a.tol} workers={a.workers}\n")

    kw = dict(maxiter=a.maxiter, popsize=a.popsize, tol=a.tol, seed=a.seed,
              polish=True, disp=False, callback=cb, mutation=(0.5, 1.0), recombination=0.7)
    if a.workers != 1:
        kw.update(workers=a.workers, updating="deferred")
    res = differential_evolution(loss, BOUNDS, **kw)
    CC.final_report(TAG, res.x, loss, unpack, TARGETS, REPORT, PARAM_NAMES)


if __name__ == "__main__":
    main()

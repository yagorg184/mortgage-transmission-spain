
# calib_B_moretargets.py
# Experiment B: over-identify. Keep the six parameters but face them with more of
# the EFF's measured moments, so there are more targets than parameters and the
# leftover fit is a test rather than a guarantee. The MPC stays out of the loss
# (reported as a check).
#
# Targets (7 > 6 params):
#   median liquid/income, median illiquid/income, median net-wealth/income,
#   poor-HtM share, wealthy-HtM share, liquid/income p10, liquid/income p90.
# The HtM split (poor vs wealthy) uses illiquid > income/24, exactly as the EFF
# code, so it also disciplines the wealthy-hand-to-mouth mass directly.
#
#   python calib_B_moretargets.py --quick
#   python calib_B_moretargets.py
#   python calib_B_moretargets.py --workers -1

import argparse
from scipy.optimize import differential_evolution
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # find _calib_common in code/model/
import _calib_common as CC

TAG = "moretargets"

BOUNDS = [(0.93, 0.975), (0.96, 0.995), (0.60, 1.00), (2.0, 8.0), (0.20, 0.45), (0.020, 0.045)]
PARAM_NAMES = ["beta_modal", "beta_patient", "pi_modal", "chi1", "sig_eps", "ra"]

# targeted moments (csv_key, weight). Core wealth levels and the HtM split at 1.0;
# the liquid dispersion tails a bit lower (noisier, and partly redundant with p50).
TARGETS = [("liquid_income_p50", 1.0),
           ("illiquid_income_p50", 0.8),
           ("net_wealth_income_p50", 0.8),
           ("htm_poor_share", 1.0),
           ("htm_wealthy_share", 1.0),
           ("liq_income_pctile_p10", 0.5),
           ("liq_income_pctile_p90", 0.5)]

REPORT = ["liquid_income_p50", "liq_income_pctile_p10", "liq_income_pctile_p25",
          "liq_income_pctile_p75", "liq_income_pctile_p90",
          "illiquid_income_p50", "net_wealth_income_p50",
          "htm_share", "htm_poor_share", "htm_wealthy_share", "mpc_annual"]


def unpack(x):
    bm, bp, pim, chi1, se, ra = x
    if bp < bm:
        return None
    return ([bm, bp], [pim, 1 - pim], chi1, se, ra)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxiter", type=int, default=150)
    ap.add_argument("--popsize", type=int, default=20)
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
    print(f"[{TAG}] targets ({len(TARGETS)}): {[k for k, _ in TARGETS]}")
    print(f"[{TAG}] 6 params vs {len(TARGETS)} targets -> over-identified; MPC untargeted")
    print(f"[{TAG}] maxiter={a.maxiter} popsize={a.popsize} tol={a.tol} workers={a.workers}\n")

    kw = dict(maxiter=a.maxiter, popsize=a.popsize, tol=a.tol, seed=a.seed,
              polish=True, disp=False, callback=cb, mutation=(0.5, 1.0), recombination=0.7)
    if a.workers != 1:
        kw.update(workers=a.workers, updating="deferred")
    res = differential_evolution(loss, BOUNDS, **kw)
    CC.final_report(TAG, res.x, loss, unpack, TARGETS, REPORT, PARAM_NAMES)


if __name__ == "__main__":
    main()

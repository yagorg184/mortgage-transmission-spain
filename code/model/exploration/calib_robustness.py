
# calib_robustness.py
# Robustness of the 4v4 calibration to the two externally-fixed parameters
# (ra, beta_patient). For each (ra, beta_patient) on a small grid, the four free
# parameters {beta_modal, pi_modal, chi1, sig_eps} are re-optimised to the same
# four targets. If the loss stays ~0 and the UNTARGETED MPC stays ~0.25 across the
# grid, the exact values of ra and beta_patient are second-order: the data
# (top-10% share) pins the mass of the rich type (pi_modal), not the fixed values.
#
#   python calib_robustness.py                 default 3x3 grid (~20-40 min)
#   python calib_robustness.py --maxfev 300    more thorough re-optimisation
#   python calib_robustness.py --ra 0.03 0.035 0.04 --bp 0.976 0.981 0.984
#
# Reads the model from _calib_common.py (same directory). Prints each combo as it
# finishes and writes output/calib_robustness.csv.

import argparse, os, csv, time
import numpy as np
from scipy.optimize import minimize
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # find _calib_common in code/model/
import _calib_common as CC

# top-10% net wealth share, EFF 2022 (same as calib_E_4v4).
CC._T_ALL.loc["net_wealth_top10_share"] = 0.5362

TARGETS = [("liquid_income_p50", 1.0), ("illiquid_income_p50", 0.8),
           ("htm_share", 1.0), ("net_wealth_top10_share", 1.0)]
# free params: {beta_modal, pi_modal, chi1, sig_eps}. Narrow box around the known
# 4v4 optimum so the local re-optimisation stays on the beta_modal knife-edge.
BOUNDS = [(0.950, 0.958), (0.92, 0.985), (5.0, 8.0), (0.28, 0.42)]
X0 = [0.9533, 0.9757, 6.84, 0.3534]          # the 4v4 solution, used as the start


def make(ra, bp):
    def unpack(x):
        bm, pim, chi1, se = x
        if bp < bm:
            return None
        return ([bm, bp], [pim, 1 - pim], chi1, se, ra)
    return CC.Objective(TARGETS, unpack), unpack


def reopt(ra, bp, maxfev):
    obj, unpack = make(ra, bp)
    r = minimize(obj, X0, method="Powell", bounds=BOUNDS,
                 options={"maxfev": maxfev, "xtol": 1e-4, "ftol": 1e-7})
    m = CC.solve_moments(*unpack(r.x))
    return float(r.fun), r.x, m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ra", type=float, nargs="+", default=[0.030, 0.035, 0.040])
    ap.add_argument("--bp", type=float, nargs="+", default=[0.976, 0.981, 0.984])
    ap.add_argument("--maxfev", type=int, default=120)   # ~few min/combo; 60 is faster, 300 finer
    a = ap.parse_args()

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "calib_robustness.csv")
    f = open(path, "w", newline="")
    wr = csv.writer(f)
    wr.writerow(["ra", "beta_patient", "loss", "beta_modal", "pi_modal", "chi1", "sig_eps",
                 "liquid_p50", "illiquid_p50", "htm", "top10", "mpc_UNtargeted"])

    print(f"grid: ra={a.ra}  beta_patient={a.bp}   maxfev={a.maxfev}")
    print("targets: liquid 0.2632 | illiquid 3.5555 | htm 0.22 | top10 0.5362   (MPC not targeted)\n")
    print(f"{'ra':>6} {'b_pat':>6} {'loss':>9} | {'b_mod':>7} {'pi_mod':>7} {'chi1':>6} {'sig_e':>6} "
          f"| {'liq':>6} {'iliq':>6} {'htm':>6} {'top10':>6} | {'MPC':>6}  {'min'}")
    losses, mpcs = [], []
    for ra in a.ra:
        for bp in a.bp:
            t0 = time.time()
            loss, x, m = reopt(ra, bp, a.maxfev)
            losses.append(loss); mpcs.append(m["mpc_annual"])
            print(f"{ra:6.3f} {bp:6.3f} {loss:9.5f} | {x[0]:7.4f} {x[1]:7.4f} {x[2]:6.2f} {x[3]:6.4f} "
                  f"| {m['liquid_income_p50']:6.3f} {m['illiquid_income_p50']:6.2f} {m['htm_share']:6.3f} "
                  f"{m['net_wealth_top10_share']:6.3f} | {m['mpc_annual']:6.3f}  {(time.time()-t0)/60:4.1f}",
                  flush=True)
            wr.writerow([ra, bp, round(loss, 6), *[round(v, 5) for v in x],
                         round(m["liquid_income_p50"], 4), round(m["illiquid_income_p50"], 4),
                         round(m["htm_share"], 4), round(m["net_wealth_top10_share"], 4),
                         round(m["mpc_annual"], 4)])
            f.flush()
    f.close()
    print(f"\nVERDICT:  max loss = {max(losses):.5f}   MPC range = [{min(mpcs):.3f}, {max(mpcs):.3f}]")
    print("robust if max loss stays ~0 and the MPC range stays tight around 0.25.")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()

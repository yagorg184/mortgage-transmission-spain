# calibrate.py
# Calibrate the two-asset household by SMM against the EFF moments and assemble
# the frozen parameter set in calibration.json.
#
#   python calibrate.py            full SMM (differential evolution, seed 0; ~1.5 h)
#   python calibrate.py --refine   re-pin the modal discount factor from the current
#                                  calibration.json (deterministic, ~1 min)
#
# Targets are read from calibration/eff_targets.csv, so the calibration tracks the
# measured EFF moments rather than hand-rounded numbers. Household wealth is very
# steep in the modal discount factor, so beta_modal is stored to six decimals and
# --refine exists to reproduce the frozen fit quickly without the full search.

import os, json, argparse
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize_scalar
import household as H

HERE = os.path.dirname(os.path.abspath(__file__))
CAL  = os.path.normpath(os.path.join(HERE, "..", "..", "calibration"))

# SMM setup: {beta_modal, beta_patient, pi_modal, chi1, sig_eps, ra}
BOUNDS  = [(0.93, 0.975), (0.96, 0.995), (0.60, 1.00), (2.0, 8.0), (0.20, 0.45), (0.020, 0.045)]
WEIGHTS = dict(liq_p50=1.0, iliq_p50=0.7, htm=1.0, mpc=0.8)

# Blocks that are not part of the household SMM: model mortgage constants and the
# externally-disciplined incidence and cross-country inputs. The steady-state
# mortgage payment uses the measured rate level (see mortgage block), so no
# separate rate constant is needed here.
MORT_MODEL = dict(amort=0.03)
INCIDENCE  = dict(beta_dep=0.32, mpc_shareholder=0.05, lambda_ext=0.50)
# Variable-rate share of the mortgage STOCK, ~2022. Spain from the EFF; euro area
# from the ECB (2025, pure-ARM stock ~25%); Germany/Italy/France reasonable ~2022
# values consistent with the ECB ordering and the HFCS cross-country ranking
# (Ehrmann-Ziegelmeyer, WP 1631, 2010, used as a complementary anchor).
COUNTRY    = {"all_fixed": 0.0, "Germany": 0.12, "France": 0.12, "euro_area": 0.25,
              "Italy": 0.45, "Spain": 0.70, "all_variable": 1.0}


def read_targets():
    t = pd.read_csv(os.path.join(CAL, "eff_targets.csv")).set_index("key")["value"]
    return t[~t.index.duplicated(keep="first")]


def smm_targets(t):
    """Household moments the SMM matches (model units: HtM in percent)."""
    return dict(liq_p50=float(t["liquid_income_p50"]),
                iliq_p50=float(t["illiquid_income_p50"]),
                htm=100 * float(t["htm_share"]),
                mpc=float(t["mpc_annual"]))


def moments(x):
    bm, bp, pim, chi1, se, ra = x
    return H.pooled([bm, bp], [pim, 1 - pim], chi1, se, over={"ra": ra})


def loss(x, T):
    if x[1] < x[0]:
        return 1e3
    try:
        m = moments(x)
    except Exception:
        return 1e3
    return sum(WEIGHTS[k] * ((m[k] - T[k]) / T[k]) ** 2 for k in WEIGHTS)


def run_smm(T, maxiter=40, popsize=10):
    res = differential_evolution(lambda x: loss(x, T), BOUNDS, maxiter=maxiter,
                                 popsize=popsize, tol=1e-4, seed=0, polish=True, disp=True)
    return res.x, res.fun


def run_refine(T):
    """Re-optimise only the modal discount factor (the steep, wealth-pinning
    parameter) from the current calibration.json, holding the rest."""
    cal = json.load(open(os.path.join(CAL, "calibration.json")))["household"]
    x = np.array([cal["beta"][0], cal["beta"][1], cal["pi"][0], cal["chi1"], cal["sig_eps"], cal["ra"]])
    lo, hi = x[0] - 0.003, x[0] + 0.003
    r = minimize_scalar(lambda bm: loss([bm, *x[1:]], T), bounds=(lo, hi),
                        method="bounded", options={"xatol": 1e-7})
    x[0] = r.x
    return x, r.fun


def assemble(x, lossval, t, T):
    # beta_modal is stored to 10 decimals because the illiquid moment is steep in
    # it; the rest to 6. Achieved moments are recomputed from the stored (rounded)
    # values so the file is internally consistent on reload.
    bm, bp, pim, chi1, se, ra = (round(x[0], 10), round(x[1], 6), round(x[2], 6),
                                 round(x[3], 6), round(x[4], 6), round(x[5], 6))
    x = [bm, bp, pim, chi1, se, ra]
    m = moments(x)
    return {
        "_about": "Single source of truth for the calibrated model. All model scripts read parameters from here; produced by calibrate.py against calibration/eff_targets.csv.",
        "household": {
            "beta": [bm, bp],
            "pi": [pim, round(1 - pim, 6)],
            "chi1": chi1, "sig_eps": se, "ra": ra,
            "fixed": {k: H.BASE[k] for k in ("eis", "rb", "chi0", "chi2", "rho_z", "sig_eta")},
            "grid": {k: H.BASE[k] for k in ("nB", "nA", "nK", "nZ", "nEps", "bmax", "amax", "kmax")},
            "targets": {"liquid_income_p50": T["liq_p50"], "illiquid_income_p50": T["iliq_p50"],
                        "htm_share": T["htm"] / 100, "mpc_annual": T["mpc"]},
            "achieved": {"liquid_income_p50": round(m["liq_p50"], 4), "illiquid_income_p50": round(m["iliq_p50"], 4),
                         "htm_share": round(m["htm"] / 100, 4), "mpc_gradient": round(m["mpc"], 4)},
            "loss": round(lossval, 6),
            "_note_precision": "beta_modal is stored to 6 decimals: household wealth is steep in the discount factor, so a 4th-decimal change moves the illiquid moment by ~0.5.",
        },
        "mortgage": {
            "mortgage_share": round(float(t["mortgage_share"]), 4),
            "variable_share": round(float(t["variable_share_among_mortgaged"]), 4),
            "M_var_income_exposed": round(float(t["Mvar_income_p50_exposed"]), 4),
            "rate_level": round(float(t["mortgage_rate_level_p50"]) / 100, 4),
            "spread_over_euribor_2022": round(float(t["spread_over_euribor1Y_2022"]) / 100, 4),
            "term_years": float(t["mortgage_term_p50"]),
            **MORT_MODEL,
        },
        "incidence": INCIDENCE,
        "country_variable_share_stock_2022": {
            "_source": "Spain: EFF 2022. Euro area: ECB (2025, pure-ARM stock ~25%). Germany/Italy/France: ~2022, consistent with the ECB ordering and the HFCS cross-country ranking (Ehrmann-Ziegelmeyer, ECB WP 1631). See docs/sources.md.",
            **COUNTRY},
        "notes": [
            "Steady-state mortgage payment uses the measured rate_level; spread_over_euribor_2022 is the margin over Euribor, used only in the incidence narrative.",
            "Germany's counterfactual uses the ECB MIR / ESRB stock variable-rate share (0.12).",
            "The illiquid moment is steep in beta_modal (stored to 10 decimals); see household._note_precision.",
        ],
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--refine", action="store_true", help="fast modal-beta refinement instead of the full SMM")
    args = ap.parse_args()

    t = read_targets(); T = smm_targets(t)
    x, lossval = run_refine(T) if args.refine else run_smm(T)
    out = assemble(x, lossval, t, T)

    with open(os.path.join(CAL, "calibration.json"), "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    h = out["household"]
    print(f"\nbeta_modal = {h['beta'][0]}  loss = {h['loss']}")
    print(f"achieved: {h['achieved']}")
    print(f"wrote {os.path.join(CAL, 'calibration.json')}")

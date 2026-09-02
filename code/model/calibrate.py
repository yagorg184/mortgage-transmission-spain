# calibrate.py
# Production calibrator (exactly-identified 4-vs-4 design) and assembler of the
# frozen parameter set in calibration.json. This is the single calibrator the
# replication driver uses; the calib_*.py scripts in this folder are the
# exploration record (mpc-free, more-targets, fewer-params, three-type, 5v5) that
# led here.
#
#   python calibrate.py            full SMM (differential evolution, seed 0)
#   python calibrate.py --refine   re-pin the modal discount factor from the
#                                  current calibration.json (deterministic, fast)
#
# Design: free {beta_modal, pi_modal, chi1, sig_eps}; external ra and beta_patient
# (see below); targets {liquid p50, illiquid p50, HtM, top-10% net-wealth share}
# read from calibration/eff_targets.csv. The MPC is not targeted; it and the
# poor/wealthy-HtM split and the top-1% net-wealth share are computed and
# reported for reference.

import os, json, argparse
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize_scalar
import _calib_common as CC

HERE = os.path.dirname(os.path.abspath(__file__))
CAL = os.path.normpath(os.path.join(HERE, "..", "..", "calibration"))
JSON = os.path.join(CAL, "calibration.json")

# ---- externally fixed (not calibrated) ---------------------------------------
# ra: real net return on the housing-heavy illiquid portfolio (BdE rental yield +
#     Jorda et al.); set externally rather than calibrated.
# beta_patient: the patient (accumulating) type, beta*(1+ra) ~ 1.016; set externally.
RA = 0.035
BETA_PAT = 0.981

# free params {beta_modal, pi_modal, chi1, sig_eps}; bounds narrowed around the
# known region (illiquid wealth is steep in beta_modal).
BOUNDS = [(0.950, 0.958), (0.92, 0.985), (5.0, 8.0), (0.28, 0.42)]
TARGETS = [("liquid_income_p50", 1.0), ("illiquid_income_p50", 0.8),
           ("htm_share", 1.0), ("net_wealth_top10_share", 1.0)]

# blocks that are not part of the household SMM
MORT_MODEL = dict(amort=0.03)
INCIDENCE = dict(beta_dep=0.32, mpc_shareholder=0.05, lambda_ext=0.50)
COUNTRY = {"all_fixed": 0.0, "Germany": 0.12, "France": 0.12, "euro_area": 0.25,
           "Italy": 0.45, "Spain": 0.70, "all_variable": 1.0}


def unpack(x):
    bm, pim, chi1, se = x
    if BETA_PAT < bm:
        return None
    return ([bm, BETA_PAT], [pim, 1 - pim], chi1, se, RA)


LOSS = CC.Objective(TARGETS, unpack)


def run_smm(maxiter=120, popsize=16):
    res = differential_evolution(LOSS, BOUNDS, maxiter=maxiter, popsize=popsize,
                                 tol=1e-6, seed=0, polish=True, disp=True,
                                 mutation=(0.5, 1.0), recombination=0.7)
    return res.x, res.fun


def run_refine():
    """Re-pin only the modal discount factor (the steep, wealth-pinning parameter)
    from the current calibration.json, holding the rest."""
    h = json.load(open(JSON))["household"]
    bm0, pim, chi1, se = h["beta"][0], h["pi"][0], h["chi1"], h["sig_eps"]
    r = minimize_scalar(lambda bm: LOSS([bm, pim, chi1, se]),
                        bounds=(bm0 - 0.003, bm0 + 0.003), method="bounded",
                        options={"xatol": 1e-8})
    return [r.x, pim, chi1, se], r.fun


def assemble(x, lossval):
    bm = round(x[0], 8)           # steep in illiquid -> kept to 8 decimals
    pim, chi1, se = round(x[1], 6), round(x[2], 6), round(x[3], 6)
    betas, pis = [bm, BETA_PAT], [pim, round(1 - pim, 6)]
    m = CC.solve_moments(betas, pis, chi1, se, RA)
    t = pd.read_csv(os.path.join(CAL, "eff_targets.csv")).drop_duplicates("key").set_index("key")["value"]
    r4 = lambda v: round(float(v), 4)
    return {
        "_about": "Single source of truth for the calibrated model. All model scripts read parameters from here; produced by calibrate.py (exactly-identified 4v4 design) against calibration/eff_targets.csv.",
        "household": {
            "beta": [bm, BETA_PAT],
            "pi": [pim, round(1 - pim, 6)],
            "chi1": chi1, "sig_eps": se, "ra": RA,
            "_external": "ra and beta_patient are set externally, not calibrated. ra = real net return on the housing-heavy illiquid portfolio; beta_patient = the patient (accumulating) type, beta*(1+ra)~1.016.",
            "fixed": {k: CC.H.BASE[k] for k in ("eis", "rb", "chi0", "chi2", "rho_z", "sig_eta")},
            "grid": {k: CC.H.BASE[k] for k in ("nB", "nA", "nK", "nZ", "nEps", "bmax", "amax", "kmax")},
            "targets": {"liquid_income_p50": r4(t["liquid_income_p50"]),
                        "illiquid_income_p50": r4(t["illiquid_income_p50"]),
                        "htm_share": r4(t["htm_share"]),
                        "net_wealth_top10_share": r4(t["net_wealth_top10_share"])},
            "achieved": {"liquid_income_p50": round(m["liquid_income_p50"], 4),
                         "illiquid_income_p50": round(m["illiquid_income_p50"], 4),
                         "htm_share": round(m["htm_share"], 4),
                         "net_wealth_top10_share": round(m["net_wealth_top10_share"], 4)},
            "untargeted": {"mpc_annual": round(m["mpc_annual"], 4),
                           "net_wealth_income_p50": round(m["net_wealth_income_p50"], 4),
                           "htm_wealthy_share": round(m["htm_wealthy_share"], 4),
                           "htm_poor_share": round(m["htm_poor_share"], 4),
                           "net_wealth_top1_share": round(m["net_wealth_top1_share"], 4),
                           "_note": "MPC is not targeted; it, the poor/wealthy-HtM split, the net-wealth median and the top-1% net-wealth share are computed and reported for reference (not targeted)."},
            "loss": round(lossval, 6),
            "_note_precision": "beta_modal is steep in the illiquid moment; stored to 8 decimals and re-pinned by --refine.",
        },
        "mortgage": {
            "mortgage_share": r4(t["mortgage_share"]),
            "variable_share": r4(t["variable_share_among_mortgaged"]),
            "M_var_income_exposed": r4(t["Mvar_income_p50_exposed"]),
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
            "Household calibration is the exactly-identified 4v4 design; see household._note_precision and exploration/calib_E_4v4.py.",
        ],
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--refine", action="store_true", help="fast modal-beta refinement instead of the full SMM")
    args = ap.parse_args()

    x, lossval = run_refine() if args.refine else run_smm()
    out = assemble(x, lossval)
    with open(JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    h = out["household"]
    print(f"\nbeta_modal = {h['beta'][0]}  loss = {h['loss']}")
    print(f"achieved:   {h['achieved']}")
    print(f"untargeted: MPC {h['untargeted']['mpc_annual']} (validation)  "
          f"poor-HtM {h['untargeted']['htm_poor_share']}  top1 {h['untargeted']['net_wealth_top1_share']}")
    print(f"wrote {JSON}")

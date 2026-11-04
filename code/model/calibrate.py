# calibrate.py
# Household calibration on the full population: households without a mortgage and
# mortgagors in three terciles of dwelling debt relative to income, each paying the
# steady-state annuity on its tercile's median balance, term and rate (EFF 2022).
#
# Parameters: the discount factors of the modal and the patient type, the share of
# patient households, the adjustment-cost scale, the dispersion of transitory
# income and, optionally, a shift of the modal discount factor in each mortgage
# group. Targets: the medians of liquid and illiquid wealth relative to income,
# the hand-to-mouth share and the top-10% share of net worth, plus, in the
# tercile calibration, the hand-to-mouth share of each tercile of mortgagors.
#
# The patient type must satisfy beta_p (1 + r_a) < 1 for its wealth distribution
# to be stationary. The bound is tightened to beta_p (1 + r_a) <= 0.999, so that
# the calibration does not rely on households arbitrarily close to an explosive
# accumulation path, and every steady state is checked to be a proper
# distribution that the asset grids do not truncate.
#
# In the baseline the five common parameters are matched to the four
# population moments, which leaves one direction unidentified (a less patient
# top type with a larger population share produces the same moments); the
# option --beta-patient fixes the patient discount factor at other points of
# that set. The option --terciles keeps the patient type of the baseline, frees
# the three group shifts and replaces the top-10% share by the three tercile
# targets: with the shifts, the top share peaks below its data value as a
# function of the patient share, so the variant reports it as untargeted.
#
# Usage:
#   python calibrate.py                        baseline, writes calibration/calibration.json
#   python calibrate.py --check [FILE]         moments at a stored calibration
#   python calibrate.py --beta-patient 0.962   alternative point, calibration/alternatives/bp0.962.json
#   python calibrate.py --terciles             tercile calibration, calibration/alternatives/terciles.json

import argparse
import json
import os
import time

import numpy as np

import config as CF
import economy as EC

POPULATION_TARGETS = ["liquid_income_p50", "illiquid_income_p50", "htm_share", "net_wealth_top10_share"]
GROUPS = ["tercile1", "tercile2", "tercile3"]
TERCILE_TARGETS = [f"htm_share_{g}" for g in GROUPS]
NAMES = ["beta_modal", "beta_patient", "patient_share", "chi1", "sig_eps"] + [f"shift_{g}" for g in GROUPS]

# Conventions and external values (sources in docs/sources.md).
HOUSEHOLD_FIXED = dict(eis=0.5, rb=0.0, ra=0.035, chi0=0.25, chi2=2.0, tax=0.0, w=1.0, N=1.0,
                       kmax=1, nB=30, nK=4, nZ=5, nEps=3, rho_z=0.98, sig_eta=0.08)
GRIDS = {"modal": dict(bmax=50, amax=200, nA=50), "patient": dict(bmax=50, amax=1500, nA=80)}
STATIONARITY_MARGIN = 0.999

START = [0.951118, 0.964506, 0.054537, 7.41633, 0.34527, 0.0, 0.0, 0.0]
TERCILE_START = [0.95108, 0.964536, 0.051191, 8.533147, 0.337328, -0.004408, -0.004775, 0.000514]
TERCILE_TOL = 5e-3      # the tercile shares have standard errors of 13-16 percent of their value
STEPS = np.array([2e-4, 2e-4, 0.005, 0.1, 0.005, 2e-4, 2e-4, 2e-4])


def data_targets():
    t = {k: float(CF.EFF[k]) for k in POPULATION_TARGETS}
    for k, g in enumerate(GROUPS, start=1):
        t[f"htm_share_{g}"] = float(CF.EXPOSURE[f"t{k}_htm_share"])
    return t


def population(z):
    beta_m, beta_p, pi_p, chi1, sig_eps = z[:5]
    shift = {g: z[5 + k] for k, g in enumerate(GROUPS)} if np.any(np.asarray(z[5:]) != 0) else None
    return CF.household_population([beta_m, beta_p, pi_p, sig_eps], HOUSEHOLD_FIXED, GRIDS, chi1=chi1,
                                   modal_beta_shift=shift)


def model_moments(z, check=True):
    hh, types, groups = population(z)
    cells = EC.solve_cells(hh, types, groups, check=check)
    m = EC.moments(cells, types, groups)
    for g, row in EC.moments_by_group(cells, types, groups).items():
        m[f"htm_share_{g}"] = row["htm_share"]
    return m, cells


def residuals(z, targets):
    data = data_targets()
    try:
        m, _ = model_moments(z)
    except RuntimeError:
        return np.full(len(targets), 10.0)
    return np.array([(m[k] - data[k]) / data[k] for k in targets])


def solve(z0, targets, free, tol=2e-3, max_iter=30):
    """Gauss-Newton iterations with forward-difference derivatives over the
    parameters flagged in `free`. When there are more free parameters than
    targets, each step is the minimum-norm solution of the linearised system.
    Steps are halved until the largest relative miss falls. Iterations stop once
    every target is matched within `tol` (0.2 percent, inside the EFF sampling
    error of each moment), or when no step along the Gauss-Newton direction
    reduces the largest miss."""
    free = np.asarray(free, bool)
    ra = HOUSEHOLD_FIXED["ra"]
    bmax = STATIONARITY_MARGIN / (1 + ra)
    lo = np.array([0.90, 0.90, 0.01, 1.0, 0.10, -0.02, -0.02, -0.02])
    hi = np.array([bmax, bmax, 0.30, 20.0, 0.80, 0.02, 0.02, 0.02])
    z = np.array(z0, float)
    r = residuals(z, targets)
    t0 = time.time()
    for it in range(max_iter):
        err = np.max(np.abs(r))
        print(f"  [{(time.time() - t0) / 60:5.1f} min] iter {it:2d}  z = {np.round(z, 6)}  res = {np.round(r, 5)}",
              flush=True)
        if err < tol:
            break
        J = np.zeros((len(r), len(z)))
        for k in np.flatnonzero(free):
            step = STEPS[k] if z[k] + STEPS[k] <= hi[k] else -STEPS[k]
            zk = z.copy()
            zk[k] += step
            J[:, k] = (residuals(zk, targets) - r) / step
        scale = STEPS * 10 * free
        dz = scale * np.linalg.lstsq(J * scale, -r, rcond=None)[0]
        lam = 1.0
        while lam > 1 / 64:
            zn = np.clip(z + lam * dz, lo, hi)
            rn = residuals(zn, targets)
            if np.max(np.abs(rn)) < err:
                break
            lam /= 2
        else:
            print("  no step reduces the largest miss; stopping", flush=True)
            break
        z, r = zn, rn
    return z, r


def import_share(year=2021):
    """Import content of household consumption from the input-output tables
    (code/empirical/io/import_content.py)."""
    import pandas as pd
    path = os.path.normpath(os.path.join(CF.HERE, "..", "empirical", "io", "output", "import_content.csv"))
    t = pd.read_csv(path).set_index("year")
    return round(float(t.loc[year, "total"]), 4)


def assemble(z, m, cells, targets):
    hh, types, groups = population(z)
    diag = {f"{types[j]['name']}/{groups[g]['name']}": c["diag"] for (j, g), c in cells.items()}
    data = data_targets()
    r4 = lambda v: round(float(v), 4)
    household = {
        "fixed": HOUSEHOLD_FIXED,
        "grids": GRIDS,
        "types": [dict(name=t["name"], beta=round(float(t["beta"]), 8), mass=round(float(t["mass"]), 8))
                  for t in types],
        "chi1": round(float(z[3]), 6),
        "sig_eps": round(float(z[4]), 6),
    }
    if np.any(np.asarray(z[5:]) != 0):
        household["modal_beta_shift"] = {g: round(float(z[5 + k]), 8) for k, g in enumerate(GROUPS)}
    household.update({
        "stationarity_margin": STATIONARITY_MARGIN,
        "targets": {k: r4(data[k]) for k in targets},
        "achieved": {k: r4(m[k]) for k in targets},
        "untargeted": {k: r4(v) for k, v in m.items() if k not in targets},
        "distribution_checks": diag,
    })
    return {
        "_about": "Single source of truth for the model. Household block produced by calibrate.py against "
                  "calibration/eff_targets.csv and eff_exposure.csv; sources for external values in docs/sources.md.",
        "household": household,
        "nominal": {"inflation_ss": 0.02},
        "banks": {"foreign_share": 0.50, "corporate_tax": 0.30,
                  "policy_spread": round(float(CF.EFF["spread_over_euribor1Y_2022"]) / 100, 4)},
        "union": {"import_share": import_share()},
        "ge": {"frisch": 1.0, "kappa_w": 0.1, "beta_w": 0.96},
        "country_variable_share_stock": {
            "_source": "Spain: EFF 2022 (share of the outstanding dwelling-debt stock at variable rates). Euro "
                       "area: ECB (2025). Germany, France, Italy: HFCS ranking (Ehrmann and Ziegelmeyer, 2017) "
                       "and ECB MIR; read as ordinal.",
            "all_fixed": 0.0, "Germany": 0.12, "France": 0.12, "euro_area": 0.25,
            "Italy": 0.45, "Spain": 0.68, "all_variable": 1.0},
    }


def stored_vector(path):
    cal, hh, types, groups = CF.load(path)
    h = cal["household"]
    shift = h.get("modal_beta_shift", {})
    return [types[0]["beta"], types[1]["beta"], types[1]["mass"], h["chi1"], h["sig_eps"]] + \
        [shift.get(g, 0.0) for g in GROUPS]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", nargs="?", const=CF.JSON, default=None)
    ap.add_argument("--beta-patient", type=float, default=None)
    ap.add_argument("--terciles", action="store_true")
    args = ap.parse_args()
    out_dir = CF.ALTERNATIVES
    if args.check:
        z = stored_vector(args.check)
        m, _ = model_moments(z)
        data = data_targets()
        for k in POPULATION_TARGETS + TERCILE_TARGETS:
            print(f"  {k:24s} model {m[k]:.4f}   EFF {data[k]:.4f}")
        return
    tol = 2e-3
    if args.terciles:
        z0 = list(TERCILE_START)
        targets = [k for k in POPULATION_TARGETS if k != "net_wealth_top10_share"] + TERCILE_TARGETS
        free = [True, False, False, True, True, True, True, True]
        path = os.path.join(out_dir, "terciles.json")
        tol = TERCILE_TOL
    elif args.beta_patient is not None:
        z0 = list(START)
        z0[1] = args.beta_patient
        targets = POPULATION_TARGETS
        free = [True, False, True, True, True, False, False, False]
        path = os.path.join(out_dir, f"bp{args.beta_patient}.json")
    else:
        z0 = list(START)
        targets = POPULATION_TARGETS
        free = [True, True, True, True, True, False, False, False]
        path = CF.JSON
    os.makedirs(out_dir, exist_ok=True)
    z, rres = solve(z0, targets, free, tol=tol)
    m, cells = model_moments(z)
    cal = assemble(z, m, cells, targets)
    with open(path, "w") as f:
        json.dump(cal, f, indent=2)
    print(f"\nresiduals = {np.round(rres, 5)}; wrote {path}")


if __name__ == "__main__":
    main()

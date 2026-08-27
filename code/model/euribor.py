# euribor.py
# Observed Euribor path (ECB Data Portal) and the partial-equilibrium consumption
# response of variable-rate borrowers to it. The 1-year Euribor is the reference
# rate for Spanish variable mortgages and matches the annual model frequency.
#
# Data: CSV exports of the 1-year Euribor and Spanish HICP under data/euribor/
# (override with the EURIBOR_DIR environment variable).

import os
import csv
import glob
import numpy as np
import mortgage as M
import config as C

DDIR = os.environ.get(
    "EURIBOR_DIR",
    os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "data", "euribor")),
)


def _load(substr):
    f = [p for p in glob.glob(os.path.join(DDIR, "*.csv")) if substr in open(p).readline()][0]
    ym = {}
    with open(f) as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            if len(row) < 3 or row[2] in ("", "-"):
                continue
            try:
                v = float(row[2])
            except ValueError:
                continue
            y = int(row[0][:4]); ym.setdefault(y, []).append(v)
    return {y: float(np.mean(v)) for y, v in ym.items()}


def annual_table(y0=2019, y1=2026):
    eur = _load("Euribor 1-year")
    hicp = _load("HICP.M.ES")
    years = [y for y in range(y0, y1 + 1) if y in eur]
    return years, {y: eur[y] for y in years}, {y: hicp.get(y, np.nan) for y in years}


def shock_path(base_years=(2021,), real=False):
    """Deviation of the (nominal or real) Euribor from its pre-hike 2021 level."""
    years, eur, hicp = annual_table()
    rate = {y: (eur[y] - hicp[y]) if real else eur[y] for y in years}
    base = float(np.mean([rate[y] for y in base_years]))
    return years, {y: rate[y] - base for y in years}, base


def var_response(y0=2022, y1=2026, T=40, tail=0.7, real=False):
    """Consumption of the variable-rate group along the observed Euribor path."""
    _, dev, _ = shock_path(real=real)
    obs = [dev[y] / 100 for y in range(y0, y1 + 1)]
    path = np.zeros(T); path[:len(obs)] = obs
    for t in range(len(obs), T):
        path[t] = path[t - 1] * tail
    num = np.zeros(T); Cden = 0.0
    for beta, pib in zip(C.BETAS, C.PIS):
        ss, _ = M.solve_type(beta, (C.RATE + C.AMORT) * C.MV)
        dC = np.asarray(M.HHm.impulse_linear(ss, {"payment": path * C.MV}, outputs=["C"])["C"])
        mass = pib * C.MORT_SHARE * C.PVAR
        num += mass * dC; Cden += mass * float(ss["C"])
    return list(range(y0, y1 + 1)), num / Cden


if __name__ == "__main__":
    years, eur, hicp = annual_table()
    dev = shock_path()[1]
    print("Euribor 1Y path and deviation from 2021:")
    for y in years:
        print(f"  {y}: {eur[y]:+.2f}%   shock {100 * dev[y]:+.0f} bp")
    yrs, r = var_response()
    print(f"\n  variable-group consumption 2024 = {100 * r[2]:+.2f}%  "
          f"peak {100 * r[np.argmax(np.abs(r[:len(yrs)]))]:+.2f}%  (BdE anchor -2.7%)")

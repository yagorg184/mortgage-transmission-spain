#!/usr/bin/env python3
# lp.py
# Local projections with an external instrument (LP-IV) for the transmission of
# ECB monetary policy, quarterly, 1999-2019. The policy indicator is the 1-year
# Euribor and the instrument is the Jarocinski-Karadi (2020) monetary policy
# shock (median rotation, purged of central-bank information effects), summed
# within the quarter. For each outcome y and horizon h,
#
#     y_{t+h} - y_{t-1} = a_h + b_h (i_t - i_{t-1}) + controls_{t-1..t-4} + e_{t+h},
#
# with i_t - i_{t-1} instrumented by the shock. Controls are four lags of the
# change in the outcome, the Euribor, euro-area log HICP and log GDP, and of the
# shock, plus quarter dummies. Standard errors are Newey-West with h+1 lags.
# The same regression with the average of y over the four quarters of year k on
# the left gives annual-average responses, which is the object the annual model
# is compared with.
#
# Outputs (output/): lp_quarterly.csv, lp_annual.csv, first_stage.txt

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "..", "..", "data", "public"))
OUT = os.path.join(HERE, "output")

START, END = "1999Q1", "2019Q4"
H = 12
YEARS = 3
NLAGS = 4
SHOCK = "MP_median"


def monthly(name):
    df = pd.read_csv(os.path.join(DATA, f"{name}.csv"))
    s = pd.Series(df["value"].values, index=pd.PeriodIndex(df["date"], freq="M"))
    return s.groupby(s.index.asfreq("Q")).mean()


def quarterly(name):
    df = pd.read_csv(os.path.join(DATA, f"{name}.csv"))
    return pd.Series(df["value"].values, index=pd.PeriodIndex(df["date"].str.replace("-", ""), freq="Q"))


def build_panel():
    jk = pd.read_csv(os.path.join(DATA, "jk_shocks_monthly.csv"))
    idx = pd.PeriodIndex(pd.to_datetime(dict(year=jk.year, month=jk.month, day=1)), freq="M")
    shocks = pd.DataFrame({c: jk[c].values for c in ("MP_median", "MP_pm")}, index=idx)
    shocks = shocks.groupby(shocks.index.asfreq("Q")).sum()

    vol_on, vol_t = monthly("dep_vol_overnight_es"), monthly("dep_vol_time_es")
    r_on, r_t = monthly("dep_rate_overnight_es"), monthly("dep_rate_time_stock_es")
    dep = (r_on * vol_on + r_t * vol_t) / (vol_on + vol_t)

    P = pd.DataFrame({
        "euribor": monthly("euribor_1y"),
        "deposit_rate_es": dep,
        "hicp_es": 100 * np.log(monthly("hicp_es")),
        "hicp_ea": 100 * np.log(monthly("hicp_ea")),
        "cons_es": 100 * np.log(quarterly("cons_es")),
        "cons_de": 100 * np.log(quarterly("cons_de")),
        "cons_ea": 100 * np.log(quarterly("cons_ea")),
        "gdp_es": 100 * np.log(quarterly("gdp_es")),
        "gdp_de": 100 * np.log(quarterly("gdp_de")),
        "gdp_ea": 100 * np.log(quarterly("gdp_ea")),
    })
    P = P.join(shocks, how="left")
    P[["MP_median", "MP_pm"]] = P[["MP_median", "MP_pm"]].fillna(0.0)
    return P.loc[START:END]


def nw_cov(Xh, u, lags):
    """Newey-West covariance of the 2SLS coefficients (Xh = projected regressors)."""
    n, k = Xh.shape
    g = Xh * u[:, None]
    S = g.T @ g
    for l in range(1, lags + 1):
        w = 1 - l / (lags + 1)
        G = g[l:].T @ g[:-l]
        S += w * (G + G.T)
    B = np.linalg.inv(Xh.T @ Xh)
    return B @ S @ B


def lp_iv(P, outcome, horizon_fn, lags_nw, shock=SHOCK):
    """2SLS of the horizon-specific left-hand side on the change in the Euribor,
    instrumented by the shock. Returns coefficient, s.e., n and first-stage F."""
    d = pd.DataFrame(index=P.index)
    d["lhs"] = horizon_fn(P[outcome])
    d["di"] = P["euribor"].diff()
    d["z"] = P[shock]
    ctrl = []
    for var in dict.fromkeys([outcome, "euribor", "hicp_ea", "gdp_ea"]):
        dv = P[var].diff()
        for l in range(1, NLAGS + 1):
            d[f"{var}_l{l}"] = dv.shift(l)
            ctrl.append(f"{var}_l{l}")
    for l in range(1, NLAGS + 1):
        d[f"z_l{l}"] = P[shock].shift(l)
        ctrl.append(f"z_l{l}")
    q = pd.get_dummies(pd.Series(P.index.quarter, index=P.index), prefix="q", drop_first=True).astype(float)
    d = d.join(q)
    ctrl += list(q.columns)
    d = d.dropna()
    n = len(d)
    W = np.column_stack([np.ones(n), d[ctrl].values])
    Z = np.column_stack([d["z"].values, W])
    X = np.column_stack([d["di"].values, W])
    y = d["lhs"].values
    Pz = Z @ np.linalg.pinv(Z.T @ Z) @ Z.T
    Xh = Pz @ X
    beta = np.linalg.lstsq(Xh, y, rcond=None)[0]
    u = y - X @ beta
    V = nw_cov(Xh, u, lags_nw)
    # first stage: di on z and controls, HAC F on z
    g = np.linalg.lstsq(Z, d["di"].values, rcond=None)[0]
    v = d["di"].values - Z @ g
    Vg = nw_cov(Z, v, lags_nw)
    F = g[0] ** 2 / Vg[0, 0]
    return beta[0], np.sqrt(V[0, 0]), n, F


def at_horizon(h):
    return lambda s: s.shift(-h) - s.shift(1)


def annual_average(k):
    return lambda s: sum(s.shift(-(4 * k + j)) for j in range(4)) / 4 - s.shift(1)


OUTCOMES = ["euribor", "deposit_rate_es", "hicp_es", "cons_es", "cons_de", "cons_ea",
            "gdp_es", "gdp_de", "gdp_ea"]


def main(shock=SHOCK):
    os.makedirs(OUT, exist_ok=True)
    P = build_panel()
    rows, arows, fs = [], [], []
    for y in OUTCOMES:
        for h in range(H + 1):
            b, se, n, F = lp_iv(P, y, at_horizon(h), h + 1, shock)
            rows.append(dict(outcome=y, h=h, coef=b, se=se, n=n, F=F))
        for k in range(YEARS):
            b, se, n, F = lp_iv(P, y, annual_average(k), 4 * k + 4, shock)
            arows.append(dict(outcome=y, year=k, coef=b, se=se, n=n, F=F))
    q = pd.DataFrame(rows)
    a = pd.DataFrame(arows)
    sfx = "" if shock == SHOCK else f"_{shock}"
    q.to_csv(os.path.join(OUT, f"lp_quarterly{sfx}.csv"), index=False)
    a.to_csv(os.path.join(OUT, f"lp_annual{sfx}.csv"), index=False)
    f0 = q[(q.outcome == "cons_es") & (q.h == 0)].iloc[0]
    with open(os.path.join(OUT, f"first_stage{sfx}.txt"), "w") as f:
        f.write(f"shock {shock}; sample {START}-{END}; impact first-stage HAC F = {f0.F:.1f} (n = {f0.n})\n")
    print(a.pivot(index="outcome", columns="year", values="coef").round(3))
    print(f"first-stage F (impact, consumption sample): {f0.F:.1f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else SHOCK)

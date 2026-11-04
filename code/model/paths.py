# paths.py
# Exogenous paths fed to the model, annual, as deviations from the steady state:
#   di   : 1-year Euribor (index of Spanish variable-rate mortgages)
#   ddep : rate paid on Spanish household deposits
#   dpi  : inflation
#
# Three experiments use them:
#   stylized_shock  a persistent rise in the Euribor (for intuition)
#   svar_shock      the paths identified by a VAR specification for an ECB shock
#   episode         the observed paths of 2022-2025 relative to 2021
#
# Households learn a path either all at once in its first year (perfect foresight)
# or year by year, expecting each period that the current deviation decays
# geometrically (a sequence of surprises). Under the second assumption a path is
# represented by the revisions of expectations it implies.

import json
import os

import numpy as np
import pandas as pd

import config as CF

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
PUBLIC = os.path.join(ROOT, "data", "public")
SVAR_OUT = os.path.join(ROOT, "code", "empirical", "macro_svar", "output")

INFLATION_TARGET = json.load(open(CF.JSON))["nominal"]["inflation_ss"]
# Deposit pass-through of the stylized and VAR experiments: the cumulative change
# in the rate on the stock of household deposits between 2021 and 2024 relative to
# the change in the 1-year Euribor (0.121; results.json, facts).
DEPOSIT_BETA = 0.12


def _monthly(name):
    df = pd.read_csv(os.path.join(PUBLIC, f"{name}.csv"))
    return pd.Series(df["value"].values, index=pd.PeriodIndex(df["date"], freq="M"))


def household_deposit_rate():
    """Volume-weighted rate on the stock of Spanish household deposits (overnight
    and with agreed maturity), monthly, percent per year."""
    r_on, r_t = _monthly("dep_rate_overnight_es"), _monthly("dep_rate_time_stock_es")
    v_on, v_t = _monthly("dep_vol_overnight_es"), _monthly("dep_vol_time_es")
    return (r_on * v_on + r_t * v_t) / (v_on + v_t)


def annual_observed(first=2019, last=2025):
    """Annual averages of the Euribor, the household deposit rate and HICP
    inflation for Spain (percent)."""
    eur = _monthly("euribor_1y")
    dep = household_deposit_rate()
    hicp = _monthly("hicp_es")
    out = pd.DataFrame({
        "euribor": eur.groupby(eur.index.year).mean(),
        "deposit_rate": dep.groupby(dep.index.year).mean(),
        "inflation": 100 * np.log(hicp.groupby(hicp.index.year).mean()).diff(),
    })
    return out.loc[first:last]


def episode(T, base_year=2021, years=range(2022, 2026), tail=0.7):
    """Observed 2022-2025 deviations of the Euribor and the deposit rate from their
    base-year level, and of inflation from the 2% target; after the last observed
    year the rate deviations decay at rate `tail` and inflation is at target."""
    obs = annual_observed()
    years = list(years)
    di, ddep, dpi = np.zeros(T), np.zeros(T), np.zeros(T)
    for k, y in enumerate(years):
        di[k] = (obs.loc[y, "euribor"] - obs.loc[base_year, "euribor"]) / 100
        ddep[k] = (obs.loc[y, "deposit_rate"] - obs.loc[base_year, "deposit_rate"]) / 100
        dpi[k] = obs.loc[y, "inflation"] / 100 - INFLATION_TARGET
    for t in range(len(years), T):
        di[t], ddep[t] = tail * di[t - 1], tail * ddep[t - 1]
    return dict(di=di, ddep=ddep, dpi=dpi, years=years)


def stylized_shock(T, size=0.01, rho=0.6, deposit_beta=DEPOSIT_BETA):
    di = size * rho ** np.arange(T)
    return dict(di=di, ddep=deposit_beta * di, dpi=np.zeros(T))


def svar_shock(T, spec="baseline", deposit_beta=DEPOSIT_BETA, years=3):
    """Annual averages of the responses to an ECB shock in one VAR specification
    (code/empirical/macro_svar): the index of variable-rate mortgages (the 1-year
    Euribor if the specification includes it, the policy rate otherwise), the
    household deposit rate (its own response if included, deposit_beta times the
    index otherwise) and Spanish inflation (the change in the annual-average log
    price level). After `years` the rates are back at baseline and the price level
    stays at its last deviation."""
    a = pd.read_csv(os.path.join(SVAR_OUT, f"svar_annual_{spec}.csv")).set_index("year")
    index = "euribor" if "euribor" in a else "rate"
    di, ddep, dpi = np.zeros(T), np.zeros(T), np.zeros(T)
    p_prev = 0.0
    for k in range(years):
        di[k] = a.loc[k, index] / 100
        ddep[k] = a.loc[k, "deposit"] / 100 if "deposit" in a else deposit_beta * di[k]
        dpi[k] = a.loc[k, "lp_es"] / 100 - p_prev
        p_prev = a.loc[k, "lp_es"] / 100
    return dict(di=di, ddep=ddep, dpi=dpi)


def expectation_revisions(x, rho):
    """Represent a realised path x as the sum of surprises: at each date s households
    revise their forecast to x_s rho^(t-s) for t >= s. Returns a list of (s, path)
    with the revision announced at s, aligned to start at s."""
    T = len(x)
    revisions = []
    forecast = np.zeros(T)
    for s in range(T):
        new = forecast.copy()
        new[s:] = x[s] * rho ** np.arange(T - s)
        rev = new - forecast
        if np.any(np.abs(rev[s:]) > 0):
            revisions.append((s, rev[s:]))
        forecast = new
    return revisions

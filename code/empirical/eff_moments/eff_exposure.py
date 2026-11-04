#!/usr/bin/env python3
# eff_exposure.py
# Mortgage exposure in the EFF 2022 for the heterogeneous-borrower version of the
# model: the distribution of dwelling debt relative to income among mortgagors,
# summarised in terciles, with the variable-rate share, the hand-to-mouth share,
# the remaining term, the current rate and the debt-service ratio in each tercile.
#
# The construction of every variable follows build_eff.R line by line (dwelling
# loans 1-4, variable if the rate-reset flag is 11, liquid wealth = deposits and
# cash net of card debt, income = 2021 household income). Point estimates average
# the five implicates; standard errors use the 1,000 bootstrap replicate weights
# within each implicate and Rubin's rules across implicates, as in eff_moments.R.
# The script first reproduces the moments of eff_targets.csv that it shares with
# the R pipeline, as a check of the port.
#
# Input : data/raw/{otras_secciones_2022_imp*.csv, databol*.csv, replicate_weights_2022.csv}
# Output: output/eff_exposure.csv
#
# Run from code/empirical/eff_moments/:  python eff_exposure.py [--no-se]

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.environ.get("EFF_RAW_DIR", os.path.join(HERE, "data", "raw"))
OUT = os.path.join(HERE, "output")
LOANS = range(1, 5)
HTM_DIV = 24
REP_R = 1000

VARS_OTRAS = (["h_2022", "facine3"] + [f"p2_{s}_{k}" for s in (9, 12, 14, 18) for k in LOANS]
              + ["p2_13_1", "p2_17_1", "p4_7_1", "p4_7_2", "p4_7_3", "p8_5a"])
VARS_DATABOL = ["h_2022", "renthog21_eur22", "riquezanet"]
NEG = [c for c in VARS_OTRAS if c not in ("h_2022", "facine3")]


def build(m):
    os_ = pd.read_csv(os.path.join(RAW, f"otras_secciones_2022_imp{m}.csv"), sep=";", usecols=VARS_OTRAS)
    db = pd.read_csv(os.path.join(RAW, f"databol{m}.csv"), sep=";", usecols=VARS_DATABOL)
    d = os_.merge(db, on="h_2022", how="left")
    for c in NEG:
        d.loc[d[c] < 0, c] = np.nan
    B = np.column_stack([d[f"p2_12_{k}"].fillna(0).values for k in LOANS])
    PY = np.column_stack([d[f"p2_18_{k}"].fillna(0).values for k in LOANS])
    RV = np.column_stack([d[f"p2_14_{k}"].values for k in LOANS])
    is_var = np.nan_to_num(RV == 11).astype(bool)
    is_fix = np.nan_to_num(RV == 22).astype(bool)
    d["M_var"] = (B * is_var).sum(1)
    d["M_fix"] = (B * is_fix).sum(1)
    d["M_tot"] = B.sum(1)
    d["mortgage_pay"] = PY.sum(1)
    d["liquid"] = d[["p4_7_1", "p4_7_2", "p4_7_3"]].fillna(0).sum(1) - d["p8_5a"].fillna(0)
    d["income"] = d["renthog21_eur22"]
    d["income_ok"] = d["income"] > 0
    inc = d["income"].where(d["income_ok"])
    d["liq_income"] = d["liquid"] / inc
    d["exp_income"] = d["M_var"] / inc
    d["tot_income"] = d["M_tot"] / inc
    d["htm"] = d["income_ok"] & (d["liquid"] <= d["income"] / HTM_DIV)
    d["dsr"] = d["mortgage_pay"] * 12 / inc
    d["var_mixed"] = d["M_var"] > 0
    return d


def w_quantile(x, w, p):
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x, w = x[ok], w[ok]
    o = np.argsort(x, kind="mergesort")
    return float(np.interp(p, np.cumsum(w[o]) / w.sum(), x[o]))


def w_share(cond, w):
    ok = ~np.isnan(cond) & np.isfinite(w) & (w > 0)
    return float((w[ok] * cond[ok]).sum() / w[ok].sum())


def rubin(imps, repw, est):
    """Point estimate (mean over implicates) and Rubin standard error."""
    ests, vars_ = [], []
    for d in imps:
        theta = est(d, d["facine3"].values)
        if repw is not None:
            W = repw.loc[d["h_2022"].values].values
            reps = np.array([est(d, W[:, r]) for r in range(W.shape[1])])
            vars_.append(np.mean((reps - theta) ** 2))
        ests.append(theta)
    M = len(ests)
    q = float(np.mean(ests))
    if repw is None:
        return q, np.nan
    B = float(np.var(ests, ddof=1)) if M > 1 else 0.0
    return q, float(np.sqrt(np.mean(vars_) + (1 + 1 / M) * B))


def tercile_cuts(imps):
    """Tercile cut points of dwelling debt / income among mortgagors, averaged
    over implicates (held fixed when computing the within-tercile moments)."""
    lo, hi = [], []
    for d in imps:
        s = (d["M_tot"] > 0) & d["income_ok"]
        x, w = d.loc[s, "tot_income"].values, d.loc[s, "facine3"].values
        lo.append(w_quantile(x, w, 1 / 3))
        hi.append(w_quantile(x, w, 2 / 3))
    return float(np.mean(lo)), float(np.mean(hi))


def main(with_se=True):
    os.makedirs(OUT, exist_ok=True)
    imps = [build(m) for m in range(1, 6)]
    repw = None
    if with_se:
        cols = ["h_2022"] + [f"wt3r_{r}" for r in range(1, REP_R + 1)]
        repw = pd.read_csv(os.path.join(RAW, "replicate_weights_2022.csv"), sep=";", usecols=cols).set_index("h_2022")
    rows = []

    def emit(key, est, note=""):
        v, se = rubin(imps, repw, est)
        rows.append(dict(key=key, value=v, se=se, note=note))
        print(f"  {key:40s} {v:9.4f}  ({se:.4f})")

    print("Check of the port against eff_targets.csv:")
    emit("check_mortgage_share", lambda d, w: w_share((d["M_tot"] > 0).values.astype(float), w))
    emit("check_variable_share_among_mortgaged",
         lambda d, w: w_share(d.loc[d["M_tot"] > 0, "var_mixed"].values.astype(float), w[(d["M_tot"] > 0).values]))
    emit("check_Mvar_income_p50_exposed",
         lambda d, w: w_quantile(d.loc[d["M_var"] > 0, "exp_income"].values, w[(d["M_var"] > 0).values], 0.5))
    emit("check_htm_share",
         lambda d, w: w_share(d.loc[d["income_ok"], "htm"].values.astype(float), w[d["income_ok"].values]))

    print("Mortgage stock and exposure:")
    emit("variable_share_of_stock",
         lambda d, w: float((w * d["M_var"].values).sum() / (w * d["M_tot"].values).sum()),
         "variable-rate balance / total dwelling balance")
    s_mort = lambda d: ((d["M_tot"] > 0) & d["income_ok"]).values
    emit("tot_income_p50_mortgagors",
         lambda d, w: w_quantile(d["tot_income"].values[s_mort(d)], w[s_mort(d)], 0.5))
    lo, hi = tercile_cuts(imps)
    rows.append(dict(key="tercile_cut_low", value=lo, se=np.nan, note="dwelling debt / income"))
    rows.append(dict(key="tercile_cut_high", value=hi, se=np.nan, note="dwelling debt / income"))
    print(f"  tercile cuts of dwelling debt / income: {lo:.4f}, {hi:.4f}")
    bounds = [(-np.inf, lo), (lo, hi), (hi, np.inf)]
    for k, (a, b) in enumerate(bounds, start=1):
        sel = lambda d, a=a, b=b: (s_mort(d) & (d["tot_income"].values > a) & (d["tot_income"].values <= b))
        emit(f"t{k}_mass_among_mortgagors", lambda d, w, sel=sel: w_share(sel(d).astype(float)[s_mort(d)], w[s_mort(d)]))
        emit(f"t{k}_tot_income_p50", lambda d, w, sel=sel: w_quantile(d["tot_income"].values[sel(d)], w[sel(d)], 0.5))
        emit(f"t{k}_variable_share", lambda d, w, sel=sel: w_share(d["var_mixed"].values[sel(d)].astype(float), w[sel(d)]))
        emit(f"t{k}_htm_share", lambda d, w, sel=sel: w_share(d["htm"].values[sel(d)].astype(float), w[sel(d)]))
        emit(f"t{k}_dsr_p50", lambda d, w, sel=sel: w_quantile(d["dsr"].values[sel(d)], w[sel(d)], 0.5))
        emit(f"t{k}_term_p50", lambda d, w, sel=sel: w_quantile(d["p2_17_1"].values[sel(d)], w[sel(d)], 0.5))
        emit(f"t{k}_rate_p50", lambda d, w, sel=sel: w_quantile(d["p2_13_1"].values[sel(d)], w[sel(d)], 0.5))
    emit("htm_share_mortgagors", lambda d, w: w_share(d["htm"].values[s_mort(d)].astype(float), w[s_mort(d)]))
    emit("dsr_p50_mortgagors", lambda d, w: w_quantile(d["dsr"].values[s_mort(d)], w[s_mort(d)], 0.5))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "eff_exposure.csv"), index=False)
    print(f"wrote {os.path.join(OUT, 'eff_exposure.csv')}")


if __name__ == "__main__":
    main(with_se="--no-se" not in sys.argv)

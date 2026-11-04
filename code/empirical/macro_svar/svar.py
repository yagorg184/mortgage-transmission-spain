#!/usr/bin/env python3
# svar.py
# Structural VARs for the response of Spanish household consumption to an ECB
# monetary policy shock, quarterly, 1999-2019. Spain is treated as a small member
# of the monetary union: euro-area variables form a block that Spanish variables
# do not enter, and the oil price depends only on its own lags.
#
# Specifications (all with two lags, a constant and a linear trend):
#
#   baseline   oil price (EUR), euro-area GDP, euro-area HICP, policy rate,
#              EUR/USD | Spanish GDP, Spanish HICP, Spanish consumption.
#              The shock is the innovation in the policy rate orthogonal to oil,
#              euro-area GDP and prices (recursive ordering), normalised to
#              +100 bp on impact.
#   euribor    adds the 1-year Euribor after the policy rate. The shock is
#              identified as in the baseline; the Euribor is the index of
#              Spanish variable-rate mortgages.
#   indicator  the 1-year Euribor replaces the policy rate as the policy
#              indicator (as the 1-year rate in Gertler and Karadi, 2015): it
#              reflects expected policy, including forward guidance, without a
#              shadow-rate model. The shock is normalised to +100 bp of the
#              Euribor on impact.
#   deposits   `indicator` with the rate on Spanish household deposits added
#              after Spanish prices. Sample 2003-2019, as the deposit series
#              start in 2003.
#   baseline_2003  the baseline on the sample of `deposits`.
#   internal   the Jarocinski-Karadi (2020) monetary policy shock, summed within
#              the quarter, ordered first and depending only on its own lags
#              (an internal instrument, Plagborg-Moller and Wolf, 2021), with
#              the variables of `euribor`. The response is normalised to
#              +100 bp of the 1-year Euribor on impact.
#
# The policy rate is a shadow short rate (Wu-Xia) spliced onto the 1-month
# Euribor; HICP indices are seasonally adjusted with STL. Bands come from a
# residual bootstrap with 1,000 draws.
#
# Outputs (output/): svar_irf_<spec>.csv (point responses and 68/90% bands by
# quarter), svar_panel.csv (baseline data).
# Run:  python svar.py [spec ...]       (all specifications by default)

import glob
import os
import sys
import warnings

import numpy as np
import openpyxl
import pandas as pd
from statsmodels.tsa.seasonal import STL

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
DATA = os.environ.get("SVAR_DATA", os.path.join(ROOT, "data", "macro"))
EURI = os.environ.get("SVAR_EURI", os.path.join(ROOT, "data", "euribor"))
PUBLIC = os.path.join(ROOT, "data", "public")
OUT = os.environ.get("SVAR_OUT", os.path.join(HERE, "output"))

LAGS, H, NBOOT = 2, 16, 1000

SPECS = {
    "baseline": dict(area=["loil", "lgdp_ea", "lp_ea", "rate", "leur"],
                     dom=["lgdp_es", "lp_es", "lcons_es"],
                     years=(1999, 2019), shock="rate", norm="rate"),
    "euribor": dict(area=["loil", "lgdp_ea", "lp_ea", "rate", "euribor", "leur"],
                    dom=["lgdp_es", "lp_es", "lcons_es"],
                    years=(1999, 2019), shock="rate", norm="rate"),
    "indicator": dict(area=["loil", "lgdp_ea", "lp_ea", "euribor", "leur"],
                      dom=["lgdp_es", "lp_es", "lcons_es"],
                      years=(1999, 2019), shock="euribor", norm="euribor"),
    "deposits": dict(area=["loil", "lgdp_ea", "lp_ea", "euribor", "leur"],
                     dom=["lgdp_es", "lp_es", "deposit", "lcons_es"],
                     years=(2003, 2019), shock="euribor", norm="euribor"),
    "baseline_2003": dict(area=["loil", "lgdp_ea", "lp_ea", "rate", "leur"],
                          dom=["lgdp_es", "lp_es", "lcons_es"],
                          years=(2003, 2019), shock="rate", norm="rate"),
    "internal": dict(area=["jk", "loil", "lgdp_ea", "lp_ea", "rate", "euribor", "leur"],
                     dom=["lgdp_es", "lp_es", "lcons_es"],
                     years=(1999, 2019), shock="jk", norm="euribor"),
}
EXOGENOUS = ("loil", "jk")       # depend only on their own lags


# ---------------------------------------------------------------- data
def _find(folder, sub):
    for p in glob.glob(os.path.join(folder, "*")):
        if sub in os.path.basename(p):
            return p
    raise FileNotFoundError(sub)


def _quarter(period_index):
    return [f"{p.year}Q{p.quarter}" for p in period_index]


def _eurostat_rows(name, col, geo):
    ws = openpyxl.load_workbook(_find(DATA, name), data_only=True)["Sheet 1"]
    out = {}
    for i in range(11, ws.max_row + 1):
        g = ws.cell(i, 1).value
        if g and str(g).startswith(geo):
            v = ws.cell(i, col).value
            if isinstance(v, (int, float)):
                out[str(ws.cell(i, 2).value).replace("-Q", "Q")] = float(v)
    return pd.Series(out)


def _ecb_portal_csv(folder, sub, valcol=2):
    for p in glob.glob(os.path.join(folder, "*.csv")):
        if sub in open(p).readline():
            out = {}
            for line in open(p).read().splitlines()[1:]:
                parts = [x.strip('"') for x in line.split('","')]
                parts[0] = parts[0].strip('"')
                if len(parts) <= valcol or parts[valcol] in ("", "-"):
                    continue
                try:
                    out[parts[0][:7]] = float(parts[valcol])
                except ValueError:
                    pass
            return pd.Series(out)
    raise FileNotFoundError(sub)


def _monthly_to_quarterly(s, how="mean"):
    """s: Series indexed by 'YYYY-MM' strings."""
    df = s.to_frame("v")
    df["q"] = df.index.str[:4] + "Q" + ((df.index.str[5:7].astype(int) - 1) // 3 + 1).astype(str)
    g = df.groupby("q")["v"]
    return g.mean() if how == "mean" else g.sum()


def _hicp_sa():
    ws = openpyxl.load_workbook(_find(DATA, "prc_hicp_minr"), data_only=True)["Sheet 1"]
    ea, es = {}, {}
    for i in range(11, ws.max_row + 1):
        t = ws.cell(i, 1).value
        if not t or "-" not in str(t):
            continue
        for col, d in ((2, ea), (4, es)):
            v = ws.cell(i, col).value
            if isinstance(v, (int, float)):
                d[str(t)] = float(v)

    def sa(d):
        s = pd.Series({pd.Period(k, "M"): np.log(v) for k, v in d.items()}).sort_index()
        s.index = s.index.to_timestamp()
        adj = s - STL(s, period=12, robust=True).fit().seasonal
        df = adj.to_frame("v")
        df["q"] = df.index.year.astype(str) + "Q" + ((df.index.month - 1) // 3 + 1).astype(str)
        return df.groupby("q")["v"].mean() * 100

    return sa(ea), sa(es)


def _oil_and_fx():
    br = pd.read_excel(_find(DATA, "RBRTE"), sheet_name="Data 1", skiprows=3, header=None, names=["date", "usd"])
    br["date"] = pd.to_datetime(br["date"], errors="coerce")
    br = br.dropna()
    fx = _ecb_portal_csv(DATA, "US dollar/Euro")
    br["ym"] = br.date.dt.strftime("%Y-%m")
    br["eur"] = br.usd / br.ym.map(fx)
    br = br.dropna(subset=["eur"])
    br["q"] = br.ym.str[:4] + "Q" + ((br.date.dt.month - 1) // 3 + 1).astype(str)
    return np.log(br.groupby("q")["eur"].mean()) * 100, np.log(_monthly_to_quarterly(fx)) * 100


def _shadow_rate():
    e1 = _ecb_portal_csv(EURI, "Euribor 1-month")
    sr = pd.read_excel(_find(DATA, "shadowrate_ECB"), header=None, names=["ym", "sr"])
    sr["ym"] = sr.ym.astype(int).astype(str).str.replace(r"(\d{4})(\d{2})", r"\1-\2", regex=True)
    sh = sr.set_index("ym")["sr"]
    lo, hi = sh.index.min(), sh.index.max()
    off_lo, off_hi = sh[lo] - e1[lo], sh[hi] - e1[hi]
    rate = {}
    for m in sorted(set(e1.index) | set(sh.index)):
        if lo <= m <= hi and m in sh.index:
            rate[m] = sh[m]
        elif m < lo and m in e1.index:
            rate[m] = e1[m] + off_lo
        elif m > hi and m in e1.index:
            rate[m] = e1[m] + off_hi
    return _monthly_to_quarterly(pd.Series(rate))


def _public_monthly(name):
    df = pd.read_csv(os.path.join(PUBLIC, f"{name}.csv"))
    return pd.Series(df["value"].values, index=df["date"].astype(str).str[:7])


def _deposit_rate():
    r_on, r_t = _public_monthly("dep_rate_overnight_es"), _public_monthly("dep_rate_time_stock_es")
    v_on, v_t = _public_monthly("dep_vol_overnight_es"), _public_monthly("dep_vol_time_es")
    return _monthly_to_quarterly(((r_on * v_on + r_t * v_t) / (v_on + v_t)).dropna())


def _jk_shock():
    jk = pd.read_csv(os.path.join(PUBLIC, "jk_shocks_monthly.csv"))
    s = pd.Series(jk["MP_median"].values, index=[f"{y}-{m:02d}" for y, m in zip(jk.year, jk.month)])
    return _monthly_to_quarterly(s, how="sum")


def build_panel():
    loil, leur = _oil_and_fx()
    lp_ea, lp_es = _hicp_sa()
    P = pd.DataFrame({
        "jk": _jk_shock(),
        "loil": loil,
        "lgdp_ea": np.log(_eurostat_rows("namq_10_gdp", 5, "Euro")) * 100,
        "lp_ea": lp_ea,
        "rate": _shadow_rate(),
        "euribor": _monthly_to_quarterly(_public_monthly("euribor_1y")),
        "leur": leur,
        "lgdp_es": np.log(_eurostat_rows("namq_10_gdp", 5, "Spain")) * 100,
        "lp_es": lp_es,
        "deposit": _deposit_rate(),
        "lcons_es": np.log(_eurostat_rows("namq_10_fcs", 3, "Spain")) * 100,
    })
    P["year"] = P.index.str[:4].astype(int)
    P["qn"] = P.index.str[5:].astype(int)
    return P.sort_values(["year", "qn"]).drop(columns="qn")


def spec_data(P, spec):
    s = SPECS[spec]
    names = s["area"] + s["dom"]
    y0, y1 = s["years"]
    D = P[(P.year >= y0) & (P.year <= y1)][names].dropna()
    return D, names


# ---------------------------------------------------------------- VAR
def estimate(Y, names, area, p=LAGS):
    """Equation-by-equation OLS with block exogeneity: exogenous variables depend
    on their own lags, euro-area variables on euro-area lags, Spanish variables on
    all lags."""
    T, n = Y.shape
    Ylag = np.hstack([Y[p - 1 - k:T - 1 - k] for k in range(p)])
    yt = Y[p:]
    const = np.ones((T - p, 1))
    trend = np.arange(p, T).reshape(-1, 1) * 1.0
    A = [np.zeros((n, n)) for _ in range(p)]
    C, G = np.zeros(n), np.zeros(n)
    resid = np.zeros((yt.shape[0], n))
    area_idx = [names.index(v) for v in area]
    for i, v in enumerate(names):
        if v in EXOGENOUS:
            allowed = [i]
        elif v in area:
            allowed = area_idx
        else:
            allowed = list(range(n))
        cols = [k * n + j for k in range(p) for j in allowed]
        X = np.hstack([const, trend, Ylag[:, cols]])
        b = np.linalg.lstsq(X, yt[:, i], rcond=None)[0]
        resid[:, i] = yt[:, i] - X @ b
        C[i], G[i] = b[0], b[1]
        for k in range(p):
            for jj, j in enumerate(allowed):
                A[k][i, j] = b[2 + k * len(allowed) + jj]
    return A, C, G, np.cov(resid, rowvar=False), resid


def irf(A, Sig, shock, norm, p=LAGS, horizon=H):
    """Responses to the Cholesky shock in position `shock`, scaled so that the
    variable in position `norm` rises by one unit on impact."""
    n = Sig.shape[0]
    L = np.linalg.cholesky(Sig)
    d = L[:, shock] / L[norm, shock]
    Phi = [np.eye(n)]
    for h in range(1, horizon + 1):
        M = np.zeros((n, n))
        for k in range(p):
            if h - 1 - k >= 0:
                M += A[k] @ Phi[h - 1 - k]
        Phi.append(M)
    return np.array([Phi[h] @ d for h in range(horizon + 1)])


def bootstrap(Y, names, area, A, C, G, resid, shock, norm, p=LAGS, horizon=H, draws=NBOOT):
    T, n = Y.shape
    out = np.full((draws, horizon + 1, n), np.nan)
    for b in range(draws):
        Yb = np.zeros((T, n))
        Yb[:p] = Y[:p]
        idx = np.random.randint(0, resid.shape[0], size=T - p)
        for t in range(p, T):
            v = C + G * t
            for k in range(p):
                v = v + A[k] @ Yb[t - 1 - k]
            Yb[t] = v + resid[idx[t - p]]
        try:
            Ab, _, _, Sb, _ = estimate(Yb, names, area, p)
            out[b] = irf(Ab, Sb, shock, norm, p, horizon)
        except np.linalg.LinAlgError:
            pass
    return out


def run(spec, P=None, seed=0):
    """Point responses and bootstrap draws of one specification."""
    P = build_panel() if P is None else P
    s = SPECS[spec]
    D, names = spec_data(P, spec)
    Y = D.values
    shock, norm = names.index(s["shock"]), names.index(s["norm"])
    A, C, G, Sig, resid = estimate(Y, names, s["area"])
    point = irf(A, Sig, shock, norm)
    np.random.seed(seed)
    draws = bootstrap(Y, names, s["area"], A, C, G, resid, shock, norm)
    return dict(spec=spec, names=names, point=point, draws=draws, sample=(D.index[0], D.index[-1]), nobs=len(D))


def to_frame(r):
    cols = {"h": list(range(H + 1))}
    for i, v in enumerate(r["names"]):
        cols[v] = r["point"][:, i]
        cols[v + "_lo68"] = np.nanpercentile(r["draws"][:, :, i], 16, axis=0)
        cols[v + "_hi68"] = np.nanpercentile(r["draws"][:, :, i], 84, axis=0)
        cols[v + "_lo90"] = np.nanpercentile(r["draws"][:, :, i], 5, axis=0)
        cols[v + "_hi90"] = np.nanpercentile(r["draws"][:, :, i], 95, axis=0)
    return pd.DataFrame(cols)


def main(specs):
    os.makedirs(OUT, exist_ok=True)
    P = build_panel()
    spec_data(P, "baseline")[0].to_csv(os.path.join(OUT, "svar_panel.csv"))
    for spec in specs:
        r = run(spec, P)
        df = to_frame(r)
        df.to_csv(os.path.join(OUT, f"svar_irf_{spec}.csv"), index=False)
        c = df["lcons_es"]
        print(f"{spec:9s} {r['sample'][0]}-{r['sample'][1]} (n={r['nobs']}): consumption impact {c.iloc[0]:+.2f}%, "
              f"trough {c.iloc[:9].min():+.2f}% at h={int(c.iloc[:9].idxmin())}")


if __name__ == "__main__":
    main(sys.argv[1:] or list(SPECS))

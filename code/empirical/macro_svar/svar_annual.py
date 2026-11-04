#!/usr/bin/env python3
# svar_annual.py
# Time-aggregates the VAR responses of svar.py to annual frequency, so that they
# can be compared with the annual model. For each specification and each year
# after the shock it reports the average of the four quarterly responses of
# consumption, with bootstrap bands computed on the same draws as svar.py, and
# the annual averages of the rates and of the Spanish price level.
#
# Output: output/svar_annual_<spec>.csv
# Run:  python svar_annual.py [spec ...]     (all specifications by default)

import os
import sys

import numpy as np
import pandas as pd

import svar as V

YEARS = 3
LEVELS = ("rate", "euribor", "deposit", "lp_es")


def annual(r, years=YEARS):
    rows = []
    names = r["names"]
    c = names.index("lcons_es")
    for k in range(years):
        q = slice(4 * k, 4 * k + 4)
        draws = np.nanmean(r["draws"][:, q, c], axis=1)
        lo90, lo68, hi68, hi90 = np.nanpercentile(draws, [5, 16, 84, 95])
        row = dict(year=k, cons=r["point"][q, c].mean(),
                   cons_lo68=lo68, cons_hi68=hi68, cons_lo90=lo90, cons_hi90=hi90)
        for v in LEVELS:
            if v in names:
                row[v] = r["point"][q, names.index(v)].mean()
        rows.append(row)
    return pd.DataFrame(rows)


def main(specs):
    os.makedirs(V.OUT, exist_ok=True)
    P = V.build_panel()
    for spec in specs:
        r = V.run(spec, P)
        out = annual(r)
        out.to_csv(os.path.join(V.OUT, f"svar_annual_{spec}.csv"), index=False)
        y = out.iloc[0]
        print(f"{spec:14s} {r['sample'][0]}-{r['sample'][1]}: year-1 consumption {y['cons']:+.2f}% "
              f"(68% [{y['cons_lo68']:+.2f}, {y['cons_hi68']:+.2f}]), year 2 {out['cons'].iloc[1]:+.2f}%")


if __name__ == "__main__":
    main(sys.argv[1:] or list(V.SPECS))

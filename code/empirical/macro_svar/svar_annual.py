#!/usr/bin/env python3
# svar_annual.py
# Time-aggregates the recursive SVAR of svar.py to annual frequency, so that its
# consumption response can be compared with the annual HANK model. For each year
# after the shock it reports the average of the four quarterly responses, with
# residual-bootstrap bands computed on the same draws as svar.py (seed 0), and the
# annual average of the identified policy-rate path.
#
# Output: output/svar_annual.csv

import os
import numpy as np
import pandas as pd
import svar as V


def main(years=3):
    P = V.build_panel()
    Y = P[V.NAMES].values
    A, C, G, Sig, resid = V.estimate(Y, V.P_LAG)
    point = V.irf(A, Sig, V.P_LAG, V.H)
    np.random.seed(0)
    boot = V.bootstrap(Y, A, C, G, resid, V.P_LAG, V.H, V.NBOOT)
    rows = []
    for k in range(years):
        q = slice(4 * k, 4 * k + 4)
        draws = np.nanmean(boot[:, q, V.CONS], axis=1)
        lo90, lo68, hi68, hi90 = np.nanpercentile(draws, [5, 16, 84, 95])
        rows.append(dict(year=k,
                         cons=point[q, V.CONS].mean(),
                         cons_lo68=lo68, cons_hi68=hi68, cons_lo90=lo90, cons_hi90=hi90,
                         rate=point[q, V.RATE].mean()))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(V.OUT, "svar_annual.csv"), index=False)
    for r in rows:
        print(f"year {r['year']}: consumption {r['cons']:+.2f}% "
              f"(68% [{r['cons_lo68']:+.2f}, {r['cons_hi68']:+.2f}]; "
              f"90% [{r['cons_lo90']:+.2f}, {r['cons_hi90']:+.2f}])   "
              f"policy rate {r['rate']:.2f} pp")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# import_content.py
# Import content of Spanish household final consumption from Eurostat's symmetric
# input-output tables at basic prices (naio_10_cp1700, product by product): the
# imported final goods households buy directly plus the imported intermediate
# inputs embodied in the domestic goods they buy,
#     m_C = [ m' (I - A_d)^{-1} f_d + 1' f_m ] / [ 1' f_d + 1' f_m ],
# where A_d are domestic input coefficients, m the imported inputs per unit of
# output, and f_d, f_m household consumption of domestic and imported products.
#
# Output: output/import_content.csv. Tables are cached under data/public/io/.

import json
import os
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "..", "..", "data", "public", "io"))
OUT = os.path.join(HERE, "output")
URL = ("https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/naio_10_cp1700?"
       "geo={geo}&time={year}&unit=MIO_EUR&stk_flow={flow}&format=JSON&lang=en")


def table(geo, year, flow):
    path = os.path.join(DATA, f"naio_10_cp1700_{geo}_{year}_{flow}.json")
    if not os.path.exists(path):
        os.makedirs(DATA, exist_ok=True)
        with urllib.request.urlopen(URL.format(geo=geo, year=year, flow=flow), timeout=180) as r, \
                open(path, "wb") as f:
            f.write(r.read())
    d = json.load(open(path))
    dims = d["id"]
    sizes = d["size"]
    cats = {k: {v: i for i, v in enumerate(sorted(d["dimension"][k]["category"]["index"],
                                                  key=d["dimension"][k]["category"]["index"].get))}
            for k in dims}
    inv = {k: {i: c for c, i in cats[k].items()} for k in dims}
    rows = []
    strides = np.cumprod([1] + sizes[::-1])[::-1][1:]
    for flat, val in d["value"].items():
        flat = int(flat)
        idx = {}
        for k, s in zip(dims, strides):
            idx[k], flat = divmod(flat, s)
        rows.append((inv["prd_ava"][idx["prd_ava"]], inv["prd_use"][idx["prd_use"]], val))
    return pd.DataFrame(rows, columns=["ava", "use", "value"]).pivot(index="ava", columns="use", values="value")


def import_content(geo, year):
    D, M = table(geo, year, "DOM"), table(geo, year, "IMP")
    prods = sorted(c for c in D.columns if c.startswith("CPA_") and c != "CPA_TOTAL" and c in D.index)
    x = D.loc["P1", prods].fillna(0).values
    Zd = D.loc[prods, prods].fillna(0).values
    Zm = M.loc[prods, prods].fillna(0).values
    fd = D.loc[prods, "P3_S14"].fillna(0).values
    fm = M.loc[prods, "P3_S14"].fillna(0).values
    keep = x > 0
    Ad = np.zeros_like(Zd)
    Ad[:, keep] = Zd[:, keep] / x[keep]
    m = np.zeros(len(prods))
    m[keep] = Zm[:, keep].sum(0) / x[keep]
    indirect = m @ np.linalg.solve(np.eye(len(prods)) - Ad, fd)
    return dict(geo=geo, year=year, direct=fm.sum() / (fd.sum() + fm.sum()),
                total=(indirect + fm.sum()) / (fd.sum() + fm.sum()), hfce=fd.sum() + fm.sum())


def main():
    os.makedirs(OUT, exist_ok=True)
    res = pd.DataFrame([import_content("ES", y) for y in (2020, 2021, 2022)])
    res.to_csv(os.path.join(OUT, "import_content.csv"), index=False)
    print(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()

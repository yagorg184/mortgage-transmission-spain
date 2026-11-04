#!/usr/bin/env python3
# download.py
# Public inputs of the local projections and of the 2022-24 episode, downloaded
# from their original sources and stored under data/public/:
#   - Jarocinski-Karadi ECB monetary policy shocks (monthly), CC BY 4.0
#   - ECB Data Portal: 1-year Euribor, HICP indices, Spanish household deposit
#     rates (MIR) and deposit volumes (BSI)
#   - Eurostat: quarterly household consumption and GDP, chain-linked volumes,
#     seasonally and calendar adjusted
# Run once; the stored files are what the estimation scripts read.

import io
import json
import os
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "..", "data", "public"))

JK_URL = "https://raw.githubusercontent.com/marekjarocinski/jkshocks_update_ecb/main/shocks_ecb_mpd_me_m.csv"
ECB_API = "https://data-api.ecb.europa.eu/service/data/{}?format=csvdata&detail=dataonly"
ESTAT_API = ("https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/{}?"
             "{}&format=JSON&lang=en")

ECB_SERIES = {
    "euribor_1y": "FM/M.U2.EUR.RT.MM.EURIBOR1YD_.HSTA",
    "hicp_es": "ICP/M.ES.N.000000.4.INX",
    "hicp_de": "ICP/M.DE.N.000000.4.INX",
    "hicp_ea": "ICP/M.U2.N.000000.4.INX",
    "dep_rate_overnight_es": "MIR/M.ES.B.L21.A.R.A.2250.EUR.N",
    "dep_rate_time_new_es": "MIR/M.ES.B.L22.F.R.A.2250.EUR.N",
    "dep_rate_time_stock_es": "MIR/M.ES.B.L22.A.R.A.2250.EUR.O",
    "dep_vol_overnight_es": "BSI/M.ES.N.A.L21.A.1.U2.2250.Z01.E",
    "dep_vol_time_es": "BSI/M.ES.N.A.L22.A.1.U2.2250.Z01.E",
}

ESTAT_SERIES = {
    "cons": ("namq_10_fcs", "na_item=P31_S14&unit=CLV10_MEUR&s_adj=SCA"),
    "gdp": ("namq_10_gdp", "na_item=B1GQ&unit=CLV10_MEUR&s_adj=SCA"),
}
GEOS = {"es": "ES", "de": "DE", "ea": "EA20"}


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def ecb(key):
    df = pd.read_csv(io.BytesIO(fetch(ECB_API.format(key))))
    return df[["TIME_PERIOD", "OBS_VALUE"]].rename(columns={"TIME_PERIOD": "date", "OBS_VALUE": "value"})


def eurostat(dataset, query, geo):
    d = json.loads(fetch(ESTAT_API.format(dataset, f"{query}&geo={geo}")))
    idx = d["dimension"]["time"]["category"]["index"]
    inv = {v: k for k, v in idx.items()}
    rows = [(inv[int(k)], v) for k, v in d["value"].items()]
    return pd.DataFrame(sorted(rows), columns=["date", "value"])


def main():
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "jk_shocks_monthly.csv"), "wb") as f:
        f.write(fetch(JK_URL))
    for name, key in ECB_SERIES.items():
        ecb(key).to_csv(os.path.join(OUT, f"{name}.csv"), index=False)
    for var, (ds, q) in ESTAT_SERIES.items():
        for tag, geo in GEOS.items():
            eurostat(ds, q, geo).to_csv(os.path.join(OUT, f"{var}_{tag}.csv"), index=False)
    print(f"stored in {OUT}")


if __name__ == "__main__":
    main()

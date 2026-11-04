# robustness.py
# Sensitivity of the headline numbers to numerical choices and to external
# parameters. Writes output/robustness.json. Run:  python robustness.py
#
#   grids        denser asset grids at the calibrated parameters (numerical accuracy)
#   import       import share of consumption 0.20 and 0.30 (baseline 0.25)
#   banks        foreign ownership 0 and 1 of bank equity; no corporate tax
#   deposits     deposit pass-through 0.32 and 1 (baseline 0.12)
#   persistence  shock persistence 0.3 and 0.85 (baseline 0.6)
#   calibrations the alternative calibrations in calibration/alternatives/, written
#                by calibrate.py: other points of the set of parameters that match
#                the four population targets (--beta-patient) and the calibration
#                that also matches the hand-to-mouth share of each tercile of
#                mortgagors (--terciles)

import copy
import glob
import json
import os

import numpy as np

import config as CF
import economy as EC
import experiments as X
import paths as PT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
DEPOSIT_BETA = PT.DEPOSIT_BETA
GERMANY = X.contract_shares(CF.load()[0])["Germany"]


def headline(econ, par, paths):
    v = X.variable_share_spain()
    es = X.aggregate(econ, par, paths, v)["total"][0]
    de = X.aggregate(econ, par, paths, GERMANY)["total"][0]
    f0 = X.aggregate(econ, par, paths, 0.0)["total"][0]
    return dict(spain=round(es, 4), germany_mix=round(de, 4), fixed_only=round(f0, 4),
                ratio_es_de=round(es / de, 3), share_contracts=round((es - f0) / es, 3),
                contracts_pp=round(es - f0, 4))


def episode_checks(econ, par):
    """Consumption gap between variable- and fixed-rate borrowers and the effect
    of the Spanish contract mix relative to the German one, 2022-2025."""
    T = econ.T
    ep = PT.episode(T)
    rates = dict(di=ep["di"], ddep=ep["ddep"], dpi=np.zeros(T))
    out = {}
    for lab, learn in (("perfect_foresight", False), ("learning", True)):
        gap = X.borrower_gap(econ, par, ep, learning=learn)[:4]
        es = X.aggregate(econ, par, rates, X.variable_share_spain(), learning=learn)["total"][:4]
        de = X.aggregate(econ, par, rates, GERMANY, learning=learn)["total"][:4]
        out[lab] = dict(borrower_gap=[round(float(x), 4) for x in gap],
                        spain_minus_germany_mix=[round(float(x), 4) for x in es - de])
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    econ, par, cal = X.build()
    T = econ.T
    base = PT.stylized_shock(T, deposit_beta=DEPOSIT_BETA)
    res = {"baseline": dict(headline(econ, par, base), episode=episode_checks(econ, par))}

    _, hh, types, groups = CF.load()
    dense = copy.deepcopy(types)
    for t in dense:
        t["grid"] = dict(t["grid"], nA=int(t["grid"]["nA"] * 1.5))
    hh_dense = dict(hh, nB=45)
    econ_d = EC.Economy(hh_dense, dense, groups, T)
    res["grids_dense"] = dict(headline(econ_d, par, base), moments={k: round(v, 4) for k, v in econ_d.moments.items()})

    for m in (0.20, 0.30):
        p = copy.deepcopy(par)
        p["union"]["import_share"] = m
        res[f"import_{m}"] = headline(econ, p, base)
    for lam in (0.0, 1.0):
        p = copy.deepcopy(par)
        p["banks"]["foreign_share"] = lam
        res[f"foreign_share_{lam}"] = headline(econ, p, base)
    p = copy.deepcopy(par)
    p["banks"]["corporate_tax"] = 0.0
    res["no_corporate_tax"] = headline(econ, p, base)
    for beta in (0.32, 1.0):
        res[f"deposit_beta_{beta}"] = headline(econ, par, PT.stylized_shock(T, deposit_beta=beta))
    for rho in (0.3, 0.85):
        res[f"persistence_{rho}"] = headline(econ, par, PT.stylized_shock(T, rho=rho, deposit_beta=DEPOSIT_BETA))
    for path in sorted(glob.glob(os.path.join(CF.ALTERNATIVES, "*.json"))):
        econ_i, par_i, cal_i = X.build(T, calibration=path)
        h = cal_i["household"]
        res["calibration_" + os.path.basename(path)[:-len(".json")]] = dict(
            headline(econ_i, par_i, base),
            episode=episode_checks(econ_i, par_i),
            parameters=dict(beta_modal=h["types"][0]["beta"], beta_patient=h["types"][1]["beta"],
                            patient_share=h["types"][1]["mass"], chi1=h["chi1"], sig_eps=h["sig_eps"],
                            modal_beta_shift=h.get("modal_beta_shift")),
            moments={k: round(v, 4) for k, v in econ_i.moments.items()},
            group_moments={g: {k: round(v, 4) for k, v in row.items()} for g, row in econ_i.group_moments.items()})
    with open(os.path.join(OUT, "robustness.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

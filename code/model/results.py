# results.py
# Every number reported in the paper that comes from the model, computed from the
# calibrated economy and written to output/results.json. Run:  python results.py

import json
import os

import numpy as np
import pandas as pd

import experiments as X
import mortgages as MG
import paths as PT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")

DEPOSIT_BETA = PT.DEPOSIT_BETA
SVAR_SPECS = ("baseline", "euribor", "indicator", "deposits", "baseline_2003")


def r(x, d=4):
    return [round(float(v), d) + 0.0 for v in np.atleast_1d(x)]     # + 0.0 turns -0.0 into 0.0


def episode_block(econ, par, paths, shares, learning):
    v_es = shares["Spain"]
    es = X.aggregate(econ, par, paths, v_es, learning=learning)
    de = X.aggregate(econ, par, paths, shares["Germany"], learning=learning)
    f0 = X.aggregate(econ, par, paths, 0.0, learning=learning)
    return {
        "aggregate_spain": r(es["total"][:4]),
        "aggregate_germany_mix": r(de["total"][:4]),
        "aggregate_fixed_only": r(f0["total"][:4]),
        "spain_minus_germany_mix": r((es["total"] - de["total"])[:4]),
        "channels_spain": {k: r(v[:4]) for k, v in es["ge"].items()},
        "borrower_gap": r(X.borrower_gap(econ, par, paths, learning=learning)[:4]),
    }


def facts():
    """Observed rates and prices quoted in the paper (Section 2)."""
    obs = PT.annual_observed()
    v_on, v_t = PT._monthly("dep_vol_overnight_es"), PT._monthly("dep_vol_time_es")
    share_on = v_on / (v_on + v_t)
    overnight = PT._monthly("dep_rate_overnight_es")
    new_time = PT._monthly("dep_rate_time_new_es")
    io = pd.read_csv(os.path.join(PT.ROOT, "code", "empirical", "io", "output", "import_content.csv")).set_index("year")
    d_eur = obs["euribor"] - obs.loc[2021, "euribor"]
    d_dep = obs["deposit_rate"] - obs.loc[2021, "deposit_rate"]
    return {
        "annual_percent": {c: {int(y): round(float(v), 3) for y, v in obs[c].items()} for c in obs},
        "deposit_pass_through_since_2021": {y: round(float(d_dep[y] / d_eur[y]), 3) for y in (2023, 2024)},
        "overnight_share_of_deposits": {y: round(float(share_on[str(y)].mean()), 3) for y in (2021, 2022, 2023, 2024)},
        "overnight_rate_max_2022_2025": round(float(overnight["2022":"2025"].max()), 3),
        "new_time_deposit_rate_2024": round(float(new_time["2024"].mean()), 2),
        "import_content_2021": {"direct": round(float(io.loc[2021, "direct"]), 4),
                                "total": round(float(io.loc[2021, "total"]), 4)},
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    econ, par, cal = X.build()
    T = econ.T
    v_es = X.variable_share_spain()
    shares = X.contract_shares(cal)
    res = {"T": T, "variable_share_spain": v_es, "deposit_beta": DEPOSIT_BETA,
           "import_share": par["union"]["import_share"], "facts": facts()}

    # Mortgage groups: steady-state payment, payment response and moments
    res["groups"] = {}
    for grp in econ.groups:
        row = {k: round(float(v), 4) for k, v in econ.group_moments[grp["name"]].items()}
        if grp["balance"] > 0:
            row.update(balance_income=round(grp["balance"], 4), rate=round(grp["rate"], 4), term=grp["term"],
                       payment_income=round(grp["payment"], 4),
                       annuity_slope=round(MG.annuity_slope(grp["rate"], grp["term"]), 4),
                       origination_share=round(MG.origination_share(grp["rate"], grp["term"],
                                                                    par["nominal"]["inflation_ss"]), 4))
        res["groups"][grp["name"]] = row
    res["steady_state"] = {"C": round(econ.ss["C"], 4), "liquid_B": round(econ.ss["B"], 4),
                           "illiquid_A": round(econ.ss["A"], 4),
                           "mortgage_stock": round(sum(g["mass"] * g["balance"] for g in econ.groups), 4)}

    # Stylized persistent shock: decomposition and contract line
    sty = PT.stylized_shock(T, deposit_beta=DEPOSIT_BETA)
    a = X.aggregate(econ, par, sty, v_es)
    direct0 = sum(v[0] for v in a["direct"].values())
    res["stylized"] = {
        "total_first3": r(a["total"][:3]),
        "ge_impact_by_channel": {k: r(v[0])[0] for k, v in a["ge"].items()},
        "direct_impact_by_channel": {k: r(v[0])[0] for k, v in a["direct"].items()},
        "direct_impact_total": round(float(direct0), 4),
        "ge_multiplier_impact": round(float(a["total"][0] / direct0), 3),
        "contract_line_impact": {k: round(v, 4) for k, v in X.contract_line(econ, par, sty, shares).items()},
        "contract_line_first3": {k: r(X.aggregate(econ, par, sty, v)["total"][:3]) for k, v in shares.items()},
        "ge_impact_by_channel_fixed_only": {k: r(v[0])[0] for k, v in X.aggregate(econ, par, sty, 0.0)["ge"].items()},
        "ge_impact_by_channel_variable_only": {k: r(v[0])[0] for k, v in X.aggregate(econ, par, sty, 1.0)["ge"].items()},
    }
    for beta in (0.32, 1.0):
        s2 = PT.stylized_shock(T, deposit_beta=beta)
        res["stylized"][f"contract_line_impact_depbeta_{beta}"] = {
            k: round(v, 4) for k, v in X.contract_line(econ, par, s2, shares).items()}
    for rho in (0.3, 0.85):
        s3 = PT.stylized_shock(T, rho=rho, deposit_beta=DEPOSIT_BETA)
        res["stylized"][f"contract_line_impact_rho_{rho}"] = {
            k: round(v, 4) for k, v in X.contract_line(econ, par, s3, shares).items()}

    # Paths identified by each VAR specification, and the VAR's own consumption
    # response (code/empirical/macro_svar)
    res["svar"] = {}
    for spec in SVAR_SPECS:
        sv = PT.svar_shock(T, spec=spec, deposit_beta=DEPOSIT_BETA)
        a = pd.read_csv(os.path.join(PT.SVAR_OUT, f"svar_annual_{spec}.csv")).set_index("year")
        block = {"paths": {k: r(100 * v[:4]) for k, v in sv.items()},
                 "var_consumption": r(a["cons"].values),
                 "var_band68_year1": r([a.loc[0, "cons_lo68"], a.loc[0, "cons_hi68"]]),
                 "var_band90_year1": r([a.loc[0, "cons_lo90"], a.loc[0, "cons_hi90"]])}
        for lab, learn in (("perfect_foresight", False), ("learning", True)):
            block[lab] = {k: r(X.aggregate(econ, par, sv, v, learning=learn)["total"][:3])
                          for k, v in shares.items()}
        res["svar"][spec] = block

    # 2022-2025 episode. The observed paths are split into the interest-rate paths
    # (Euribor and deposit rate) and the inflation path; by linearity the
    # responses add up to the response to all paths together.
    ep = PT.episode(T)
    zero = np.zeros(T)
    parts = {"rates": dict(di=ep["di"], ddep=ep["ddep"], dpi=zero),
             "inflation": dict(di=zero, ddep=zero, dpi=ep["dpi"]),
             "all": dict(di=ep["di"], ddep=ep["ddep"], dpi=ep["dpi"])}
    res["episode"] = {"years": ep["years"],
                      "paths": {k: r(100 * np.asarray(ep[k])[:5]) for k in ("di", "ddep", "dpi")}}
    for part, p in parts.items():
        res["episode"][part] = {lab: episode_block(econ, par, p, shares, learn)
                                for lab, learn in (("perfect_foresight", False), ("learning", True))}
    res["episode"]["payment_change_median_borrower"] = r(X.payment_change(econ, par, parts["all"])[:4])
    res["episode"]["payment_change_median_borrower_rates_only"] = r(X.payment_change(econ, par, parts["rates"])[:4])

    # Accuracy of the linear mortgage block: largest gap between the linear and
    # the exact payment path in the first 30 years, relative to the largest
    # exact payment deviation, by group, contract and experiment.
    H = 30
    tests = {"stylized": (sty["di"][:H], np.zeros(H)),
             "episode": (np.asarray(ep["di"])[:H], np.asarray(ep["dpi"])[:H]),
             "400bp_rho0.85": (0.04 * 0.85 ** np.arange(H), np.zeros(H))}
    acc = {}
    for name, (di, dpi) in tests.items():
        worst = 0.0
        for grp in econ.groups:
            if grp["balance"] == 0:
                continue
            for k in ("V", "F"):
                args = (grp["rate"], grp["term"], par["nominal"]["inflation_ss"])
                lin = MG.unit_paths(di, dpi, k, *args)["payment"]
                ex = MG.exact_paths(di, dpi, k, *args)["payment"]
                worst = max(worst, float(np.max(np.abs(lin - ex)) / np.max(np.abs(ex))))
        acc[name] = round(worst, 4)
    res["mortgage_linearization_error"] = acc

    # Payment MPC and model moments
    res["mpc_payment"] = {f"rho_{rho}": round(X.implied_mpc_payment(econ, par, rho), 4) for rho in (0.0, 0.6, 0.85)}
    res["moments"] = {k: round(float(v), 4) for k, v in econ.moments.items()}
    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

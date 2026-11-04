# experiments.py
# The model's experiments, built on the calibrated economy:
#   build()            the household sector and the parameter dictionary
#   aggregate()        response of aggregate consumption to exogenous paths
#   borrower_gap()     consumption of variable- relative to fixed-rate borrowers
#   contract_line()    impact response as a function of the variable-rate share

import numpy as np

import config as CF
import economy as EC
import equilibrium as EQ
import mortgages as MG
import paths as PT

T_DEFAULT = 200
RHO_EXPECT = 0.7


def build(T=T_DEFAULT, use_cache=True, calibration=CF.JSON):
    cal, hh, types, groups = CF.load(calibration)
    econ = EC.Economy(hh, types, groups, T, use_cache=use_cache)
    banks = dict(cal["banks"], discount_rate=cal["household"]["fixed"]["ra"])
    par = dict(payment_rule="annuity", redistribution="banks",
               nominal=cal["nominal"], banks=banks, union=cal["union"], ge=cal["ge"])
    return econ, par, cal


def variable_share_spain():
    return float(CF.EFF["variable_share_among_mortgaged"])


def contract_shares(cal):
    """Variable-rate shares of the contract structures compared in the paper:
    Spain's share of mortgagors and the approximate shares of other countries in
    calibration.json."""
    c = cal["country_variable_share_stock"]
    return {"fixed_only": c["all_fixed"], "Germany": c["Germany"], "euro_area": c["euro_area"],
            "Italy": c["Italy"], "Spain": variable_share_spain(), "variable_only": c["all_variable"]}


def aggregate(econ, par, paths, v, learning=False, rho_e=RHO_EXPECT):
    if learning:
        return EQ.respond_learning(econ, par, paths, v, rho_e)
    return EQ.respond(econ, par, paths, v)


def _gap_direct(econ, par, paths):
    """Per-capita consumption difference between variable- and fixed-rate borrowers,
    in percent of mortgagors' consumption, under perfect foresight."""
    num, den = np.zeros(econ.T), 0.0
    for g, grp in enumerate(econ.groups):
        if grp["balance"] == 0:
            continue
        pv = sum(EQ.mortgage_paths(par, paths, "V", grp)[i]["payment"] for i in (0, 1))
        pf = sum(EQ.mortgage_paths(par, paths, "F", grp)[i]["payment"] for i in (0, 1))
        num += grp["mass"] * (econ.J_group("C", "payment", g) @ (grp["balance"] * (pv - pf)))
        den += grp["mass"] * econ.group_ss("C", g)
    return 100 * num / den


def borrower_gap(econ, par, paths, learning=False, rho_e=RHO_EXPECT):
    """Consumption of variable- relative to fixed-rate borrowers. Both groups face
    the same wages and returns, so the general-equilibrium terms cancel and only
    the difference in payments matters."""
    if not learning:
        return _gap_direct(econ, par, paths)
    T = econ.T
    revs = {k: dict(PT.expectation_revisions(np.asarray(paths[k], float), rho_e)) for k in ("di", "ddep", "dpi")}
    out = np.zeros(T)
    for s in sorted(set().union(*[r.keys() for r in revs.values()])):
        p = {k: np.zeros(T) for k in ("di", "ddep", "dpi")}
        for k in p:
            if s in revs[k]:
                p[k][:T - s] = revs[k][s][:T - s]
        out[s:] += _gap_direct(econ, par, p)[:T - s]
    return out


def payment_change(econ, par, paths, contract="V"):
    """Change in the annual mortgage payment of a median borrower, in percent of
    income, along the given paths."""
    grp = [g for g in econ.groups if g["name"] == "tercile2"][0]
    p = sum(EQ.mortgage_paths(par, paths, contract, grp)[i]["payment"] for i in (0, 1))
    return 100 * grp["balance"] * p


def contract_line(econ, par, paths, shares):
    """Impact response of aggregate consumption for each variable-rate share."""
    return {k: float(aggregate(econ, par, paths, v)["total"][0]) for k, v in shares.items()}


def implied_mpc_payment(econ, par, rho):
    """First-year consumption response of mortgagors per euro of a persistent
    increase in their payment (decaying at rate rho)."""
    T = econ.T
    num, den = 0.0, 0.0
    for g, grp in enumerate(econ.groups):
        if grp["balance"] == 0:
            continue
        dp = rho ** np.arange(T) * 0.01
        num += grp["mass"] * (econ.J_group("C", "payment", g) @ dp)[0]
        den += grp["mass"] * dp[0]
    return -num / den

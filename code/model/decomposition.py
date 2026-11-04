# decomposition.py
# Decomposes the aggregate consumption response to a +100 bp shock into the three
# components of the direct impulse in ge.py: the liquid-return (saving) channel
# common to every mortgage structure, the repricing of variable-rate payments, and
# the spending of the extra interest by its recipients. The model is linear in the
# shock, so each component can be pushed through the general-equilibrium closure
# on its own and the pieces add up to the total.

import numpy as np
import config as C
import ge as GE


def components(dbp=100, rho=0.6, vshare=None, kappaw=0.1):
    jc = GE.build_jacs()
    spend = GE.endog_spend()
    T = GE.T
    massV = jc["massV"] if vshare is None else C.MORT_SHARE * vshare
    dr = (dbp / 1e4) * rho ** np.arange(T)
    dpayV = C.MV * dr
    dT = spend * massV * C.MV * dr
    Jpay_agg = jc["massN"] * jc["JCpay0"] + C.MORT_SHARE * jc["JCpayV"]
    direct = {
        "liquid return": jc["JCrb"] @ dr,
        "repricing": massV * (jc["JCpayV"] @ dpayV),
        "interest spending": -Jpay_agg @ dT,
    }
    I = np.eye(T)
    W = GE.wage_op(kappaw)
    A = np.block([[I - jc["JCw"], -jc["JCw"]], [-(1 / GE.NU) * W, I]])   # GHH closure (no wealth effect)
    out = {}
    for name, d in direct.items():
        sol = np.linalg.solve(A, np.concatenate([d, np.zeros(T)]))
        out[name] = (100 * d[0] / jc["Cagg"], 100 * sol[0] / jc["Cagg"])
    return out


if __name__ == "__main__":
    for v, label in [(None, "Spain (EFF share)"), (C.VSHARE["Spain"], "Spain (v = 0.70)")]:
        res = components(vshare=v)
        print(f"{label}: impact response by component, % of steady-state consumption")
        print(f"  {'component':18s} {'direct':>8s} {'GE':>8s}")
        for name, (pe, ge) in res.items():
            print(f"  {name:18s} {pe:+8.3f} {ge:+8.3f}")
        print(f"  {'total':18s} {sum(p for p, _ in res.values()):+8.3f} {sum(g for _, g in res.values()):+8.3f}")

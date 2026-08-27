# ge.py
# General-equilibrium closure (frozen housing) in sequence space, on the
# calibrated household. Blocks: goods/production (Y = C + G, Y = N), the household
# Jacobians, a labour-supply block (the wealth effect, i.e. the explicit GHH
# switch) closed by a wage Phillips curve, and an intermediary/redistribution
# block. Delivers the aggregate consumption impulse response, the derived
# redistribution parameter, and the mortgage-structure counterfactual.

import numpy as np
from sequence_jacobian.examples import two_asset as ta
import household as H
import config as C

NU = 1.0        # Frisch elasticity of labour supply
T = 25          # truncation horizon
BETA_W = 0.96   # discount factor in the wage Phillips curve


def income_pay(e_grid, tax, w, N, payment):
    z_grid = (1 - tax) * w * N * e_grid - payment
    return z_grid


HHp = ta.hh.add_hetinputs([H.make_grids_es, income_pay])


def ss_and_J(payment):
    cal = dict(H.BASE)
    cal.update(dict(beta=C.BETAS[0], chi1=C.CHI1, sig_eps=C.SIG_EPS, ra=C.RA, payment=payment))
    ss = HHp.steady_state(cal)
    J = HHp.jacobian(ss, inputs=["rb", "w", "payment"], T=T)
    return ss, J


def build_jacs():
    """Aggregate household Jacobians, mixing non-mortgaged and mortgaged (variable
    plus fixed) households by their EFF masses. The steady-state mortgage payment
    uses the measured mortgage rate."""
    massN = 1 - C.MORT_SHARE
    massV = C.MORT_SHARE * C.PVAR
    massF = C.MORT_SHARE * (1 - C.PVAR)
    ss0, J0 = ss_and_J(0.0)
    ssV, JV = ss_and_J((C.RATE + C.AMORT) * C.MV)
    agg = lambda o, x: massN * J0[o][x] + (massV + massF) * JV[o][x]
    return dict(JCrb=agg("C", "rb"), JCw=agg("C", "w"),
                JUrb=agg("UCE", "rb"), JUw=agg("UCE", "w"),
                JCpayV=JV["C"]["payment"], JUpayV=JV["UCE"]["payment"],
                JCpay0=J0["C"]["payment"], JUpay0=J0["UCE"]["payment"],
                UCEss=massN * float(ss0["UCE"]) + (massV + massF) * float(ssV["UCE"]),
                Cagg=massN * float(ss0["C"]) + (massV + massF) * float(ssV["C"]),
                massN=massN, massV=massV, massF=massF)


def wage_op(kappaw):
    """Wage Phillips curve operator: dw = W . dMRS. kappaw=None -> flexible wages."""
    if kappaw is None:
        return np.eye(T)
    A = (1 + kappaw + BETA_W) * np.eye(T) - np.eye(T, k=-1) - BETA_W * np.eye(T, k=1)
    return kappaw * np.linalg.inv(A)


def ge_solve(jc, wealth, dr, arm, kappaw=0.1, spend=1.0, vshare=None):
    """Solve the joint {consumption, wage} system for a monetary path dr.

    wealth : 1 standard preferences (wealth effect), 0 GHH (no wealth effect).
    spend  : fraction of the extra mortgage interest that returns to demand
             (0 = savers hoard it, 1 = average recipient).
    vshare : variable-rate share among mortgaged households (None = Spain).
    """
    I = np.eye(T); Wop = wage_op(kappaw)
    massV = jc["massV"] if vshare is None else C.MORT_SHARE * vshare
    dpayV = C.MV * dr * arm
    dTreb = spend * massV * C.MV * dr * arm
    Jpay_agg = jc["massN"] * jc["JCpay0"] + C.MORT_SHARE * jc["JCpayV"]
    JUpay_agg = jc["massN"] * jc["JUpay0"] + C.MORT_SHARE * jc["JUpayV"]
    dCdir = jc["JCrb"] @ dr + massV * (jc["JCpayV"] @ dpayV) - Jpay_agg @ dTreb
    dUdir = (jc["JUrb"] @ dr + massV * (jc["JUpayV"] @ dpayV) - JUpay_agg @ dTreb) / jc["UCEss"]
    JwC = jc["JCw"]; JwU = jc["JUw"] / jc["UCEss"]
    A = np.block([[I - JwC, -JwC],
                  [Wop @ (wealth * JwU) - (1 / NU) * Wop, I + Wop @ (wealth * JwU)]])
    rhs = np.concatenate([dCdir, -Wop @ (wealth * dUdir)])
    sol = np.linalg.solve(A, rhs)
    return 100 * sol[:T] / jc["Cagg"]


def household_mpcs():
    """Average MPC and the deposit-weighted MPC of the calibrated household."""
    cal = dict(H.BASE)
    cal.update(dict(beta=C.BETAS[0], chi1=C.CHI1, sig_eps=C.SIG_EPS, ra=C.RA))
    ss = H.HH.steady_state(cal); I = ss.internals["hh"]
    D = np.array(I["D"]); b = np.array(I["b_grid"]); c = np.array(I["c"])
    mpc = np.clip(np.gradient(c, b, axis=1), 0, 1)
    Bpos = np.maximum(b, 0)[None, :, None] * np.ones_like(c)
    return float((D * mpc).sum()), float((D * mpc * Bpos).sum() / (D * Bpos).sum())


def endog_spend(beta_dep=None, mpc_share=None, lam=None):
    """Redistribution parameter derived from the interest incidence: the share of
    the extra interest that returns to demand = [beta_dep * MPC_deposit +
    (1 - beta_dep)(1 - lambda) * MPC_shareholder] / MPC_average."""
    beta_dep = C.BETA_DEP if beta_dep is None else beta_dep
    mpc_share = C.MPC_SHARE if mpc_share is None else mpc_share
    lam = C.LAMBDA_EXT if lam is None else lam
    avg, mpc_dep = household_mpcs()
    return (beta_dep * mpc_dep + (1 - beta_dep) * (1 - lam) * mpc_share) / avg


def consumption_irf(dbp=100, rho=0.6, wealth=0.0, kappaw=0.1, spend=None):
    """Aggregate consumption IRF of the Spanish (variable-rate) economy to +dbp."""
    if spend is None:
        spend = endog_spend()
    jc = build_jacs()
    dr = (dbp / 1e4) * rho ** np.arange(T)
    return ge_solve(jc, wealth, dr, arm=1.0, kappaw=kappaw, spend=spend)


def counterfactual(spend=None, dbp=100, rho=0.6, wealth=0.0, kappaw=0.1):
    """Same Spain, with the fixed/variable mix of each country. The aggregate
    response is linear in the variable-rate share."""
    if spend is None:
        spend = endog_spend()
    jc = build_jacs()
    dr = (dbp / 1e4) * rho ** np.arange(T)
    resp = lambda f: ge_solve(jc, wealth, dr, arm=1.0, kappaw=kappaw, spend=spend, vshare=f)[0]
    tab = {k: resp(f) for k, f in C.VSHARE.items()}
    es = resp(C.VSHARE["Spain"]); de = resp(C.VSHARE["Germany"]); f0 = resp(0.0)
    return dict(tab=tab, es=es, de=de, ratio=es / de, frac=100 * (es - f0) / es, spend=spend)


if __name__ == "__main__":
    irf = consumption_irf()
    m = counterfactual()
    print(f"spend (derived) = {m['spend']:.3f}")
    print(f"aggregate consumption IRF, +100 bp: impact {irf[0]:+.4f}%")
    print(f"Spain ({C.VSHARE['Spain']}) {m['es']:+.4f}%  |  Germany ({C.VSHARE['Germany']}) {m['de']:+.4f}%  "
          f"-> ratio {m['ratio']:.3f}x  ({m['frac']:.0f}% attributable to the variable structure)")

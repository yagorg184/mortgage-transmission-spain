# mortgage.py
# Frozen mortgage block on the calibrated household. Each household carries a
# fixed mortgage position (type and principal from the EFF) whose payment drains
# income; positions are not re-optimised. Provides the partial-equilibrium
# cash-flow response of variable-rate borrowers to a persistent rate path (the
# object that validates the -2.7% BdE anchor) and the level/covariance
# decomposition of the aggregate response.

import numpy as np
from sequence_jacobian.examples import two_asset as ta
import household as H
import config as C


def income_pay(e_grid, tax, w, N, payment):
    z_grid = (1 - tax) * w * N * e_grid - payment
    return z_grid


HHm = ta.hh.add_hetinputs([H.make_grids_es, income_pay])


def solve_type(beta, payment, over=None):
    cal = dict(H.BASE)
    cal.update(dict(beta=beta, chi1=C.CHI1, sig_eps=C.SIG_EPS, ra=C.RA, payment=payment))
    if over:
        cal.update(over)
    ss = HHm.steady_state(cal)
    return ss, ss.internals["hh"]


def dynamic_response(dbp=400, rho=0.85, T=40):
    """Consumption of the variable-rate group along a persistent mortgage-rate
    path. This is the object comparable to the BdE's -2.7%: a persistent shock,
    not a transitory one. Uses the household's linear impulse response."""
    drm = (dbp / 1e4) * rho ** np.arange(T)
    num = np.zeros(T); Cden = 0.0
    for beta, pib in zip(C.BETAS, C.PIS):
        ss, _ = solve_type(beta, (C.RATE + C.AMORT) * C.MV)
        dC = np.asarray(HHm.impulse_linear(ss, {"payment": drm * C.MV}, outputs=["C"])["C"])
        mass = pib * C.MORT_SHARE * C.PVAR
        num += mass * dC; Cden += mass * float(ss["C"])
    return num / Cden


def cashflow_object(dbp=100):
    """Partial-equilibrium consumption response to +dbp on the variable rate, and
    its decomposition into the level term E[MPC]*E[exposure] and the covariance
    term Cov(MPC, exposure) that would amplify it."""
    types = [("none", 0.0, 1 - C.MORT_SHARE),
             ("V", C.MV, C.MORT_SHARE * C.PVAR),
             ("F", C.MV, C.MORT_SHARE * (1 - C.PVAR))]
    Wl, MPCl, EXPl = [], [], []
    Dvl, MPCvl, Cvl = [], [], []
    b = None
    for beta, pib in zip(C.BETAS, C.PIS):
        for lab, M, mass in types:
            pay = (C.RATE + C.AMORT) * M if lab != "none" else 0.0
            ss, I = solve_type(beta, pay)
            D = np.array(I["D"]); c = np.array(I["c"])
            if b is None:
                b = np.array(I["b_grid"])
            mpc = np.clip(np.gradient(c, b, axis=1), 0, 1)
            expo = M if lab == "V" else 0.0
            Wl.append((pib * mass * D).ravel()); MPCl.append(mpc.ravel())
            EXPl.append(np.full(D.size, expo))
            if lab == "V":
                Dvl.append((pib * mass * D).ravel()); MPCvl.append(mpc.ravel()); Cvl.append(c.ravel())
    W = np.concatenate(Wl); W = W / W.sum()
    MPC = np.concatenate(MPCl); EXP = np.concatenate(EXPl)
    EMPC = float((W * MPC).sum()); EEXP = float((W * EXP).sum())
    cov = float((W * MPC * EXP).sum()) - EMPC * EEXP
    Dv = np.concatenate(Dvl); MPCv = np.concatenate(MPCvl); Cv = np.concatenate(Cvl)
    drm = dbp / 1e4
    resp_var = -drm * C.MV * float((Dv * MPCv).sum()) / float((Dv * Cv).sum())
    return dict(EMPC=EMPC, EEXP=EEXP, level=EMPC * EEXP, cov=cov, resp_var_400=resp_var * 4)


if __name__ == "__main__":
    r = dynamic_response(dbp=400, rho=0.85)
    peak = 100 * r[np.argmax(np.abs(r))]
    print(f"Variable-rate group, +400 bp persistent: impact {100 * r[0]:+.2f}%  peak {peak:+.2f}%  (BdE anchor -2.7%)")
    m = cashflow_object(dbp=100)
    print(f"E[MPC]={m['EMPC']:.3f}  E[exposure]={m['EEXP']:.3f}  Cov={m['cov']:+.4f}  "
          f"variable group +400bp = {100 * m['resp_var_400']:+.2f}%")

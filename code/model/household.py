# household.py
# Two-asset household with discount-factor heterogeneity, annual frequency,
# partial equilibrium. Liquid asset b and illiquid asset a (with an adjustment
# cost), CRRA utility, a persistent-plus-transitory income process, and an
# optional permanent fixed effect. The transitory shock generates the
# precautionary demand for liquidity; the distribution of permanent discount
# types delivers the wealthy-hand-to-mouth tail.
#
# The model structure lives here; the calibrated parameters (beta, chi1,
# sig_eps, ra, and the type masses) are supplied by the caller, which reads
# them from calibration.json. Built on the sequence-space Jacobian toolkit.

import numpy as np
from sequence_jacobian.examples import two_asset as ta
from sequence_jacobian import grids as G


def statdist(Pi, tol=1e-14):
    """Stationary distribution of a Markov transition matrix."""
    d = np.ones(Pi.shape[0]) / Pi.shape[0]
    for _ in range(200000):
        d2 = d @ Pi
        if np.max(np.abs(d2 - d)) < tol:
            return d2
        d = d2
    return d


def make_grids_es(bmax, amax, kmax, nB, nA, nK, nZ, nEps, rho_z, sig_eta, sig_eps):
    """Asset grids and the income process: persistent (Rouwenhorst) x transitory
    (Gauss-Hermite), normalised to unit mean."""
    b_grid = G.agrid(amax=bmax, n=nB)
    a_grid = G.agrid(amax=amax, n=nA)
    k_grid = G.agrid(amax=kmax, n=nK)[::-1].copy()
    zp, _, Pi_p = G.markov_rouwenhorst(rho=rho_z, sigma=sig_eta / np.sqrt(1 - rho_z**2), N=nZ)
    if nEps > 1 and sig_eps > 0:
        nod, wt = np.polynomial.hermite.hermgauss(nEps)
        eps = np.exp(np.sqrt(2) * sig_eps * nod); pe = wt / np.sqrt(np.pi)
    else:
        eps = np.array([1.0]); pe = np.array([1.0])
    e_grid = np.kron(zp, eps)
    Pi = np.kron(Pi_p, np.tile(pe, (len(eps), 1)))
    e_grid = e_grid / (statdist(Pi) @ e_grid)
    return b_grid, a_grid, k_grid, e_grid, Pi


def income_es(e_grid, tax, w, N, fe_mult, payment):
    """After-tax labour income scaled by the permanent fixed effect, net of the
    mortgage payment. The returned name (z_grid) is the het-input output."""
    z_grid = (1 - tax) * w * N * e_grid * fe_mult - payment
    return z_grid


HH = ta.hh.add_hetinputs([make_grids_es, income_es])

# Fixed structural parameters and solution grid. The calibrated parameters
# (beta, chi1, sig_eps) are passed in by solve_type / pooled.
BASE = dict(eis=0.5, rb=0.0, ra=0.04, chi0=0.25, chi2=2.0, tax=0.0, w=1.0, N=1.0,
            bmax=50, amax=80, kmax=1, nB=30, nA=40, nK=4, nZ=5, nEps=3,
            rho_z=0.9831, sig_eta=0.072, fe_mult=1.0, payment=0.0)


def fe_types(nFE=1, sig_fe=0.30):
    """Discretised permanent fixed effect (Gauss-Hermite), levels normalised to
    unit mean. nFE=1 recovers the no-fixed-effect case."""
    if nFE <= 1 or sig_fe <= 0:
        return np.array([1.0]), np.array([1.0])
    nod, wt = np.polynomial.hermite.hermgauss(nFE)
    lev = np.exp(np.sqrt(2) * sig_fe * nod); p = wt / np.sqrt(np.pi)
    lev = lev / (p @ lev)
    return lev, p


def solve_type(beta, chi1, sig_eps, fe_mult=1.0, payment=0.0, over=None):
    """Steady state of one (beta, fixed-effect) household type."""
    cal = dict(BASE)
    cal.update(dict(beta=beta, chi1=chi1, sig_eps=sig_eps, fe_mult=fe_mult, payment=payment))
    if over:
        cal.update(over)
    ss = HH.steady_state(cal)
    return ss, ss.internals['hh']


def pooled(betas, pis, chi1, sig_eps, fes=None, pfes=None, over=None):
    """Mixture over permanent discount types (x fixed-effect types). Returns the
    aggregate moments used for calibration: median liquid/income, median
    illiquid/income, hand-to-mouth share, and the annual MPC."""
    assert abs(sum(pis) - 1) < 1e-9
    if fes is None:
        fes, pfes = np.array([1.0]), np.array([1.0])
    b = a = e = None
    A = B = C = 0.0
    LR, AR, MP, WT = [], [], [], []
    for beta, pib in zip(betas, pis):
        for fe, pf in zip(fes, pfes):
            ss, I = solve_type(beta, chi1, sig_eps, fe_mult=fe, over=over)
            D = np.array(I['D']); c = np.array(I['c']); w = pib * pf
            if b is None:
                b = np.array(I['b_grid']); a = np.array(I['a_grid']); e = np.array(I['e_grid'])
            inc = (e * fe)[:, None, None] * np.ones_like(c)
            Bm = b[None, :, None] * np.ones_like(inc)
            Am = a[None, None, :] * np.ones_like(inc)
            mpc = np.clip(np.gradient(c, b, axis=1), 0, 1)
            LR.append((Bm / inc).ravel()); AR.append((Am / inc).ravel())
            MP.append(mpc.ravel()); WT.append((w * D).ravel())
            A += w * float(ss['A']); B += w * float(ss['B']); C += w * float(ss['C'])
    LR = np.concatenate(LR); AR = np.concatenate(AR)
    MP = np.concatenate(MP); WT = np.concatenate(WT); WT /= WT.sum()

    def wq(x, p):
        o = np.argsort(x)
        return np.interp(p, np.cumsum(WT[o]) / WT.sum(), x[o])

    return dict(A=A, B=B, C=C,
                liq_p50=wq(LR, .5), liq_p10=wq(LR, .1), liq_p90=wq(LR, .9),
                iliq_p50=wq(AR, .5), htm=100 * float(WT[LR <= 1 / 24].sum()),
                mpc=float((WT * MP).sum()))


if __name__ == '__main__':
    import json, os
    cal = json.load(open(os.path.join(os.path.dirname(__file__), '..', '..', 'calibration', 'calibration.json')))
    h = cal['household']
    m = pooled(h['beta'], h['pi'], h['chi1'], h['sig_eps'], over={'ra': h['ra']})
    t = h['targets']
    print("Household fit at the frozen calibration:")
    print(f"  liquid/income  {m['liq_p50']:.4f}  (target {t['liquid_income_p50']})")
    print(f"  illiquid/income {m['iliq_p50']:.4f}  (target {t['illiquid_income_p50']})")
    print(f"  hand-to-mouth  {m['htm']:.2f}%  (target {100 * t['htm_share']:.0f}%)")
    print(f"  annual MPC     {m['mpc']:.4f}  (target {t['mpc_annual']})")

# household.py
# Two-asset household block at annual frequency, built on the sequence-space
# Jacobian toolkit (Auclert, Bardoczy, Rognlie and Straub, 2021). Households hold
# a liquid asset b and an illiquid asset a subject to a convex adjustment cost,
# have CRRA utility, and face a persistent-plus-transitory income process. A
# mortgage enters only through its real payment, which is subtracted from labour
# income. Discount-factor heterogeneity is handled by solving one block per type.

import numpy as np
from sequence_jacobian import grids as G
from sequence_jacobian.hetblocks.hh_twoasset import hh as _hh_twoasset


def statdist(Pi, tol=1e-14):
    """Stationary distribution of a Markov transition matrix."""
    d = np.ones(Pi.shape[0]) / Pi.shape[0]
    for _ in range(200000):
        d2 = d @ Pi
        if np.max(np.abs(d2 - d)) < tol:
            return d2
        d = d2
    return d


def make_grids(bmax, amax, kmax, nB, nA, nK, nZ, nEps, rho_z, sig_eta, sig_eps):
    """Asset grids and the income process: persistent AR(1) (Rouwenhorst) times an
    i.i.d. transitory shock (Gauss-Hermite), normalised to unit mean."""
    b_grid = G.agrid(amax=bmax, n=nB)
    a_grid = G.agrid(amax=amax, n=nA)
    k_grid = G.agrid(amax=kmax, n=nK)[::-1].copy()
    zp, _, Pi_p = G.markov_rouwenhorst(rho=rho_z, sigma=sig_eta / np.sqrt(1 - rho_z**2), N=nZ)
    if nEps > 1 and sig_eps > 0:
        nod, wt = np.polynomial.hermite.hermgauss(nEps)
        eps = np.exp(np.sqrt(2) * sig_eps * nod)
        pe = wt / np.sqrt(np.pi)
    else:
        eps, pe = np.array([1.0]), np.array([1.0])
    e_grid = np.kron(zp, eps)
    Pi = np.kron(Pi_p, np.tile(pe, (len(eps), 1)))
    e_grid = e_grid / (statdist(Pi) @ e_grid)
    return b_grid, a_grid, k_grid, e_grid, Pi


def income(e_grid, tax, w, N, payment):
    """After-tax labour income net of the real mortgage payment."""
    z_grid = (1 - tax) * w * N * e_grid - payment
    return z_grid


HH = _hh_twoasset.add_hetinputs([make_grids, income])

INPUTS = ("rb", "ra", "w", "payment")


def solve_ss(params):
    """Steady state of one household block; `params` must contain every input of
    HH. Returns the steady-state dict and its internals."""
    ss = HH.steady_state(params)
    return ss, ss.internals["hh"]


def jacobian(ss, T, inputs=INPUTS):
    """Jacobians of consumption and of the two asset stocks with respect to the
    model's household inputs."""
    return HH.jacobian(ss, inputs=list(inputs), outputs=["C", "A", "B"], T=T)


def check_distribution(internals, tol_neg=1e-10, tol_top=1e-6):
    """A steady-state distribution is accepted only if it is a proper probability
    distribution and essentially no mass sits on the last node of either asset
    grid (otherwise the grid truncates the wealth distribution)."""
    D = np.asarray(internals["D"])
    total = float(D.sum())
    neg = float(D.min())
    top_a = float(D[:, :, -1].sum())
    top_b = float(D[:, -1, :].sum())
    ok = abs(total - 1) < 1e-8 and neg > -tol_neg and top_a < tol_top and top_b < tol_top
    return ok, dict(total=total, min=neg, top_a=top_a, top_b=top_b)

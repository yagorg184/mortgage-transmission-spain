"""Checks of the mortgage block against exact (nonlinear) loan arithmetic."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mortgages as MG  # noqa: E402

RATE, TERM, PI = 0.0206, 14.0, 0.02


def test_annuity_slope_matches_finite_difference():
    h = 1e-7
    fd = (MG.annuity(RATE + h, TERM) - MG.annuity(RATE - h, TERM)) / (2 * h)
    assert abs(MG.annuity_slope(RATE, TERM) - fd) < 1e-6


def test_annuity_identity():
    """A level payment repays the loan exactly over its term."""
    a = MG.annuity(RATE, TERM)
    b = 1.0
    for _ in range(int(TERM)):
        b = b * (1 + RATE) - a
    assert abs(b) < 1e-12


def test_linear_paths_match_exact_dynamics_for_small_shocks():
    T = 30
    eps = 1e-5
    di = eps * 0.8 ** np.arange(T)
    dpi = eps * 0.5 * 0.7 ** np.arange(T)
    for contract in ("V", "F"):
        lin = MG.unit_paths(di, dpi, contract, RATE, TERM, PI)
        ex = MG.exact_paths(di, dpi, contract, RATE, TERM, PI)
        pay, bal = ex["payment"], ex["balance"]
        assert np.max(np.abs(lin["payment"] - pay)) < 1e-3 * eps
        assert np.max(np.abs(lin["balance"] - bal)) < 1e-3 * eps


def test_variable_payment_rises_by_annuity_slope_on_impact():
    di = np.r_[0.01, np.zeros(19)]
    p = MG.unit_paths(di, np.zeros(20), "V", RATE, TERM, PI)["payment"]
    assert abs(p[0] - MG.annuity_slope(RATE, TERM) * 0.01) < 1e-12
    assert 0.55 < p[0] / 0.01 < 0.65


def test_fixed_rate_payment_unaffected_on_impact_by_rates():
    di = 0.01 * np.ones(20)
    p = MG.unit_paths(di, np.zeros(20), "F", RATE, TERM, PI)["payment"]
    phi = MG.origination_share(RATE, TERM, PI)
    assert abs(p[0] - MG.annuity_slope(RATE, TERM) * phi * 0.01) < 1e-12


def test_inflation_erodes_real_payments_of_both_contracts():
    dpi = np.r_[0.01, np.zeros(19)]
    for c in ("V", "F"):
        p = MG.unit_paths(np.zeros(20), dpi, c, RATE, TERM, PI)["payment"]
        assert p[0] < 0 and np.all(p <= 1e-15)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("passed:", name)

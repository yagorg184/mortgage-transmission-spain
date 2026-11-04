"""Properties of the equilibrium code that must hold by construction."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import equilibrium as EQ  # noqa: E402
import experiments as X  # noqa: E402
import paths as PT  # noqa: E402


def test_expectation_revisions_sum_to_path():
    x = np.array([0.0158, 0.0436, 0.0376, 0.0271, 0.019, 0.0133, 0.0, 0.0])
    T = len(x)
    total = np.zeros(T)
    for s, rev in PT.expectation_revisions(x, 0.7):
        total[s:] += rev[:T - s]
    assert np.max(np.abs(total - x)) < 1e-15


def test_response_is_affine_in_variable_share_and_sums_over_channels():
    econ, par, _ = X.build()
    p = PT.stylized_shock(econ.T, deposit_beta=0.12)
    r = {v: EQ.respond(econ, par, p, v) for v in (0.0, 0.5, 1.0)}
    mid = 0.5 * (r[0.0]["total"] + r[1.0]["total"])
    assert np.max(np.abs(r[0.5]["total"] - mid)) < 1e-12
    whole = EQ.respond(econ, par, p, 0.7, by_channel=False)["total"]
    parts = EQ.respond(econ, par, p, 0.7)["total"]
    assert np.max(np.abs(whole - parts)) < 1e-12


def test_general_equilibrium_does_not_depend_on_truncation():
    out = []
    for T in (200, 400):
        econ, par, _ = X.build(T=T)
        out.append(EQ.respond(econ, par, PT.stylized_shock(T, deposit_beta=0.12), 0.7)["total"][:5])
    assert np.max(np.abs(out[0] - out[1])) < 1e-4


def test_no_bank_profit_change_without_rate_change():
    econ, par, _ = X.build()
    T = econ.T
    p = dict(di=np.zeros(T), ddep=np.zeros(T), dpi=np.zeros(T))
    d = EQ.direct_effects(econ, par, p, 0.7)
    assert all(np.max(np.abs(v)) < 1e-15 for v in d.values())


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("passed:", name)

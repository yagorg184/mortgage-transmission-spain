"""The current code, with every new feature switched off, reproduces the previous
version of the model (git tag v2) to machine precision."""

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import economy as EC  # noqa: E402
import equilibrium as EQ  # noqa: E402

REF = json.load(open(os.path.join(HERE, "v2_reference.json")))


def _v2_economy():
    mu, M = REF["mortgage_share"], REF["balance"]
    types = [dict(beta=REF["beta_modal"], mass=1.0)]
    groups = [dict(name="none", mass=1 - mu, balance=0.0, payment=0.0),
              dict(name="mortgage", mass=mu, balance=M, payment=REF["payment_factor"] * M)]
    # The first version's illiquid grid leaves 2e-5 of mass on its last node, so
    # the distribution check is relaxed here and only here.
    return EC.Economy(dict(REF["household"]), types, groups, REF["T"], use_cache=False, check=False)


def _v2_par():
    return dict(payment_rule="level", redistribution="s_transfer", s=REF["s"], closure="wage_phillips",
                banks=dict(policy_spread=0.0), union=dict(import_share=0.0),
                ge=dict(kappa_w=0.1, beta_w=0.96, frisch=1.0))


def test_reproduces_v2():
    econ = _v2_economy()
    par = _v2_par()
    paths = EQ.stylized_paths(REF["T"], rho=REF["rho"], deposit_beta=1.0)
    out = REF["outputs"]
    cases = [(REF["variable_share"], out["impact_eff_share"]), (0.70, out["impact_spain_070"]),
             (0.12, out["impact_germany_012"]), (0.0, out["impact_fixed"])]
    for v, target in cases:
        got = EQ.respond(econ, par, paths, v, by_channel=False)["total"][0]
        assert abs(got - target) < 1e-10, (v, got, target)


if __name__ == "__main__":
    test_reproduces_v2()
    print("nesting test passed")

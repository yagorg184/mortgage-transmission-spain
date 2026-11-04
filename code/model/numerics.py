# numerics.py
# Numerical checks reported in Appendix B, written to output/numerics.json:
#   - the response to the stylized shock computed with horizons of 200 and 400 years;
#   - the closed-economy case (no import leakage): the present value of the first
#     column of the income Jacobian, the conditioning of the equilibrium system and
#     the first-year response, at both horizons;
#   - the largest mass on the last illiquid-asset node across household blocks.
# Run:  python numerics.py   (builds the T = 400 economy, a few minutes)

import copy
import json
import os

import numpy as np

import experiments as X
import paths as PT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
HORIZONS = (200, 400)


def main():
    os.makedirs(OUT, exist_ok=True)
    out = {"horizons": list(HORIZONS)}
    first5 = {}
    for T in HORIZONS:
        econ, par, cal = X.build(T=T)
        v = X.variable_share_spain()
        shock = PT.stylized_shock(T)
        first5[T] = X.aggregate(econ, par, shock, v)["total"][:5]
        Jy = econ.J_agg("C", "w")
        disc = (1 + par["banks"]["discount_rate"]) ** -np.arange(T)
        closed = copy.deepcopy(par)
        closed["union"]["import_share"] = 0.0
        m = par["union"]["import_share"]
        out[f"T{T}"] = {
            "pv_first_column_income_jacobian": round(float(disc @ Jy[:, 0]), 3),
            "condition_number_closed": round(float(np.linalg.cond(np.eye(T) - Jy)), 1),
            "condition_number_union": round(float(np.linalg.cond(np.eye(T) - (1 - m) * Jy)), 2),
            "first_year_response_closed": round(float(X.aggregate(econ, closed, shock, v)["total"][0]), 4),
            "first_year_response_union": round(float(first5[T][0]), 4),
        }
    gap = np.max(np.abs(first5[HORIZONS[0]] - first5[HORIZONS[1]]))
    out["max_abs_change_first5_years_when_horizon_doubles"] = float(f"{gap:.2e}")
    checks = cal["household"]["distribution_checks"]
    out["max_mass_last_illiquid_node"] = float(max(c["top_a"] for c in checks.values()))
    out["max_mass_last_liquid_node"] = float(max(c["top_b"] for c in checks.values()))
    out["min_mass"] = float(min(c["min"] for c in checks.values()))
    with open(os.path.join(OUT, "numerics.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

# results.py
# Reproduce the paper's headline numbers from the frozen calibration, reading
# every parameter from calibration.json. Run:  python results.py

import numpy as np
import config as C
import mortgage as MG
import euribor as EU
import ge as GE


def main():
    print("=" * 66)
    print("Headline results at the frozen calibration")
    print(f"  beta = {C.BETAS}   mortgage rate = {C.RATE:.4f}   exposure MV = {C.MV}")
    print("=" * 66)

    # Micro validation: variable-rate group vs the BdE -2.7% anchor
    r = MG.dynamic_response(dbp=400, rho=0.85)
    peak = 100 * r[np.argmax(np.abs(r))]
    yrs, re = EU.var_response()
    print("\nMicro (out-of-sample vs BdE -2.7%):")
    print(f"  +400 bp persistent (AR):        impact {100 * r[0]:+.2f}%")
    print(f"  observed Euribor path, 2024:    {100 * re[2]:+.2f}%   (peak {peak:+.2f}%)")

    # Amplification object: extensive margin, not fragility
    m = MG.cashflow_object(dbp=100)
    print("\nCash-flow decomposition (+100 bp on the variable rate):")
    print(f"  E[MPC] = {m['EMPC']:.3f}   Cov(MPC, exposure) = {m['cov']:+.4f}  (~0 -> extensive margin)")

    # General equilibrium and the counterfactual
    irf = GE.consumption_irf()
    cf = GE.counterfactual()
    print("\nGeneral equilibrium (+100 bp):")
    print(f"  derived redistribution spend = {cf['spend']:.3f}")
    print(f"  aggregate consumption impact = {irf[0]:+.4f}%")
    print(f"  Spain ({C.VSHARE['Spain']}) {cf['es']:+.4f}%  vs Germany ({C.VSHARE['Germany']}) {cf['de']:+.4f}%")
    print(f"  ratio ES/DE = {cf['ratio']:.3f}x   ({cf['frac']:.0f}% of the response attributable to ARM structure)")


if __name__ == "__main__":
    main()

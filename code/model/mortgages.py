# mortgages.py
# Nominal, long-term, level-payment mortgages (French amortisation). A mortgage
# group is described by its steady-state balance relative to mean income and by
# its contract: variable rate (the rate resets every year with the Euribor and the
# payment is recomputed over the remaining term) or fixed rate (the rate of the
# outstanding stock moves only as new loans replace amortised balances).
#
# All paths are first-order deviations from the steady state, per unit of the
# steady-state balance, in real terms (deflated by the price level). The aggregate
# stock follows
#     s_{t+1} (1 + pi_{t+1}) = s_t (1 + i_t - a_t) + o,
# where s is the real balance, i the contract rate, a the annuity factor, pi
# inflation and o real new lending, which is held at its steady-state value.

import numpy as np


def annuity(r, n):
    """Payment per unit of outstanding balance of a level-payment loan at rate r
    with n years to maturity."""
    return r / (1 - (1 + r) ** (-n))


def annuity_slope(r, n):
    """Derivative of the annuity factor with respect to the rate."""
    q = (1 + r) ** (-n)
    return ((1 - q) - r * n * q / (1 + r)) / (1 - q) ** 2


def origination_share(rate, term, infl_ss):
    """Real new lending per unit of the stock that keeps the real stock constant
    in the steady state: o/s = pi + a - i."""
    return infl_ss + annuity(rate, term) - rate


def unit_paths(di, dpi, contract, rate, term, infl_ss):
    """Deviations of the contract rate, the real balance, the real payment and
    the real interest accrued, per unit of the steady-state balance.

    di  : path of the Euribor (decimal deviation from steady state)
    dpi : path of inflation (decimal deviation from steady state)
    contract : "V" (rate resets every year) or "F" (only new lending reprices)
    """
    di = np.asarray(di, dtype=float)
    dpi = np.asarray(dpi, dtype=float)
    T = len(di)
    a = annuity(rate, term)
    ar = annuity_slope(rate, term)
    if contract == "V":
        drate = di.copy()
    elif contract == "F":
        phi = origination_share(rate, term, infl_ss)
        drate = np.zeros(T)
        prev = 0.0
        for t in range(T):
            prev = (1 - phi) * prev + phi * di[t]
            drate[t] = prev
    else:
        raise ValueError(f"unknown contract {contract!r}")
    da = ar * drate
    dx = np.zeros(T)
    dx[0] = -dpi[0] / (1 + infl_ss)
    for t in range(T - 1):
        dx[t + 1] = ((1 + rate - a) * dx[t] + drate[t] - da[t] - dpi[t + 1]) / (1 + infl_ss)
    return dict(rate=drate, balance=dx, payment=a * dx + da, interest=rate * dx + drate)


def exact_paths(di, dpi, contract, rate, term, infl_ss):
    """Nonlinear counterpart of unit_paths: the same stock dynamics and contract
    rules without linearisation. Returns the deviations of the real payment and
    the real balance from the steady state, per unit of the steady-state balance."""
    di = np.asarray(di, dtype=float)
    dpi = np.asarray(dpi, dtype=float)
    T = len(di)
    a0 = annuity(rate, term)
    o = origination_share(rate, term, infl_ss)
    rate_f = rate
    s = (1 + infl_ss) / (1 + infl_ss + dpi[0])
    pay, bal = np.zeros(T), np.zeros(T)
    for t in range(T):
        if contract == "V":
            r = rate + di[t]
        elif contract == "F":
            rate_f = (1 - o) * rate_f + o * (rate + di[t])
            r = rate_f
        else:
            raise ValueError(f"unknown contract {contract!r}")
        a = annuity(r, term)
        pay[t] = a * s - a0
        bal[t] = s - 1.0
        if t + 1 < T:
            s = (s * (1 + r - a) + o) / (1 + infl_ss + dpi[t + 1])
    return dict(payment=pay, balance=bal)


def level_payment_paths(di, contract):
    """Reduced-form mapping used in the first version of the model: the payment of
    a variable-rate borrower rises one for one with the rate on the whole balance
    and nothing else moves. Kept only to reproduce earlier results in the tests."""
    di = np.asarray(di, dtype=float)
    z = np.zeros_like(di)
    if contract == "V":
        return dict(rate=di.copy(), balance=z, payment=di.copy(), interest=di.copy())
    return dict(rate=z, balance=z, payment=z, interest=z)

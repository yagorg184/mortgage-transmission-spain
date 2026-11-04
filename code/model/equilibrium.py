# equilibrium.py
# Transmission of exogenous paths of the Euribor, the deposit rate and inflation
# to household consumption in a small member of a monetary union. Spain takes
# the ECB's policy rate and euro-area inflation as given; domestic output is
# demand-determined and part of household spending is imported. The model is
# linear around its steady state, so the response is the sum of the responses to
# each channel, and the response to any contract mix v is affine in v.

import numpy as np

import mortgages as MG

CHANNELS = ("liquid_rate", "liquid_inflation", "repricing_variable", "repricing_fixed",
            "mortgage_inflation", "bank_mortgage", "bank_deposits")


def wage_operator(T, kappa_w, beta_w):
    """Linear map from the path of the marginal rate of substitution to the path of
    the real wage implied by a wage Phillips curve."""
    A = (1 + kappa_w + beta_w) * np.eye(T) - np.eye(T, k=-1) - beta_w * np.eye(T, k=1)
    return kappa_w * np.linalg.inv(A)


def mortgage_paths(par, paths, contract, group):
    """Per-unit-balance paths of one contract held by one mortgage group, split
    into the part driven by the Euribor and the part driven by inflation."""
    T = len(paths["di"])
    zero = np.zeros(T)
    if par.get("payment_rule", "annuity") == "level":
        return MG.level_payment_paths(paths["di"], contract), \
            MG.level_payment_paths(zero, contract)
    args = (group["rate"], group["term"], par["nominal"]["inflation_ss"])
    return (MG.unit_paths(paths["di"], zero, contract, *args),
            MG.unit_paths(zero, paths["dpi"], contract, *args))


def direct_effects(econ, par, paths, v):
    """Impact of the exogenous paths on aggregate consumption before any change in
    labour income, by channel, in levels of per-capita consumption."""
    T = econ.T
    di = np.asarray(paths["di"], float)
    dpi = np.asarray(paths["dpi"], float)
    ddep = np.asarray(paths["ddep"], float)
    out = {k: np.zeros(T) for k in CHANNELS}

    Jrb = econ.J_agg("C", "rb")
    out["liquid_rate"] = Jrb @ ddep
    out["liquid_inflation"] = Jrb @ (-dpi)

    kappa = par["banks"]["policy_spread"]
    dprofit_mortgage = np.zeros(T)
    for g, grp in enumerate(econ.groups):
        if grp["balance"] == 0:
            continue
        Jp = econ.J_group("C", "payment", g)
        for contract, share in (("V", v), ("F", 1 - v)):
            mass = grp["mass"] * share
            if mass == 0:
                continue
            p_rate, p_infl = mortgage_paths(par, paths, contract, grp)
            key = "repricing_variable" if contract == "V" else "repricing_fixed"
            out[key] += mass * (Jp @ (grp["balance"] * p_rate["payment"]))
            out["mortgage_inflation"] += mass * (Jp @ (grp["balance"] * p_infl["payment"]))
            total_rate = p_rate["rate"] + p_infl["rate"]
            total_balance = p_rate["balance"] + p_infl["balance"]
            dprofit_mortgage += mass * grp["balance"] * (total_rate - di + kappa * total_balance)

    mode = par.get("redistribution", "banks")
    if mode == "banks":
        # Changes in bank profits reach domestic shareholders after tax. News about
        # future profits is capitalised in the value of bank equity when it
        # arrives, so households receive the present value of the change as a
        # one-off revaluation of illiquid wealth; future returns are unchanged.
        keep = (1 - par["banks"]["corporate_tax"]) * (1 - par["banks"]["foreign_share"])
        disc = (1 + par["banks"]["discount_rate"]) ** -np.arange(T)
        Jra = econ.J_agg("C", "ra")
        for key, dprofit in (("bank_mortgage", dprofit_mortgage), ("bank_deposits", econ.ss["B"] * (di - ddep))):
            dra = np.zeros(T)
            dra[0] = keep * float(disc @ dprofit) / econ.ss["A"]
            out[key] = Jra @ dra
    elif mode == "s_transfer":
        # Reduced-form rule of the first version: a share s of the extra interest
        # paid by variable-rate borrowers returns as a lump-sum transfer.
        Jp_all = sum(grp["mass"] * econ.J_group("C", "payment", g) for g, grp in enumerate(econ.groups))
        windfall = sum(grp["mass"] * v * grp["balance"] for grp in econ.groups) * di
        out["bank_mortgage"] = -Jp_all @ (par["s"] * windfall)
    elif mode != "none":
        raise ValueError(f"unknown redistribution mode {mode!r}")
    return out


def general_equilibrium(econ, par, dC_dir):
    """Consumption path that solves the demand-income fixed point.

    Baseline ("union"): output is demand-determined, a share m of consumption is
    imported, and household income moves one for one with domestic output in
    proportion to labour productivity (firms' revenues accrue to households with
    uniform incidence):
        dC = dC_dir + J_w dY,   dY = (1 - m) dC.
    Leakage abroad makes this fixed point well defined (the response does not
    depend on the truncation horizon), unlike the closed-economy version in which
    the present value of the columns of J_w equals one.

    "wage_phillips": the closure of the first version of the model, kept only to
    reproduce it in the tests."""
    T = econ.T
    Jw = econ.J_agg("C", "w")
    m = par["union"]["import_share"]
    if par.get("closure", "union") == "union":
        A = np.eye(T) - (1 - m) * Jw
    elif par["closure"] == "wage_phillips":
        g = par["ge"]
        W = wage_operator(T, g["kappa_w"], g["beta_w"])
        A = np.eye(T) - (1 - m) * Jw @ (np.eye(T) + W / g["frisch"])
    else:
        raise ValueError(f"unknown closure {par['closure']!r}")
    return np.linalg.solve(A, dC_dir)


def respond(econ, par, paths, v, by_channel=True):
    """Consumption response in percent of steady-state consumption: the total and,
    because the model is linear, the contribution of each channel."""
    direct = direct_effects(econ, par, paths, v)
    scale = 100 / econ.ss["C"]
    res = {"direct": {k: scale * x for k, x in direct.items()}}
    if by_channel:
        res["ge"] = {k: scale * general_equilibrium(econ, par, x) for k, x in direct.items()}
        res["total"] = sum(res["ge"].values())
    else:
        res["total"] = scale * general_equilibrium(econ, par, sum(direct.values()))
    return res


def respond_learning(econ, par, paths, v, rho_e):
    """Response when households learn the paths year by year: each year they
    observe the current deviations and expect them to decay at rate rho_e. The
    realised path is the sum of the revisions of these forecasts; by linearity and
    time invariance, the response to a revision announced in year s is the
    response to the same path announced in year 0, shifted by s."""
    import paths as PT
    T = econ.T
    revs = {k: dict(PT.expectation_revisions(np.asarray(paths[k], float), rho_e)) for k in ("di", "ddep", "dpi")}
    dates = sorted(set().union(*[r.keys() for r in revs.values()]))
    total = np.zeros(T)
    ge = {k: np.zeros(T) for k in CHANNELS}
    for s in dates:
        p = {k: np.zeros(T) for k in ("di", "ddep", "dpi")}
        for k in p:
            if s in revs[k]:
                p[k][:T - s] = revs[k][s][:T - s]
        r = respond(econ, par, p, v)
        total[s:] += r["total"][:T - s]
        for k in CHANNELS:
            ge[k][s:] += r["ge"][k][:T - s]
    return {"total": total, "ge": ge}


def stylized_paths(T, size=0.01, rho=0.6, deposit_beta=1.0):
    """A persistent increase in the Euribor with no inflation response."""
    di = size * rho ** np.arange(T)
    return dict(di=di, ddep=deposit_beta * di, dpi=np.zeros(T))

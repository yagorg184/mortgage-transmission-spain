# config.py
# Load the frozen calibration (calibration/calibration.json) and expose it to the
# model modules. This is the single source of truth: no parameter is hard-coded
# anywhere else.

import os
import json

_HERE = os.path.dirname(os.path.abspath(__file__))
_PATH = os.path.normpath(os.path.join(_HERE, "..", "..", "calibration", "calibration.json"))

with open(_PATH) as _f:
    CAL = json.load(_f)

# Household
_h = CAL["household"]
BETAS   = _h["beta"]
PIS     = _h["pi"]
CHI1    = _h["chi1"]
SIG_EPS = _h["sig_eps"]
RA      = _h["ra"]

# Mortgage block (frozen from the EFF, plus model amortisation)
_m = CAL["mortgage"]
MORT_SHARE = _m["mortgage_share"]          # share of households with dwelling debt
PVAR       = _m["variable_share"]          # variable-rate share of the mortgaged stock
MV         = _m["M_var_income_exposed"]    # variable balance / income among exposed
RATE       = _m["rate_level"]              # steady-state mortgage rate (measured, ~2%)
AMORT      = _m["amort"]                   # amortisation rate

# Interest-incidence inputs (externally disciplined)
_i = CAL["incidence"]
BETA_DEP   = _i["beta_dep"]
MPC_SHARE  = _i["mpc_shareholder"]
LAMBDA_EXT = _i["lambda_ext"]

# Variable-rate share of the mortgage stock by country (~2022), for the counterfactual
VSHARE = {k: v for k, v in CAL["country_variable_share_stock_2022"].items()
          if not k.startswith("_")}

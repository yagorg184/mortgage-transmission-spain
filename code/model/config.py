# config.py
# Single entry point to the model's parameters. Everything is read from
# calibration/: calibration.json (household parameters produced by calibrate.py
# and external values with their sources) and the two EFF tables produced by the
# survey pipeline (eff_targets.csv, eff_exposure.csv).

import json
import os

import pandas as pd

import mortgages as MG

HERE = os.path.dirname(os.path.abspath(__file__))
CAL_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "calibration"))
JSON = os.path.join(CAL_DIR, "calibration.json")
ALTERNATIVES = os.path.join(CAL_DIR, "alternatives")


def eff_table(name):
    t = pd.read_csv(os.path.join(CAL_DIR, name)).drop_duplicates("key").set_index("key")
    return t["value"].astype(float)


EFF = eff_table("eff_targets.csv")
EXPOSURE = eff_table("eff_exposure.csv")


def mortgage_groups(exposure=EXPOSURE, eff=EFF):
    """Mortgagors split into terciles of dwelling debt relative to income, each
    with the median balance, remaining term and rate of its tercile."""
    share = float(eff["mortgage_share"])
    groups = [dict(name="none", mass=1 - share, balance=0.0, payment=0.0, rate=0.0, term=0.0)]
    masses = [float(exposure[f"t{k}_mass_among_mortgagors"]) for k in (1, 2, 3)]
    tot = sum(masses)
    for k in (1, 2, 3):
        M = float(exposure[f"t{k}_tot_income_p50"])
        rate = float(exposure[f"t{k}_rate_p50"]) / 100
        term = float(exposure[f"t{k}_term_p50"])
        groups.append(dict(name=f"tercile{k}", mass=share * masses[k - 1] / tot, balance=M,
                           payment=MG.annuity(rate, term) * M, rate=rate, term=term))
    return groups


def household_population(x, fixed, grids, chi1=None, modal_beta_shift=None):
    """Household parameters, discount types and mortgage groups for a vector
    x = (beta_modal, beta_patient, patient_share, sig_eps). If given,
    `modal_beta_shift` maps mortgage groups to the amount by which the discount
    factor of their modal-type households differs from beta_modal."""
    beta_m, beta_p, pi_p, sig_eps = x
    hh = dict(fixed, chi1=chi1, sig_eps=sig_eps)
    types = [dict(name="modal", beta=beta_m, mass=1 - pi_p, grid=grids["modal"]),
             dict(name="patient", beta=beta_p, mass=pi_p, grid=grids["patient"])]
    groups = mortgage_groups()
    if modal_beta_shift:
        for g in groups:
            if g["name"] in modal_beta_shift:
                g["modal_beta_shift"] = float(modal_beta_shift[g["name"]])
    return hh, types, groups


def load(path=JSON):
    """Calibrated parameters and the population they define."""
    cal = json.load(open(path))
    h = cal["household"]
    x = [h["types"][0]["beta"], h["types"][1]["beta"], h["types"][1]["mass"], h["sig_eps"]]
    hh, types, groups = household_population(x, h["fixed"], h["grids"], chi1=h["chi1"],
                                              modal_beta_shift=h.get("modal_beta_shift"))
    return cal, hh, types, groups

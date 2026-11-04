# make_figures.py
# Figures of the paper, drawn from the outputs of the other scripts:
#   data/public/*.csv                          observed rates and prices
#   code/model/output/results.json             model experiments (results.py)
# Writes PDF files to paper/figures/. Run:  python make_figures.py

import json
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "code", "model"))
import paths as PT  # noqa: E402

OUT = os.path.join(ROOT, "paper", "figures")
RESULTS = os.path.join(ROOT, "code", "model", "output", "results.json")

# Categorical slots in fixed order; line styles and markers duplicate the colour
# coding so that the figures survive grey-scale printing.
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, MUTED = "#0b0b0b", "#52514e"
mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10, "legend.fontsize": 9,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.labelcolor": INK, "text.color": INK, "lines.linewidth": 1.6, "lines.markersize": 5,
    "legend.frameon": False, "savefig.bbox": "tight", "figure.dpi": 150,
})

# Calendar-time estimates of CaixaBank Research (Mestres Domenech et al., 2026):
# spending of variable- relative to fixed-rate mortgagors, percent.
CAIXABANK = [(2024 + 1.5 / 12, -2.7, "2024Q1"), (2025 + 11.5 / 12, -1.0, "Dec. 2025")]


def zero_line(ax):
    ax.axhline(0, color=MUTED, lw=0.6)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print("wrote", name)


def fig_rates():
    """Euribor, household deposit rate and inflation in Spain, 2019-2025."""
    eur = PT._monthly("euribor_1y")
    dep = PT.household_deposit_rate()
    hicp = PT._monthly("hicp_es")
    infl = 100 * (hicp / hicp.shift(12) - 1)
    s = slice(pd.Period("2019-01", "M"), pd.Period("2025-12", "M"))
    x = lambda ser: ser.loc[s].index.to_timestamp()
    fig, ax = plt.subplots(figsize=(6.3, 3.0))
    zero_line(ax)
    ax.plot(x(eur), eur.loc[s], color=C1, label="12-month Euribor")
    ax.plot(x(dep), dep.loc[s], color=C2, ls=(0, (5, 2)), label="Rate on household deposits")
    ax.plot(x(infl), infl.loc[s], color=C3, ls=(0, (1, 1.5)), lw=1.8, label="HICP inflation, Spain (12-month)")
    ax.set_ylabel("percent per year")
    ax.legend(loc="upper left")
    save(fig, "fig_rates.pdf")


def fig_stylized(res):
    """Response to the stylized tightening by contract structure, and its
    decomposition by channel on impact."""
    st = res["stylized"]
    years = np.arange(3)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 2.9), gridspec_kw=dict(width_ratios=[1, 1.15]))
    zero_line(a1)
    for key, lab, col, ls, mk in (("Spain", "Spanish mix", C1, "-", "o"),
                                  ("Germany", "German mix", C2, "--", "s"),
                                  ("fixed_only", "Fixed rates only", C3, ":", "^")):
        a1.plot(years, st["contract_line_first3"][key], color=col, ls=ls, marker=mk, label=lab)
    a1.set_xticks(years, ["1", "2", "3"])
    a1.set_xlabel("year after the shock")
    a1.set_ylabel("consumption, % deviation")
    a1.set_title("(a) Response by contract structure", loc="left")
    a1.set_ylim(-0.17, 0.025)
    a1.legend(loc="upper left", fontsize=8)

    ge = st["ge_impact_by_channel"]
    rows = [("Repricing, variable rate", ge["repricing_variable"]),
            ("Repricing, fixed rate", ge["repricing_fixed"]),
            ("Deposit rate", ge["liquid_rate"]),
            ("Bank profits, deposits", ge["bank_deposits"]),
            ("Bank profits, mortgages", ge["bank_mortgage"])]
    y = np.arange(len(rows))[::-1]
    vals = [v for _, v in rows]
    a2.barh(y, vals, height=0.6, color=[C1 if v < 0 else C3 for v in vals], edgecolor="white", lw=1)
    a2.axvline(0, color=MUTED, lw=0.6)
    a2.set_yticks(y, [n for n, _ in rows])
    for yi, v in zip(y, vals):
        a2.text(v + (0.004 if v >= 0 else -0.004), yi, f"{v:+.3f}", va="center",
                ha="left" if v >= 0 else "right", fontsize=8.5, color=INK)
    a2.set_xlim(-0.24, 0.07)
    a2.set_xlabel("contribution on impact, pp")
    a2.set_title("(b) Channels, Spanish mix, year 1", loc="left")
    a2.spines["left"].set_visible(False)
    a2.tick_params(axis="y", length=0)
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig_stylized.pdf")


def fig_contract_line(res):
    """Impact response as a function of the variable-rate share."""
    st = res["stylized"]
    line = st["contract_line_impact"]
    shares = {"fixed_only": 0.0, "Germany": 0.12, "euro_area": 0.25, "Italy": 0.45,
              "Spain": res["variable_share_spain"], "variable_only": 1.0}
    labels = {"Germany": ("Germany, France", (2, -13)), "euro_area": ("Euro area", (3, 6)),
              "Italy": ("Italy", (3, 6)), "Spain": ("Spain", (3, 6))}
    v = np.array([0.0, 1.0])
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    ax.plot(v, [line["fixed_only"], line["variable_only"]], color=C1, label="deposit pass-through 0.12")
    for beta, ls in ((0.32, "--"), (1.0, ":")):
        alt = st[f"contract_line_impact_depbeta_{beta}"]
        ax.plot(v, [alt["fixed_only"], alt["variable_only"]], color=MUTED, ls=ls, lw=1.1,
                label=f"deposit pass-through {beta:g}")
    for k, (lab, off) in labels.items():
        ax.plot(shares[k], line[k], "o", color=C1, ms=5, zorder=3)
        ax.annotate(lab, (shares[k], line[k]), xytext=off, textcoords="offset points", fontsize=8.5)
    zero_line(ax)
    ax.set_xlabel("variable-rate share of mortgagors")
    ax.set_ylabel("consumption on impact, % deviation")
    ax.set_xlim(-0.02, 1.02)
    ax.legend(loc="lower left")
    save(fig, "fig_contract_line.pdf")


def fig_episode(res):
    """2022-2025: aggregate response to the observed interest-rate paths by
    contract structure, and its decomposition by channel."""
    ep = res["episode"]["rates"]
    years = np.array(res["episode"]["years"])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.0))
    zero_line(a1)
    for key, lab, col, mk in (("aggregate_spain", "Spanish mix", C1, "o"),
                              ("aggregate_germany_mix", "German mix", C2, "s")):
        a1.plot(years, ep["perfect_foresight"][key], color=col, marker=mk, label=f"{lab}, perfect foresight")
        a1.plot(years, ep["learning"][key], color=col, marker=mk, ls="--", mfc="white",
                label=f"{lab}, learning")
    a1.set_xticks(years)
    a1.set_ylabel("consumption, % deviation")
    a1.set_title("(a) Aggregate consumption", loc="left")
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=8, ncol=1)
    a1.set_ylim(-1.3, 0.1)

    ch = ep["perfect_foresight"]["channels_spain"]
    groups = [("Repricing, variable rate", np.array(ch["repricing_variable"]), C1),
              ("Repricing, fixed rate", np.array(ch["repricing_fixed"]), C2),
              ("Deposit rate", np.array(ch["liquid_rate"]), C3),
              ("Bank profits", np.array(ch["bank_deposits"]) + np.array(ch["bank_mortgage"]), C4)]
    pos, neg = np.zeros(len(years)), np.zeros(len(years))
    for lab, vals, col in groups:
        base = np.where(vals >= 0, pos, neg)
        a2.bar(years, vals, bottom=base, width=0.6, color=col, edgecolor="white", lw=1, label=lab)
        pos += np.where(vals >= 0, vals, 0)
        neg += np.where(vals < 0, vals, 0)
    total = np.array(ep["perfect_foresight"]["aggregate_spain"])
    a2.plot(years, total, "D", color=INK, ms=4.5, label="Total")
    zero_line(a2)
    a2.set_xticks(years)
    a2.set_title("(b) Channels, Spanish mix, perfect foresight", loc="left")
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=8, ncol=2)
    a2.set_ylim(-1.3, 0.35)
    fig.tight_layout(w_pad=1.2)
    save(fig, "fig_episode.pdf")


def fig_borrowers(res):
    """Consumption of variable- relative to fixed-rate borrowers, 2022-2025."""
    years = np.array(res["episode"]["years"])
    ep = res["episode"]["all"]
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    zero_line(ax)
    ax.plot(years, ep["perfect_foresight"]["borrower_gap"], color=C1, marker="o", label="Model, perfect foresight")
    ax.plot(years, ep["learning"]["borrower_gap"], color=C1, marker="o", ls="--", mfc="white",
            label="Model, learning")
    for x, y, lab in CAIXABANK:
        ax.plot(x, y, "s", color=C2, ms=6, zorder=3)
        ax.annotate(lab, (x, y), xytext=(6, -2), textcoords="offset points", fontsize=8.5)
    ax.plot([], [], "s", color=C2, label="CaixaBank Research (2026)")
    ax.set_xticks(years)
    ax.set_xlim(years[0] - 0.3, years[-1] + 1.3)
    ax.set_ylim(-3.1, 0.2)
    ax.set_ylabel("consumption gap, %")
    ax.legend(loc="lower left")
    save(fig, "fig_borrowers.pdf")


def main():
    res = json.load(open(RESULTS))
    fig_rates()
    fig_stylized(res)
    fig_contract_line(res)
    fig_episode(res)
    fig_borrowers(res)


if __name__ == "__main__":
    main()

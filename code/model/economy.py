# economy.py
# The household sector: discount-factor types crossed with mortgage groups. Each
# (type, group) cell is one household block with its own steady state; mortgage
# groups differ only in the steady-state payment, and variable- and fixed-rate
# borrowers with the same balance share a steady state (they differ only in how
# their payment responds to shocks). Aggregates are mass-weighted sums of cells.

import hashlib
import json
import os
import pickle

import numpy as np

import household as H

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(HERE, "output", "cache")


def _cell_params(hh, typ, group):
    p = dict(hh)
    p.update(typ.get("grid", {}))
    p["beta"] = typ["beta"] + (group.get("modal_beta_shift", 0.0) if typ.get("name") == "modal" else 0.0)
    p["payment"] = group["payment"]
    return p


def solve_cells(hh, types, groups, check=True):
    """Steady states of every (type, group) cell."""
    cells = {}
    for j, typ in enumerate(types):
        for g, grp in enumerate(groups):
            ss, I = H.solve_ss(_cell_params(hh, typ, grp))
            ok, diag = H.check_distribution(I)
            if check and not ok:
                raise RuntimeError(f"invalid stationary distribution in cell "
                                   f"(type {j}, group {grp['name']}): {diag}")
            cells[(j, g)] = dict(ss=ss, I=I, diag=diag)
    return cells


def _topshare(x, w, top):
    """Weighted share of total x held by the top `top` fraction of households."""
    o = np.argsort(x)
    xs, ws = x[o], w[o]
    cw = np.cumsum(ws)
    k = min(int(np.searchsorted(cw, (1 - top) * cw[-1])), len(xs) - 1)
    sel = x >= xs[k]
    tot = float((w * x).sum())
    return float((w[sel] * x[sel]).sum() / tot)


def _spacing(grid):
    """Width of the interval each grid node represents (half-way to its neighbours)."""
    g = np.asarray(grid, float)
    mid = np.r_[g[0], 0.5 * (g[1:] + g[:-1]), g[-1]]
    return np.diff(mid)


class SmoothDistribution:
    """Distribution of a variable whose mass at each grid node is spread uniformly
    over the interval the node represents. Quantiles and threshold shares are then
    continuous in the model's parameters, which the calibration needs."""

    def __init__(self, value, width, weight):
        self.lo = value - width / 2
        self.hi = value + width / 2
        self.w = weight / weight.sum()

    def cdf(self, x):
        span = np.maximum(self.hi - self.lo, 1e-14)
        return float((self.w * np.clip((x - self.lo) / span, 0.0, 1.0)).sum())

    def quantile(self, p, tol=1e-12):
        a, b = float(self.lo.min()), float(self.hi.max())
        for _ in range(200):
            m = 0.5 * (a + b)
            if self.cdf(m) < p:
                a = m
            else:
                b = m
            if b - a < tol:
                break
        return 0.5 * (a + b)


def moments(cells, types, groups):
    """Cross-sectional moments of the pooled economy, in the units of the EFF:
    wealth relative to gross annual labour income, hand-to-mouth if liquid wealth
    is at most half a month of income, net worth = liquid + illiquid."""
    cols = {k: [] for k in ("lv", "lw", "av", "aw", "nv", "nw", "nl", "mpc", "wt")}
    for (j, g), cell in cells.items():
        I = cell["I"]
        D = np.asarray(I["D"])
        c = np.asarray(I["c"])
        b = np.asarray(I["b_grid"])
        a = np.asarray(I["a_grid"])
        e = np.asarray(I["e_grid"])
        ones = np.ones_like(D)
        inc = e[:, None, None] * ones
        Bm = b[None, :, None] * ones
        Am = a[None, None, :] * ones
        db = _spacing(b)[None, :, None] * ones
        da = _spacing(a)[None, None, :] * ones
        cols["lv"].append((Bm / inc).ravel())
        cols["lw"].append((db / inc).ravel())
        cols["av"].append((Am / inc).ravel())
        cols["aw"].append((da / inc).ravel())
        cols["nv"].append(((Bm + Am) / inc).ravel())
        cols["nw"].append(((db + da) / inc).ravel())
        cols["nl"].append((Bm + Am).ravel())
        cols["mpc"].append(np.clip(np.gradient(c, b, axis=1), 0, 1).ravel())
        cols["wt"].append((types[j]["mass"] * groups[g]["mass"] * D).ravel())
    X = {k: np.concatenate(v) for k, v in cols.items()}
    W = X["wt"] / X["wt"].sum()
    liq = SmoothDistribution(X["lv"], X["lw"], W)
    ill = SmoothDistribution(X["av"], X["aw"], W)
    nw = SmoothDistribution(X["nv"], X["nw"], W)
    thr = 1 / 24
    span = np.maximum(liq.hi - liq.lo, 1e-14)
    p_htm = np.clip((thr - liq.lo) / span, 0.0, 1.0)          # probability of being HtM, by node
    span_a = np.maximum(ill.hi - ill.lo, 1e-14)
    p_low_a = np.clip((thr - ill.lo) / span_a, 0.0, 1.0)
    return dict(
        liquid_income_p50=liq.quantile(0.5),
        illiquid_income_p50=ill.quantile(0.5),
        net_wealth_income_p50=nw.quantile(0.5),
        htm_share=float((W * p_htm).sum()),
        htm_poor_share=float((W * p_htm * p_low_a).sum()),
        htm_wealthy_share=float((W * p_htm * (1 - p_low_a)).sum()),
        net_wealth_top10_share=_topshare(X["nl"], W, 0.10),
        net_wealth_top1_share=_topshare(X["nl"], W, 0.01),
        mpc_annual=float((W * X["mpc"]).sum()),
    )


def moments_by_group(cells, types, groups):
    """Moments of each mortgage group, pooling the discount types."""
    out = {}
    for g, grp in enumerate(groups):
        sub = {(j, 0): cells[(j, g)] for j in range(len(types))}
        m = moments(sub, types, [dict(grp, mass=1.0)])
        out[grp["name"]] = dict(mass=grp["mass"], htm_share=m["htm_share"],
                                liquid_income_p50=m["liquid_income_p50"], mpc_annual=m["mpc_annual"])
    return out


CACHE_VERSION = 2


class Economy:
    """Steady states and Jacobians of the household sector."""

    def __init__(self, hh, types, groups, T, use_cache=True, check=True):
        self.hh, self.types, self.groups, self.T = dict(hh), types, groups, T
        key = hashlib.sha1(json.dumps([CACHE_VERSION, hh, types, groups, T], sort_keys=True,
                                      default=float).encode()).hexdigest()[:16]
        path = os.path.join(CACHE_DIR, f"economy_{key}.pkl")
        if use_cache and os.path.exists(path):
            with open(path, "rb") as f:
                stored = pickle.load(f)
        else:
            cells = solve_cells(hh, types, groups, check=check)
            stored = dict(J={}, agg={}, diag={}, moments=moments(cells, types, groups),
                          group_moments=moments_by_group(cells, types, groups))
            for k, cell in cells.items():
                J = H.jacobian(cell["ss"], T)
                stored["J"][k] = {o: {i: np.asarray(J[o][i]) for i in H.INPUTS} for o in ("C", "A", "B")}
                stored["agg"][k] = {o: float(cell["ss"][o]) for o in ("C", "A", "B")}
                stored["diag"][k] = cell["diag"]
            if use_cache:
                os.makedirs(CACHE_DIR, exist_ok=True)
                with open(path, "wb") as f:
                    pickle.dump(stored, f)
        self.J, self.cell_agg, self.diag = stored["J"], stored["agg"], stored["diag"]
        self.moments = stored["moments"]
        self.group_moments = stored["group_moments"]
        self.ss = {o: sum(self._w(k) * v[o] for k, v in self.cell_agg.items()) for o in ("C", "A", "B")}

    def _w(self, k):
        j, g = k
        return self.types[j]["mass"] * self.groups[g]["mass"]

    def group_index(self, name):
        return [g for g, grp in enumerate(self.groups) if grp["name"] == name][0]

    def J_agg(self, output, inp):
        """Per-capita aggregate Jacobian of `output` with respect to an input that
        moves for every household."""
        return sum(self._w(k) * J[output][inp] for k, J in self.J.items())

    def J_group(self, output, inp, g):
        """Jacobian of the per-capita `output` of mortgage group g (averaged over
        discount types) with respect to an input that moves only for that group."""
        return sum(self.types[j]["mass"] * self.J[(j, gg)][output][inp]
                   for (j, gg) in self.J if gg == g)

    def group_ss(self, output, g):
        return sum(self.types[j]["mass"] * self.cell_agg[(j, gg)][output]
                   for (j, gg) in self.cell_agg if gg == g)

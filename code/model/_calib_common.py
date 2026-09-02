
# _calib_common.py
# Shared machinery for the three calibration experiments (mpc-free, more-targets,
# fewer-params). Kept in one place so the three scripts stay short and identical
# where they should be; each experiment script imports from here and only defines
# its own free parameters, bounds, and target list.
#
# Nothing here touches the rest of the model: these scripts recalibrate and report,
# they do not write calibration.json.

import os, sys, json, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                      # so `import household` works from here
import household as H                          # noqa: E402

CAL_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "calibration"))
OUT_DIR = os.path.join(HERE, "output")
os.makedirs(OUT_DIR, exist_ok=True)

# ---- EFF targets (measured, with their csv keys) -------------------------------
_T_ALL = (pd.read_csv(os.path.join(CAL_DIR, "eff_targets.csv"))
          .drop_duplicates("key").set_index("key")["value"])

def eff(key):
    return float(_T_ALL[key])

# ---- model moments -------------------------------------------------------------
# One solve per discount type; aggregate against the stationary distribution.
# Returns every moment the three experiments might target or report, in model
# units matching the EFF definitions (shares in [0,1]; ratios in income units;
# the poor/wealthy-HtM split uses illiquid > income/24, exactly as eff_moments.R).

def _topshare(x, w, top):
    """Weighted share of total x held by the top `top` fraction (e.g. 0.10).
    Threshold at the (1-top) weighted quantile; matches w_topshare in the R."""
    o = np.argsort(x)
    xs = x[o]; ws = w[o]
    cw = np.cumsum(ws)
    q = xs[min(int(np.searchsorted(cw, (1 - top) * cw[-1])), len(xs) - 1)]
    sel = x >= q
    tot = float((w * x).sum())
    return float((w[sel] * x[sel]).sum() / tot) if tot > 0 else float("nan")


def solve_moments(betas, pis, chi1, sig_eps, ra):
  with np.errstate(all="ignore"):          # some parameter regions warn on power/log; guarded in the loss
    b = None
    LRl, ARl, NWl, NWLl, MPl, WTl = [], [], [], [], [], []
    for beta, pib in zip(betas, pis):
        ss, I = H.solve_type(beta, chi1, sig_eps, over={"ra": ra})
        D = np.array(I["D"]); c = np.array(I["c"])
        if b is None:
            b = np.array(I["b_grid"]); a = np.array(I["a_grid"]); e = np.array(I["e_grid"])
        inc = e[:, None, None] * np.ones_like(c)
        Bm = b[None, :, None] * np.ones_like(inc)
        Am = a[None, None, :] * np.ones_like(inc)
        mpc = np.clip(np.gradient(c, b, axis=1), 0, 1)
        LRl.append((Bm / inc).ravel()); ARl.append((Am / inc).ravel())
        NWl.append(((Bm + Am) / inc).ravel())
        NWLl.append((Bm + Am).ravel())
        MPl.append(mpc.ravel()); WTl.append((pib * D).ravel())
    LR = np.concatenate(LRl); AR = np.concatenate(ARl); NW = np.concatenate(NWl)
    NWL = np.concatenate(NWLl)
    MP = np.concatenate(MPl); WT = np.concatenate(WTl); WT = WT / WT.sum()

    def wq(x, p):
        o = np.argsort(x)
        return float(np.interp(p, np.cumsum(WT[o]), x[o]))

    htm = LR <= 1 / 24
    poor = htm & (AR <= 1 / 24)
    wealthy = htm & (AR > 1 / 24)
    return dict(
        liquid_income_p50=wq(LR, .5),
        liq_income_pctile_p10=wq(LR, .10), liq_income_pctile_p25=wq(LR, .25),
        liq_income_pctile_p75=wq(LR, .75), liq_income_pctile_p90=wq(LR, .90),
        illiquid_income_p50=wq(AR, .5),
        net_wealth_income_p50=wq(NW, .5),
        htm_share=float(WT[htm].sum()),
        htm_poor_share=float(WT[poor].sum()),
        htm_wealthy_share=float(WT[wealthy].sum()),
        net_wealth_top10_share=_topshare(NWL, WT, 0.10),
        net_wealth_top1_share=_topshare(NWL, WT, 0.01),
        mpc_annual=float((WT * MP).sum()),
    )

# ---- loss ----------------------------------------------------------------------
# Weighted sum of squared percentage deviations from the EFF targets. `targets`
# is a list of (csv_key, weight); `unpack` is a module-level function mapping a
# raw parameter vector to (betas, pis, chi1, sig_eps, ra), returning None for an
# invalid point. Implemented as a picklable callable so differential_evolution
# can run with workers>1 (multi-core).

class Objective:
    def __init__(self, targets, unpack):
        self.targets = targets
        self.unpack = unpack
        self.T = {k: eff(k) for k, _ in targets}
    PENALTY = 1e6

    def __call__(self, x):
        up = self.unpack(x)
        if up is None:
            return self.PENALTY
        try:
            m = solve_moments(*up)
        except Exception:
            return self.PENALTY
        # A divergent solve returns non-finite moments without raising. Penalise
        # that region so the optimiser avoids it instead of propagating a nan
        # (a nan loss silently poisons differential_evolution's incumbent).
        if not all(np.isfinite(m[k]) for k, _ in self.targets):
            return self.PENALTY
        L = sum(w * ((m[k] - self.T[k]) / self.T[k]) ** 2 for k, w in self.targets)
        return float(L) if np.isfinite(L) else self.PENALTY

def make_loss(targets, unpack):
    obj = Objective(targets, unpack)
    return obj, obj.T

# ---- checkpoint callback -------------------------------------------------------

def make_callback(tag, loss, unpack, targets, report_keys, every=15):
    state = {"gen": 0, "t0": time.time()}
    hist_path = os.path.join(OUT_DIR, f"calib_{tag}_history.jsonl")
    ckpt_path = os.path.join(OUT_DIR, f"calib_{tag}_checkpoint.json")
    open(hist_path, "w").close()          # fresh history each run (don't append to a stale one)

    def snapshot(xk):
        up = unpack(xk)
        m = solve_moments(*up)
        L = float(loss(xk))
        rec = {"tag": tag, "gen": state["gen"], "elapsed_min": round((time.time() - state["t0"]) / 60, 1),
               "loss": L, "params": [float(v) for v in xk],
               "moments": {k: round(m[k], 5) for k in report_keys}}
        with open(hist_path, "a") as f:
            f.write(json.dumps(rec) + "\n")
        with open(ckpt_path, "w") as f:
            json.dump(rec, f, indent=2)
        print(f"\n[gen {state['gen']:3d} | {rec['elapsed_min']:.1f} min] loss = {L:.5f}")
        print("   params:", "  ".join(f"{v:.5f}" for v in xk))
        Tt = {k: eff(k) for k, _ in targets}
        for k in report_keys:
            tag_k = "T" if k in Tt else " "
            tv = Tt.get(k, eff(k) if k in _T_ALL.index else float("nan"))
            dev = 100 * (m[k] - tv) / tv if tv == tv and tv != 0 else float("nan")
            print(f"     [{tag_k}] {k:<26} model {m[k]:8.4f}   eff {tv:8.4f}   {dev:+6.1f}%")

    def cb(*args, **kwargs):
        xk = args[0] if args else kwargs.get("xk")
        if hasattr(xk, "x"):
            xk = xk.x
        xk = np.asarray(xk, dtype=float)
        state["gen"] += 1
        if state["gen"] % every == 0:
            snapshot(xk)
        return False
    return cb

# ---- final report --------------------------------------------------------------

def final_report(tag, xk, loss, unpack, targets, report_keys, params_names):
    up = unpack(xk)
    m = solve_moments(*up)
    Tt = {k: eff(k) for k, _ in targets}
    print("\n" + "=" * 74)
    print(f"FINAL — {tag}    loss = {float(loss(xk)):.6f}")
    print("=" * 74)
    print("free parameters:")
    for name, v in zip(params_names, xk):
        print(f"   {name:<14} {v:.6f}")
    print("\nmoments  ([T]=targeted, [ ]=untargeted check):")
    for k in report_keys:
        is_t = k in Tt
        tv = Tt.get(k, eff(k) if k in _T_ALL.index else float("nan"))
        dev = 100 * (m[k] - tv) / tv if tv == tv and tv != 0 else float("nan")
        print(f"   [{'T' if is_t else ' '}] {k:<26} model {m[k]:9.4f}   eff {tv:9.4f}   {dev:+7.1f}%")
    out = {"tag": tag, "loss": float(loss(xk)),
           "params": dict(zip(params_names, [float(v) for v in xk])),
           "moments_model": {k: round(m[k], 5) for k in report_keys},
           "moments_eff": {k: (round(Tt[k], 5) if k in Tt else
                               (round(eff(k), 5) if k in _T_ALL.index else None)) for k in report_keys},
           "targeted": [k for k, _ in targets]}
    path = os.path.join(OUT_DIR, f"calib_{tag}_result.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {path}")
    return out

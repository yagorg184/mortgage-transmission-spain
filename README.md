# Mortgage rate structure and the transmission of ECB monetary policy to Spanish consumption

Replication package for the working paper *(title TBD)*.

**Question.** How much does Spain's near-universal **adjustable-rate mortgage (ARM)**
structure amplify the transmission of ECB monetary policy to household consumption,
relative to a fixed-rate (German-style) counterfactual?

**Method.** Stylized facts disciplined with the Spanish Survey of Household Finances
(EFF 2022) + a two-asset HANK model with discount-factor heterogeneity and a frozen
ARM block, solved in sequence space (SSJ), plus an own empirical macro validation
(recursive block-exogenous SVAR).

**Headline (all reproduced by execution).** Spain transmits **~1.19×** what it would
under the German fixed-rate structure; **~20%** of the Spanish consumption response to
a +100 bp ECB shock is attributable to the variable-rate structure. The mechanism is a
redistribution channel (Auclert 2019): the rate rise moves income from higher-MPC
borrowers to lower-MPC creditors (banks → shareholders/foreign investors). The aggregate
and distributional magnitudes, and their interpretation, are developed in the paper.

## Results

| Object | Value | Source |
|---|---|---|
| Household fit (loss) | 0.000004 | `calibrate.py` |
| Variable group, +400 bp (impact) | −2.8% | `results.py` |
| Variable group, observed Euribor path, 2024 | −2.4% | `results.py` |
| Aggregate consumption IRF, +100 bp | −0.37% | `results.py` |
| Counterfactual Spain / Germany | 1.19× (20% attributable) | `results.py` |
| SVAR consumption, +100 bp (impact) | −0.27% | `svar.py` |

The micro anchor (−2.4% vs the BdE's −2.7%) is an **out-of-sample** validation — that
number was never targeted. The SVAR (−0.27% impact, ≈−1% trough) confirms the GE
multiplier is not inflated.

## Reproduce

You need Python 3 with the packages in `requirements.txt`. Install them and run the
whole thing (works on Windows, macOS and Linux):

```
python -m pip install -r requirements.txt
python run_all.py
```

On Windows PowerShell, if `python`/`pip` are not found, use the launcher that ships with
the python.org installer, or the environment you normally use for Python:

```
py -m pip install -r requirements.txt
py run_all.py
```

`run_all.py` takes a few minutes and prints, in order: the household calibration
(`calibrate.py --refine`, ~1 min), the structural headlines (`results.py`:
−2.8% / −2.4% / −0.37% / 1.19×), and the SVAR (`svar.py`: −0.27%, with figures written
to `output/`). A Unix shell version, `run_all.sh`, is also provided.

### Full replication from scratch (optional)

The committed `calibration/eff_targets.csv` lets everything above run without the survey
microdata. To regenerate those moments from the raw EFF instead:

1. Request the EFF 2022 microdata from the Banco de España (see `DATA_AVAILABILITY.md`),
   unzip it under `code/empirical/eff_moments/data/raw/`, and install R with `data.table`.
2. Run the pipeline (one step — the cleaning is integrated; set `EFF_RAW_DIR` to read the
   raw files from another location):
   ```
   cd code/empirical/eff_moments
   Rscript eff_moments.R    # raw EFF -> targets, with replicate-weight + Rubin SEs
   ```
   Then copy `output/eff_targets.csv` over `calibration/eff_targets.csv`.
3. For the full SMM search (~1.5 h) instead of the fast refine:
   ```
   cd code/model
   python calibrate.py
   ```

## Layout

- `calibration/` — single source of truth: `calibration.json` (all frozen parameters)
  and `eff_targets.csv` (the EFF moments, each mapped to the object it disciplines).
- `code/model/` — `config.py` (loads the JSON), `household.py`, `mortgage.py`,
  `euribor.py`, `ge.py`, `results.py`, `calibrate.py`, and `_calib_common.py` (shared
  calibration helpers). No parameter is hard-coded outside the JSON.
- `code/model/exploration/` — calibration-design experiments (alternative target and
  parameter sets, and robustness checks) that informed the final calibration. Not part of
  the reproduction; kept for transparency.
- `code/empirical/eff_moments/` — the R pipeline: `eff_moments.R` runs raw EFF → all
  targets and facts (replicate-weight + Rubin SEs) in one step; `build_eff.R` provides the
  cleaning it calls, and can be run alone to cache `data/clean/`.
- `code/empirical/macro_svar/` — `svar.py`, the recursive SVAR.
- `docs/number_to_script.md` — every paper number → the script that produces it.
- `data/`, `output/` — see `DATA_AVAILABILITY.md`.

## Method notes and limitations

- The household calibration matches four EFF moments (median liquid/income, median
  illiquid/income, hand-to-mouth share, and the top-10% net-wealth share); the illiquid
  return and the patient discount factor are set externally, and the annual MPC is not
  targeted (reported alongside the fit). Household wealth is steep in the modal discount
  factor, so `beta_modal` is stored to several decimals and re-pinned by `calibrate.py --refine`.
- The steady-state mortgage payment uses the measured rate level (~2%); the ~1 pp margin
  over Euribor enters only the interest-incidence narrative.
- The German counterfactual uses the ECB MIR / ESRB stock variable-rate share (0.12).
- The SVAR shows a mild euro-area-GDP output puzzle (an information-effect residual) that
  does not affect the Spanish consumption response. The sign-restricted VAR and the
  1-SD local-projection exercise are not part of this package (a different normalisation).
- The transitory income component is calibrated to the MPC rather than measured; the debt
  wedge and loan-to-value blocks are not identified here. See the paper's limitations.

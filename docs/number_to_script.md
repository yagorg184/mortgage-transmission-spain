# Where each number comes from

Every figure and headline number in the paper, mapped to the script that produces
it and the file it lands in. Run `run_all.sh` to regenerate them.

## Empirical facts (EFF 2022)

| Number | Value | Script | Output |
|---|---|---|---|
| Liquid / income (median) | 0.26 | `code/empirical/eff_moments/eff_moments.R` | `calibration/eff_targets.csv` |
| Illiquid / income (median) | 3.56 | same | same |
| Hand-to-mouth share | 22% | same | same |
| Poor / wealthy HtM split | 9.1% / 12.9% | same | same |
| Variable-rate share (stock) | 70% | same | same |
| Exposure M_var / income (exposed) | 1.44 | same | same |
| Mortgage rate level (median) | 2.06% | same | same |
| Cov(MPC, exposure), partial | −0.03 (≈0) | same | same |
| Debt-service tail DSR>40% (pre / post +400 bp) | 4.4% / 8.4% | same | same |

## Calibration

| Number | Value | Script | Output |
|---|---|---|---|
| Household parameters (beta, pi, chi1, sig_eps, ra) | see file | `code/model/calibrate.py` | `calibration/calibration.json` |
| Fit (loss) | 0.0002 | same | same |

## Structural results (`code/model/results.py`)

| Number | Value | Script |
|---|---|---|
| Variable group, +400 bp AR (impact) | −2.8% | `mortgage.py` |
| Variable group, observed Euribor path, 2024 | −2.4% | `euribor.py` |
| Derived redistribution `spend` | 0.28 | `ge.py` |
| Aggregate consumption IRF, +100 bp (impact) | −0.38% | `ge.py` |
| Counterfactual ratio Spain / Germany | 1.27× | `ge.py` |
| Share attributable to the ARM structure | 26% | `ge.py` |

## Macro validation (`code/empirical/macro_svar/svar.py`)

| Number | Value | Output |
|---|---|---|
| Spanish consumption, SVAR, +100 bp (impact) | −0.27% | `svar_irf.csv` |
| Trough (within a year) | ≈ −1.1% | same |

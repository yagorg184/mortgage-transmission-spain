# Adjustable-Rate Mortgages and the Transmission of ECB Monetary Policy to Household Consumption in Spain

Code and data for the paper of the same title. About two-thirds of the mortgage debt of Spanish households
carries a rate indexed to the Euribor. The paper measures how much this contract structure adds to the
response of household consumption to ECB monetary policy. It uses a two-asset heterogeneous-agent model
calibrated to the Spanish Survey of Household Finances (EFF 2022), with:

- nominal long-term annuity mortgages for three groups of borrowers;
- banks that receive the extra interest;
- a small member of a currency union.

## Main results

| Result | Value | Output |
|---|---|---|
| First-year response of consumption to a persistent 100 bp rise in the Euribor: Spanish contract mix / German mix / fixed rates only | -0.158% / -0.097% / -0.084% | `results.json`: `stylized.contract_line_impact` |
| Contribution of the Spanish contracts, relative to fixed rates only | -0.074 pp (range -0.053 to -0.120 across shock persistence) | `robustness.json` |
| Effect of the Spanish contract mix relative to the German one, 2023-2025 | about -0.3 pp a year | `results.json`: `episode.rates.*.spain_minus_germany_mix` |
| Payment of the median borrower in 2023 | +3.0% of income (bank data: about +3 pp) | `results.json`: `episode.payment_change_median_borrower` |
| Consumption of variable- relative to fixed-rate borrowers, 2023-2024 | -1.2% to -1.4% (bank data: -2.1% to -2.7%) | `results.json`: `episode.all.*.borrower_gap` |

`docs/number_to_script.md` maps every number in the paper to the file and key that contains it.

## Reproducing the results

Python 3.10 or later with the packages in `requirements.txt`:

```
python -m pip install -r requirements.txt
python run_all.py
```

`run_all.py` runs, in order:

- the import content of consumption;
- the VAR specifications and the local projections;
- a check of the stored calibration;
- the model experiments (`results.py`), the sensitivity analysis (`robustness.py`) and the numerical checks (`numerics.py`);
- the figures.

It takes about ten minutes, most of it computing the Jacobians of the household blocks. These are cached in `code/model/output/cache/`. Options:

- `--calibrate` re-solves the baseline and the three alternative calibrations, about 30 minutes more;
- `--download` downloads the public data again;
- `--tests` runs the tests.

The outputs `code/model/output/results.json` and `robustness.json` are identical bit for bit across the two environments tested (Python 3.10 with numpy 2.2, Python 3.13 with numpy 2.5). The VAR outputs differ by less than 1e-9.

The survey moments in `calibration/` come from the restricted EFF microdata. `DATA_AVAILABILITY.md` explains how to request the microdata and regenerate them (R for `eff_moments.R`, Python for `eff_exposure.py`).

## Layout

| Path | Contents |
|---|---|
| `calibration/` | `eff_targets.csv` and `eff_exposure.csv` (EFF moments with standard errors), `calibration.json` (every model parameter), `alternatives/` (alternative calibrations used in the sensitivity analysis) |
| `code/empirical/eff_moments/` | EFF pipeline: `build_eff.R`, `eff_functions.R`, `eff_moments.R`, `eff_exposure.py` |
| `code/empirical/public_data/` | `download.py`: ECB, Eurostat and Jarociński-Karadi series |
| `code/empirical/io/` | `import_content.py`: import content of household consumption from input-output tables |
| `code/empirical/macro_svar/` | `svar.py` and `svar_annual.py`: VAR specifications of Appendix C |
| `code/empirical/lp/` | `lp.py`: local projections with high-frequency shocks (reported as uninformative) |
| `code/model/` | the model: `household.py` (two-asset household block), `mortgages.py` (annuity mortgages), `economy.py` (household population, steady states, Jacobians, moments), `equilibrium.py` (channels, banks, currency-union closure, expectations), `paths.py` (observed and identified paths), `config.py`, `calibrate.py`, `experiments.py`, `results.py`, `robustness.py`, `numerics.py`, `tests/` |
| `code/figures/` | `make_figures.py` |
| `data/` | public data (see `docs/sources.md`) |
| `docs/` | `sources.md` (every input and its source), `number_to_script.md` |

## Notes

- Household, mortgage, bank and currency-union parameters are read from `calibration/calibration.json`. The settings of the experiments (paths of the shocks, the deposit pass-through of the stylized shock and the persistence of expectations) are in `code/model/paths.py` and `code/model/experiments.py`.
- The test `tests/test_nesting.py` checks that the code reproduces the previous version of the model (git tag `v2`) to machine precision when its new features are switched off.
- The model is linear around its steady state and shocks are unanticipated paths. The paper discusses what it leaves out: house prices, refinancing, contract choice and income incidence.
- The code is released under the MIT License (`LICENSE`). The data in `data/` and the EFF moments in `calibration/` remain subject to the terms of their sources (`docs/sources.md`, `DATA_AVAILABILITY.md`).

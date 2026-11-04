# Where each number in the paper comes from

Every number in the paper is either an EFF moment (`calibration/eff_targets.csv`,
`calibration/eff_exposure.csv`), a calibrated parameter (`calibration/calibration.json`) or an output
of a script in this repository. `python run_all.py` regenerates all outputs. Keys below are paths in the
JSON files (`a.b.c`) or rows of the CSV files. In `results.json`, `Spain` is the variable-rate share of
Spanish mortgagors (0.704), `Germany` 0.12, `fixed_only` 0 and `variable_only` 1; lists of years start
in the first year after the shock, or in 2022 for the episode.

## Section 2 and Table 1: facts

| Number in the paper | Source |
|---|---|
| Mortgage share 0.277; variable-rate share of mortgagors 0.704 | `eff_targets.csv`: `mortgage_share`, `variable_share_among_mortgaged` |
| Variable-rate share of the balance 0.680 | `eff_exposure.csv`: `variable_share_of_stock` |
| Tercile cut points 1.01 and 2.09; tercile medians of debt/income, rate, term, DSR; variable-rate and hand-to-mouth shares by tercile; all mortgagors | `eff_exposure.csv`: `tercile_cut_*`, `t{1,2,3}_*`, `tot_income_p50_mortgagors`, `dsr_p50_mortgagors`, `htm_share_mortgagors`, `check_variable_share_among_mortgaged` |
| Liquid, illiquid and net worth medians; top shares; hand-to-mouth shares; H0 and H1; DSR above 40% | `eff_targets.csv`: `liquid_income_p50`, `illiquid_income_p50`, `net_wealth_income_p50`, `net_wealth_top{10,1}_share`, `htm_*`, `dsr_over40_*` |
| Euribor −0.49 (2021) and 3.87 (2023); deposit rate 0.02, 0.22, 0.48; inflation 8.0 and 3.3 | `results.json`: `facts.annual_percent` |
| Deposit pass-through 5% (2023) and 12% (2024) | `results.json`: `facts.deposit_pass_through_since_2021` |
| Overnight share of deposits 90% and 85%; overnight rate below 0.2%; new time deposits 2.4% | `results.json`: `facts.overnight_share_of_deposits`, `facts.overnight_rate_max_2022_2025`, `facts.new_time_deposit_rate_2024` |
| Figure 1 | `code/figures/make_figures.py`, `fig_rates` |

## Sections 3 and 4: model and calibration

| Number in the paper | Source |
|---|---|
| α_i ≈ 0.59 (0.591, 0.582, 0.592) | `results.json`: `groups.tercile*.annuity_slope` |
| New lending 6 to 16% of the stock (0.156, 0.077, 0.060) | `results.json`: `groups.tercile*.origination_share` |
| Table 2: calibrated household parameters | `calibration.json`: `household.types`, `household.chi1`, `household.sig_eps` |
| Mortgage spread 0.97 pp | `eff_targets.csv`: `spread_over_euribor1Y_2022` |
| Import content 0.252, of which 14% direct and 11% indirect | `results.json`: `facts.import_content_2021`; `code/empirical/io/import_content.py` |
| Table 3: targeted and untargeted moments | `calibration.json`: `household.achieved`, `household.untargeted`; `results.json`: `moments`, `groups.*.htm_share`, `groups.*.payment_income` |
| Top-1% share between 0.12 and 0.16 | `results.json`: `moments.net_wealth_top1_share`; `robustness.json`: `grids_dense.moments`, `calibration_bp*.moments` |
| MPC 0.247; MPC by group 0.24, 0.23, 0.22 and 0.25 | `results.json`: `moments.mpc_annual`, `groups.*.mpc_annual` |
| MPC out of payments 0.23, 0.39, 0.59 | `results.json`: `mpc_payment` |
| Alternative calibrations (β_p = 0.962, 0.965; liquidity by tercile) | `calibration/alternatives/*.json` (`calibrate.py --beta-patient`, `--terciles`) |

## Section 5: stylized tightening

| Number in the paper | Source |
|---|---|
| −0.158, −0.123, −0.101 | `results.json`: `stylized.total_first3` |
| Table 4 | `results.json`: `stylized.direct_impact_by_channel`, `stylized.ge_impact_by_channel`, `stylized.direct_impact_total` |
| Multiplier 1.7 | `results.json`: `stylized.ge_multiplier_impact` |
| Deposit margin offsets 15% of repricing; revaluation effect tripled | ratios of entries of `stylized.ge_impact_by_channel` and `stylized.direct_impact_by_channel` |
| Equation (12): −0.084 and −0.106 | `results.json`: `stylized.contract_line_impact.fixed_only`, and `variable_only` minus `fixed_only` |
| Mortgage-margin offset about a tenth of the slope | `results.json`: `stylized.ge_impact_by_channel_fixed_only.bank_mortgage` and `ge_impact_by_channel_variable_only.bank_mortgage` |
| Table 5 | `results.json`: `stylized.contract_line_first3` |
| Ratio 1.64; 0.062 and 0.074 pp; 47%; 1.7 in year 3; 1.67 and 1.51 for a German share of 0.10 or 0.20 | computed from `stylized.contract_line_impact` and `stylized.contract_line_first3` |
| Table 6 | `robustness.json` (one entry per row) |
| Figures 2 and 3 | `make_figures.py`, `fig_stylized`, `fig_contract_line` |

## Section 6: the 2022-2025 tightening

| Number in the paper | Source |
|---|---|
| Euribor and deposit-rate deviations from 2021 | `results.json`: `episode.paths` |
| Aggregate response, perfect foresight and learning | `results.json`: `episode.rates.{perfect_foresight,learning}.aggregate_spain` |
| Repricing about 0.8% a year; deposit margins about 0.2%; deposit-rate channel | `results.json`: `episode.rates.perfect_foresight.channels_spain` |
| Spanish minus German mix: 0.32, 0.33, 0.31, 0.27 and 0.11, 0.33, 0.34, 0.30 | `results.json`: `episode.rates.*.spain_minus_germany_mix` |
| Inflation above target: 6.0, 1.3, 0.8, 0.7 | `results.json`: `episode.paths.dpi` |
| Inflation channels: 0.4% a year on mortgages; 0.04% and 0.8% (perfect foresight) and 1.4% (learning) on deposits | `results.json`: `episode.inflation.*.channels_spain` |
| Payment of the median borrower 3.0 and 2.7 (3.7 and 3.4 without inflation) | `results.json`: `episode.payment_change_median_borrower`, `episode.payment_change_median_borrower_rates_only` |
| Gap between variable- and fixed-rate borrowers | `results.json`: `episode.all.*.borrower_gap` |
| Gap with liquidity matched by tercile, 1.3% in 2023 | `robustness.json`: `calibration_terciles.episode.perfect_foresight.borrower_gap` |
| Figures 4 and 5 | `make_figures.py`, `fig_episode`, `fig_borrowers` |

## Appendices

| Number in the paper | Source |
|---|---|
| Linearization errors 2.7%, 6.0%, 9.5% | `results.json`: `mortgage_linearization_error` |
| Grids, income states, horizon | `calibration.json`: `household.fixed`, `household.grids`; `results.json`: `T` |
| Mass on the last node below 1e-11; horizon doubling changes the response by less than 1e-4 | `numerics.json`: `max_mass_last_illiquid_node`, `max_abs_change_first5_years_when_horizon_doubles` |
| Closed economy: present value 0.98, condition numbers 119 and 276, responses −0.815 and −0.823 | `numerics.json`: `T200`, `T400` |
| Final calibration misses below 0.04% | `calibration.json`: `household.targets`, `household.achieved` |
| Table 7 (VAR specifications and the model) | `results.json`: `svar.*`; `code/empirical/macro_svar/output/svar_annual_*.csv` |
| Internal instrument: +2.93 in year 1 | `svar_annual_internal.csv` |
| First-stage F of the local projections, about 10 | `code/empirical/lp/output/first_stage.txt` |

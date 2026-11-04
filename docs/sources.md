# Sources

Where every input of the model and every external number in the paper comes from. Values measured in
the EFF are produced by the scripts in `code/empirical/eff_moments/` and stored in
`calibration/eff_targets.csv` and `calibration/eff_exposure.csv`; the series downloaded from public
sources are stored in `data/`.

## Household survey: EFF 2022 (Banco de España)

Restricted microdata (see `DATA_AVAILABILITY.md`). Five multiply imputed implicates and 1,000 bootstrap
replicate weights; standard errors combine both by Rubin's rules.

| Object | Value | Variable(s) | Script |
|---|---|---|---|
| Households with dwelling debt | 0.277 | loans 1-4 (`p2_12_*`) | `eff_moments.R` |
| Variable-rate share of mortgagors | 0.704 | rate type `p2_14_*` = 11 | `eff_moments.R` |
| Variable-rate share of the outstanding balance | 0.680 | idem, balance-weighted | `eff_exposure.py` |
| Terciles of dwelling debt / income: median balance, rate, remaining term, DSR, variable-rate share, hand-to-mouth share | Table 1, panel B | `p2_12_*`, `p2_13_1`, `p2_17_1`, `p2_18_*`, `renthog21_eur22` | `eff_exposure.py` |
| Median liquid and illiquid wealth / income, net worth / income | 0.263, 3.56, 4.18 | `p4_7_*`, `p8_5a`, `riquezanet`, `renthog21_eur22` | `eff_moments.R` |
| Hand-to-mouth share (liquid wealth at most half a month of income), poor and wealthy | 0.220, 0.091, 0.129 | idem | `eff_moments.R` |
| Top-10% and top-1% shares of net worth | 0.536, 0.211 | `riquezanet` | `eff_moments.R` |
| Hand-to-mouth share of variable-rate borrowers, by household and by balance | 0.202, 0.224 | idem | `eff_moments.R` |
| DSR above 40% before and after a +400 bp reset of the variable-rate balance | 0.044, 0.084 | idem | `eff_moments.R` |
| Mortgage spread over the Euribor | 0.97 pp | median rate of variable-rate loans minus the 2022 average 1-year Euribor (1.09%) | `eff_moments.R` |

## External parameters

| Parameter | Value | Source |
|---|---|---|
| Real return on illiquid assets, r^a | 0.035 | Banco de España, *Informe Anual 2023*, chapter 4: gross rental yield of about 5.5%, of which about two points are running costs and risk. |
| Persistent income, ρ_z and σ_η | 0.98, 0.08 | Calibrated: Spanish earnings close to a random walk (Pijoan-Mas and Sánchez-Marcos, 2010); annual persistence of 0.98-0.99 in Storesletten, Telmer and Yaron (2004) and Guvenen (2009); σ_η set for a cross-sectional variance of persistent log income of about 0.16. |
| Elasticity of intertemporal substitution | 0.5 | Convention. |
| Adjustment-cost curvature χ0, χ2 | 0.25, 2 | Kaplan, Moll and Violante (2018). |
| Steady-state inflation | 0.02 | ECB target. |
| Deposit pass-through, β^d | 0.12 | Rate on the stock of Spanish household deposits (overnight and time deposits, ECB MFI interest-rate statistics, weighted by BSI volumes), change 2021-2024 relative to the change in the 1-year Euribor. Computed in `code/model/paths.py`; reported in `results.json` (`facts`). |
| Corporate tax on banks, τ | 0.30 | Spanish corporate income tax rate for credit institutions (Ley 27/2014, art. 29.6). |
| Foreign share of bank equity, λ | 0.50 | Bolsas y Mercados Españoles, report on the ownership of listed Spanish shares (2023): non-residents held 50.3% of the market. Used for bank shareholders in the absence of a bank-specific figure. |
| Import content of household consumption, m | 0.2519 | Eurostat symmetric input-output table for Spain, 2021 (`naio_10_cp1700`, product by product): imported final consumption plus the imported inputs embodied in domestic consumption. `code/empirical/io/import_content.py`. |
| Variable-rate shares of other countries | DE/FR 0.12, EA 0.25, IT 0.45 | Ordinal placement only. Germany and France: between a tenth and a fifth in household-survey and bank data (Ehrmann and Ziegelmeyer, 2017; Albertazzi, Fringuellotti and Ongena, 2024). Euro area: ECB (2025), about a quarter of the stock at pure variable rates. Because the response is linear in the share, any other value can be read off Table 5. |

## Validation evidence

| Object | Value | Source |
|---|---|---|
| Spending of variable- relative to fixed-rate borrowers | -1.2 / -2.1 / -2.5% at 6 / 12 / 18 months after the first payment revision; -2.7% in 2024Q1; -1.0% in December 2025 | Mestres Domènech, García Arenas, Graziano and Alcobé Garcia, CaixaBank Research, 10 April 2026. Staggered difference-in-differences on CaixaBank customers with a single mortgage; spending measured from card payments, direct debits and cash withdrawals. |
| Change in the median debt-service ratio of variable-rate borrowers | close to +3 pp, March 2022 to July 2024 | idem |
| Spending response to the payment cuts of late 2025 | about three-quarters of the fall in payments (median payment -40 euros, spending +30 euros) | idem |
| Response of consumption to adjustable-rate payment increases | about half of the increase in interest payments for high-debt borrowers | Flodén, Kilström, Sigurdsson and Vestman (2021). |

## Time series

All stored in `data/`; `code/empirical/public_data/download.py` downloads the series in `data/public/`.

| Series | Source | File |
|---|---|---|
| 1-year Euribor, monthly | ECB Data Portal, `FM/M.U2.EUR.RT.MM.EURIBOR1YD_.HSTA` | `data/public/euribor_1y.csv` |
| HICP, Spain, euro area and Germany, monthly | ECB Data Portal, `ICP/M.{ES,U2,DE}.N.000000.4.INX` | `data/public/hicp_*.csv` |
| Rates on Spanish household deposits: overnight, time deposits (stock and new business) | ECB MFI interest-rate statistics, `MIR/M.ES.B.L21.A.R.A.2250.EUR.N`, `MIR/M.ES.B.L22.A.R.A.2250.EUR.O`, `MIR/M.ES.B.L22.F.R.A.2250.EUR.N` | `data/public/dep_rate_*.csv` |
| Volumes of Spanish household deposits | ECB BSI, `BSI/M.ES.N.A.L21.A.1.U2.2250.Z01.E`, `BSI/M.ES.N.A.L22.A.1.U2.2250.Z01.E` | `data/public/dep_vol_*.csv` |
| Household consumption and GDP, quarterly, chain-linked volumes, seasonally and calendar adjusted | Eurostat `namq_10_fcs`, `namq_10_gdp` | `data/public/cons_*.csv`, `gdp_*.csv`; `data/macro/*.xlsx` |
| Jarociński-Karadi monetary policy shocks, monthly (CC BY 4.0) | github.com/marekjarocinski/jkshocks_update_ecb | `data/public/jk_shocks_monthly.csv` |
| Input-output tables, Spain 2020-2022 | Eurostat `naio_10_cp1700` | `data/public/io/*.json` |
| Monthly HICP indices used in the VAR | Eurostat `prc_hicp_minr` | `data/macro/prc_hicp_minr_*.xlsx` |
| Brent crude oil price (USD) and the USD/EUR exchange rate | U.S. Energy Information Administration (`RBRTEm.xls`); ECB Data Portal | `data/macro/` |
| Euro-area shadow short rate | Wu and Xia (`shadowrate_ECB.xls`), spliced onto the 1-month Euribor (ECB Data Portal) | `data/macro/`, `data/euribor/` |

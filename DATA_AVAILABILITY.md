# Data Availability Statement

This repository contains **code only**. Data fall into two groups.

## 1. Restricted — NOT included in this repository
**Spanish Survey of Household Finances (EFF 2022), Banco de España.**
The EFF microdata are distributed by the Banco de España under terms that do
**not** permit redistribution. They are therefore excluded from this repo (see
`.gitignore`). To reproduce the household calibration and the stylized facts:

1. Request the EFF 2022 microdata from the Banco de España
   (https://www.bde.es/ → Statistics → Survey of Household Finances / EFF).
2. Place the files as described in `code/empirical/eff_moments/README.md`.
3. Run the R pipeline; it regenerates `calibration/calibration_targets*.csv`.

The **aggregate calibration targets** derived from the EFF (medians, shares,
Hand-to-Mouth rate, MPC target, percentiles) are publishable summary statistics
and **are** included in `calibration/` so the model can be run without the
microdata.

## 2. Public — obtainable from source (not committed to keep the repo light)
- **ECB EA-MPD** (Altavilla et al. 2019) — high-frequency monetary surprises.
- **ECB Data Portal / SDW** — Euribor 1Y, HICP (ES, EA), MIR mortgage rates, 2Y rate.
- **Eurostat** — national accounts (consumption, GDP), unemployment, industrial production, retail.
- **Wu-Xia / ECB shadow rate**; **oil price (€)**.

Exact series IDs and download instructions: `data/public/README.md`.

# Data availability

## Restricted: EFF 2022 microdata (not included)

The Spanish Survey of Household Finances (*Encuesta Financiera de las Familias*, EFF) 2022 is produced
by the Banco de España, which distributes the microdata to researchers under a user agreement that does
not allow redistribution. The microdata are therefore not in this repository, and `.gitignore` excludes
the raw and cleaned files by name wherever they are placed.

To reproduce the survey moments:

1. Request the EFF 2022 microdata from the Banco de España (www.bde.es, Statistics, Survey of Household
   Finances). The files used are `otras_secciones_2022_imp{1..5}.csv`, `databol{1..5}.csv` and
   `replicate_weights_2022.csv`.
2. Place them under `code/empirical/eff_moments/data/raw/`. The environment variable `EFF_RAW_DIR`
   points the scripts to another folder (`eff_moments.R` reads the replicate weights from `data/raw/` in
   any case).
3. Run, from `code/empirical/eff_moments/`:
   - `Rscript eff_moments.R` (R 4.3 or later with `data.table`), which writes `output/eff_targets.csv`;
   - `python eff_exposure.py`, which writes `output/eff_exposure.csv`.
4. Copy both files to `calibration/`.

Both scripts report standard errors that combine the 1,000 replicate weights with the variation across
the five implicates (Rubin's rules). Without the replicate weights they report point estimates only.

The survey moments themselves are aggregate statistics and are included in `calibration/`, so that the
model runs without the microdata.

## Public: included

All other data are public and included in `data/`:

- `data/public/`: series downloaded from the ECB Data Portal and Eurostat, the Jarociński-Karadi
  monetary policy shocks (CC BY 4.0) and the Eurostat input-output tables for Spain. Running
  `code/empirical/public_data/download.py` downloads them again.
- `data/macro/` and `data/euribor/`: the original downloads used by the VAR (Eurostat national accounts
  and HICP, Brent price, exchange rate, Wu-Xia shadow rate, 1-month Euribor).

`docs/sources.md` lists every series with its source and code.

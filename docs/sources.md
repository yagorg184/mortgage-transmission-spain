# Sources for the calibrated parameters

Every number the model uses, and where it comes from. Parameters measured from the
EFF 2022 are reproduced by the R pipeline (`calibration/eff_targets.csv`); the rest are
external and pinned to a citation here. Status: **confirmed** (checked against the
source), **cited** (source named, figure not re-checked line by line), **open** (needs
a better anchor).

## Measured from the EFF 2022 (own, reproducible)

| Parameter | Value | Source |
|---|---|---|
| Liquid / income (median) | 0.26 | EFF 2022 → `eff_moments.R` |
| Illiquid / income (median) | 3.56 | EFF 2022 |
| Hand-to-mouth share | 22% | EFF 2022 |
| Mortgage prevalence | 28% | EFF 2022 |
| **Variable-rate share, Spain (STOCK)** | **70%** | EFF 2022 — the outstanding stock, which is what the 2022–23 shock hit |
| Exposure M_var / income (exposed) | 1.44 | EFF 2022 (p2_13) |
| Mortgage rate level (median) | 2.06% | EFF 2022 (p2_13) |
| Residual term (median) | 14 y | EFF 2022 (p2_17) |

## External parameters

| Parameter | Value | Source | Status |
|---|---|---|---|
| Deposit beta | 0.32 | **Banco de España, DO 2312 (2023):** household time-deposit remuneration rose 136 bp, "just **32%** of the change in the 12-month Euribor"; vs 53–57% elsewhere in the euro area. | **confirmed (exact)** |
| Foreign ownership of Spanish equity | 50.3% | **BME (2023 record):** non-resident investors hold 50.3% of Spanish listed shares. | **confirmed (exact)** |
| Shareholder MPC | 0.05 | Financial-wealth MPC. For Spain, Bover / panel evidence finds the MPC out of **financial** wealth is ≈0 (housing ≈3 c/€); 0.05 is a conservative-high ceiling. | confirmed (order of magnitude) |
| Income persistence ρ / s.d. σ_η | 0.9831 / 0.072 | Banco de España, DT 2043 (Spanish earnings process) | cited |
| Annual MPC (aggregate) | 0.25 | Disciplined by the EFF wealth distribution; consistent with the euro-area HANK benchmark (Slacalek-Tristani-Violante 2020). See note below. | **resolved** |
| −2.7% consumption drop, variable borrowers 2024 | −2.7% | Banco de España / CaixaBank Research (2024) | cited |

## Cross-country variable-rate shares (the counterfactual) — STOCK, not flow

The counterfactual compares Spain's **outstanding stock** (70% variable, EFF) with the stock of
another country — **not** the flow of new loans. The two diverge a lot in Spain because the
new-business mix swings with the rate cycle while the stock moves slowly (see the Spain-flow note
below). The stock is what the 2022–23 shock repriced.

| Country | Value in the model | Best available anchor |
|---|---|---|
| Spain | 0.70 | EFF 2022, stock — solid |
| Germany | 0.12 | Fixed-rate country on every source (ECB; Albertazzi et al. new-business ARM ≈17%; HFCS 2010 ≈19.5%, falling since); headline ratio robust to 0.10–0.17 |
| Italy | 0.45 | ARM-heavy but shifting to fixed (HFCS 2010 ≈51.6%; Albertazzi new-business ≈66%); 0.45 is a reasonable ~2022 stock |
| Euro area | 0.25 | ECB (2025): euro-area *pure-ARM* stock ≈25% |

Sources for a clean cross-country series (in order of consistency with the EFF-based 0.70):
- **HFCS** (ECB) — same survey family as the EFF (the EFF *is* its Spanish component), same
  mortgage-type question for every euro-area country → an apples-to-apples **stock** share. The
  cleanest published figures are in **Ehrmann & Ziegelmeyer, ECB WP 1631, Table 3** (HFCS Wave 1,
  ~2010): variable-rate share of outstanding mortgages — Spain 82.9%, Portugal 84.5%, Netherlands
  82.7%, Austria 66.7%, Italy 51.6%, Greece 48.2%, Belgium 31.6%, **Germany 19.5%**, France 12.8%,
  euro area 45.5%. **Caveat — these are ~2010 levels**: Spain's own figure fell from 82.9% (2010)
  to 70.3% (EFF 2022) as fixed-rate lending accumulated, and every country shifted the same way, so
  the 2010 levels overstate 2022. Use the **latest HFCS wave (2021)** for a time-consistent set;
  the 2010 table is the right cite for the cross-country **ordering and method**.
- **ECB Data Portal / MIR — "outstanding amounts" by initial rate fixation** (not "new business"):
  the authoritative statistical **stock** series; build the share as (floating + ≤1y)/total.
- **ECB blog (May 2025), *Monetary policy transmission: from mortgage rates to consumption*** —
  euro-area pure-ARM stock ≈25%; ARMs prevalent in Spain/Italy, FRMs in Germany/France (qualitative
  by country).
- **Albertazzi, Fringuellotti & Ongena, ECB WP 2322**, *Fixed rate versus adjustable rate mortgages*
  — new-business ARM shares by country (ES ≈81%, IT ≈66%, DE ≈17%, FR ≈16%, 2007–15). Flow, but
  fixes the ordering and magnitudes.
- **EMF Hypostat** — "outstanding residential loans by interest-rate type" per country (industry).

Headline safety: Germany is low-ARM on every source, so the Spain-vs-Germany ratio (≈1.27×) is
robust regardless. Italy and the euro-area rows are illustrative (the response is linear in the
share); if cited, move them to the anchors above.

## Spain new business (flow) — the "de-risking" evidence

The **stock** stayed ~70% variable, but the **flow** of new mortgages swung with the cycle:
Spain issued mostly **fixed** in 2021–22 (variable new business fell to ~26% in 2H-2022; ECB), then
partly back to variable as the Euribor rose — **INE**, December 2023: **54.2% fixed / 45.8% variable**
of newly constituted dwelling mortgages. This flow story supports the "de-risking of new lending,
but the standing stock is what the shock hit" narrative. Source: INE, *Estadística de Hipotecas*.

## The annual MPC — resolved: it is disciplined by the wealth data

The 0.25 is a **calibration target, not an EFF number** (a wealth survey doesn't measure MPCs),
but it is **not a free choice** — it is the aggregate MPC the Spanish wealth distribution implies.
Two facts pin it down:

1. **The model can't produce a higher aggregate MPC while matching Spanish wealth.** Raising the
   transitory-income risk (the parameter that governs the MPC) *lowers* the aggregate MPC — more
   risk means more precautionary liquid saving, so households are better buffered — and it blows up
   the liquid and illiquid targets. Verified by execution (σ_ε sweep): at the calibrated point
   MPC ≈ 0.25 with liquid 0.26 / illiquid 3.5; at higher risk the MPC falls toward 0.20 while
   liquid/illiquid explode. So ~0.25 is the aggregate MPC consistent with the observed distribution.

2. **It agrees with the euro-area HANK benchmark.** Slacalek-Tristani-Violante (2020) use **annual**
   MPCs of 0.50 (hand-to-mouth) and 0.05 (non-hand-to-mouth); with a ~22% HtM share that aggregates
   to ≈0.15–0.20. The model's 0.25 sits just above that — a central, not a low, aggregate.

The much higher **direct-survey** MPCs (euro-area HFCS ≈0.46–0.50) are (i) self-reported/hypothetical
(they overstate) and (ii) the MPC of *constrained* households, not the population aggregate; matching
them would require a far more liquidity-poor distribution than Spain's. So the statement for the paper
is: **the aggregate annual MPC (0.25) is disciplined by the EFF wealth distribution and consistent
with the euro-area HANK benchmark — it is not a free parameter.** No robustness at a higher MPC is
needed or feasible.

Source: Slacalek, Tristani & Violante (2020), *Household balance sheet channels of monetary policy:
a back-of-the-envelope calculation for the euro area*, J. Econ. Dynamics & Control.

## Note for the paper's motivation

The euro-area household-finance literature (Christelis-Georgarakos-Jappelli et al., NBER w25082;
Hintermaier-Koeniger, *Quantitative Economics* 2024) finds monetary transmission to consumption
is **stronger in Spain and Italy than in Germany and France** — the same asymmetry this paper
quantifies, with independent empirical support from outside the model.

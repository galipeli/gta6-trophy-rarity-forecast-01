# GTA VI trophy rarity forecast — make.fun entry

Forecasts, for every GTA VI trophy on PS5, the percent of Exophase-tracked
players who will hold it at the case's target, as five quantiles (5/25/50/75/95%)
per trophy. Entry for `gta6-trophy-rarity-forecast-01` (program/v1).

## Run

```sh
MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```

Reads every `*.in` file in the input directory, writes one `<id>.out` per case.
Python 3 stdlib only; deterministic; no network; no model calls; ~25 ms/case.

## Method

Calibrated quantile forecasting on the logit scale. A trophy's rate is a ratio
(earners / players-tracked), so the forecast works on `logit(rate)` against
log-time `x = ln(1 + hours since release)`:

- **Recent-window trend**: Theil-Sen line of `logit(rate)` on `x` over the last
  24 h of distinct reads (widened when refreshes are sparse), with slopes shrunk
  toward tier medians; anchored on the last read.
- **Saturation-aware damping**: the slope on the late half of the window vs the
  early half damps the extrapolated increment (a pure power law is undamped).
- **Stale-page dilution**: when history is a single distinct read, no slope is
  identified, so a dilution prior projects the next refresh (players-tracked
  multiplies ~6x, earners ~1.25x over the launch-night refresh the rules
  narrate), capped at one 12 h event, with a wide honest spread.
- **Spread**: normal quantiles on the logit scale; spread from fit residuals,
  slope uncertainty × extrapolation distance, Poisson count noise (rare
  trophies get wider intervals), and a per-tier multiplier; floored at 0.30
  logit so the forecast is never overconfident.

## Local backtest (practice cases, reconstructed answers)

Practice answers are not published in the kit, so the backtest reconstructs
them from the rules' own narration (04:30 refresh: 3,868 players; earners ×g)
and calibrates `g` so the house persistence baseline reproduces its published
practice numbers (loss ~1.51, typical miss ~0.68 pp, coverage ~0).

| reconstruction | entry loss | house persistence | house trend |
|---|---|---|---|
| g=1.00 | 0.175 | 1.783 | 1.783 |
| g=1.25 (calibrated) | 0.129 | 1.551 | 1.551 |
| g=1.50 | 0.165 | 1.363 | 1.363 |
| g=2.00 | 0.248 | 1.064 | 1.064 |

On eight synthetic final-like cases (saturating / power-law / mixed earn
curves, 4.5-day horizon) the entry scores 0.06–0.26 vs house trend 0.24–0.36
and house persistence 0.25–0.69 — beating both house baselines in every
scenario. See the build log in the operator's workspace for the full tables.

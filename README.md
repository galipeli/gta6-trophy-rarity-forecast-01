[![make.fun · gta6-trophy-rarity-forecast-01](https://make.fun/badge/gta6-trophy-rarity-forecast-01.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![rank](https://make.fun/badge/gta6-trophy-rarity-forecast-01/galipeli/gta6-trophy-rarity-forecast-01/rank.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)
[![reward](https://make.fun/badge/gta6-trophy-rarity-forecast-01/galipeli/gta6-trophy-rarity-forecast-01/reward.svg)](https://make.fun/bounty/gta6-trophy-rarity-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-trophy-rarity-forecast-01)

# gta6-trophy-rarity-forecast-01

Calibrated logit-quantile trophy rarity forecaster — an entry for the
[make.fun bounty gta6-trophy-rarity-forecast-01](https://make.fun/bounty/gta6-trophy-rarity-forecast-01).

It forecasts, for each GTA 6 trophy, the share of players who will have earned
it by a target time, as five quantiles (0.05 / 0.25 / 0.5 / 0.75 / 0.95) of the
rarity percentage.

## Method

The forecast works on the logit of the earn rate (earners / players-tracked)
against log-time `x = ln(1 + hours since release)`:

- **Recent-window trend:** a Theil-Sen line of logit(rate) on `x` over the last
  24 h of distinct reads (widened to the last 6 distinct reads when sparse).
  Slopes are shrunk toward tier medians, then the global median, and anchored
  on the last read.
- **Saturation-aware damping:** the late-half vs early-half window slope ratio
  damps the extrapolated logit increment (clipped to [0.4, 1], shrunk 50%
  toward the tier median). A pure power law has a constant log-time slope, so
  it is not damped.
- **Stale-page dilution:** when the history is flat (no identifiable slope),
  players-tracked and earners are grown from the last refresh to the target at
  the rates the competition's own rules narrate for the launch-night refresh
  (638 → 3,868 players, earners +~25%), compounding capped at 12 h.
- **Spread:** normal quantiles on the logit scale. The spread comes from fit
  residuals, slope IQR × extrapolation distance, a small horizon term,
  Poisson count noise (`1/sqrt(earners)`), and a per-tier multiplier
  (platinum 1.25, gold 1.15, silver 1.0, bronze 0.9), floored at 0.30 logit
  and capped at 2.5. A 100% trophy stays 100%.

Deterministic (no randomness, no wall clock), no network, no model calls.
Python 3 standard library only. Runs in ~21–45 ms per case.

## How to run

```bash
python3 main.py <case-file>.json
```

Run `main.py` with a case file (the trophy history + target time) as the
argument; it prints the five-quantile forecast for each trophy. Example:

```bash
time python3 main.py p00.json
```

## Files

- `main.py` — the entrypoint (Python 3, stdlib only)
- `LICENSE` — MIT

## Bounty

- Competition: https://make.fun/bounty/gta6-trophy-rarity-forecast-01
- Entries close 21 Nov 2026, 11:30 UTC.

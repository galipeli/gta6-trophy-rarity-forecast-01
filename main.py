#!/usr/bin/env python3
"""Forecast GTA VI trophy earn rates at a target instant.

Reads one case per *.in file in the input directory (schema
makefun/gta6-trophies-case/v1): a list of trophies and a history of reads of
Exophase's trophy page, each read giving every trophy's earn rate (percent of
tracked players holding it) and the players-tracked count. Writes one
<id>.out per case: five quantiles (5/25/50/75/95%) of each trophy's earn rate
at the case's target, as {"trophies": [{"id": n, "q": [q05, q25, q50, q75, q95]}]}.

Method (calibrated quantile forecasting on the logit scale)
-----------------------------------------------------------
A trophy's rate is a ratio: earners / players-tracked. Both move, and the
page the house reads is stale between Exophase refreshes, so the forecast
works on the rate's logit against log-time, x = ln(1 + hours since release):

* Recent-window trend. Trophy curves saturate, so the slope is fit on the
  last 24 hours of distinct reads only (widened when refreshes are sparse),
  with a robust Theil-Sen line of logit(rate) on x per trophy. Slopes a
  trophy cannot identify shrink toward its tier's median slope, then the
  global median. The forecast anchors on the trophy's last read rate and
  extrapolates the fitted slope to the target.
* Saturation-aware damping. The slope on the late half of the window is
  compared with the slope on the early half: a declining slope means the
  curve is flattening, so the extrapolated logit increment is damped by
  that ratio (clipped to [0.4, 1], shrunk halfway to the tier median). A
  pure power-law curve has a constant log-time slope, so its ratio is ~1
  and it is not damped.
* Stale-page dilution. When the whole history is a single distinct read
  (the page never refreshed inside the window), no slope is identified.
  The forecast then uses a dilution prior: the tracked-player base
  multiplies at ~0.386 log-points per hour and earners at ~0.048 per hour
  from the last refresh to the target (the launch-night refresh the
  competition's own rules narrate took players 638 -> 3,868 while earners
  grew only ~25%), so rates fall by that ratio. The horizon the prior
  compounds over is capped at 12 hours: it describes one refresh event,
  not indefinite compounding. The spread is wide and honest there.
* Spread. Quantiles are normal on the logit scale around the median, with
  a spread from the fit residuals, the slope uncertainty (interquartile
  range of the pairwise slopes) times the extrapolation distance, a small
  horizon term, and Poisson count noise (1/sqrt(earners)) so rare trophies
  get wider intervals; a per-tier multiplier (platinum widest). The spread
  is floored at 0.30 logit and capped at 2.5, so the forecast is never
  overconfident. A trophy every tracked player holds stays at 100%.

Stdlib only, deterministic (no randomness, no wall clock: every time comes
from the case file), no network, no model calls.
"""

import json
import math
import os
from datetime import datetime, timezone

Z = (-1.6448536269514722, -0.6744897501960817, 0.0,
     0.6744897501960817, 1.6448536269514722)

# Dilution prior (log-points per hour), used only when history is flat:
# players-tracked multiplies ~6.06x and earners ~1.25x over the 4.67 h from
# the last stale read to the first refresh the rules narrate.
PI_P = math.log(3868.0 / 638.0) / 4.6667
PI_E = math.log(1.25) / 4.6667
# The prior compounds over at most this many hours: one refresh event.
DT_CAP_HOURS = 12.0

S_DILUTION = 0.55   # logit spread of the flat-history dilution prior
S_FLOOR = 0.30      # never narrower than this (logit)
S_CAP = 2.5         # never wider than this (logit)
TIER_SPREAD = {"platinum": 1.25, "gold": 1.15, "silver": 1.0, "bronze": 0.9}


def clamp(x, lo, hi):
    return lo if x < lo else hi if x > hi else x


def logit(pct):
    p = clamp(pct, 0.01, 99.99) / 100.0
    return math.log(p / (1.0 - p))


def inv_logit(l):
    if l >= 0:
        return 100.0 / (1.0 + math.exp(-l))
    e = math.exp(l)
    return 100.0 * e / (1.0 + e)


def parse_time(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def median(v):
    v = sorted(v)
    n = len(v)
    if n == 0:
        return None
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def theil_sen(xs, ys):
    """Robust line fit. Returns (intercept, slope, slope_iqr) or Nones."""
    n = len(xs)
    if n < 2:
        return None, None, None
    slopes = []
    for i in range(n):
        for j in range(i + 1, n):
            if xs[j] != xs[i]:
                slopes.append((ys[j] - ys[i]) / (xs[j] - xs[i]))
    if not slopes:
        return None, None, None
    slopes.sort()
    m = len(slopes)
    b = slopes[m // 2] if m % 2 else 0.5 * (slopes[m // 2 - 1] + slopes[m // 2])
    q1 = slopes[m // 4]
    q3 = slopes[(3 * m) // 4]
    inter = sorted(y - b * x for x, y in zip(xs, ys))
    a = inter[len(inter) // 2] if len(inter) % 2 else 0.5 * (
        inter[len(inter) // 2 - 1] + inter[len(inter) // 2])
    return a, b, q3 - q1


def distinct_reads(history):
    """Compress history to reads that change players_tracked or any rate."""
    out = []
    prev = None
    for h in history:
        key = (h.get("players_tracked"), tuple(h.get("rates") or []))
        if key != prev:
            out.append(h)
            prev = key
    return out


def last_known_rate(history, k):
    """Most recent non-null rate for trophy index k (None if never seen)."""
    for h in reversed(history):
        rates = h.get("rates") or []
        if k < len(rates) and rates[k] is not None:
            return rates[k]
    return None


def dilution_multiplier(first_at, target):
    """Rate multiplier of the stale-page dilution prior, capped at one event."""
    dt = (target - parse_time(first_at)).total_seconds() / 3600.0
    dt = min(max(dt, 0.0), DT_CAP_HOURS)
    return math.exp((PI_E - PI_P) * dt)


def quantiles(median_rate, spread):
    base = logit(median_rate)
    qs = [round(clamp(inv_logit(base + z * spread), 0.0, 100.0), 4)
          for z in Z]
    for i in range(1, 5):
        if qs[i] < qs[i - 1]:
            qs[i] = qs[i - 1]
    return qs


def spread_with_count_noise(base_spread, earners_last, tier):
    s = base_spread
    if earners_last is not None and earners_last > 0:
        s = math.sqrt(s * s + 1.0 / max(earners_last, 1.0))
    s *= TIER_SPREAD.get(tier, 1.0)
    return min(max(s, S_FLOOR), S_CAP)


def default_rate(history, trophies):
    """Prior rate for a trophy never observed: the cross-trophy median of the
    last-known rates (1.0% if nothing was ever observed)."""
    seen = []
    for k in range(len(trophies)):
        r = last_known_rate(history, k)
        if r is not None and r < 100.0:
            seen.append(r)
    m = median(seen)
    return m if m is not None else 1.0


def forecast_flat(case, trophies, history, target, rel0):
    """Dilution prior: no slope is identifiable, so project the refresh."""
    mult = dilution_multiplier(history[0]["at"], target)
    last = history[-1]
    last_p = max(last.get("players_tracked") or 1, 1)
    fallback = default_rate(history, trophies)
    out = {}
    for k, t in enumerate(trophies):
        r = last_known_rate(history, k)
        if r is None:
            r = fallback
        if r >= 100.0:
            out[t["id"]] = [100.0, 100.0, 100.0, 100.0, 100.0]
            continue
        m = clamp(r * mult, 0.01, 99.99)
        s = spread_with_count_noise(S_DILUTION, r * last_p / 100.0, t["tier"])
        out[t["id"]] = quantiles(m, s)
    return out


def forecast_trend(case, trophies, history, target, rel0):
    """Recent-window logit trend with saturation-aware damping."""
    hist = distinct_reads(history)
    xs_all = [math.log1p(max((parse_time(h["at"]) - rel0).total_seconds()
                             / 3600.0, 0.0)) for h in hist]
    x_tgt = math.log1p(max((target - rel0).total_seconds() / 3600.0, 0.0))
    x_last = xs_all[-1]
    dx = x_tgt - x_last

    # Recent window: last 24 h of distinct reads, widened when sparse.
    last_dt = parse_time(hist[-1]["at"])
    idx = [i for i, h in enumerate(hist)
           if (last_dt - parse_time(h["at"])).total_seconds() <= 24 * 3600]
    if len(idx) < 4:
        idx = list(range(max(0, len(hist) - 6), len(hist)))
    if len(idx) < 3:
        idx = list(range(len(hist)))
    win = [hist[i] for i in idx]
    xs_w = [xs_all[i] for i in idx]

    last = history[-1]
    last_p = max(last.get("players_tracked") or 1, 1)
    fallback = default_rate(history, trophies)

    # Per-trophy fits on the window, with tier/global slope shrinkage.
    fits = {}
    tier_slopes = {}
    all_slopes = []
    for k, t in enumerate(trophies):
        ex, ey = [], []
        for j, h in enumerate(win):
            rates = h.get("rates") or []
            if k >= len(rates) or rates[k] is None:
                continue
            ex.append(xs_w[j])
            ey.append(logit(rates[k]))
        a, b, iqr = theil_sen(ex, ey)
        resid = 0.0
        if a is not None and len(ex) > 2:
            rr = [y - (a + b * x) for x, y in zip(ex, ey)]
            resid = math.sqrt(sum(v * v for v in rr) / max(len(rr) - 2, 1))
        fits[k] = (a, b, resid, iqr if iqr is not None else 0.0, ex, ey)
        if b is not None:
            all_slopes.append(b)
            tier_slopes.setdefault(t["tier"], []).append(b)
    global_slope = median(all_slopes)
    tier_med = {k: median(v) for k, v in tier_slopes.items()}

    # Saturation damping: late-half slope vs early-half slope per trophy.
    damps = {}
    tier_damps = {}
    for k, t in enumerate(trophies):
        _, _, _, _, ex, ey = fits[k]
        d = 1.0
        if len(ex) >= 4:
            mid = (ex[0] + ex[-1]) / 2.0
            e_early = [(x, y) for x, y in zip(ex, ey) if x <= mid]
            e_late = [(x, y) for x, y in zip(ex, ey) if x > mid]
            if len(e_early) >= 2 and len(e_late) >= 2:
                _, b_e, _ = theil_sen([p[0] for p in e_early],
                                      [p[1] for p in e_early])
                _, b_l, _ = theil_sen([p[0] for p in e_late],
                                      [p[1] for p in e_late])
                if b_e is not None and b_l is not None and b_e > 0.05:
                    d = max(0.4, min(1.0, b_l / b_e))
        damps[k] = d
        tier_damps.setdefault(t["tier"], []).append(d)
    tier_med_damp = {k: median(v) for k, v in tier_damps.items()}

    out = {}
    for k, t in enumerate(trophies):
        r_last = last_known_rate(history, k)
        if r_last is None:
            r_last = fallback
        if r_last >= 100.0:
            out[t["id"]] = [100.0, 100.0, 100.0, 100.0, 100.0]
            continue
        _, b, resid, iqr, _, _ = fits[k]
        if b is None:
            b = tier_med.get(t["tier"], global_slope)
            iqr = 0.0
        if b is None:
            # No slope anywhere: fall back to the dilution prior.
            m = clamp(r_last * dilution_multiplier(history[0]["at"], target),
                      0.01, 99.99)
            s = spread_with_count_noise(S_DILUTION, r_last * last_p / 100.0,
                                        t["tier"])
            out[t["id"]] = quantiles(m, s)
            continue
        d = 0.5 * damps[k] + 0.5 * tier_med_damp.get(t["tier"], 1.0)
        y_tgt = logit(r_last) + d * b * dx
        m = clamp(inv_logit(y_tgt), 0.01, 99.99)
        s = math.sqrt(resid ** 2 + ((iqr / 2.0) * dx) ** 2 + (0.05 * dx) ** 2)
        s = spread_with_count_noise(s, r_last * last_p / 100.0, t["tier"])
        out[t["id"]] = quantiles(m, s)
    return out


def forecast_case(case):
    trophies = case["trophies"]
    history = case.get("history") or []
    target = parse_time(case["target"])
    rel_date = (case.get("game") or {}).get("release_date") or "2026-11-19"
    rel0 = datetime.fromisoformat(rel_date + "T00:00:00+00:00")
    if not history:
        return {t["id"]: quantiles(1.0, 1.0) for t in trophies}
    if len(distinct_reads(history)) < 2:
        return forecast_flat(case, trophies, history, target, rel0)
    return forecast_trend(case, trophies, history, target, rel0)


def main():
    in_dir = os.environ.get("MAKEFUN_INPUT_DIR") or "/input"
    out_dir = os.environ.get("MAKEFUN_OUTPUT_DIR") or "/output"
    os.makedirs(out_dir, exist_ok=True)
    names = sorted(n for n in os.listdir(in_dir) if n.endswith(".in"))
    for name in names:
        stem = name[:-3]
        with open(os.path.join(in_dir, name), "r", encoding="utf-8") as f:
            case = json.load(f)
        fc = forecast_case(case)
        rows = [{"id": t["id"], "q": fc[t["id"]]} for t in case["trophies"]]
        with open(os.path.join(out_dir, stem + ".out"), "w",
                  encoding="utf-8") as f:
            json.dump({"trophies": rows}, f, separators=(",", ":"))


if __name__ == "__main__":
    main()

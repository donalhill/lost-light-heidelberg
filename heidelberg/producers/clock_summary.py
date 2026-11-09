"""The asteroid clock's summary (requirement N2): from clock_batch.py's rows (one per segment near a predicted asteroid
track; the brightest flagged as the trail), the trails that pass the stated rule, each plate's offset (the median over
its trails of the trail's offset along the track from its predicted middle, as time), and over plates: the number, the
median offset with its bootstrap standard error, the robust spread (1.4826 times the median absolute deviation), the
plates off by 10 to 25 minutes and by more than an hour, and the years. The rule (checked against the trails of the
clock's first run): the astrometric fit of the cut-out uses at least
MIN_STARS Gaia stars; the trail lies within ACROSS of the predicted track; its measured length along the track is between
LENGTH[0] and LENGTH[1] times the predicted motion over the exposure (shorter or longer segments are not the asteroid).
python -E -s -m heidelberg.producers.clock_summary OUT.json CLOCK.csv [CLOCK.csv ...]"""
import json
import sys

import numpy as np
import pandas as pd

from heidelberg import config

K = config.A.clock
MIN_STARS, ACROSS, LENGTH = K.clock_min_stars, K.clock_across_arcsec, K.clock_length   # stars; arcsec; fraction of the predicted length
NBOOT, SEED = K.clock_boot, K.clock_seed


def trails(A):
    """The trails that pass the rule."""
    T = A[(A.get("trail") == True) & A.get("error", pd.Series(np.nan, index=A.index)).isna()].copy()  # noqa: E712
    ratio = (T.end_hi - T.end_lo) / T.length_pred
    return T[(T.fit_stars >= MIN_STARS) & (T.across.abs() <= ACROSS) & ratio.between(*LENGTH)]


def repeat_tests(median, se, n_far):
    """Whether the median offset is beyond clock_offset_sigmas standard errors of zero, or any plate beyond clock_far_min
    (reported; whether the times are shifted is stages/clock_correct.py's rule)."""
    return bool(abs(median) > K.clock_offset_sigmas * se or n_far > 0)


def summary(A):
    T = trails(A); S = T.groupby("plate").dt_min.median()
    rng = np.random.default_rng(SEED); boot = [np.median(rng.choice(S.values, len(S))) for _ in range(NBOOT)]
    years = T.groupby("plate").epoch.first().astype(int)
    return dict(rule=dict(min_stars=MIN_STARS, across_arcsec=ACROSS, length=list(LENGTH)), rows=int(len(A)), trails=int(len(T)),
                plates=int(len(S)), median_min=float(S.median()), median_se_min=float(np.std(boot)),
                robust_spread_min=float(1.4826 * np.median(np.abs(S - S.median()))),
                n_10_25=int(((S.abs() >= K.clock_band_min[0]) & (S.abs() <= K.clock_band_min[1])).sum()),
                n_25_60=int(((S.abs() > K.clock_band_min[1]) & (S.abs() <= K.clock_far_min)).sum()), n_over_60=int((S.abs() > K.clock_far_min).sum()),
                year_first=int(years.min()), year_last=int(years.max()), by_plate={k: round(float(v), 2) for k, v in S.items()},
                repeat_tests=repeat_tests(float(S.median()), float(np.std(boot)), int((S.abs() > K.clock_far_min).sum())))


def read_parts(files):
    """The batches' rows; a batch with none (no plate of it with a kept scan or a trail) adds none."""
    parts = [pd.read_csv(f) for f in files if open(f).read().strip()]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


if __name__ == "__main__":
    A = read_parts(sys.argv[2:])
    out = summary(A); json.dump(out, open(sys.argv[1], "w"), indent=1)
    print(f"{out['plates']} plates ({out['year_first']}-{out['year_last']}), {out['trails']} trails: median {out['median_min']:.2f} min "
          f"(se {out['median_se_min']:.2f}), robust spread {out['robust_spread_min']:.2f} min, {out['n_10_25']} off by 10-25 min, {out['n_25_60']} by 25-60 min, {out['n_over_60']} by more than an hour")

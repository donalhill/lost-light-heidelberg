"""The Bruce twin pairs from GAVO's plate metadata (requirement S1): the one rule for a twin pair, shared by the record
set (heidelberg/stages/record_set.py) and the shadow geometry (heidelberg/stages/record_geometry.py). A pair is an a
plate and the b plate that starts within pair_start_min of it with the nearest centre, within centre_sep_deg, both
exposed for at least min_exposure_min, from first_date to before end_before (config/analysis.json, section sample)."""
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

from heidelberg.config import A

S = A.sample


def plates(L):
    """The plates the pairing can use: a start time, a centre, an exposure of at least min_exposure_min, in the dates."""
    L = L.assign(t0=pd.to_datetime(L.startTime, errors="coerce")).dropna(subset=["t0", "centerAlpha", "exposure"])
    return L[(L.exposure >= S.min_exposure_min * 60) & (L.t0 >= S.first_date) & (L.t0 < S.end_before)]


def twin_pairs(L):
    """Every twin pair among the plates of L (GAVO's table): a, b, their GAVO URLs, plate a's mid-exposure (MJD),
    exposure, centre, year, and month."""
    L = plates(L); A_, B_ = L[L.plateID.str.endswith("a")], L[L.plateID.str.endswith("b")]
    rows = []
    for _, a in A_.iterrows():
        d = (B_.t0 - a.t0).abs().dt.total_seconds(); c = B_[d < S.pair_start_min * 60]
        if not len(c):
            continue
        sep = SkyCoord(a.centerAlpha, a.centerDelta, unit="deg").separation(SkyCoord(c.centerAlpha.values, c.centerDelta.values, unit="deg")).deg
        if sep.min() < S.centre_sep_deg:
            b = c.iloc[int(np.argmin(sep))]
            rows.append(dict(a=a.plateID, b=b.plateID, ua=a.accref, ub=b.accref, mjd=a.dateObs, exp=a.exposure, ra=a.centerAlpha,
                             dec=a.centerDelta, year=a.t0.year, month=a.t0.month))
    return pd.DataFrame(rows)


def match_sunlit(shadow, pool, S, n=None):
    """Sunlit pairs matched to shadow pairs, nearest in declination, season, and year, without replacement: the shadow pairs
    taken in a fixed random order (S.match_seed), n of them (default S.sunlit_per_shadow per shadow pair), so that a
    smaller draw is the start of a larger one. Beyond one per shadow pair the order repeats: a second pass matches every
    shadow pair again, in the same order. Returns the chosen pool indices in order."""
    n = round(S.sunlit_per_shadow * len(shadow)) if n is None else n
    targets = shadow.iloc[np.resize(np.random.default_rng(S.match_seed).permutation(len(shadow)), n)]
    feat = lambda d: np.c_[d.dec / S.match_dec_deg, np.cos(2 * np.pi * d.month / 12) * S.match_season,
                           np.sin(2 * np.pi * d.month / 12) * S.match_season, (d.year - 1900) / S.match_year]   # months on a circle
    chosen, free = [], pool.copy()
    for _, t in targets.iterrows():
        dd = np.linalg.norm(feat(free) - feat(t.to_frame().T.astype({"dec": float, "month": float, "year": float})), axis=1)
        k = int(np.argmin(dd)); chosen.append(free.index[k]); free = free.drop(free.index[k])
    return chosen

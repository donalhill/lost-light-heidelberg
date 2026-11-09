"""Plate times and eligibility from the plate catalogue (requirements S6, S7). The catalogue (res/plattenkatalog.csv,
GAVO's SVN, revision 9476) gives per plate Date (the civil date of the UT start), TM_START and TM_END (UT, written past
24 h when the exposure runs into the next day), EXPTIME (s, net), and remarks. Checked against the scanned logbooks for
twelve plates: Date with TM_START and TM_END agree with them within 30 s. GAVO's dateObs is 12 h early when TM_END is
written as 00:xx after a start before 0 h UT (it halves a negative span), and the catalogue's UT and JD are a day late on
469 plates (UT written past 24 h on a Date already moved on), so neither is used.

Mid-exposure = Date + the midpoint of TM_START and TM_END, a time past 24 h taken into the next day and TM_END moved a
day on when earlier than TM_START. A plate is eligible when its exposure is known (span within 5 min of EXPTIME; no
remark of a multiple exposure, an exposure over two nights, or an uncertain time) and its time and centre are plausible
(at mid-exposure the Sun below -12 deg and the plate centre above 10 deg altitude at Koenigstuhl). Plates missing from
the catalogue keep GAVO's times and are not eligible.

Writes the GAVO plate table (and the coverage table, if given) with dateObs, startTime, and endTime replaced, and the
columns time_source, span_s, remark, exposure_known, night_ok, eligible added.
python -E -s -m heidelberg.stages.plate_times CATALOGUE.csv LSW.csv OUT_LSW.csv [COVERAGE.csv OUT_COVERAGE.csv]"""
import sys

import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, SkyCoord, get_sun
from astropy.time import Time

from heidelberg.config import A
from heidelberg.lib.shadow_geom import KOENIGSTUHL

GAP_MIN = A.times.gap_min         # the largest difference between span and net exposure for one continuous exposure (min)
SUN_MAX, ALT_MIN = A.times.sun_max, A.times.alt_min    # deg at mid-exposure
REMARKS = r"multiple|night|uncertain"
MJD0 = pd.Timestamp("1858-11-17")


def hours(s):
    """'HH:MM:SS' (HH may exceed 23) to hours; NaN where unreadable."""
    x = s.astype(str).str.extract(r"^(\d+):(\d+):(\d+(?:\.\d*)?)$").astype(float)
    return x[0] + x[1] / 60 + x[2] / 3600


def times(K):
    """Per catalogue row: start, mid, and end (MJD) and the span (s) by the rule of S6."""
    s, e = hours(K.TM_Start), hours(K.TM_End)
    e = e.where(e >= s, e + 24)
    day = (pd.to_datetime(K.Date, errors="coerce") - MJD0).dt.days
    return pd.DataFrame(dict(start=day + s / 24, mid=day + (s + e) / 48, end=day + e / 24, span_s=(e - s) * 3600), index=K.index)


def exposure_known(span_s, exptime_s, remark):
    return ((span_s - exptime_s).abs() <= 60 * GAP_MIN) & ~remark.fillna("").str.lower().str.contains(REMARKS)


def night_ok(mjd, ra, dec):
    t = Time(np.asarray(mjd, float), format="mjd"); fr = AltAz(obstime=t, location=KOENIGSTUHL)
    sun = get_sun(t).transform_to(fr).alt.deg
    alt = SkyCoord(np.asarray(ra, float), np.asarray(dec, float), unit="deg").transform_to(fr).alt.deg
    return (sun < SUN_MAX) & (alt > ALT_MIN)


def isot(mjd):
    return [Time(m, format="mjd").isot[:19] if np.isfinite(m) else "" for m in mjd]


def retime(L, K):
    """The GAVO table L with the catalogue's times and the eligibility columns (K indexed by plate)."""
    L = L.copy(); k = K.reindex(L.plateID); T = times(k); ok = T.mid.notna().to_numpy()
    L["time_source"] = np.where(ok, "catalogue", "gavo")
    L.loc[ok, "dateObs"] = T.mid.to_numpy()[ok]
    for col, src in (("startTime", "start"), ("endTime", "end")):
        if col in L:
            L.loc[ok, col] = np.array(isot(T[src].to_numpy()), dtype=object)[ok]
    L["span_s"] = T.span_s.to_numpy(); L["remark"] = k.Remark.to_numpy()
    L["exposure_known"] = ok & exposure_known(T.span_s, pd.to_numeric(k.Exptime, errors="coerce"), k.Remark).to_numpy()
    good = ok & L.centerAlpha.notna().to_numpy() & L.centerDelta.notna().to_numpy()
    L["night_ok"] = False
    L.loc[good, "night_ok"] = night_ok(L.dateObs[good], L.centerAlpha[good], L.centerDelta[good])
    L["eligible"] = L.exposure_known & L.night_ok
    return L


def catalogue(path):
    K = pd.read_csv(path, sep=";", dtype=str).iloc[:, :16]
    return K.drop_duplicates("Plate").set_index("Plate")     # 11 duplicated IDs, none a Bruce plate


if __name__ == "__main__":
    CAT, LSW, OUT = sys.argv[1:4]
    K = catalogue(CAT); L = retime(pd.read_csv(LSW), K); L.to_csv(OUT, index=False)
    old = pd.read_csv(LSW).dateObs.to_numpy(); d = (L.dateObs.to_numpy() - old) * 1440
    print(f"{len(L)} plates; times from the catalogue: {(L.time_source == 'catalogue').sum()}; moved from GAVO's dateObs by more than "
          f"1 min: {int((np.abs(d) > 1).sum())} (12 h: {int((np.abs(np.abs(d) - 720) <= 1).sum())})")
    print(f"exposure not known: {int((~L.exposure_known).sum())}; time or centre implausible: {int((~L.night_ok).sum())}; eligible: {int(L.eligible.sum())}")
    if len(sys.argv) > 5:
        C = pd.read_csv(sys.argv[4]); T = times(K.reindex(C.plateID)); ok = T.mid.notna().to_numpy()
        C.loc[ok, "dateObs"] = T.mid.to_numpy()[ok]; C.to_csv(sys.argv[5], index=False)
        print(f"coverage: {len(C)} plates, {ok.sum()} retimed")

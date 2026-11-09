"""The plate times part one uses: the catalogue's (stages/plate_times.py), shifted by the asteroid clock's median offset
on the main sample's plates when that offset is clear, that is, larger than clock_offset_sigmas standard errors and larger
than clock_correct_min minutes; otherwise unchanged. The offset is the time at which the asteroids stood where their
trails are, less the catalogue's time (producers/clock_summary.py), so the corrected time is the catalogue's plus the
offset. Every plate is shifted alike (the two plates of a twin pair were exposed together). A clock that failed, or found
no trail, changes nothing. The coverage map and the shadow geometry read the times written here.
python -E -s -m heidelberg.stages.clock_correct CLOCK_SUMMARY.json LSW_PLATES.csv OUT_PLATES.csv OUT.json"""
import json
import shutil
import sys

import pandas as pd

from heidelberg import config

K = config.A.clock


def correction(S):
    """The minutes to add to every plate's time, and why."""
    if S.get("failed") or "median_min" not in S:
        return 0.0, "the clock failed or found no trail"
    m, se = float(S["median_min"]), float(S["median_se_min"])
    if abs(m) > K.clock_offset_sigmas * se and abs(m) > K.clock_correct_min:
        return m, f"median {m:.2f} min, beyond {K.clock_offset_sigmas} standard errors ({se:.2f}) and {K.clock_correct_min} min"
    return 0.0, f"median {m:.2f} min (standard error {se:.2f}): within {K.clock_offset_sigmas} standard errors or {K.clock_correct_min} min"


if __name__ == "__main__":
    SUM, LSW, OUT, OUT_J = sys.argv[1:5]
    off, why = correction(json.load(open(SUM)))
    if not off:
        shutil.copyfile(LSW, OUT)
    else:
        L = pd.read_csv(LSW)
        L["dateObs"] = L.dateObs + off / 1440
        for c in ("startTime", "endTime"):
            if c in L:
                t = pd.to_datetime(L[c], errors="coerce") + pd.to_timedelta(off, unit="min")
                L[c] = t.dt.strftime("%Y-%m-%dT%H:%M:%S").where(t.notna(), L[c])
        L.to_csv(OUT, index=False)
    json.dump(dict(applied=bool(off), offset_min=off, reason=why), open(OUT_J, "w"), indent=1)
    print(f"plate times {'shifted by ' + format(off, '.2f') + ' min' if off else 'unchanged'}: {why}")

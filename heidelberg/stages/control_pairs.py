"""Different-night control pairs for the record run: for record pairs chosen at random, plate a with a plate b of the
same field (centres within control_sep_deg) from another night (more than control_min_days and less than control_max_days away; exposure
at least the sample's minimum); then those not eligible (requirement S7) removed (config/analysis.json, section controls). No flash can
coincide across nights, so what passes the selection is background, including the static sources the shifted pairs
cannot see. Plates of the tuning pairs are never used (requirement S2).
python -E -s -m heidelberg.stages.control_pairs RECORD_SET.csv LSW_PLATES.csv TUNING_PAIRS.csv OUT.csv"""
import sys

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

from heidelberg.config import A

K = A.controls
RS, LSW, TUN, OUT = sys.argv[1:5]; N = K.control_n_drawn
R = pd.read_csv(RS); L = pd.read_csv(LSW); T = pd.read_csv(TUN)
L = L[(L.exposure >= A.sample.min_exposure_min * 60) & L.plateID.str.endswith("b") & ~L.plateID.isin(set(T.a) | set(T.b))]
cB = SkyCoord(L.centerAlpha.values, L.centerDelta.values, unit="deg"); rows = []
for _, p in R.iloc[np.random.default_rng(K.control_seed).permutation(len(R))].iterrows():
    if len(rows) >= N:
        break
    dt = L.dateObs.values - p.mjd; sep = SkyCoord(p.ra, p.dec, unit="deg").separation(cB).deg
    k = np.flatnonzero((np.abs(dt) > K.control_min_days) & (np.abs(dt) < K.control_max_days) & (sep < K.control_sep_deg))
    if len(k):
        j = k[np.argmin(np.abs(dt[k]))]
        rows.append(dict(a=p.a, b=L.plateID.values[j], ua=p.ua, ub=L.accref.values[j], gb=np.nan, days=round(float(dt[j]), 2), sep_deg=round(float(sep[j]), 3), role=p.role))
C = pd.DataFrame(rows)
E = set(pd.read_csv(LSW).query("eligible").plateID) if "eligible" in L else None
if E is not None:                                               # requirement S7, after the choice, so that it is unchanged
    C = C[C.a.isin(E) & C.b.isin(E)]
if "examined" in R:                                             # pairs examined while building the pipeline are not measured
    C = C[~C.a.isin(R[R.examined].a)]
C.to_csv(OUT, index=False)
print(f"{len(C)} control pairs ({(C.role == 'shadow').sum()} on shadow pairs' fields); days apart median {C.days.abs().median():.1f}; offset median {C.sep_deg.median():.3f} deg")

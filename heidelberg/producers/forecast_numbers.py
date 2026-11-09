"""What the main sample should hold if no source is a glint, at the tuning pairs' background per pair: selected sources in
all, and in Earth's shadow (weighted by each pair's mean shadow share), both allowing for the share of pairs the selection
refused on the tuning pairs. Rounded to two figures, since they are forecasts. Stated in the draft's Introduction. Also how
far the shadow's axis moves on the sky, seen from the site, over each shadow pair's common exposure (median, degrees).
python -E -s -m heidelberg.producers.forecast_numbers RECORD_SET.csv J3_NUMBERS.json OUT.json"""
import json
import sys

import numpy as np
import pandas as pd
from astropy.coordinates import get_sun
from astropy.time import Time

from heidelberg.lib import shadow_model as sh

RS, J3, OUT = sys.argv[1:4]
R = pd.read_csv(RS); J = json.load(open(J3))
B, keep = J["tuning_background"], 1 - J["refused_share"]
equivalent = float(R[R.role == "shadow"].s_mean.sum())                 # pairs wholly in shadow with the same shadowed sky
sources, in_shadow = B * len(R) * keep, B * equivalent * keep
two = lambda x: float(f"{x:.2g}")


def axis(t):
    """The direction from the site to the point of the shadow's axis at geostationary radius (GCRS unit vectors)."""
    obs = sh.SITE.get_gcrs(t).cartesian.xyz.to("km").value.T
    a = -get_sun(t).cartesian.xyz.value.T; a = a / np.linalg.norm(a, axis=1, keepdims=True)
    d = sh.R_GEO * a - obs
    return d / np.linalg.norm(d, axis=1, keepdims=True)


S = R[R.role == "shadow"]; half = S.exp.values / 86400 / 2            # exp: the common exposure, in seconds
moved = np.degrees(np.arccos(np.clip((axis(Time(S.mjd.values - half, format="mjd")) * axis(Time(S.mjd.values + half, format="mjd"))).sum(1), -1, 1)))
out = dict(shadow_equivalent_pairs=round(equivalent), sources=two(sources), sources_in_shadow=two(in_shadow),
           sources_exact=sources, sources_in_shadow_exact=in_shadow, exposure_median_min=float(np.median(S.exp.values)) / 60,
           shadow_motion_median_deg=float(np.median(moved)), shadow_motion_range_deg=[float(x) for x in np.percentile(moved, [10, 90])])
json.dump(out, open(OUT, "w"), indent=1); print(out)

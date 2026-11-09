"""The record set, from GAVO metadata alone (no measurement): every Bruce twin pair from 1900 to 1957-10-03 (before
Sputnik) touched by Earth's shadow at GEO over the 6 x 6 deg box, and sunlit pairs matched to a random subset of them in
declination, season, and year (nearest neighbour, without replacement), the tuning pairs excluded (already seen). The
pairs not eligible (requirement S7: exposure not known, or time or centre implausible) are then removed, after the
matching, so that the matching is unchanged. Per pair the covariates that could differ by role for reasons other than
glints: galactic and ecliptic latitude, distance from the anti-solar point, the Moon's distance and phase, altitude at
mid-exposure, exposure, year.
Writes OUT.csv (one row per pair: role shadow | sunlit, the GAVO URLs, the covariates, examined) and prints the comparison by role.
The twin-pair rule is heidelberg/lib/sample.py's; the matching's sunlit pairs per shadow pair, seed, and scales, and the examined pair, are
in config/analysis.json (section sample).
python -E -s -m heidelberg.stages.record_set LSW_PLATES.csv GEOMETRY.csv TUNING_PAIRS.csv OUT.csv"""
import sys

import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, GeocentricTrueEcliptic, SkyCoord, get_body, get_sun
from astropy.time import Time

from heidelberg.config import A
from heidelberg.lib.sample import match_sunlit, plates, twin_pairs
from heidelberg.lib.shadow_geom import KOENIGSTUHL

S = A.sample
EXAMINED = tuple(S.examined)                  # studied while building the pipeline; flagged
LSW, GEOM, TUN, OUT = sys.argv[1:5]
L = pd.read_csv(LSW); Lp = plates(L)
E = set(Lp.plateID[Lp.eligible]) if "eligible" in Lp else set(Lp.plateID)   # requirement S7 (stages/plate_times.py)
P = twin_pairs(L)
tuning = set(pd.read_csv(TUN).a)
G = pd.read_csv(GEOM).set_index("a")
P["role"] = np.where(P.a.isin(G.index), "shadow", "sunlit"); n_tuning = int(P.a.isin(tuning).sum()); P = P[~P.a.isin(tuning)].reset_index(drop=True)
shadow = P[P.role == "shadow"]; pool = P[P.role == "sunlit"].copy()
chosen = match_sunlit(shadow, pool, S)                         # sample.sunlit_per_shadow per shadow pair
R = pd.concat([shadow, P.loc[chosen]]).reset_index(drop=True)
R["s_mean"] = R.a.map(G.s_mean).fillna(0.0); R["partial"] = R.a.map(G.partial).fillna(0.0)
t = Time(R.mjd.values, format="mjd"); c = SkyCoord(R.ra.values, R.dec.values, unit="deg")
R["gal_b"] = c.galactic.b.deg; R["ecl_lat"] = c.transform_to(GeocentricTrueEcliptic(equinox=t)).lat.deg
sun = get_sun(t); anti = SkyCoord((sun.ra.deg + 180) % 360, -sun.dec.deg, unit="deg"); R["antisun_deg"] = c.separation(anti).deg
moon = get_body("moon", t); R["moon_deg"] = moon.separation(c).deg                   # in the Moon's (geocentric) frame: from the barycentre it lies along the Earth's direction
R["moon_phase"] = 0.5 * (1 - np.cos(np.radians(sun.separation(moon).deg)))     # illuminated fraction
R["altitude"] = c.transform_to(AltAz(obstime=t, location=KOENIGSTUHL)).alt.deg
R["examined"] = R.a.isin(EXAMINED)                              # flagged after matching, so the matching is unchanged
el = R.a.isin(E) & R.b.isin(E); n_inel = R[~el].role.value_counts().to_dict(); R = R[el].reset_index(drop=True)   # S7, after matching
R.to_csv(OUT, index=False)
print(f"examined while building the pipeline, flagged and not measured: {', '.join(R[R.examined].a)}")
print(f"record set: {len(R)} pairs ({(R.role == 'shadow').sum()} shadow, {(R.role == 'sunlit').sum()} sunlit); tuning pairs excluded: {n_tuning}; pairs not eligible (S7): {n_inel}")
cov = ["dec", "year", "exp", "gal_b", "ecl_lat", "antisun_deg", "moon_deg", "moon_phase", "altitude"]
S_, U_ = R[R.role == "shadow"], R[R.role == "sunlit"]
print(f"{'covariate':12s} {'shadow mean':>12s} {'sunlit mean':>12s} {'std diff':>9s}")
for v in cov:
    sd = np.sqrt(0.5 * (S_[v].var() + U_[v].var()))
    print(f"{v:12s} {S_[v].mean():12.2f} {U_[v].mean():12.2f} {(S_[v].mean() - U_[v].mean()) / sd if sd else 0:9.2f}")

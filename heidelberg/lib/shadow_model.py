"""Earth's shadow at an orbit's radius, seen from the observing site: how much of an exposure each line of sight spends
in shadow. Used by heidelberg/lib/shadow_geom.py; checked against an implementation that shares no code with it by
heidelberg/tests/test_shadow.py. The site, radii, and instants are in config/analysis.json (section physics)."""
import warnings
from dataclasses import dataclass, replace
from functools import lru_cache

import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, get_sun
from erfa import ErfaWarning
from astropy.utils.exceptions import AstropyWarning
from astropy.utils import iers

from heidelberg.config import A, C
from heidelberg.lib.geometry import deproject, unit                             # noqa: F401  (used here and re-exported)

iers.conf.auto_download = False     # the Sun's position needs no fresh Earth-rotation tables; downloading them stalls runs

warnings.simplefilter("ignore", ErfaWarning)
warnings.simplefilter("ignore", AstropyWarning)

PH = A.physics
SITE = EarthLocation.from_geodetic(lon=PH.site_lon_deg * u.deg, lat=PH.site_lat_deg * u.deg, height=PH.site_height_m * u.m)   # Koenigstuhl
R_GEO = PH.geo_radius_km


def distinct(t):
    """The distinct times in t, and where each element of t sits among them (the Sun is then found once per time)."""
    if t.isscalar:
        return t, None
    _, i, inv = np.unique(np.c_[t.jd1, t.jd2], axis=0, return_index=True, return_inverse=True)
    return t[i], inv.ravel()


def antisun_axis(t):
    """Unit vector of the anti-solar direction (GCRS: ICRS-aligned, J2000 axes)."""
    tu, inv = distinct(t)
    s = get_sun(tu)
    a = -unit(s.ra.deg, s.dec.deg)
    return a if inv is None else a[inv]


# ------------------------------------------------------------------ the shadow model
# A glint is a flash at one unknown instant, equally likely at any moment of the exposure. At that instant it is dark if
# the point where the site's line of sight meets the orbit is in Earth's shadow; a glint's brightness follows how much of
# the Sun's disc that point sees, so its chance of being recorded is taken to follow it too. A point on the plate is
# therefore in shadow for a share of the exposure: the average over the exposure of the share of the Sun's disc the Earth
# hides from the orbit there. Seen from the site the shadow slides across the field during an exposure (the site moves
# about 1,000 km in an hour at Heidelberg's latitude); from the centre of the Earth it barely moves. Under the hypothesis that r is the rate of
# transients in shadow over the rate in sunlight, a point in shadow for a share s of the exposure has rate b (1 - (1 - r) s),
# b the rate with no effect; with s 0 or 1 everywhere this is the plain count in shadow against sunlight.
# Edge: "disc" (the hidden share of the Sun's disc, from umbra, 8.42 deg, to penumbra, 8.96 deg, at the geostationary
# radius, half hidden at 8.69); "penumbra" (dark wherever any of it is hidden), with earth > 1 enlarging the Earth for its
# atmosphere (Danjon's 1/85), gives the shadow at its widest, for deciding which plates it can reach.

R_EARTH, R_SUN = PH.earth_radius_km, PH.sun_radius_km
DANJON = 1 + 1 / PH.danjon_parts
NT = PH.shadow_instants                                                          # instants per exposure
CHUNK = C.memory.shadow_chunk                                                    # points evaluated at once


@dataclass(frozen=True)
class Shadow:
    frame: str = "topo"
    edge: str = "disc"
    earth: float = 1.0
    r_orbit: float = R_GEO
    nt: int = NT

    def but(self, **k):
        return replace(self, **k)


MAIN = Shadow("topo")


def wide(frame, r_orbit=R_GEO):
    """The shadow at its widest (penumbra's outer edge, the Earth enlarged for its atmosphere), over the whole exposure:
    plates it misses are clear of the shadow."""
    return Shadow(frame, edge="penumbra", earth=DANJON, r_orbit=r_orbit)


def instants(p, S):
    """The times (datetime64) at which plate p's exposure is sampled: S.nt instants spread evenly over it."""
    mid = np.datetime64(np.asarray(p.ut, dtype="datetime64[ns]"))
    k = (np.arange(S.nt) + 0.5) / S.nt - 0.5
    return mid + (k * float(p.exp) * 60e9).astype("timedelta64[ns]")


@lru_cache(maxsize=C.memory.sun_cache)
def _at(t):
    """The site's and the Sun's positions (km, GCRS) and the anti-sun axis at time t (datetime64, as an int of ns)."""
    from astropy.time import Time
    tt = Time(np.datetime64(t, "ns"))
    return (SITE.get_gcrs(tt).cartesian.xyz.to(u.km).value.reshape(3), get_sun(tt).cartesian.xyz.to(u.km).value.reshape(3),
            np.asarray(antisun_axis(tt)).reshape(3))


def orbit_point(d, t, S):
    """Where each line of sight d (unit vectors, n x 3) meets the orbit's sphere at time t (km, GCRS): from the centre of
    the Earth (geo) or from the site (topo)."""
    if S.frame == "geo":
        return S.r_orbit * d
    obs = _at(int(np.datetime64(t, "ns").astype(np.int64)))[0]
    a = d @ obs; rng = -a + np.sqrt(a ** 2 - (obs @ obs - S.r_orbit ** 2))
    return obs + rng[:, None] * d


def hidden(P, sun, earth=1.0):
    """Share of the Sun's disc hidden by the Earth, seen from points P (km, n x 3); sun: the Sun's position (km). The Sun
    is small (0.27 deg), so the overlap is the area shared by a small disc and the Earth's limb near it, the limb kept at
    its true angular distance and curvature."""
    dE = np.linalg.norm(P, axis=1); vS = sun[None, :] - P; dS = np.linalg.norm(vS, axis=1)
    rE = np.arcsin(np.clip(earth * R_EARTH / dE, -1, 1)); rS = np.arcsin(R_SUN / dS)
    th = np.arccos(np.clip(np.sum(-P * vS, 1) / (dE * dS), -1, 1))
    b = np.tan(rE); d = b + (th - rE); a = rS                                    # planar: limb's curvature, distance kept
    with np.errstate(invalid="ignore", divide="ignore"):
        lens = (a ** 2 * np.arccos(np.clip((d ** 2 + a ** 2 - b ** 2) / (2 * d * a), -1, 1))
                + b ** 2 * np.arccos(np.clip((d ** 2 + b ** 2 - a ** 2) / (2 * d * b), -1, 1))
                - 0.5 * np.sqrt(np.clip((-d + a + b) * (d + a - b) * (d - a + b) * (d + a + b), 0, None)))
    return np.where(d >= a + b, 0.0, np.where(d <= b - a, 1.0, lens / (np.pi * a ** 2)))


def darkness(d, t, S):
    """How dark the orbit is along lines of sight d (unit vectors) at one time t, under shadow model S (0 sunlit, 1 dark)."""
    sun = _at(int(np.datetime64(t, "ns").astype(np.int64)))[1]
    P = orbit_point(d, t, S)
    h = hidden(P, sun, S.earth)
    return h if S.edge == "disc" else (h > 0).astype(float)


def dark_dirs(d, p, S=MAIN):
    """Share of plate p's exposure in shadow along lines of sight d (unit vectors, n x 3), under model S. p: a plate row
    (ut mid-exposure, exp in min). Taken in chunks of points, the instants summed as they come, so memory stays small."""
    ts = instants(p, S); out = np.empty(len(d))
    for a in range(0, len(d), CHUNK):
        acc = np.zeros(min(CHUNK, len(d) - a))
        for t in ts:
            acc += darkness(d[a:a + CHUNK], t, S)
        out[a:a + CHUNK] = acc / len(ts)
    return out


def dark(xi, eta, p, S=MAIN):
    """dark_dirs at plate-frame offsets (xi, eta) (deg) about p's centre (ra_c, dec_c)."""
    d = unit(*deproject(np.ravel(xi), np.ravel(eta), p.ra_c, p.dec_c)).reshape(-1, 3)
    return dark_dirs(d, p, S).reshape(np.shape(xi))

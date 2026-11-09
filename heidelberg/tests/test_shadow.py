"""Independent check of the shadow geometry used for the Heidelberg pairs (heidelberg/lib/shadow_model.py, Koenigstuhl as the
site): darkness along random lines of sight near the anti-sun point, at random times 1900-1957, against a separate
implementation that shares no code with it (WGS84 site, GMST rotation, IAU 1976 precession, Meeus's low-precision Sun,
its own ray-sphere intersection and umbra test). Points within EDGE deg of the umbra edge, as seen from the orbit, are
skipped (the two Sun positions differ by ~0.01 deg and the site by precession residuals of a few km).
python -E -s -m heidelberg.tests.test_shadow [N]"""
import sys

import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation

from heidelberg.lib import shadow_model as sh

LAT, LON, H = 49.3978, 8.7206, 0.560
sh.SITE = EarthLocation.from_geodetic(lon=LON * u.deg, lat=LAT * u.deg, height=H * 1000 * u.m)
RE, RS, AU, RGEO, EDGE = 6371.0, 696000.0, 1.495978707e8, 42164.0, 0.3


def jd(t):
    return t.astype("datetime64[ns]").astype(np.int64) / 86400e9 + 2440587.5


def precess_to_j2000(jdate):
    """IAU 1976 precession matrix, mean equator of date -> J2000."""
    T = (jdate - 2451545.0) / 36525; a = np.radians(1 / 3600)
    z_, z, th = (2306.2181 * T + 0.30188 * T ** 2) * a, (2306.2181 * T + 1.09468 * T ** 2) * a, (2004.3109 * T - 0.42665 * T ** 2) * a
    Rz = lambda g: np.array([[np.cos(g), -np.sin(g), 0], [np.sin(g), np.cos(g), 0], [0, 0, 1]])
    Ry = lambda g: np.array([[np.cos(g), 0, np.sin(g)], [0, 1, 0], [-np.sin(g), 0, np.cos(g)]])
    return (Rz(z) @ Ry(-th) @ Rz(z_)).T


def site(jdate):
    """Koenigstuhl, km, J2000 axes: WGS84 -> rotated by GMST (mean equator of date) -> precessed."""
    a, f = 6378.137, 1 / 298.257223563; e2 = f * (2 - f); la = np.radians(LAT)
    N = a / np.sqrt(1 - e2 * np.sin(la) ** 2)
    rho, zc = (N + H) * np.cos(la), (N * (1 - e2) + H) * np.sin(la)
    T = (jdate - 2451545.0) / 36525
    gmst = np.radians((280.46061837 + 360.98564736629 * (jdate - 2451545.0) + 0.000387933 * T ** 2) % 360)
    th = gmst + np.radians(LON)
    return precess_to_j2000(jdate) @ np.array([rho * np.cos(th), rho * np.sin(th), zc])


def sun(jdate):
    """Meeus (low precision, ~0.01 deg): geocentric Sun, km, J2000 axes."""
    T = (jdate - 2451545.0) / 36525
    L0, M = 280.46646 + 36000.76983 * T, np.radians(357.52911 + 35999.05029 * T)
    e = 0.016708634 - 0.000042037 * T
    C = (1.914602 - 0.004817 * T) * np.sin(M) + (0.019993 - 0.000101 * T) * np.sin(2 * M) + 0.000289 * np.sin(3 * M)
    lon = np.radians(L0 + C); nu = M + np.radians(C); R = 1.000001018 * (1 - e ** 2) / (1 + e * np.cos(nu)) * AU
    eps = np.radians(23.4392911 - 0.0130042 * T)                     # mean obliquity of date
    v = R * np.array([np.cos(lon), np.cos(eps) * np.sin(lon), np.sin(eps) * np.sin(lon)])
    return precess_to_j2000(jdate) @ v


def independent(d, jdate):
    """1 if the orbit point along d is in full umbra, 0 if fully sunlit, nan near the edge."""
    o = site(jdate); s = sun(jdate)
    b = d @ o; rng = -b + np.sqrt(b ** 2 - (o @ o - RGEO ** 2)); P = o + rng[:, None] * d
    toE, toS = -P, s - P
    sep = np.degrees(np.arccos(np.clip(np.sum(toE * toS, 1) / (np.linalg.norm(toE, axis=1) * np.linalg.norm(toS, axis=1)), -1, 1)))
    rE = np.degrees(np.arcsin(RE / np.linalg.norm(P, axis=1))); rS = np.degrees(np.arcsin(RS / np.linalg.norm(toS, axis=1)))
    out = np.full(len(d), np.nan)
    out[sep < rE - rS - EDGE] = 1.0; out[sep > rE + rS + EDGE] = 0.0
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    S = sh.Shadow("topo", edge="disc", r_orbit=RGEO) if hasattr(sh, "Shadow") else sh.MAIN
    g = np.random.default_rng(7); bad = tested = 0
    t0, t1 = np.datetime64("1900-01-01"), np.datetime64("1957-10-01")
    for k in range(n):
        t = t0 + (g.uniform() * (t1 - t0).astype("timedelta64[s]").astype(int)).astype("timedelta64[s]")
        jdate = jd(np.array([t]))[0]; s = sun(jdate); a = -s / np.linalg.norm(s)
        d = a + g.normal(0, 0.18, (400, 3)); d /= np.linalg.norm(d, axis=1)[:, None]       # ~10 deg around anti-sun
        ref = independent(d, jdate); ours = sh.darkness(d, np.datetime64(t, "ns"), S)
        m = np.isfinite(ref); tested += m.sum(); bad += int((np.round(ours[m]) != ref[m]).sum())
    print(f"{tested} lines of sight away from the edge at {n} times: {bad} disagree")
    assert bad == 0

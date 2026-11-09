"""Earth's shadow at geostationary radius seen from Koenigstuhl: the share of a pair's exposure each line of sight spends
in shadow, from the shadow model (lib/shadow_model.py, independently checked by heidelberg/tests/test_shadow.py). Frames: "topo" (seen from the site: the test) and "geo" (from the Earth's centre: an anti-solar disc of
the same size, the discriminator)."""
from types import SimpleNamespace

import numpy as np
from astropy.time import Time

from heidelberg.config import A
from heidelberg.lib import shadow_model as sh

KOENIGSTUHL = sh.SITE                                         # the site (config/analysis.json, section physics)
FRAMES = {"topo": sh.MAIN, "geo": sh.MAIN.but(frame="geo")}


def unit(ra, dec):
    r, d = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.c_[np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]


def plate(mjd_mid, exp_s, shift_days=0.0):
    """The model's plate record: mid-exposure time (MJD, moved by shift_days for a replica) and exposure (min)."""
    t = Time(mjd_mid + shift_days, format="mjd").datetime64
    return SimpleNamespace(ut=np.datetime64(t, "ns"), exp=exp_s / 60.0)


def dark(ra, dec, mjd_mid, exp_s, frame="topo", shift_days=0.0):
    """Share of the exposure each position (deg) spends in shadow (0 sunlit, 1 dark throughout)."""
    if not len(np.atleast_1d(ra)):
        return np.zeros(0)
    return sh.dark_dirs(unit(ra, dec), plate(mjd_mid, exp_s, shift_days), FRAMES[frame])


POSITION = ("gaia_bright", "gaia_any", "blend", "halo", "fixed_vector", "scanner_mark", "validity", "asteroid", "ps1",
            "variable", "meteor", "galaxy", "sn_host")       # the vetoes that depend only on where a detection lies


def midpoint(ra, dec, rb, db):
    """The midpoint of two positions (deg)."""
    u = unit(ra, dec) + unit(rb, db); u /= np.linalg.norm(u, axis=1)[:, None]
    return np.degrees(np.arctan2(u[:, 1], u[:, 0])) % 360, np.degrees(np.arcsin(u[:, 2]))


def replica_shifts(n, npair, seed=0):
    """Whole-day moves to other nights, avoiding the month around the true date (the shadow moves ~1 deg per day): an
    independent move for every pair in every replica (n x npair)."""
    lo, hi = A.checks.replica_min_days, A.checks.replica_max_days
    k = np.r_[np.arange(-hi, -lo + 1), np.arange(lo, hi + 1)]
    return np.random.default_rng(seed).choice(k, (n, npair)).astype(float)

"""The twin flux-ratio reference for the twin-ratio cut (selection.select). A source's ratio is log10 of b's aperture significance over a's, measured at a's detection (b forced there)
and at b's detection (a forced there) and averaged, so neither plate is favoured by being the one the source was found
on. It is compared with Gaia stars selected as candidates are (own and forced significance >= Z_MIN on both twins, detected on both)
of the same geometric-mean significance (within 0.05 dex) within 816" (800 px; each camera's vignetting differs), or
over the pair when fewer than 15 are there. Compared with all stars instead, any source near threshold sits above the
stars' ratio: b is the shallower twin in most pairs, so only sources that came out bright on b pass (twin_asymmetry.py)."""
import numpy as np
from scipy.spatial import cKDTree

from heidelberg.config import A

Z_MIN = A.selection.z_min                     # aperture significance required on both twins
MARGIN = A.selection.ratio_margin_dex         # dex allowed between a source's twin ratio and stars'
K, FLOOR = A.selection.ratio_k, A.selection.ratio_floor_dex          # noise multiple; dex scatter of stars' ratio at high z
NEAR, BIN = A.selection.ratio_near_arcsec, A.selection.ratio_bin_dex
N_LOCAL, N_MIN, N_REF = A.selection.ratio_n_local, A.selection.ratio_n_min, A.selection.ratio_n_ref
PAIR_R, GAIA_SEP, Z_REF_MAX = A.pairing.pair_radius_arcsec, A.selection.ref_gaia_arcsec, A.selection.ref_zmax


def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec); return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


class TwinRatio:
    def __init__(self, za, zb):
        sep, jb = cKDTree(unit(zb["ra"], zb["dec"])).query(unit(za["ra"], za["dec"])); sep = np.degrees(sep) * 3600
        zs, zt = za["zself"], np.nan_to_num(za["ztwin"]); Bs, Bt = zb["zself"][jb], np.nan_to_num(zb["ztwin"][jb])
        ok = (sep < PAIR_R) & (zs > 0) & (zt > 0) & (Bs > 0) & (Bt > 0)
        one = lambda x: np.where(ok, x, 1.0)
        d = np.where(ok, 0.5 * (np.log10(one(zt) / one(zs)) + np.log10(one(Bs) / one(Bt))), np.nan)
        # reference stars chosen with a and b alike: each other's nearest detection, matched to Gaia
        # on both plates, inside both boxes, own and forced >= Z_MIN on both (the candidates' rule), placed at the midpoint
        _, ia_ = cKDTree(unit(za["ra"], za["dec"])).query(unit(zb["ra"], zb["dec"])); mutual = ia_[jb] == np.arange(len(jb))
        k = np.flatnonzero(za["inside"] & zb["inside"][jb] & za["gaia"] & (za["gsep"] < GAIA_SEP) & zb["gaia"][jb] & (zb["gsep"][jb] < GAIA_SEP) & mutual
                           & (zs >= Z_MIN) & (zt >= Z_MIN) & (Bs >= Z_MIN) & (Bt >= Z_MIN) & (zs < Z_REF_MAX) & (Bs < Z_REF_MAX) & ok)
        self.lz, self.v = 0.5 * np.log10(zs[k] * Bs[k]), d[k]
        u = unit(za["ra"][k], za["dec"][k]) + unit(zb["ra"][jb[k]], zb["dec"][jb[k]])
        self.tree = cKDTree(u / np.linalg.norm(u, axis=1)[:, None] if len(k) else np.zeros((0, 3))); self.n = len(k)

    def offset(self, ra, dec, d, z1, z2):
        """A source's ratio minus stars' (dex) and its noise (dex); nan where undefined."""
        ra, dec, d, z1, z2 = (np.atleast_1d(np.asarray(q, float)) for q in (ra, dec, d, z1, z2))
        dev, s = np.full(len(ra), np.nan), 0.434 * np.sqrt(1 / z1 ** 2 + 1 / z2 ** 2)
        if self.n < N_REF:
            return dev, s
        U = unit(np.nan_to_num(ra), np.nan_to_num(dec))
        for n in range(len(ra)):
            if not (np.isfinite(ra[n]) and np.isfinite(d[n]) and z1[n] > 0 and z2[n] > 0):
                continue
            l = 0.5 * np.log10(z1[n] * z2[n]); nb = np.asarray(self.tree.query_ball_point(U[n], np.radians(NEAR / 3600)), int)
            nb = nb[np.abs(self.lz[nb] - l) < BIN]
            if len(nb) < N_LOCAL:
                nb = np.flatnonzero(np.abs(self.lz - l) < BIN)
            if len(nb) >= N_MIN:
                dev[n] = d[n] - np.median(self.v[nb])
        return dev, s

    def pull(self, ra, dec, d, z1, z2):
        dev, s = self.offset(ra, dec, d, z1, z2); return dev / np.hypot(s, FLOOR)

    def ok(self, ra, dec, d, z1, z2):
        """True unless the ratio differs from stars' beyond MARGIN + K x noise (undefined ratios pass, as before)."""
        dev, s = self.offset(ra, dec, d, z1, z2)
        return ~(np.abs(dev) > MARGIN + K * np.hypot(s, FLOOR))

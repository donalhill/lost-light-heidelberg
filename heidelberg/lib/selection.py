"""Selection of Heidelberg twin-plate sources: the catalogue cross-matches (cached, frozen for the record run), the vetoes,
the twin agreement and twin ratio tests, location validity, and the one selection function, select(), that judges every
population (candidates, chance pairs, injected glints) given as a pair table. No side effects at import; settings in
CFG, set by configure()."""
import os
import warnings

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from heidelberg.config import A
from heidelberg.lib.twin_ratio import TwinRatio

warnings.filterwarnings("ignore")
CFG = dict(source="primary", xcache="xmatch_cache", asteroids="asteroids", scanner_marks="scanner_marks.csv", offline=False, rule="symmetric", snapshots=None)
USED = set()                                  # cache files read or written (their hashes go in the stage manifest)


def configure(**kw):
    unknown = set(kw) - set(CFG)
    if unknown:
        raise KeyError(unknown)
    CFG.update(kw)
    from heidelberg.lib import catalogues
    catalogues.SNAP["dir"] = CFG["snapshots"]


def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec); return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]
def xmatch(ra, dec, cat, radius):
    """The cross-match, cached on disk by (source, catalogue, radius, positions to 1e-7 deg): a rerun of the same
    selection reads the same answer, offline; a changed position list is a new query. With CFG["offline"] a query not
    already cached is an error (the record run uses frozen answers only); every cache file read is listed in USED."""
    import hashlib
    from astropy.table import Table
    if not len(ra) or CFG["source"] == "none":             # "none": synthetic tests only (no catalogue has these fields)
        return None
    if CFG["source"] == "primary":
        from heidelberg.lib import catalogues
        if cat in catalogues.LOCAL:                          # frozen local snapshots (hashed in their manifest): matched here, not cached
            USED.update(catalogues.snapshot_files()); return catalogues.xmatch(ra, dec, cat, radius)
    h = hashlib.sha256(np.round(np.c_[np.asarray(ra, float), np.asarray(dec, float)], 7).tobytes()).hexdigest()[:24]
    f = f"{CFG['xcache']}/{CFG['source']}_{cat.replace('/', '_')}_{radius:g}_{h}.ecsv"
    if os.path.exists(f):
        USED.add(f); t = Table.read(f, format="ascii.ecsv"); return t if len(t) else None
    if CFG["offline"]:
        raise RuntimeError(f"offline: {cat} at {len(ra)} positions not in the frozen cache")
    from heidelberg.lib.catalogues import xmatch as live
    m = live(ra, dec, cat, radius)
    os.makedirs(CFG["xcache"], exist_ok=True); USED.add(f)
    (m if m is not None else Table({"i": np.zeros(0, int)})).write(f + ".tmp", format="ascii.ecsv", overwrite=True); os.replace(f + ".tmp", f)
    return m


def gaia_near(ra, dec, epoch, rmax):
    """For each position, every Gaia DR3 star within rmax (arcsec) moved to the plate epoch: (index, separation, G)."""
    m = xmatch(ra, dec, "I/355/gaiadr3", rmax)
    if m is None or not len(m):
        return np.array([], int), np.array([]), np.array([])
    dt = epoch - A.physics.gaia_epoch_jyear; i = np.asarray(m["i"], int)
    gd = np.asarray(m["DEdeg"], float) + np.nan_to_num(np.asarray(m["pmDE"], float)) * dt / 3.6e6
    gr = np.asarray(m["RAdeg"], float) + np.nan_to_num(np.asarray(m["pmRA"], float)) * dt / 3.6e6 / np.cos(np.radians(gd))
    sep = np.degrees(np.linalg.norm(unit(gr, gd) - unit(ra[i], dec[i]), axis=1)) * 3600
    return i, sep, np.asarray(m["Gmag"], float)


def gaia_local(ra, dec, g, rmax=None):
    """Gaia stars brighter than G 17.5 within rmax (arcsec; the widest star veto by default) of each position, from the pair's own Gaia list at the plate
    epoch (the _gaia.npz the pair's run saved): (index, separation, G). Used for the halo and ghost veto and for blends,
    which only stars bright enough to make or mimic a plate source can cause."""
    rmax = max(BRIGHT_R, GHOST_R[1], max(RG[1])) if rmax is None else rmax
    T = cKDTree(unit(g["ra"], g["dec"])); P = unit(ra, dec)
    lst = T.query_ball_point(P, np.radians(rmax / 3600))
    gi = np.array([i for i, js in enumerate(lst) for _ in js], int); gj = np.array([j for js in lst for j in js], int)
    if not len(gi):
        return gi, np.array([]), np.array([])
    sep = np.degrees(np.linalg.norm(P[gi] - unit(g["ra"][gj], g["dec"][gj]), axis=1)) * 3600
    return gi, sep, g["G"][gj]


def any_within(ra, dec, cat, r, mag=None, lim=None):
    """Whether each position has a source of cat within r arcsec (with mag column brighter than that position's lim, when
    given: a static object bright enough to have made the plate source)."""
    m = xmatch(ra, dec, cat, r); hit = np.zeros(len(ra), bool)
    if m is not None and len(m):
        i = np.asarray(m["i"], int)
        mv = None if mag is None else np.nanmin(np.c_[[np.nan_to_num(np.asarray(m[c], float), nan=99) for c in np.atleast_1d(mag)]], 0)
        keep = np.ones(len(i), bool) if mag is None else mv < np.asarray(lim)[i]
        hit[i[keep]] = True
    return hit


RG = (A.selection.halo_g, A.selection.halo_r)          # star G -> veto radius (arcsec)
GHOST_G, GHOST_R, BRIGHT_G, BRIGHT_R = A.selection.ghost_g, A.selection.ghost_r, A.selection.bright_g, A.selection.bright_r
STREAK_MAX = A.selection.streak_max_arcsec    # arcsec: the longest glint streak accepted (about 2 s at geostationary radius)
Z_MIN = A.selection.z_min                     # aperture significance required on both twins (chance pairs of noise below)


def star_veto(n, gi, gs, gg):
    """Positions near a Gaia star's halo, spikes, or ghost (measured on the tuning pairs: the excess of two-plate sources
    around stars by magnitude): within RG's radius for the star's G (none fainter than the table), in the ghost ring
    GHOST_R for stars brighter than GHOST_G, or within BRIGHT_R of a star brighter than BRIGHT_G."""
    hit = (gs < np.interp(gg, *RG, right=0.0)) | ((gg < GHOST_G) & (gs > GHOST_R[0]) & (gs < GHOST_R[1])) | ((gg < BRIGHT_G) & (gs < BRIGHT_R))
    v = np.zeros(n, bool); v[gi[hit]] = True
    return v


AST_R = A.selection.ast_arcsec                # arcsec around a known asteroid's track (orbits propagated back to 1900-57)


def asteroid_veto(plate, ra, dec, vmax=A.selection.ast_vmax):
    """Positions within AST_R of a known asteroid's track over the exposure (start to end; our ASSIST ephemeris from
    Koenigstuhl, ASTEROIDS/*.parquet with plate, when, desig, ra, dec, V), for asteroids brighter than vmax. None when
    no ephemeris covers the plate (reported, so a missing ephemeris cannot pass silently)."""
    import glob as _g, os as _o
    d_ = CFG["asteroids"]
    fs = _g.glob(f"{d_}/*.parquet")
    if not fs:
        return None
    Eph = pd.concat([pd.read_parquet(f) for f in fs]); Eph = Eph[(Eph.plate == plate) & (Eph.V < vmax)]
    if Eph.empty:                                     # computed and empty, or never computed: ASTEROIDS/plates.csv lists the computed plates
        pl = f"{d_}/plates.csv"
        covered = set(pd.read_csv(pl).plate) if _o.path.exists(pl) else set(pd.concat([pd.read_parquet(f, columns=["plate"]) for f in fs]).plate)
        return np.zeros(len(ra), bool) if plate in covered else None
    v = np.zeros(len(ra), bool); P = unit(ra, dec)
    for _, g in Eph.groupby("desig"):
        g = g.set_index("when")
        if not {"start", "end"} <= set(g.index):
            continue
        s0, e0 = unit(g.loc["start", "ra"], g.loc["start", "dec"])[0], unit(g.loc["end", "ra"], g.loc["end", "dec"])[0]
        seg = e0 - s0; t = np.clip(((P - s0) @ seg) / max(seg @ seg, 1e-30), 0, 1)
        dist = np.degrees(np.linalg.norm(P - (s0 + t[:, None] * seg), axis=1)) * 3600
        v |= dist < AST_R
    return v


def fixed_vectors(za, g):
    """Fixed-offset copies on this pair (a double exposure or guiding jump copies every star at one vector on both
    twins): two-plate sources with no Gaia star (significance at least vector_zmin), stacked by offset (east, north;
    arcsec) from Gaia stars brighter than vector_star_gmax within vector_rmax_arcsec, against the same offsets rotated at
    random (radial profile kept); bins of vector_bin_arcsec more than vector_cut_sigma standard deviations above it.
    Returns the vectors (n x 2, arcsec)."""
    S = A.selection; zmin, rmax, cut, bw = S.vector_zmin, S.vector_rmax_arcsec, S.vector_cut_sigma, S.vector_bin_arcsec
    rng = np.random.default_rng(S.vector_seed); src = za["inside"] & za["twin"] & ~za["gaia"] & (za["zself"] >= zmin)
    st = g["inbox"] & (g["G"] < S.vector_star_gmax); c0 = np.median(za["dec"][za["inside"]])
    S_ = np.c_[za["ra"][src] * np.cos(np.radians(c0)), za["dec"][src]] * 3600; G_ = np.c_[g["ra"][st] * np.cos(np.radians(c0)), g["dec"][st]] * 3600
    pr = cKDTree(G_).query_ball_point(S_, rmax)
    off = np.array([S_[i] - G_[j] for i, js in enumerate(pr) for j in js]) if any(len(x) for x in pr) else np.zeros((0, 2))
    if len(off) < S.vector_min_offsets:
        return np.zeros((0, 2))
    r = np.hypot(*off.T); rng_ = [-rmax, rmax]
    H, xe, ye = np.histogram2d(off[:, 0], off[:, 1], bins=int(2 * rmax / bw), range=[rng_, rng_])
    Hr = np.mean([np.histogram2d(r * np.cos(t), r * np.sin(t), bins=int(2 * rmax / bw), range=[rng_, rng_])[0]
                  for t in [rng.uniform(0, 2 * np.pi, len(r)) for _ in range(S.vector_rotations)]], 0)
    ks = np.argwhere((H - Hr) / np.sqrt(Hr + 1) > cut)
    return np.array([[(xe[a] + xe[a + 1]) / 2, (ye[b] + ye[b + 1]) / 2] for a, b in ks]) if len(ks) else np.zeros((0, 2))


def fixed_vectors_both(za, zb, g):
    """Fixed vectors found on either plate's two-plate sources (a mount movement copies every star on both twins)."""
    V = [v for v in list(fixed_vectors(za, g)) + list(fixed_vectors(zb, g))]
    out = []
    for v in V:                                    # one entry per vector (the two plates find the same one within a bin)
        if all(np.hypot(*(np.asarray(v) - np.asarray(w))) > A.selection.vector_bin_arcsec for w in out):
            out.append(v)
    return out


def vector_veto(ra, dec, g, V, r=A.selection.vector_arcsec):
    """Positions within r of (a Gaia star brighter than G 17.5 + one of the fixed vectors V)."""
    if not len(V) or not len(ra):
        return np.zeros(len(ra), bool)
    c0 = np.median(dec); G_ = np.c_[g["ra"] * np.cos(np.radians(c0)), g["dec"]] * 3600
    T = cKDTree(np.concatenate([G_ + v for v in V])); P = np.c_[np.asarray(ra) * np.cos(np.radians(c0)), np.asarray(dec)] * 3600
    return T.query(P)[0] < r


SHOWERS = [tuple(tuple(x) if isinstance(x, (list, tuple)) else x for x in sh) for sh in A.selection.showers]   # IMO calendar (config)
RADIANT_R = A.selection.radiant_deg          # deg around an active radiant (a meteor seen head-on is a brief point)


def meteor_veto(ra, dec, epoch):
    """Positions within RADIANT_R of the radiant of a major shower active on the plate's date (epoch in Julian years)."""
    import datetime as _dt
    from astropy.time import Time
    d0 = Time(epoch, format="jyear").datetime.date(); y = d0.year; v = np.zeros(len(ra), bool)   # the plate's calendar date
    for name, a, b, pk, rr, dd, dr, ddc in SHOWERS:
        for yy in (y - 1, y, y + 1):
            try:
                s0, e0, p0 = _dt.date(yy, *a), _dt.date(yy if b >= a else yy + 1, *b), _dt.date(yy if pk >= a else yy + 1, *pk)
            except ValueError:
                continue
            if s0 <= d0 <= e0:
                n = (d0 - p0).days; r0, de0 = (rr + dr * n) % 360, dd + ddc * n
                v |= np.degrees(np.linalg.norm(unit(ra, dec) - unit(np.full(len(ra), r0), np.full(len(ra), de0)), axis=1)) < RADIANT_R
    return v


VARCATS = (("B/vsx/vsx", A.selection.var_arcsec), ("B/cb/cbdata", A.selection.var_arcsec), ("B/sn/sncat", A.selection.sn_arcsec))   # variables and novae, cataclysmic binaries, supernovae (historical positions are rough)


def variable_veto(ra, dec):
    """Positions with a known variable star, nova, cataclysmic binary, or supernova within the catalogue's radius."""
    hit = np.zeros(len(ra), bool)
    for cat, r in VARCATS:
        hit |= any_within(ra, dec, cat, r)
    return hit



SN_KMAX, SN_RFAC, SN_RMIN = A.selection.host_kmax, A.selection.host_rfac, A.selection.host_rmin_arcsec   # hosts: 2MASS K < SN_KMAX; radius


def sn_hosts():
    """Galaxies that could host a supernova reaching the plates' limit: 2MASS extended sources with K < SN_KMAX and the
    2MASS Large Galaxy Atlas (the largest nearby galaxies, absent from the extended-source list), from the frozen all-sky
    lists (heidelberg/stages/snapshots.py). Returns (ra, dec, veto radius in arcsec = max(SN_RMIN, SN_RFAC x K-band
    isophotal radius))."""
    if CFG["source"] == "none":                           # synthetic tests only
        return np.zeros(0), np.zeros(0), np.zeros(0)
    from astropy.table import Table, vstack
    from heidelberg.lib.catalogues import host_snapshots
    fs = [f"{CFG['snapshots']}/{f}" for f in host_snapshots(SN_KMAX)]
    for f in fs:
        if not os.path.exists(f):
            raise FileNotFoundError(f"{f} (heidelberg/stages/snapshots.py)")
    USED.update(fs); T = vstack([Table.read(f, format="ascii.ecsv") for f in fs]); r = np.nan_to_num(np.asarray(T["r_k20fe"], float), nan=0.0)
    return np.asarray(T["ra"], float), np.asarray(T["dec"], float), np.maximum(SN_RMIN, SN_RFAC * r)


def sn_host_veto(ra, dec, hosts):
    """Positions within a host galaxy's veto radius (supernovae, which history did not always catalogue)."""
    hra, hde, hr = hosts
    if not len(ra) or not len(hra):
        return np.zeros(len(ra), bool)
    t = cKDTree(unit(hra, hde)); H = unit(hra, hde); P = unit(ra, np.asarray(dec))
    nb = t.query_ball_point(P, np.radians(hr.max() / 3600))
    i = np.repeat(np.arange(len(nb)), [len(x) for x in nb]); j = np.fromiter((x for js in nb for x in js), int, len(i))
    hit = np.zeros(len(ra), bool)
    if len(i):
        sep = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", P[i], H[j]), -1, 1))) * 3600
        hit[i[sep < hr[j]]] = True
    return hit

REMOTE = ("gaia_any", "ps1", "variable", "galaxy")   # the cuts answered by remote catalogues (xmatch)
RING_MIN = A.selection.ring_min              # dark-ring significance below which a source is a development defect
ASTRO_MIN, ASTRO_K = A.selection.validity_share, A.selection.validity_k   # share of the nearest bright detections matched to Gaia below which a location is invalid


def astro_ok(za, zb):
    """Where the plate solutions and the Gaia list both hold: at least ASTRO_MIN of the ASTRO_K nearest bright detections
    (significance above validity_zmin) on each twin are matched to Gaia. Elsewhere real stars pass as sources with no Gaia counterpart (a solution
    good over part of the plate; a Gaia query truncated in a dense field). Returns ok(ra, dec) -> bool, applied to
    candidates and injections alike."""
    refs = []
    for z in (za, zb):
        br = np.flatnonzero(z["inside"] & (z["zself"] > A.selection.validity_zmin)); refs.append((cKDTree(unit(z["ra"][br], z["dec"][br])), z["gaia"][br], len(br)))

    def ok(ra, dec):
        out = np.ones(len(ra), bool); U = unit(np.nan_to_num(np.asarray(ra, float)), np.nan_to_num(np.asarray(dec, float)))
        for t, gm, n in refs:
            if n < ASTRO_K:
                return np.zeros(len(ra), bool)
            out &= gm[t.query(U, k=ASTRO_K)[1]].mean(1) >= ASTRO_MIN
        return out
    return ok


def ps1_galaxy(ra, dec, r=A.selection.galaxy_arcsec):
    """Positions with an extended Pan-STARRS source (Kron g at least galaxy_dmag brighter than PSF g, Kron g below
    galaxy_kron_g) within r arcsec, at any brightness. Read from the ps1 cut's match at the same positions when it
    carries separations (angDist), so Pan-STARRS is asked once."""
    R5 = A.selection.ps1_arcsec
    m = xmatch(ra, dec, "II/349/ps1", R5) if r <= R5 else None
    if m is not None and len(m) and "angDist" in m.colnames:
        m = m[np.asarray(m["angDist"], float) <= r]
    elif r > R5 or (m is not None and len(m)):            # no separations in this answer (an older cache entry): ask at r
        m = xmatch(ra, dec, "II/349/ps1", r)
    hit = np.zeros(len(ra), bool)
    if m is not None and len(m):
        g, gk = np.nan_to_num(np.asarray(m["gmag"], float), nan=99), np.nan_to_num(np.asarray(m["gKmag"], float), nan=99)
        k = (g - gk > A.selection.galaxy_dmag) & (gk < A.selection.galaxy_kron_g); hit[np.asarray(m["i"], int)[k]] = True
    return hit


def g_equiv(g, tag):
    """The plate's flux -> Gaia G of stars (median G per geq_bin_dex of log flux, Gaia stars detected on plate tag): a
    function."""
    bw = A.selection.geq_bin_dex; nb = round(1 / bw); f, G_ = g[f"flux_{tag}"], g["G"]; ok = np.isfinite(f) & (f > 0)
    lf = np.log10(f[ok]); b = np.arange(np.floor(lf.min() * nb) / nb, lf.max() + bw, bw)          # bins from a multiple of bw
    mid, med = [], []
    for lo in b:
        k = (lf >= lo) & (lf < lo + bw)
        if k.sum() >= A.selection.geq_min_stars:
            mid.append(lo + bw / 2); med.append(np.median(G_[ok][k]))
    return lambda flux: np.interp(np.log10(np.clip(flux, 1, None)), mid, np.minimum.accumulate(med))


def consistency(za, zb):
    """How alike the two images of a source are, as the two images of a star are. Each image is first described against
    the stars around it on its own plate (stars of that flux within 500 px: each camera's shape changes across the
    field), then the two descriptions are differenced, so the features measure only disagreement between the twins:
    twin separation over the median at that flux; log flux ratio b / a minus its median at that flux; the difference
    between the twins of log width, of log axis ratio, and of orientation (each relative to the local stars; the
    orientation term weighted by the mean elongation). A real source of any shape (sharper than stars, or streaked along
    RA, as a sub-second glint from geostationary radius is) looks the same on both twins; a chance pair of defects does
    not. Calibrated on Gaia stars detected on both twins: robust centre, covariance, and the 99th percentile of their
    Mahalanobis distance. Also returns the 99th percentiles, on each plate, of log minor axis relative to local stars:
    the one-sided "not broader than a point source" test, which a streak passes and an extended object fails."""
    # calibration stars chosen with a and b alike, so that exchanging the plates changes no verdict:
    # each other's nearest detection within agree_pair_arcsec, Gaia on both, inside both; up to agree_max_stars by a hash
    # of the pair's midpoint (the same stars whichever plate is called a); each star's flux the mean over the twins
    S = A.selection
    tb_ = cKDTree(unit(zb["ra"], zb["dec"])); d, j = tb_.query(unit(za["ra"], za["dec"])); d = np.degrees(d) * 3600
    _, i_ = cKDTree(unit(za["ra"], za["dec"])).query(unit(zb["ra"], zb["dec"])); mutual = i_[j] == np.arange(len(j))
    k = np.flatnonzero(za["inside"] & zb["inside"][j] & za["gaia"] & (d < S.agree_pair_arcsec) & zb["gaia"][j] & mutual)
    um = unit(za["ra"][k], za["dec"][k]) + unit(zb["ra"][j[k]], zb["dec"][j[k]]); hk = np.round(um / np.linalg.norm(um, axis=1)[:, None] * 1e6).astype(np.int64)
    k = k[np.argsort((hk[:, 0] * 73856093 ^ hk[:, 1] * 19349663 ^ hk[:, 2] * 83492791) % 1000003, kind="stable")[:S.agree_max_stars]]
    la, lb = np.log10(np.clip(za["flux"], 1, None)), np.log10(np.clip(zb["flux"], 1, None))
    has_t = "theta" in za.files and "theta" in zb.files
    def plate_ref(z, l):
        s_ = np.flatnonzero(z["inside"] & z["gaia"]); q = z["a"] / np.maximum(z["b"], 1e-3)
        return dict(t=cKDTree(np.c_[z["x"][s_], z["y"][s_]]), l=l[s_], w=z["fwhm"][s_], q=q[s_], bm=z["b"][s_], am=z["a"][s_],
                    th=np.radians(z["theta"][s_]) if has_t else np.zeros(len(s_)))
    RA_, RB_ = plate_ref(za, la), plate_ref(zb, lb)
    def local(Rf, x, y, lf):
        for rad in (*S.agree_px, 1e6):                        # 1e6 px: the whole plate
            if np.isfinite(x) and np.isfinite(y) and rad < 1e6:
                nb = np.asarray(Rf["t"].query_ball_point([x, y], rad), int)
            else:                                               # no position (older files) or too few nearby: all stars
                nb = np.arange(len(Rf["l"]))
            nb = nb[np.abs(Rf["l"][nb] - lf) < S.agree_flux_dex] if len(nb) else nb
            if len(nb) >= S.agree_min_stars:
                th = Rf["th"][nb]; mth = 0.5 * np.arctan2(np.median(np.sin(2 * th)), np.median(np.cos(2 * th)))
                return np.median(Rf["w"][nb]), np.median(Rf["q"][nb]), mth, np.median(Rf["bm"][nb]), np.median(Rf["am"][nb])
        return np.nan, np.nan, 0.0, np.nan, np.nan
    lfk = 0.5 * (la[k] + lb[j[k]]); dk = d[k]; rk = lb[j[k]] - la[k]
    def feats(fa, fb, wa, wb, q_a, q_b, sep, xa, ya, xb, yb, tha=None, thb=None, ba=None, bb=None, aa=None, ab_=None):
        out, mins = [], []
        for i in range(len(fa)):
            l1, l2 = np.log10(max(fa[i], 1)), np.log10(max(fb[i], 1)); m = np.abs(lfk - 0.5 * (l1 + l2)) < S.agree_feat_dex
            if m.sum() < S.agree_feat_min_stars:
                out.append([np.nan] * 5); mins.append([np.nan] * 4); continue
            w1, q1, t1, b1, a1 = local(RA_, xa[i], ya[i], l1); w2, q2, t2, b2, a2 = local(RB_, xb[i], yb[i], l2)
            dth = 0.0 if (tha is None or thb is None) else np.sin(2 * ((np.radians(tha[i]) - t1) - (np.radians(thb[i]) - t2))) * (0.5 * (q_a[i] + q_b[i]) - 1)
            out.append([sep[i] / np.median(dk[m]), l2 - l1 - np.median(rk[m]), np.log(wa[i] / w1) - np.log(wb[i] / w2),
                        np.log(q_a[i] / q1) - np.log(q_b[i] / q2), dth])
            smear = (STREAK_MAX / A.plates.pixel_arcsec) ** 2 / 12   # the longest glint streak considered, as a variance (px^2; a uniform streak's variance is L^2/12)
            mins.append([np.log(ba[i] / b1) if ba is not None else np.nan, np.log(bb[i] / b2) if bb is not None else np.nan,
                         aa[i] / np.sqrt(a1 ** 2 + smear) if aa is not None else np.nan, ab_[i] / np.sqrt(a2 ** 2 + smear) if ab_ is not None else np.nan])
        return np.array(out), np.array(mins)
    qa_all, qb_all = za["a"] / np.maximum(za["b"], 1e-3), zb["a"] / np.maximum(zb["b"], 1e-3)
    F, M = feats(za["flux"][k], zb["flux"][j[k]], za["fwhm"][k], zb["fwhm"][j[k]], qa_all[k], qb_all[j[k]], d[k],
                 za["x"][k], za["y"][k], zb["x"][j[k]], zb["y"][j[k]],
                 za["theta"][k] if has_t else None, zb["theta"][j[k]] if has_t else None, za["b"][k], zb["b"][j[k]], za["a"][k], zb["a"][j[k]])
    ok = np.isfinite(F).all(1); F, M = F[ok], M[ok]
    if len(F) < S.agree_min_cal:                           # too few stars on both twins to calibrate: nothing passes on this pair
        return feats, (lambda X: np.full(len(X), np.inf)), 0.0, (-np.inf, -np.inf)
    use = [c for c in range(F.shape[1]) if np.std(F[:, c]) > 0]            # orientation term only when saved
    G_ = F[:, use]; mu = np.median(G_, 0); keep = np.ones(len(G_), bool)
    for _ in range(S.agree_robust_iter):                        # robust: centre and covariance from the inner agree_robust_pct
        Ci = np.linalg.inv(np.cov((G_[keep] - mu).T)); dd = np.sqrt(np.einsum("ij,jk,ik->i", G_ - mu, Ci, G_ - mu))
        keep = dd <= np.percentile(dd, S.agree_robust_pct); mu = G_[keep].mean(0)
    Ci = np.linalg.inv(np.cov((G_[keep] - mu).T))
    score = lambda X: np.sqrt(np.einsum("ij,jk,ik->i", X[:, use] - mu, Ci, X[:, use] - mu))
    lim = float(np.percentile(score(F), S.agree_pct))
    bmax = (float(np.nanpercentile(M[:, 0], S.minor_axis_pct)), float(np.nanpercentile(M[:, 1], S.minor_axis_pct)))
    return feats, score, lim, bmax


GAIA_R = A.selection.gaia_bright_arcsec       # arcsec: a detection with a G < 17.5 star this near is that star (the solve's match radius)
CUTS = ("pairing", "gaia_bright", "gaia_any", "blend", "halo", "fixed_vector", "scanner_mark", "significance", "agreement", "point_like", "twin_ratio",
        "validity", "asteroid", "ps1", "variable", "meteor", "dark_ring", "galaxy", "sn_host")
DECISION = CUTS                               # every cut decides (the vetoes are optical only: no infrared catalogue)
BLEND_FWHM, BLEND_DMAG, SCANNER_PX = A.selection.blend_fwhm, A.selection.blend_dmag, A.selection.scanner_px


ELONG_Q = A.selection.elong_pct               # point-like: elongation no larger than this percentile of the plate's own Gaia stars at that flux


def elongation_limit(z, q=ELONG_Q):
    """One plate's point-like limit (Solano et al.'s elongation criterion in tuned form): the q-th
    percentile of the major over minor axis of that plate's own reference stars (inside, matched to Gaia within
    ref_gaia_arcsec, significance from z_min to ref_zmax) in bins of elong_bin_dex in log flux, as a function of flux.
    Built from the plate alone. On the tuning set q = 98 costs 4.5% of star-like injections."""
    S = A.selection; bw = S.elong_bin_dex; nb = round(1 / bw)
    s_ = z["inside"] & z["gaia"] & (z["gsep"] < S.ref_gaia_arcsec) & (np.nan_to_num(z["zself"]) >= Z_MIN) & (np.nan_to_num(z["zself"]) < S.ref_zmax) & (z["flux"] > 0)
    lf, el = np.log10(z["flux"][s_]), z["a"][s_] / np.maximum(z["b"][s_], 1e-3); mid, qq = [], []
    if not len(lf):
        return lambda f: np.full(len(np.atleast_1d(f)), np.inf)
    for lo in np.arange(np.floor(lf.min() * nb) / nb, lf.max() + bw, bw):                 # bins from a multiple of bw
        k = (lf >= lo) & (lf < lo + bw)
        if k.sum() >= S.elong_min_stars:
            mid.append(lo + bw / 2); qq.append(np.percentile(el[k], q))
    if not mid:
        return lambda f: np.full(len(np.atleast_1d(f)), np.inf)
    return lambda f: np.interp(np.log10(np.clip(np.asarray(f, float), 1, None)), mid, qq)


def context(base, o, za, zb, g):
    """What every cut needs for one pair, built once: the agreement calibration, the twin-ratio reference, the
    location-validity map, fixed vectors, scanner marks, host galaxies, the G scale, the asteroid ephemeris."""
    pair = o["plates"][0].split("/")[-1][:-5]
    feats, score, CLIM, WMAX = consistency(za, zb)
    SMf = CFG["scanner_marks"]
    if SMf not in ("", None, "none") and not os.path.exists(SMf):       # off only when asked ("none"), never by a missing file
        raise FileNotFoundError(f"scanner marks {SMf} (heidelberg/stages/backgrounds.py PRODUCTS_DIR; HD_SCANNER_MARKS=none to switch the cut off)")
    return dict(pair=pair, epoch=o["epoch"], za=za, zb=zb, g=g, feats=feats, score=score, CLIM=CLIM, WMAX=WMAX, TR=TwinRatio(za, zb),
                AOK=astro_ok(za, zb), FV=fixed_vectors_both(za, zb, g), ELa=elongation_limit(za), ELb=elongation_limit(zb), geq=g_equiv(g, "a"), geqb=g_equiv(g, "b"),
                SM=cKDTree(pd.read_csv(SMf)[["x", "y"]].to_numpy()) if SMf not in ("", None, "none") else None,
                HOSTS=sn_hosts(), tg=cKDTree(unit(g["ra"], g["dec"])))


def select(T, C):
    """Every cut, on any population of two-plate pairs (candidates, shifted chance pairs, injected glints) given as a
    table with, per pair: ra, dec (a's detection), ra_b, dec_b, x_a, y_a, x_b, y_b, flux_a, flux_b, fwhm_a, fwhm_b,
    a_a, b_a, a_b, b_b, theta_a, theta_b, ring_a, ring_b, and the four significances: zo_a, zo_b (each plate at its own
    detection) and zf_a, zf_b (each forced at the twin's detection). One function, so no population can be judged by
    different rules. Returns a table of booleans, True where the pair fails that cut (CUTS order), and the column
    remote_evaluated. The cuts that need a remote catalogue (REMOTE) are evaluated on the rows that fail at most one of
    the others, in every population alike: a row failing two others fails whatever the catalogues say and can be no
    cut's only failure, so passes and N-1 counts are exact; only the total failing a remote cut is among evaluated rows.
    (Every row would be 23,000 positions per pair, most of them Gaia stars, against per-position archive calls.)"""
    n = len(T); F = pd.DataFrame(False, index=T.index, columns=list(CUTS))
    if "population" in T:                         # injections paired as candidates are (candidates and chance pairs are by construction)
        inj = (T["population"] == "injection").to_numpy()
        F.loc[inj, "pairing"] = ~injection_pairing(T.ra[inj], T.dec[inj], T.ra_bo[inj], T.dec_bo[inj], C["za"], C["zb"])
    if not n:
        return F.assign(remote_evaluated=False)
    ra, de, rb, db = (T[c].to_numpy(float) for c in ("ra", "dec", "ra_b", "dec_b"))
    # every veto that depends on where a detection lies is judged at BOTH detections' own positions, each with its own
    # plate's brightness and width, and fails the pair if either fails (judged at a's position only, a shifted chance
    # pair's b detection, 20" away, would go unvetted, and stars and halo debris recorded on b alone would enter the
    # chance estimate). For candidates and injections the two positions agree to
    # within the pairing radius. ra_bo, dec_bo: b's detection at its own position (ra_b, dec_b when absent).
    rbo, dbo = (T[c].to_numpy(float) for c in ("ra_bo", "dec_bo")) if "ra_bo" in T else (rb, db)
    Gp = C["geq"](T.flux_a.to_numpy(float)); Gpb = C["geqb"](T.flux_b.to_numpy(float))
    sides = ((ra, de, Gp, T.fwhm_a.to_numpy(float)), (rbo, dbo, Gpb, T.fwhm_b.to_numpy(float)))
    for r_, d_, G_, fw in sides:
        F["gaia_bright"] |= np.degrees(C["tg"].query(unit(r_, d_))[0]) * 3600 < GAIA_R      # a G < 17.5 star: the source is that star
        gi, gs, gg = gaia_local(r_, d_, C["g"])
        nbl = np.zeros(n, int); np.add.at(nbl, gi, (gs < BLEND_FWHM * fw[gi]) & (gg < G_[gi] + BLEND_DMAG)); F["blend"] |= nbl >= 2   # two or more stars
        F["halo"] |= star_veto(n, gi, gs, gg); F["fixed_vector"] |= vector_veto(r_, d_, C["g"], C["FV"])
    if C["SM"] is not None:
        F["scanner_mark"] = (C["SM"].query(T[["x_a", "y_a"]].to_numpy(float))[0] < SCANNER_PX) | (C["SM"].query(T[["x_b", "y_b"]].to_numpy(float))[0] < SCANNER_PX)
    zo_a, zo_b, zf_a, zf_b = (np.nan_to_num(T[c].to_numpy(float)) for c in ("zo_a", "zo_b", "zf_a", "zf_b"))
    F["significance"] = ~((zo_a >= Z_MIN) & (zo_b >= Z_MIN) & (zf_a >= Z_MIN) & (zf_b >= Z_MIN))   # own and forced, both plates
    sep = np.degrees(np.linalg.norm(unit(ra, de) - unit(rb, db), axis=1)) * 3600
    X, Mn = C["feats"](T.flux_a.to_numpy(float), T.flux_b.to_numpy(float), T.fwhm_a.to_numpy(float), T.fwhm_b.to_numpy(float),
                       T.a_a.to_numpy(float) / np.maximum(T.b_a.to_numpy(float), 1e-3), T.a_b.to_numpy(float) / np.maximum(T.b_b.to_numpy(float), 1e-3), sep,
                       T.x_a.to_numpy(float), T.y_a.to_numpy(float), T.x_b.to_numpy(float), T.y_b.to_numpy(float), T.theta_a.to_numpy(float), T.theta_b.to_numpy(float),
                       T.b_a.to_numpy(float), T.b_b.to_numpy(float), T.a_a.to_numpy(float), T.a_b.to_numpy(float))
    okx = np.isfinite(X).all(1); sc = np.full(n, np.inf); sc[okx] = C["score"](X[okx])
    F["agreement"] = ~((sc <= C["CLIM"]) & (np.nan_to_num(Mn[:, 0], nan=9) <= C["WMAX"][0]) & (np.nan_to_num(Mn[:, 1], nan=9) <= C["WMAX"][1])
                       & (np.nan_to_num(Mn[:, 2], nan=9) <= 1) & (np.nan_to_num(Mn[:, 3], nan=9) <= 1))
    if CFG["rule"] != "streaks":                    # point-like on each plate by its own stars (the streak secondary keeps streaks)
        el_a = T.a_a.to_numpy(float) / np.maximum(T.b_a.to_numpy(float), 1e-3); el_b = T.a_b.to_numpy(float) / np.maximum(T.b_b.to_numpy(float), 1e-3)
        F["point_like"] = (el_a > C["ELa"](T.flux_a.to_numpy(float))) | (el_b > C["ELb"](T.flux_b.to_numpy(float)))
    with np.errstate(divide="ignore", invalid="ignore"):      # position-symmetric twin ratio (twin_ratio.py)
        d = 0.5 * (np.log10(zf_b / zo_a) + np.log10(zo_b / zf_a))
    d = np.where((zo_a > 0) & (zo_b > 0) & (zf_a > 0) & (zf_b > 0), d, np.nan)
    um = unit(ra, de) + unit(rbo, dbo); um /= np.linalg.norm(um, axis=1)[:, None]                 # the reference looked up at the midpoint
    F["twin_ratio"] = ~C["TR"].ok(np.degrees(np.arctan2(um[:, 1], um[:, 0])) % 360, np.degrees(np.arcsin(um[:, 2])), d, zo_a, zo_b)
    for r_, d_, G_, fw in sides:
        F["validity"] |= ~C["AOK"](r_, d_)
        av = asteroid_veto(C["pair"], r_, d_)
        if av is not None:
            F["asteroid"] |= av
        F["meteor"] |= meteor_veto(r_, d_, C["epoch"]); F["sn_host"] |= sn_host_veto(r_, d_, C["HOSTS"])
    F["dark_ring"] = (np.nan_to_num(T.ring_a.to_numpy(float), nan=0) < RING_MIN) | (np.nan_to_num(T.ring_b.to_numpy(float), nan=0) < RING_MIN)
    ev = F[[c for c in CUTS if c not in REMOTE]].sum(axis=1).to_numpy() <= 1
    if "population" in T:                         # injections: every row, so each one's position vetoes are known (the shadow
        ev |= (T["population"] == "injection").to_numpy()   # stage's area controls); a remote cut depends on position only, so no verdict changes
    e = np.flatnonzero(ev)
    if len(e):
        ne = len(e); re_, dee, Ge = np.r_[ra[e], rbo[e]], np.r_[de[e], dbo[e]], np.r_[Gp[e], Gpb[e]]     # both positions, one query
        xi, xs, _ = gaia_near(re_, dee, C["epoch"], A.selection.gaia_query_arcsec); dmin = np.full(2 * ne, np.inf); np.minimum.at(dmin, xi, xs)
        either = lambda v: np.asarray(v[:ne]) | np.asarray(v[ne:])
        for c, v in (("gaia_any", ~(dmin > A.selection.gaia_any_arcsec)),
                     ("ps1", any_within(re_, dee, "II/349/ps1", A.selection.ps1_arcsec, mag=["gmag", "gKmag"], lim=Ge + A.selection.ps1_dmag)),
                     ("variable", variable_veto(re_, dee)), ("galaxy", ps1_galaxy(re_, dee))):
            F[c] = np.isin(np.arange(n), e[either(v)])
    return F.assign(remote_evaluated=ev)


def det_table(za, zb, ia, jb, zf_a=None, zf_b=None, sh_b=None):
    """The pair table for detections ia on a paired with jb on b (candidates, or chance pairs with b moved by sh_b
    degrees (dRA, dDec) and the forced significances measured at the moved positions)."""
    rbo, dbo = zb["ra"][jb].astype(float), zb["dec"][jb].astype(float); rb, db = rbo, dbo
    if sh_b is not None:
        rb, db = rb + sh_b[0], db + sh_b[1]
    return pd.DataFrame(dict(ra=za["ra"][ia], dec=za["dec"][ia], ra_b=rb, dec_b=db, ra_bo=rbo, dec_bo=dbo, x_a=za["x"][ia], y_a=za["y"][ia], x_b=zb["x"][jb], y_b=zb["y"][jb],
                             flux_a=za["flux"][ia], flux_b=zb["flux"][jb], fwhm_a=za["fwhm"][ia], fwhm_b=zb["fwhm"][jb], a_a=za["a"][ia], b_a=za["b"][ia],
                             a_b=zb["a"][jb], b_b=zb["b"][jb], theta_a=za["theta"][ia], theta_b=zb["theta"][jb], ring_a=za["ring"][ia], ring_b=zb["ring"][jb],
                             zo_a=za["zself"][ia], zo_b=zb["zself"][jb], zf_a=zb["ztwin"][jb] if zf_a is None else zf_a, zf_b=za["ztwin"][ia] if zf_b is None else zf_b))


def candidates(za, zb):
    """Two-plate sources: detections on a and b that are each other's nearest within 3" (the double scan: the pairing
    does not depend on which plate is listed first), inside both boxes. Returns (ia, jb)."""
    R = np.radians(PAIR_R / 3600); dab, j = cKDTree(unit(zb["ra"], zb["dec"])).query(unit(za["ra"], za["dec"]), distance_upper_bound=R)
    _, i = cKDTree(unit(za["ra"], za["dec"])).query(unit(zb["ra"], zb["dec"]))
    ia = np.flatnonzero(np.isfinite(dab)); jb = j[ia]; m = (i[jb] == ia) & za["inside"][ia] & zb["inside"][jb]
    return ia[m], jb[m]


def chance_tables(za, zb, radius=None):
    """Chance pairs: single-plate detections on a (no detection within 3" on b) against single-plate detections on b
    moved by radius (20" for the estimate; the further radii in zshift_x test its dependence on the lag) in eight
    directions, each other's nearest within 3" after the move, with the forced significances measured at the moved
    positions (measure_pair.py zshift, zshift_x). One table per direction."""
    radius = A.background.chance_shift_arcsec if radius is None else radius
    if radius == A.background.chance_shift_arcsec:
        if "zshift" not in za.files or "zshift" not in zb.files:
            return None
        Za, Zb = za["zshift"], zb["zshift"]
    else:
        if "zshift_x" not in za.files or "zshift_x" not in zb.files or radius not in list(za["xshift_radii"]):
            return None
        q = list(za["xshift_radii"]).index(radius); Za, Zb = za["zshift_x"][:, :, q], zb["zshift_x"][:, :, q]
    R = np.radians(PAIR_R / 3600); sa = np.flatnonzero(~za["twin"] & za["inside"]); sb = np.flatnonzero(~zb["twin"] & zb["inside"]); out = []
    for k in range(Za.shape[1]):
        t = 2 * np.pi * k / Za.shape[1]; dd = radius * np.sin(t) / 3600; dr = radius * np.cos(t) / 3600 / np.cos(np.radians(zb["dec"][sb]))
        Ub = unit(zb["ra"][sb] + dr, zb["dec"][sb] + dd); dab, j = cKDTree(Ub).query(unit(za["ra"][sa], za["dec"][sa]), distance_upper_bound=R)
        _, i = cKDTree(unit(za["ra"][sa], za["dec"][sa])).query(Ub); ia = np.flatnonzero(np.isfinite(dab)); jj = j[ia]; ia = ia[i[jj] == ia]; jj = j[ia]
        A_, B_ = sa[ia], sb[jj]
        out.append(det_table(za, zb, A_, B_, zf_a=Zb[B_, k], zf_b=Za[A_, k], sh_b=(dr[jj], np.full(len(jj), dd))))
    return out


SAME = A.pairing.same_arcsec                  # arcsec: an original detection this near an injection's is the same object
PAIR_R = A.pairing.pair_radius_arcsec         # arcsec: a two-plate source's detections are each other's nearest within this


def injection_pairing(ra_a, dec_a, ra_b, dec_b, za, zb):
    """Whether each injection (found on both plates: a detection within 3 px of where it was placed on each) is paired as
    a candidate would be (candidates(): each other's nearest within 3"): its two detections within 3" of each other, and
    neither with a nearer detection on the other plate among that plate's sources (za, zb; a source within SAME of the
    glint's own detection there is the same object, merged with the glint in the injected scan). The cut "pairing"."""
    if not len(ra_a):
        return np.zeros(0, bool)
    ua, ub = unit(np.asarray(ra_a, float), np.asarray(dec_a, float)), unit(np.asarray(ra_b, float), np.asarray(dec_b, float))
    sep = np.degrees(np.linalg.norm(ua - ub, axis=1)) * 3600
    def nearer(z, u_from, u_own):                          # a source of z nearer to u_from than the glint's twin, not the twin itself
        t = cKDTree(unit(z["ra"], z["dec"])); out = np.zeros(len(sep), bool)
        for n, js in enumerate(t.query_ball_point(u_from, np.radians(np.maximum(sep, 1e-6) / 3600))):
            if js:
                out[n] = bool((np.degrees(np.linalg.norm(unit(z["ra"][js], z["dec"][js]) - u_own[n], axis=1)) * 3600 > SAME).any())
        return out
    return (sep <= PAIR_R) & ~nearer(zb, ua, ub) & ~nearer(za, ub, ua)


def injection_table(J):
    """Injected glints found on both twins, as a pair table (own and forced significance per plate); whether each is
    paired as a candidate would be is the cut "pairing" (injection_pairing)."""
    J = J[J.found_a.fillna(False).astype(bool) & J.found_b.fillna(False).astype(bool)]
    fw = lambda a, b: 2.355 * np.sqrt((a ** 2 + b ** 2) / 2)
    return J, pd.DataFrame(dict(ra=J.ra_f_a, dec=J.dec_f_a, ra_b=J.ra_f_b, dec_b=J.dec_f_b, ra_bo=J.ra_f_b, dec_bo=J.dec_f_b, x_a=J.x_a, y_a=J.y_a, x_b=J.x_b, y_b=J.y_b, flux_a=J.flux_a, flux_b=J.flux_b,
                                fwhm_a=fw(J.a_a, J.b_a), fwhm_b=fw(J.a_b, J.b_b), a_a=J.a_a, b_a=J.b_a, a_b=J.a_b, b_b=J.b_b, theta_a=J.theta_a, theta_b=J.theta_b,
                                ring_a=J.ring_a, ring_b=J.ring_b, zo_a=J.z_a, zo_b=J.z_b,
                                # files made before the forced measurements: own-centroid z stands in (flagged in the row: inj_forced False)
                                zf_a=J.zf_a if "zf_a" in J else J.z_a, zf_b=J.zf_b if "zf_b" in J else J.z_b), index=J.index)


def efficiency_table(J, passed, base):
    """Efficiency by target G, duration, and sharpening (injected = all placed; found on both; passing every cut)."""
    J = J.assign(passed=passed.reindex(J.index).fillna(False).astype(bool), found_both=J.found_a.fillna(False).astype(bool) & J.found_b.fillna(False).astype(bool))
    g0, g1, gw = A.injections.eff_g_bins
    J["Gbin"] = pd.cut(J.G, np.arange(g0, g1 + gw / 2, gw)); keys = [c for c in ("dur", "sharpen") if c in J] + ["Gbin"]
    E = J.groupby(keys, observed=True).agg(n=("passed", "size"), found_both=("found_both", "mean"), eff=("passed", "mean")).round(3)
    E.to_csv(f"{base}_eff.csv")
    E0 = E.xs((A.injections.durations[0], 1.0), level=("dur", "sharpen")) if "dur" in J else E   # a star-like flash (no duration, not sharpened): the headline
    e = E0.eff.values; mids = np.array([iv.mid for iv in E0.index]); k = np.flatnonzero(e >= A.injections.eff_half * np.nanmax(e)) if len(e) and np.nanmax(e) > 0 else []
    return float(mids[k[-1]]) if len(k) else np.nan


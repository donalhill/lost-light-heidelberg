"""Measurement of one Heidelberg twin pair: detection, the astrometric solve against Gaia (header pre-correction, blind
search, completeness grid), forced aperture photometry, the dark ring, and the injected glints. Functions only: no side
effects at import (the pair's flow is heidelberg/stages/measure_pair.py). Every
setting is in config/analysis.json (sections detection, solve, injections, and background)."""
import os
import sys
import warnings

import astropy.units as u
import numpy as np
import pandas as pd
from astropy import log
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, convolve_fft
from scipy.ndimage import convolve1d
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SourceCatalog, detect_sources
from scipy.spatial import cKDTree

from heidelberg import config
from heidelberg.lib import meridian

warnings.filterwarnings("ignore"); log.setLevel("ERROR")
DET, SOL, INJ, PHY = config.A.detection, config.A.solve, config.A.injections, config.A.physics
SCALE = config.A.plates.pixel_arcsec
BOX = int(config.A.sample.test_area_deg * 3600 / SCALE)     # the central test area (6 x 6 deg: 21,176 px at 1.02"/px)
TILE, MARGIN = DET.tile_px, DET.tile_margin_px       # processed in tiles; sources kept in each tile's core
THRESH = DET.detect_sigma
N_THREADS = max(1, int(os.environ.get("HD_THREADS", 1)))   # threads for the background interpolation (results identical)
ORDER = SOL.solve_order                                      # polynomial order of the final plate solution (6 deg box)
RAP = DET.aperture_px
CH_SHIFT, CH_DIRS = config.A.background.chance_shift_arcsec, config.A.background.shifts   # one twin moved in each direction
CH_XRADII = config.A.background.lag_radii_arcsec             # further shift radii (zshift_x): the estimate's dependence on the lag
DURATIONS = INJ.durations                           # glint durations injected (s)
SHARPEN = (1.0, INJ.sharpen)                        # image size relative to the local stars (a flash samples an instant)
G_INJ_MIN, G_INJ_MAX = INJ.g_lo, INJ.g_hi           # injection targets stop 0.2 mag inside the Gaia list (G < 17.5): donors within 0.2


def detect(f, extra=None):
    """Sources in the central BOX px, tile by tile: background (128 px boxes), the image smoothed by a 1.5 px Gaussian,
    sources at THRESH times the smoothed image's own noise in that tile (the grain is correlated), each kept in the tile
    whose core holds its centroid. Keeps the background-subtracted image of the whole box for forced photometry. With
    extra (a box-sized layer), it is added to the scan's values before anything else, clipped at the scan's saturation
    level: injected glints then pass through exactly the detection and measurement that real sources do."""
    from astropy.stats import sigma_clipped_stats
    k = Gaussian1DKernel(DET.smooth_sigma_px, x_size=DET.smooth_size_px); k.normalize(); k1 = k.array   # the 2-D Gaussian, as two 1-D passes
    with fits.open(f, memmap=True, do_not_scale_image_data=True) as h:
        hd = h[0].header; ny, nx = hd["NAXIS2"], hd["NAXIS1"]; W = min(BOX, nx, ny)
        y0, x0 = (ny - W) // 2, (nx - W) // 2
        img = np.zeros((W, W), np.float32); rows, rmss, csds = [], [], []
        edges = np.linspace(0, W, int(np.ceil(W / TILE)) + 1).astype(int)
        for ty0, ty1 in zip(edges[:-1], edges[1:]):
            for tx0, tx1 in zip(edges[:-1], edges[1:]):
                py0, py1, px0, px1 = max(ty0 - MARGIN, 0), min(ty1 + MARGIN, W), max(tx0 - MARGIN, 0), min(tx1 + MARGIN, W)
                d = np.asarray(h[0].data[y0 + py0:y0 + py1, x0 + px0:x0 + px1], dtype=np.float32) * hd.get("BSCALE", 1.0) + hd.get("BZERO", 0.0)
                if extra is not None:
                    d = np.minimum(d + extra[py0:py1, px0:px1], saturation(hd))
                bk = Background2D(d, (DET.bkg_box_px,) * 2, filter_size=(DET.bkg_filter_boxes,) * 2, bkg_estimator=MedianBackground(), n_threads=N_THREADS)   # identical for any thread count
                r = (d - bk.background).astype(np.float32); rmss.append(float(np.median(bk.background_rms_mesh))); del d, bk
                cv = (convolve1d(convolve1d(r, k1, axis=0, mode="constant", cval=0.0), k1, axis=1, mode="constant", cval=0.0) if np.isfinite(r).all()
                      else convolve_fft(r, Gaussian2DKernel(DET.smooth_sigma_px, x_size=DET.smooth_size_px, y_size=DET.smooth_size_px), allow_huge=True)).astype(np.float32)   # NaN: astropy interpolates them
                nst = DET.noise_step_px; csd = sigma_clipped_stats(cv[::nst, ::nst], sigma=DET.clip_sigma, maxiters=DET.noise_clip_iters)[2]
                segm = detect_sources(cv, THRESH * csd, n_pixels=DET.detect_npix); csds.append((ty0, ty1, tx0, tx1, float(csd)))
                img[ty0:ty1, tx0:tx1] = r[ty0 - py0:ty1 - py0, tx0 - px0:tx1 - px0]
                if segm is None:
                    continue
                T = SourceCatalog(r, segm, convolved_data=cv).to_table(columns=["x_centroid", "y_centroid", "segment_flux", "semimajor_axis", "semiminor_axis", "max_value", "area", "orientation"])
                xc, yc = np.asarray(T["x_centroid"], float) + px0, np.asarray(T["y_centroid"], float) + py0
                keep = np.isfinite(xc) & (xc >= tx0) & (xc < tx1) & (yc >= ty0) & (yc < ty1)
                rows.append(np.c_[xc, yc, np.asarray(T["segment_flux"], float), np.asarray(T["semimajor_axis"].value, float),
                                  np.asarray(T["semiminor_axis"].value, float), np.asarray(T["max_value"], float), np.asarray(T["area"].value, float),
                                  np.asarray(T["orientation"].value, float)][keep])
                del r, cv, segm
    R = np.concatenate(rows)
    x, y, flux, a_, b_, peak, npix, theta = (R[:, i] for i in range(8))
    h2 = hd.copy(); h2["CTYPE1"], h2["CTYPE2"] = "RA---TAN-SIP", "DEC--TAN-SIP"
    w = WCS(h2)
    t = Time(hd["JD"], format="jd") if "JD" in hd else Time(hd["DATE-OBS"])
    return dict(hd=hd, x=x + x0, y=y + y0, flux=flux, a=a_, b=b_, w=w, epoch=t.jyear, shape=(ny, nx), rms=float(np.median(rmss)),
                img=img, x0=x0, y0=y0, W=W, csd=csds, theta=theta, fwhm=2.355 * np.sqrt((a_ ** 2 + b_ ** 2) / 2), npix=npix, peak=peak)


def saturation(hd):
    """The largest value the scan can hold (16-bit: 32767 x BSCALE + BZERO)."""
    bits = int(hd.get("BITPIX", -32)); top = {8: 255, 16: 32767}.get(bits)        # FITS: 8-bit unsigned, 16-bit signed
    return np.float32(top * hd.get("BSCALE", 1.0) + hd.get("BZERO", 0.0)) if top else np.float32(np.inf)


def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


def gaia_tiles(ra0, dec0, width, epoch, gmax=config.A.selection.gaia_bright_g, tile=config.C.network.gaia_tile_deg, tries=config.C.network.async_tries):
    """Gaia DR3 G < gmax over a width x width (deg) box about (ra0, dec0), in tiles of about tile deg by async jobs at
    ESA's archive, else IRSA's copy of the same table (catalogues.gaia_async; no row cap; a dense field in one query came
    back truncated from VizieR). Each tile is bounded in RA and Dec without overlap (RA wrapped), so no star is counted
    twice. Positions moved to the plate epoch by proper motion."""
    from heidelberg.lib.catalogues import gaia_async
    h = width / 2; d_edges = np.linspace(max(dec0 - h, -90), min(dec0 + h, 90), int(np.ceil(width / tile)) + 1)
    out = []
    for d1, d2 in zip(d_edges[:-1], d_edges[1:]):
        cosd = max(np.cos(np.radians(max(abs(d1), abs(d2)))), 0.05); rw = min(width / cosd, 360.0)
        r_edges = np.linspace(ra0 - rw / 2, ra0 + rw / 2, int(np.ceil(rw / tile)) + 1)
        for r1, r2 in zip(r_edges[:-1], r_edges[1:]):
            a, b = r1 % 360, r2 % 360
            rcond = f"ra >= {a} AND ra < {b}" if a < b else f"(ra >= {a} OR ra < {b})"
            out.append(gaia_async("ra, dec, pmra, pmdec, phot_g_mean_mag, bp_rp",
                                  f"phot_g_mean_mag < {gmax} AND {rcond} AND dec >= {d1} AND dec < {d2}", tries=tries))
    from astropy.table import vstack
    g = vstack(out); dt = epoch - PHY.gaia_epoch_jyear
    pmra, pmde = (np.nan_to_num(np.asarray(g[c], float)) for c in ("pmra", "pmdec"))
    dec = np.asarray(g["dec"], float) + pmde * dt / 3.6e6
    ra = np.asarray(g["ra"], float) + pmra * dt / 3.6e6 / np.cos(np.radians(dec))
    return ra, dec, np.asarray(g["phot_g_mean_mag"], float), np.asarray(np.ma.filled(np.ma.asarray(g["bp_rp"], float), np.nan), float)


def gaia(c, width, epoch):
    return gaia_tiles(c.ra.deg, c.dec.deg, width.to(u.deg).value, epoch)


def displacement_field(x, y, D, flux, GX, gm, shape, n=SOL.field_grid, half=SOL.field_half_arcsec, binw=SOL.field_bin_arcsec, nmin=SOL.field_min_peak):
    """The header's error as a smooth field: in an n x n grid of patches, the offset (Gaia minus detection, tangent plane,
    arcsec) at the peak of the histogram of all pairs between the patch's 200 brightest detections and Gaia stars
    (G < 15) near it; peaks with fewer than nmin pairs, or not 3x the next-best peak away from it, are dropped; a
    quadratic in pixels is fitted to the rest with 3-sigma clipping. Returns f(x, y) -> (n, 2) offsets, and the patches."""
    ex = np.linspace(0, shape[1], n + 1); ey = np.linspace(0, shape[0], n + 1); pts = []
    gb = np.flatnonzero(gm < SOL.field_gaia_g); tg = cKDTree(GX[gb])
    for a in range(n):
        for b in range(n):
            k = np.flatnonzero((x >= ex[a]) & (x < ex[a + 1]) & (y >= ey[b]) & (y < ey[b + 1]))
            if len(k) < SOL.field_min_det:
                continue
            k = k[np.argsort(-flux[k])[:SOL.field_n_det]]; cen = D[k].mean(0)
            gg = gb[tg.query_ball_point(cen, np.hypot(*(D[k] - cen).T).max() + half)]
            if len(gg) < SOL.field_min_gaia:
                continue
            off = (GX[gg][None] - D[k][:, None]).reshape(-1, 2); off = off[(np.abs(off) < half).all(1)]
            nb = int(2 * half / binw); H, xe, ye = np.histogram2d(off[:, 0], off[:, 1], bins=nb, range=[[-half, half], [-half, half]])
            from scipy.ndimage import uniform_filter
            Hs = uniform_filter(H, 2, mode="constant") * 4                          # 2 x 2 sums: a peak split across bins
            kk = np.unravel_index(Hs.argmax(), Hs.shape); top = Hs[kk]; mb = SOL.field_peak_mask_bins
            mask = np.ones_like(Hs, bool); mask[max(kk[0] - mb, 0):kk[0] + mb + 1, max(kk[1] - mb, 0):kk[1] + mb + 1] = False
            if top < nmin or top < SOL.field_peak_ratio * Hs[mask].max():
                continue
            c = np.array([xe[kk[0]], ye[kk[1]]])                                 # refine: median of pair offsets near the peak
            near = off[np.hypot(*(off - c).T) < SOL.field_refine_bins * binw]
            pts.append((x[k].mean(), y[k].mean(), *np.median(near, 0), top))
    P = np.array(pts)
    if len(P) < SOL.field_min_patches:
        return None, P
    B = lambda xx, yy: np.c_[np.ones(len(xx)), xx, yy, xx * xx, xx * yy, yy * yy] if np.ndim(xx) else None
    sx, sy = shape[1], shape[0]; X = B(P[:, 0] / sx - 0.5, P[:, 1] / sy - 0.5); keep = np.ones(len(P), bool)
    for _ in range(SOL.field_clip_iters):
        cx = np.linalg.lstsq(X[keep], P[keep, 2], rcond=None)[0]; cy = np.linalg.lstsq(X[keep], P[keep, 3], rcond=None)[0]
        r = np.hypot(X @ cx - P[:, 2], X @ cy - P[:, 3]); s = 1.4826 * np.median(r[keep]) + SOL.field_floor_arcsec
        keep = r < SOL.field_clip_sigma * s
    f = lambda xx, yy: np.c_[B(np.asarray(xx) / sx - 0.5, np.asarray(yy) / sy - 0.5) @ cx, B(np.asarray(xx) / sx - 0.5, np.asarray(yy) / sy - 0.5) @ cy]
    return f, np.c_[P, keep]


def similarity(D, flux, GX, gm, scales=np.arange(*SOL.sim_scales), rots=np.arange(*SOL.offset_rot_deg), half=SOL.offset_half_arcsec, binw=SOL.offset_bin_arcsec):
    """The best global scale, rotation, and shift taking header tangent-plane positions D to Gaia: the peak of the
    histogram of offsets between the 400 brightest detections and the 1,500 brightest Gaia stars, over a grid of scale and
    rotation (some headers are 0.9% off in scale: ~90" at the plate edge). Returns (scale, rot_deg, shift, peak)."""
    b = np.argsort(-flux)[:SOL.offset_n_det]; g = GX[np.argsort(gm)[:SOL.offset_n_gaia]]; best = (0, 1.0, 0.0, np.zeros(2))
    nb = int(2 * half / binw)
    for s in scales:
        for r in rots:
            c, s_ = np.cos(np.radians(r)), np.sin(np.radians(r)); Dr = s * D[b] @ np.array([[c, s_], [-s_, c]])
            off = (g[None] - Dr[:, None]).reshape(-1, 2); off = off[(np.abs(off) < half).all(1)]
            H, xe, ye = np.histogram2d(off[:, 0], off[:, 1], bins=nb, range=[[-half, half], [-half, half]])
            k = np.unravel_index(H.argmax(), H.shape)
            if H[k] > best[0]:
                best = (H[k], s, r, np.array([(xe[k[0]] + xe[k[0] + 1]) / 2, (ye[k[1]] + ye[k[1] + 1]) / 2]))
    pk, s, r, t = best; c, s_ = np.cos(np.radians(r)), np.sin(np.radians(r))
    return (lambda X: s * X @ np.array([[c, s_], [-s_, c]]) + t), dict(scale=round(float(s), 4), rot=round(float(r), 2), shift=[round(float(v), 1) for v in t], peak=int(pk))


def precorrect(P, G):
    """GAVO's header WCS can be wrong by a smooth field reaching ~200" across a plate (found on 4 of 45 tuning pairs: the
    staged fit then locks onto one corner and extrapolates). The field is measured in a 5 x 5 grid of patches and removed,
    and a TAN-SIP WCS is refitted to the corrected positions, before the staged fit. Plates with a good header are
    unchanged (offsets within the fit's first radius). Records the patches in P["precorrect"]."""
    from astropy.wcs.utils import fit_wcs_from_points
    gra, gde, gm = G; c0 = P["w"].pixel_to_world(P["shape"][1] / 2, P["shape"][0] / 2); fr = c0.skyoffset_frame()
    s = P["w"].pixel_to_world(P["x"], P["y"]).transform_to(fr); D = np.c_[s.lon.arcsec, s.lat.arcsec]
    g = SkyCoord(gra, gde, unit="deg").transform_to(fr); GX = np.c_[g.lon.arcsec, g.lat.arcsec]
    P["w_header"] = P.get("w_header", P["w"])
    sim, info = similarity(D, P["flux"], GX, gm); D = sim(D)                  # global scale, rotation, shift first
    f, pts = displacement_field(P["x"], P["y"], D, P["flux"], GX, gm, P["shape"])
    if f is None:
        P["precorrect"] = dict(applied=False, patches=int(len(pts)), similarity=info); return
    D2 = D + f(P["x"], P["y"]); back = SkyCoord(D2[:, 0] * u.arcsec, D2[:, 1] * u.arcsec, frame=fr).transform_to("icrs")
    sel = np.random.default_rng(SOL.field_fit_seed).choice(len(P["x"]), min(SOL.field_fit_n, len(P["x"])), replace=False)
    P["w"] = fit_wcs_from_points((P["x"][sel], P["y"][sel]), back[sel], proj_point="center", sip_degree=SOL.field_sip_degree)
    P["precorrect"] = dict(applied=True, patches=int(len(pts)), used=int(pts[:, -1].sum()), max_offset_arcsec=round(float(np.abs(pts[:, 2:4]).max()), 1), similarity=info)


def match_grid(P, ra, de, G, n=SOL.grid_n):
    """Astrometric completeness over the box, independent of depth: in an n x n grid on the sky, the share of Gaia stars
    from 2.5 to 1 mag brighter than the plate's 50% G (and G >= 11, below saturation) with a detection within 3" (~0.85
    on a good plate; a solution good over part of the plate, or a plate of another field, gives ~0). Returns (G50, grid)."""
    gra, gde, gm = G; inb = (P["x"] >= P["x0"]) & (P["x"] < P["x0"] + P["W"]) & (P["y"] >= P["y0"]) & (P["y"] < P["y0"] + P["W"])
    t_ = cKDTree(unit(ra[inb], de[inb])); r0, r1, d0, d1 = ra[inb].min(), ra[inb].max(), de[inb].min(), de[inb].max()
    k = (gra > r0) & (gra < r1) & (gde > d0) & (gde < d1); d = np.degrees(t_.query(unit(gra[k], gde[k]))[0]) * 3600; g = gm[k]
    g0, g1, gw = SOL.grid_g_bins; rm = SOL.grid_match_arcsec; e = np.arange(g0, g1, gw)
    fr = np.array([np.mean(d[(g >= a) & (g < a + gw)] < rm) if ((g >= a) & (g < a + gw)).sum() > SOL.grid_min_bin else np.nan for a in e])
    ok = np.isfinite(fr) & (fr >= SOL.grid_depth_share); g50 = float(e[ok].max() + gw) if ok.any() else SOL.grid_default_g50
    band = (g >= max(SOL.grid_gmin, g50 - SOL.grid_band[1])) & (g < g50 - SOL.grid_band[0]); r_, d_, hit = gra[k][band], gde[k][band], d[band] < rm
    er, ed = np.linspace(r0, r1, n + 1), np.linspace(d0, d1, n + 1); grid = []
    for b in range(n):
        row = []
        for a in range(n):
            m = (r_ >= er[a]) & (r_ <= er[a + 1]) & (d_ >= ed[b]) & (d_ <= ed[b + 1]); row.append(round(float(hit[m].mean()), 2) if m.sum() >= SOL.grid_min_cell else None)
        grid.append(row)
    return round(g50, 2), grid


def grid_ok(grid, lim=SOL.grid_share):
    v = [x for r in grid for x in r if x is not None]; return bool(v) and float(np.median(v)) >= lim


def basis(x, y, s, order=ORDER):
    """Polynomial terms x^i y^j (i + j <= order) of pixel coordinates scaled to -0.5..0.5 over shape s."""
    x, y = (x - s[1] / 2) / s[1], (y - s[0] / 2) / s[0]
    return np.column_stack([x ** (n - j) * y ** j for n in range(order + 1) for j in range(n + 1)])


def solve(P, G):
    """The plate's sky positions: header WCS shifted to Gaia, then a cubic in pixels for (xi, eta) about the field centre."""
    gra, gde, gm = G; c0 = P["w"].pixel_to_world(P["shape"][1] / 2, P["shape"][0] / 2)
    sky = P["w"].pixel_to_world(P["x"], P["y"]); ra, de = sky.ra.deg, sky.dec.deg
    # blind offset: offsets between the brightest detections and the brightest Gaia stars (tangent plane about the
    # field centre), over small rotations; the peak of the 2-D histogram, refined on a 1" grid
    tp = lambda r, d: SkyCoord(r, d, unit="deg").transform_to(c0.skyoffset_frame())
    bright = np.argsort(-P["flux"])[:SOL.offset_n_det]; sd = tp(ra[bright], de[bright]); D0 = np.c_[sd.lon.arcsec, sd.lat.arcsec]
    gb = np.argsort(gm)[:SOL.offset_n_gaia]; sg = tp(gra[gb], gde[gb]); GB = np.c_[sg.lon.arcsec, sg.lat.arcsec]
    top = (0, 0.0, 0.0, 0.0); half = SOL.offset_half_arcsec; nb = int(2 * half / SOL.offset_bin_arcsec)
    for rot in np.arange(*SOL.offset_rot_deg):
        c_, s_ = np.cos(np.radians(rot)), np.sin(np.radians(rot)); Dr = D0 @ np.array([[c_, s_], [-s_, c_]])
        off = (GB[None, :, :] - Dr[:, None, :]).reshape(-1, 2); off = off[(np.abs(off) < half).all(1)]
        H, xe, ye = np.histogram2d(off[:, 0], off[:, 1], bins=nb, range=[[-half, half], [-half, half]]); k = np.unravel_index(H.argmax(), H.shape)
        if H.max() > top[0]:
            top = (H.max(), rot, (xe[k[0]] + xe[k[0] + 1]) / 2, (ye[k[1]] + ye[k[1] + 1]) / 2)
    c_, s_ = np.cos(np.radians(top[1])), np.sin(np.radians(top[1]))
    Dall = tp(ra, de); Dall = np.c_[Dall.lon.arcsec, Dall.lat.arcsec] @ np.array([[c_, s_], [-s_, c_]]) + [top[2], top[3]]
    back = SkyCoord(Dall[:, 0] * u.arcsec, Dall[:, 1] * u.arcsec, frame=c0.skyoffset_frame()).transform_to("icrs"); ra, de = back.ra.deg, back.dec.deg
    T = cKDTree(np.c_[gra * np.cos(np.radians(c0.dec.deg)), gde])
    best = (0, 0.0, 0.0); rh, rs, rmt = SOL.refine_half_arcsec, SOL.refine_step_arcsec, SOL.refine_match_arcsec
    for dx in np.arange(-rh, rh + rs, rs) / 3600:
        for dy in np.arange(-rh, rh + rs, rs) / 3600:
            q = np.c_[(ra[bright] + dx / np.cos(np.radians(c0.dec.deg))) * np.cos(np.radians(c0.dec.deg)), de[bright] + dy]
            n = int(np.sum(T.query(q, distance_upper_bound=rmt / 3600)[0] < rmt / 3600))
            if n > best[0]:
                best = (n, dx, dy)
    ra = ra + best[1] / np.cos(np.radians(c0.dec.deg)); de = de + best[2]
    tan = lambda r, d: SkyCoord(r, d, unit="deg").transform_to(c0.skyoffset_frame())
    gx = tan(gra, gde); G2 = cKDTree(np.c_[gx.lon.arcsec, gx.lat.arcsec])
    GBt = cKDTree(np.c_[gx.lon.arcsec, gx.lat.arcsec][gm < SOL.stage_gaia_g]); gb_idx = np.flatnonzero(gm < SOL.stage_gaia_g)
    br = np.zeros(len(P["x"]), bool); br[np.argsort(-P["flux"])[:SOL.stage_n_bright]] = True
    # wide radii on bright stars only (sparse, so unconfused) with a low order over the whole box, then all stars, a
    # shrinking radius, and the full order: the outer box is matched before the polynomial is free to bend there
    for rad, order, bright in SOL.solve_stages:
        X = basis(P["x"], P["y"], P["shape"], order)
        s = tan(ra, de); q = np.c_[s.lon.arcsec, s.lat.arcsec]
        if bright:
            d, jb_ = GBt.query(q, distance_upper_bound=rad); m = np.isfinite(d) & br; j = np.where(m, gb_idx[np.minimum(jb_, len(gb_idx) - 1)], 0)
        else:
            d, j = G2.query(q, distance_upper_bound=rad); m = np.isfinite(d)
        cx, *_ = np.linalg.lstsq(X[m], gx.lon.arcsec[j[m]], rcond=None); cy, *_ = np.linalg.lstsq(X[m], gx.lat.arcsec[j[m]], rcond=None)
        lon, lat = X @ cx, X @ cy
        back = SkyCoord(lon * u.arcsec, lat * u.arcsec, frame=c0.skyoffset_frame()).transform_to("icrs")
        ra, de = back.ra.deg, back.dec.deg
    s = tan(ra, de); d, j = G2.query(np.c_[s.lon.arcsec, s.lat.arcsec], distance_upper_bound=rad)
    L = np.c_[s.lon.arcsec, s.lat.arcsec]; sp = 3600.0 * SOL.inverse_span_deg
    Q = lambda l: basis(l[:, 0] + sp / 2, l[:, 1] + sp / 2, (sp, sp))
    ix, *_ = np.linalg.lstsq(Q(L), P["x"], rcond=None); iy, *_ = np.linalg.lstsq(Q(L), P["y"], rcond=None)
    def to_pix(r, d_):
        t_ = tan(r, d_); l = np.c_[t_.lon.arcsec, t_.lat.arcsec]; return Q(l) @ ix, Q(l) @ iy
    P["to_pix"] = to_pix
    def to_sky(x, y):
        Xb = basis(np.atleast_1d(x), np.atleast_1d(y), P["shape"])
        b_ = SkyCoord((Xb @ cx) * u.arcsec, (Xb @ cy) * u.arcsec, frame=c0.skyoffset_frame()).transform_to("icrs"); return b_.ra.deg, b_.dec.deg
    P["to_sky"] = to_sky
    return ra, de, d, j, dict(blind=[float(top[0]), float(top[1]), round(float(top[2]), 1), round(float(top[3]), 1)], shift_arcsec=[round(best[1] * 3600, 1), round(best[2] * 3600, 1)], shift_matches=best[0],
                              gaia_matched=int(np.isfinite(d).sum()), resid_median_arcsec=round(float(np.median(d[np.isfinite(d)])), 3),
                              match_radius=round(rad, 2))


def blind(P, cache):
    """Replace P's header WCS by the meridian search's (see meridian.py)."""
    r = meridian.search(P["x"], P["y"], P["flux"], P["shape"][1], P["shape"][0], meridian.mid_time(P["hd"]), cache)
    P["w"] = meridian.wcs_from(r["ra"], r["dec"], r["par"], r["rot"], P["shape"][1], P["shape"][0]); P["wcs"] = dict(source="meridian", **r)


def in_box(P, r, d, m=DET.box_margin_px):
    px, py = P["to_pix"](r, d)
    return (px > P["x0"] + m) & (py > P["y0"] + m) & (px < P["x0"] + P["W"] - m) & (py < P["y0"] + P["W"] - m)


def aperture_noise(P, n=DET.random_n, seed=DET.apnoise_seed):
    """The spread of aperture sums at random positions (sigma-clipped): the noise a forced flux is judged against."""
    from astropy.stats import sigma_clipped_stats
    g = np.random.default_rng(seed); m = DET.apnoise_margin_px; xy = g.uniform(m, P["W"] - m, (n, 2))
    f = np.asarray(aperture_photometry(P["img"], CircularAperture(xy, RAP))["aperture_sum"], float)
    return sigma_clipped_stats(f, sigma=DET.clip_sigma, maxiters=DET.clip_iters)[2]


def ring(P, x, y, fwhm, n_null=0, seed=config.A.selection.ring_seed):
    """The dark ring: the mean of the background-subtracted image in an annulus around each source (inner radius 1.5 FWHM,
    at least 6 px; width 6 px), in units of the spread of the same annulus mean at random positions. Emulsion defects
    from development often carry a dark halo (negative); stars and glints do not. With n_null, also returns the random
    positions' values (the null distribution the cut is set from)."""
    from photutils.aperture import CircularAnnulus
    SE = config.A.selection; rm = SE.ring_margin_px
    xb, yb = np.asarray(x) - P["x0"], np.asarray(y) - P["y0"]; r1 = np.maximum(SE.ring_inner_fwhm * np.asarray(fwhm), SE.ring_inner_min_px)
    out = np.full(len(xb), np.nan); ok = (xb > rm) & (yb > rm) & (xb < P["W"] - rm) & (yb < P["W"] - rm)
    for r in np.unique(np.round(r1[ok])):                  # one annulus size at a time (vectorised per size)
        k = ok & (np.round(r1) == r)
        ap = CircularAnnulus(np.c_[xb[k], yb[k]], r, r + SE.ring_width_px); out[k] = np.asarray(aperture_photometry(P["img"], ap)["aperture_sum"], float) / ap.area
    if "ringnoise" not in P:                                # spread of the annulus mean on blank sky (the random annulus)
        g = np.random.default_rng(seed); m = SE.ringnoise_margin_px; xy = g.uniform(m, P["W"] - m, (DET.random_n, 2))
        ap = CircularAnnulus(xy, *SE.ring_random_px); v = np.asarray(aperture_photometry(P["img"], ap)["aperture_sum"], float) / ap.area
        from astropy.stats import sigma_clipped_stats
        P["ringnoise"] = sigma_clipped_stats(v, sigma=DET.clip_sigma, maxiters=DET.clip_iters)[2]
    return out / P["ringnoise"]


def forced(P, x, y):
    """Aperture flux and its significance at pixel (x, y) of P's residual image (inside its box), else nan."""
    xb, yb = x - P["x0"], y - P["y0"]; m = DET.forced_margin_px; ok = (xb > m) & (yb > m) & (xb < P["W"] - m) & (yb < P["W"] - m)
    f = np.full(len(x), np.nan); z = np.full(len(x), np.nan)
    if ok.any():
        f[ok] = np.asarray(aperture_photometry(P["img"], CircularAperture(np.c_[xb[ok], yb[ok]], RAP))["aperture_sum"], float)
        z[ok] = f[ok] / P["apnoise"]
    return f, z


STAMP = DET.stamp_half_px                                   # stamps: (2 STAMP + 1) px square, about 25" (images kept for offline tests)


def stamps(img, x, y, x0=0, y0=0, half=STAMP):
    """Stamps of img (box pixel = full pixel - (x0, y0)) centred on the nearest pixel to each (x, y), zero outside the
    image; float32 (n, 2 half + 1, 2 half + 1), and the sub-pixel offsets of each centre (dx, dy)."""
    x, y = np.atleast_1d(np.asarray(x, float)) - x0, np.atleast_1d(np.asarray(y, float)) - y0; n = len(x); S = 2 * half + 1
    out = np.zeros((n, S, S), np.float32); X, Y = np.round(x).astype(int), np.round(y).astype(int)
    H, W = img.shape
    inner = (Y >= half) & (Y + half < H) & (X >= half) & (X + half < W); k = np.flatnonzero(inner); g = np.arange(-half, half + 1)
    out[k] = img[Y[k, None, None] + g[None, :, None], X[k, None, None] + g[None, None, :]]   # wholly inside: one gather
    for i in np.flatnonzero(~inner):                                                      # at an edge: the part inside
        ya, yb, xa, xb = Y[i] - half, Y[i] + half + 1, X[i] - half, X[i] + half + 1
        sy0, sy1, sx0, sx1 = max(0, ya), min(H, yb), max(0, xa), min(W, xb)
        if sy1 > sy0 and sx1 > sx0:
            out[i, sy0 - ya:sy1 - ya, sx0 - xa:sx1 - xa] = img[sy0:sy1, sx0:sx1]
    return out, np.c_[x - X, y - Y].astype(np.float32)


def shifted_forced(P_other, ra, dec, z, sign, radius=CH_SHIFT):
    """Forced significance on the twin at each detection's position moved by sign x the k-th shift (k = 0..7), for the
    chance pairs (aggregate.py moves b's sources by +shift: a pseudo pair's b measurement lies at a's position minus the
    shift, its a measurement at b's position plus it). Only detections with z >= 6; others nan."""
    out = np.full((len(ra), CH_DIRS), np.nan); k_ = np.flatnonzero(np.nan_to_num(z) >= DET.chance_zmin)
    for k in range(CH_DIRS):
        t = 2 * np.pi * k / CH_DIRS; dd = sign * radius * np.sin(t) / 3600; dr = sign * radius * np.cos(t) / 3600 / np.cos(np.radians(dec[k_]))
        out[k_, k] = forced(P_other, *P_other["to_pix"](ra[k_] + dr, dec[k_] + dd))[1]
    return out


def csd_at(P, x, y):
    """The detection noise of the tile holding box pixel (x, y)."""
    for ty0, ty1, tx0, tx1, c in P["csd"]:
        if ty0 <= y < ty1 and tx0 <= x < tx1:
            return c
    return P["csd"][0][4]


def fshift(a, dy, dx):
    """a moved by (dy, dx) px with a Fourier shift (exact for a band-limited image; no interpolation smoothing)."""
    from scipy import ndimage
    return np.fft.ifft2(ndimage.fourier_shift(np.fft.fft2(a), (dy, dx))).real


def stamp_centroid(s_, sig=DET.centroid_sigma_px, it=DET.centroid_iters):
    """A stamp's light centre relative to its middle pixel (x, y): Gaussian-weighted, iterated; nan if no light."""
    n = s_.shape[0]; g = np.arange(n) - (n - 1) / 2; cx = cy = 0.0
    for _ in range(it):
        w = np.exp(-((g[None] - cx) ** 2 + (g[:, None] - cy) ** 2) / (2 * sig ** 2)) * np.clip(s_, 0, None)
        if w.sum() <= 0:
            return np.nan, np.nan
        cx, cy = float((w.sum(0) @ g) / w.sum()), float((w.sum(1) @ g) / w.sum())
    return cx, cy


def glint_image(stamp, scale, length, ang):
    """A glint's image from the local stars' average image: shrunk about its centre by scale (flux kept), then smeared
    by a line of length px at angle ang (the plate's RA direction: a geostationary object moves 15.04" cos dec per second
    on a sidereally guided plate)."""
    from scipy import ndimage
    n = stamp.shape[0]; c = (n - 1) / 2
    if scale != 1.0:
        m = np.eye(2) / scale; off = c - m @ np.array([c, c])
        stamp = ndimage.affine_transform(stamp, m, offset=off, order=INJ.glint_spline_order) / scale ** 2
    if length > 1.0:
        L = int(np.ceil(length)) | 1; ker = np.zeros((L, L)); t = np.linspace(-length / 2, length / 2, INJ.glint_line_samples * L)
        np.add.at(ker, (np.clip(np.round(L // 2 + t * np.sin(ang)).astype(int), 0, L - 1), np.clip(np.round(L // 2 + t * np.cos(ang)).astype(int), 0, L - 1)), 1.0)
        stamp = ndimage.convolve(stamp, ker / ker.sum(), mode="constant")
    return stamp


N_DONOR, DONOR_DG, DONOR_R = INJ.donor_n, INJ.donor_dmag, INJ.donor_r_px   # local image model: nearest donors within DONOR_DG and DONOR_R


def _nanmedian0(A):
    """np.nanmedian(A, 0) for a 3-D float64 array (compiled; identical values)."""
    global _NANMED0
    if _NANMED0 is None:
        from numba import njit
        @njit(cache=True)
        def f(A):
            n, ny, nx = A.shape; out = np.empty((ny, nx)); buf = np.empty(n)
            for i in range(ny):
                for j in range(nx):
                    k = 0
                    for t in range(n):                              # finite values, kept sorted (insertion; n <= 40)
                        v = A[t, i, j]
                        if v == v:
                            q = k
                            while q > 0 and buf[q - 1] > v:
                                buf[q] = buf[q - 1]; q -= 1
                            buf[q] = v; k += 1
                    if k == 0:
                        out[i, j] = np.nan
                    else:
                        h = k // 2
                        out[i, j] = buf[h] if k % 2 else (buf[h - 1] + buf[h]) / 2.0
            return out
        _NANMED0 = f
    return _NANMED0(np.ascontiguousarray(A, dtype=np.float64))


_NANMED0 = None


def local_model(P, pg, gm, donors, x, y, Gt, D, tdet):
    """A low-noise image of a star of G = Gt at plate pixel (x, y): the nearest N_DONOR donors (stars of the Sun's
    colour within DONOR_DG mag and DONOR_R px), each stamp ((2D + 1) px) centred on its own light with a Fourier shift
    (of the continuous stamp), with every other detection masked where it lies after the shift, its local pedestal (the
    median of its unmasked outer annulus) removed, and scaled to Gt; then combined pixel by pixel with a sigma-clipped
    mean (3 sigma about the median, by the spread across donors), so that light near any one donor (faint neighbours,
    halos of other stars) is not carried into every glint; then tapered to zero beyond the radius at which its azimuthal
    profile falls to the stack's own residual noise. Returns (model centred on the middle pixel, donors used, residual
    noise of the stack, mean donor G); None when no donor is usable."""
    c = donors[np.abs(gm[donors] - Gt) < DONOR_DG]; dist = np.hypot(pg[0][c] - x, pg[1][c] - y); c = c[dist < DONOR_R]
    c = c[np.argsort(np.hypot(pg[0][c] - x, pg[1][c] - y))][:INJ.donor_pool * N_DONOR]
    S = 2 * D + 1; g = np.arange(S) - D; R = np.hypot(g[None], g[:, None]); outer = R > INJ.donor_outer * D; stack, Gu = [], []
    cache = P.setdefault("_donor_cache", {})                     # a donor's prepared stamp depends only on (donor, D)
    for dd in c:
        if len(stack) >= N_DONOR:
            break
        if (dd, D) in cache:
            prep = cache[(dd, D)]
            if prep is not None:
                stack.append(np.where(prep[1], prep[0] * 10 ** (-0.4 * (Gt - gm[dd])), np.nan)); Gu.append(gm[dd])
            continue
        cache[(dd, D)] = None
        dx, dy = pg[0][dd], pg[1][dd]; dq_, jq_ = tdet.query([dx, dy])
        if dq_ < INJ.donor_snap_px:
            dx, dy = P["x"][jq_], P["y"][jq_]
        X, Y = int(np.round(dx)) - P["x0"], int(np.round(dy)) - P["y0"]
        if not (D + 1 <= X < P["W"] - D - 1 and D + 1 <= Y < P["W"] - D - 1):
            continue
        st = P["img"][Y - D:Y + D + 1, X - D:X + D + 1].astype(float)
        others = [(P["x"][k] - P["x0"] - X, P["y"][k] - P["y0"] - Y, max(INJ.donor_mask_a * P["a"][k], INJ.donor_mask_min_px))
                  for k in tdet.query_ball_point([X + P["x0"], Y + P["y0"]], D * INJ.donor_mask_reach)]
        others = [o for o in others if np.hypot(o[0], o[1]) >= INJ.donor_self_px]   # every other detection in the stamp
        def mask(dx_=0.0, dy_=0.0):
            w = np.ones((S, S), bool)
            for ox, oy, rr in others:
                w[(g[None] - ox - dx_) ** 2 + (g[:, None] - oy - dy_) ** 2 < rr ** 2] = False
            return w
        cx_, cy_ = stamp_centroid(st * mask())
        if not np.isfinite(cx_):
            continue
        w = mask(-cx_, -cy_); sh = fshift(st, -cy_, -cx_)                 # shift the continuous stamp, then mask
        if (w & outer).sum() < INJ.donor_min_outer:
            continue
        sh0 = sh - np.median(sh[w & outer]); cache[(dd, D)] = (sh0, w)        # local pedestal off
        sh = sh0 * 10 ** (-0.4 * (Gt - gm[dd]))                                   # to the target's brightness
        stack.append(np.where(w, sh, np.nan)); Gu.append(gm[dd])
    if not stack:
        return None
    A = np.array(stack); med = _nanmedian0(A); sp = 1.4826 * _nanmedian0(np.abs(A - med))
    A = np.where(np.abs(A - med) <= INJ.donor_clip * np.maximum(sp, 1e-9), A, np.nan)   # clipped: one donor's neighbour is rejected
    with np.errstate(invalid="ignore"):
        m = np.nan_to_num(np.nanmean(A, 0))
    nr = D - INJ.model_noise_ring_px
    noise = float(1.4826 * np.median(np.abs(m[R > nr] - np.median(m[R > nr])))) if (R > nr).any() else 0.0
    prof = np.array([m[(R >= r) & (R < r + 1)].mean() for r in range(D)]); low = np.flatnonzero(prof < max(noise, 1e-9))
    rc = (low[0] if len(low) else nr) + INJ.taper_pad_px
    taper = np.clip(1 - (R - rc) / INJ.taper_px, 0, 1)                # unity inside rc, to zero over taper_px
    return m * taper, len(stack), noise, float(np.mean(Gu))


def inject(PA, PB, G, BPRP, FA, FB, n=INJ.n_inj, seed=INJ.inj_seed):
    """Artificial glints through the same detection and measurement as real sources: n positions drawn uniformly in a's
    box pixels over the candidates' region (30 px from both boxes' edges), at least 2 x 70 px apart; at each, on each
    twin, the plate's local image model of a star of the target G (local_model: uniform in G_INJ_MIN to G_INJ_MAX), shrunk by
    the sharpening and smeared along that plate's RA direction for the duration (glint_image), placed to a fraction of a
    pixel (Fourier shift); all glints of a plate added to its scan (detect(extra), clipped at saturation) and the whole
    detection rerun; each glint is the new detection nearest its target within 3 px, measured with the functions used
    for real detections (forced at its own centroid and at the twin's, ring). Returns (table, stamps {"a", "b"}: each
    glint's stamp on each plate, from the injected scan; zero where not found)."""
    rng = np.random.default_rng(seed); gra, gde, gm = G
    pa_g, pb_g = PA["to_pix"](gra, gde), PB["to_pix"](gra, gde)
    tg = cKDTree(np.c_[pa_g]); lone = np.array([len(v) for v in tg.query_ball_point(np.c_[pa_g], INJ.donor_lone_px)]) == 1   # no other Gaia star near
    (c0, c1), dm = INJ.donor_bprp, INJ.donor_margin_px
    donors = np.flatnonzero(lone & (BPRP > c0) & (BPRP < c1) & in_box(PA, gra, gde, dm) & in_box(PB, gra, gde, dm))
    S = 2 * STAMP + 1; empty = (pd.DataFrame(), {t: np.zeros((0, S, S), np.float32) for t in "ab"})
    if len(donors) < INJ.donor_min:
        print(f"inject: only {len(donors)} donors; skipped", file=sys.stderr); return empty
    pts = []; tries = 0; bm = DET.box_margin_px
    while len(pts) < n and tries < INJ.inj_tries * n:                    # uniform in a's box pixels, kept apart
        tries += 1; xa = rng.uniform(PA["x0"] + bm, PA["x0"] + PA["W"] - bm); ya = rng.uniform(PA["y0"] + bm, PA["y0"] + PA["W"] - bm)
        if not pts or np.min((np.array(pts)[:, 0] - xa) ** 2 + (np.array(pts)[:, 1] - ya) ** 2) > INJ.inj_sep_px ** 2:
            pts.append((xa, ya))
    pts = np.array(pts); r_, d_ = PA["to_sky"](pts[:, 0], pts[:, 1]); r_, d_ = np.asarray(r_, float), np.asarray(d_, float)
    keep = in_box(PA, r_, d_, bm) & in_box(PB, r_, d_, bm); r_, d_ = r_[keep], d_[keep]
    rows = [dict(ra=r_[i], dec=d_[i], G=float(rng.uniform(G_INJ_MIN, G_INJ_MAX)), dur=float(rng.choice(DURATIONS)), sharpen=float(rng.choice(SHARPEN))) for i in range(len(r_))]
    (g1, g2), (h1, h2, h3) = INJ.inj_model_g, INJ.inj_model_half_px
    layers, tdets = {}, {}
    for tag, P, pg, f in (("a", PA, pa_g, FA), ("b", PB, pb_g, FB)):
        L = np.zeros((P["W"], P["W"]), np.float32); tdets[tag] = tdet = cKDTree(np.c_[P["x"], P["y"]])
        x, y = P["to_pix"](r_, d_)
        xd, yd = P["to_pix"](r_ + INJ.inj_dir_arcsec / 3600 / np.cos(np.radians(d_)), d_)   # each target moved along RA: the plate's RA direction
        for i, row in enumerate(rows):
            L_arcsec = PHY.geo_trail_arcsec_s * np.cos(np.radians(row["dec"])) * row["dur"]
            D = (h1 if row["G"] < g1 else h2 if row["G"] < g2 else h3) + int(L_arcsec / SCALE / 2) + INJ.inj_model_pad_px
            mod = local_model(P, pg, gm, donors, x[i], y[i], row["G"], D, tdet)
            if mod is None:
                row[f"placed_{tag}"] = False; continue
            img, nd, nz, gd = mod
            v = (np.array([xd[i], yd[i]], float) - np.array([x[i], y[i]], float)) / INJ.inj_dir_arcsec; ang = np.arctan2(v[1], v[0])
            X, Y = int(np.floor(x[i])) - P["x0"], int(np.floor(y[i])) - P["y0"]
            img = fshift(glint_image(img, row["sharpen"], L_arcsec * np.hypot(*v), ang), y[i] - np.floor(y[i]), x[i] - np.floor(x[i]))
            ys0, ys1, xs0, xs1 = max(0, Y - D), min(P["W"], Y + D + 1), max(0, X - D), min(P["W"], X + D + 1)
            L[ys0:ys1, xs0:xs1] += img[ys0 - (Y - D):ys1 - (Y - D), xs0 - (X - D):xs1 - (X - D)].astype(np.float32)
            row.update({f"placed_{tag}": True, f"n_donor_{tag}": nd, f"model_noise_{tag}": nz, f"donor_G_{tag}": gd, f"x_t_{tag}": float(x[i]), f"y_t_{tag}": float(y[i])})
        layers[tag] = L
    Q = {}
    for tag, P, f in (("a", PA, FA), ("b", PB, FB)):                  # the whole detection again, on the scan with the glints
        Qi = detect(f, extra=layers[tag]); del layers[tag]
        Qi.update(apnoise=P["apnoise"], ringnoise=P.get("ringnoise", np.nan), to_sky=P["to_sky"], to_pix=P["to_pix"]); Q[tag] = Qi
    st = {}
    for tag in "ab":
        Qi = Q[tag]; t_ = cKDTree(np.c_[Qi["x"], Qi["y"]]); xt = np.array([r.get(f"x_t_{tag}", np.nan) for r in rows]); yt = np.array([r.get(f"y_t_{tag}", np.nan) for r in rows])
        ok = np.isfinite(xt); d, j = np.full(len(rows), np.inf), np.zeros(len(rows), int)
        d[ok], j[ok] = t_.query(np.c_[xt[ok], yt[ok]]); found = d < INJ.inj_found_px
        z = np.full(len(rows), np.nan); z[found] = forced(Qi, Qi["x"][j[found]], Qi["y"][j[found]])[1]
        rg = np.full(len(rows), np.nan); rg[found] = ring(Qi, Qi["x"][j[found]], Qi["y"][j[found]], Qi["fwhm"][j[found]])
        rs, ds = (np.full(len(rows), np.nan) for _ in range(2))
        if found.any():
            a_, b_ = Qi["to_sky"](Qi["x"][j[found]], Qi["y"][j[found]]); rs[found], ds[found] = a_, b_
        for i, row in enumerate(rows):
            row[f"found_{tag}"] = bool(found[i])
            if found[i]:
                k = j[i]; row.update({f"flux_{tag}": float(Qi["flux"][k]), f"a_{tag}": float(Qi["a"][k]), f"b_{tag}": float(Qi["b"][k]), f"peak_{tag}": float(Qi["peak"][k]),
                                      f"theta_{tag}": float(Qi["theta"][k]), f"x_{tag}": float(Qi["x"][k]), f"y_{tag}": float(Qi["y"][k]), f"z_{tag}": float(z[i]),
                                      f"ring_{tag}": float(rg[i]), f"ra_f_{tag}": float(rs[i]), f"dec_f_{tag}": float(ds[i])})
        st[tag] = np.zeros((len(rows), S, S), np.float32)
        if found.any():
            st[tag][found] = stamps(Qi["img"], Qi["x"][j[found]], Qi["y"][j[found]], Qi["x0"], Qi["y0"])[0]
    for tag, other in (("a", "b"), ("b", "a")):                       # each plate forced at the twin's detected position
        Qi = Q[tag]; ii = [i for i, r in enumerate(rows) if r.get(f"found_{other}")]
        zf = np.full(len(rows), np.nan)
        if ii:
            px, py = Qi["to_pix"](np.array([rows[i][f"ra_f_{other}"] for i in ii]), np.array([rows[i][f"dec_f_{other}"] for i in ii])); zf[ii] = forced(Qi, px, py)[1]
        for i, row in enumerate(rows):
            row[f"zf_{tag}"] = float(zf[i])
    del Q
    return pd.DataFrame(rows), st

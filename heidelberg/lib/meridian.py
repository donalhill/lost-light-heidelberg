"""Blind centre for a plate whose header astrometry is wrong (GAVO's centre below the horizon at the plate time, or a
header solve that fails): the field must have been above the horizon at mid-exposure (tiles above 8 deg altitude, hour
angle within 9 h, so comet plates far from the meridian are found), at the instrument's scale (1.02"/px) and scan
orientation (rotated 176-184 deg, either parity). Gaia DR3 G < 11 all sky (CACHE/allsky_g11.npz, pinned), 6-deg
tangent-plane tiles, FFT cross-correlation of 40" rasters of the 150 brightest detections against Gaia. Its settings are
in config/analysis.json (section solve, meridian_*; the site in section physics)."""
import os
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u

from heidelberg import config

SOL, PH = config.A.solve, config.A.physics
SITE = EarthLocation(lon=PH.site_lon_deg * u.deg, lat=PH.site_lat_deg * u.deg, height=PH.site_height_m * u.m)
SCALE = config.A.plates.pixel_arcsec


def mid_time(hd):
    t = Time(f"{hd['DATE-OBS']}T{hd['TM_START']}", scale="utc") + hd.get("EXPTIME", 0) / 2 * u.s
    return t


def header_ok(hd, w):
    t = mid_time(hd); c = w.pixel_to_world(hd["NAXIS1"] / 2, hd["NAXIS2"] / 2)
    from astropy.coordinates import AltAz
    alt = c.transform_to(AltAz(obstime=t, location=SITE)).alt.deg
    return alt > SOL.header_alt_min, alt


def strip(lst, cache):
    """Gaia DR3 G < 11, all sky (cached in CACHE/allsky_g11.npz, pinned; else from ESA's archive or IRSA's copy, 1.25M
    stars, about 2 min to fetch)."""
    f = os.path.join(cache, "allsky_g11.npz")
    if not os.path.exists(f):
        from heidelberg.lib.catalogues import gaia_async
        t = gaia_async("ra, dec, phot_g_mean_mag", f"phot_g_mean_mag < {SOL.meridian_gmax}", maxrec=config.C.network.allsky_maxrec, wait=config.C.network.allsky_wait_s)
        os.makedirs(cache, exist_ok=True)
        np.savez(f, ra=np.asarray(t["ra"], float), dec=np.asarray(t["dec"], float), g=np.asarray(t["phot_g_mean_mag"], float))
    z = np.load(f); return z["ra"], z["dec"]


def raster(xy, half, pix):
    n = int(2 * half / pix); i = np.floor((xy + half) / pix).astype(int); ok = ((i >= 0) & (i < n)).all(1)
    im = np.zeros((n, n), np.float32); np.add.at(im, (i[ok, 0], i[ok, 1]), 1.0)
    return np.minimum(im, 1)


def search(x, y, flux, nx, ny, t, cache, nbright=SOL.meridian_n_det, ha_max=SOL.meridian_ha_max_h):
    lst = t.sidereal_time("apparent", longitude=SITE.lon).deg
    gra, gde = strip(lst, cache)
    br = np.argsort(-flux)[:nbright]; D0 = np.c_[x[br] - nx / 2, y[br] - ny / 2] * SCALE
    half, pix = SOL.meridian_tile_deg / 2 * 3600, SOL.meridian_raster_arcsec
    from scipy.ndimage import maximum_filter
    FD = {}
    for par in (1, -1):
        for rot in np.arange(*SOL.meridian_rot_deg):
            c_, s_ = np.cos(np.radians(rot)), np.sin(np.radians(rot))
            FD[par, float(rot)] = np.conj(np.fft.rfft2(maximum_filter(raster((D0 * [par, 1]) @ np.array([[c_, s_], [-s_, c_]]), half, pix), SOL.meridian_dilate_px)))
    best = (0, None)
    lat = np.radians(SITE.lat.deg)
    d0, d1, dstep = SOL.meridian_dec_deg
    for dc in np.arange(d0, d1, dstep):
        step = dstep / max(np.cos(np.radians(dc)), SOL.meridian_cos_floor)
        for rc in np.arange(lst - ha_max * 15, lst + ha_max * 15 + step, step):
            alt = np.degrees(np.arcsin(np.sin(np.radians(dc)) * np.sin(lat) + np.cos(np.radians(dc)) * np.cos(lat) * np.cos(np.radians(lst - rc))))
            if alt < SOL.meridian_alt_min:
                continue
            nd = SOL.meridian_near_deg; dra = ((gra - rc + 180) % 360 - 180) * np.cos(np.radians(dc)); near = (np.abs(dra) < nd) & (np.abs(gde - dc) < nd)
            if near.sum() < SOL.meridian_min_stars:
                continue
            c0 = SkyCoord(rc % 360, dc, unit="deg")
            s = SkyCoord(gra[near], gde[near], unit="deg").transform_to(c0.skyoffset_frame())
            G = raster(np.c_[s.lon.arcsec, s.lat.arcsec], half, pix); FG = np.fft.rfft2(G); n = G.shape[0]
            for key, F in FD.items():
                cc = np.fft.irfft2(FG * F, s=G.shape); k = np.unravel_index(cc.argmax(), cc.shape)
                if cc[k] > best[0]:
                    off = np.array([(k[0] + n // 2) % n - n // 2, (k[1] + n // 2) % n - n // 2]) * pix
                    best = (float(cc[k]), dict(rc=rc % 360, dc=float(dc), par=key[0], rot=key[1], off=off.tolist()))
    v, b = best
    c0 = SkyCoord(b["rc"], b["dc"], unit="deg")
    cen = SkyCoord(b["off"][0] * u.arcsec, b["off"][1] * u.arcsec, frame=c0.skyoffset_frame()).transform_to("icrs")
    return dict(votes=v, ra=float(cen.ra.deg), dec=float(cen.dec.deg), par=b["par"], rot=b["rot"], lst=float(lst))


def wcs_from(ra, dec, par, rot, nx, ny):
    from astropy.wcs import WCS
    w = WCS(naxis=2); w.wcs.ctype = ["RA---TAN", "DEC--TAN"]; w.wcs.crval = [ra, dec]; w.wcs.crpix = [nx / 2 + 1, ny / 2 + 1]
    c_, s_ = np.cos(np.radians(rot)), np.sin(np.radians(rot))
    w.wcs.cd = SCALE / 3600 * np.array([[par * c_, -s_], [par * s_, c_]])
    return w

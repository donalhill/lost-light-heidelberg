"""The asteroid clock over many plates, reading the compressed scans on S3 by section (no whole-scan download): for
each bright predicted asteroid, a cut-out around its mid-exposure position, local astrometry against Gaia (shift
search, then affine), the trail on the predicted track (segments within 6" across), and the offsets of its centroid
and ends along the track from the predicted middle, start, and end, as times (offset / rate). One row per segment near
the track, the brightest flagged as the trail.
Gaia comes from the pairs' saved queries (the measurement's PAIR_gq_RA_DEC.npz: G < 17.5 over both plates' boxes, at
the plate epoch), the one centred nearest the cut-out, so no catalogue service is asked.
python -E -s -m heidelberg.stages.clock_batch PRED_DIR PLATES.csv S3_FITS_PREFIX GQ_DIR OUT.csv (settings: config/analysis.json, section clock)
  PRED_DIR: our ephemeris (asteroid_ephem_024.py) parts; PLATES.csv: plate, ra_c, dec_c, wcs (header or meridian);
  S3_FITS_PREFIX: bucket/path holding PLATE.fits.fz (credentials from the environment); GQ_DIR: the saved queries."""
import glob, sys, warnings
import numpy as np, pandas as pd
import s3fs
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_sources, SourceCatalog
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.stats import sigma_clipped_stats
from scipy.spatial import cKDTree

from heidelberg import config

warnings.filterwarnings("ignore")
K, D_ = config.A.clock, config.A.detection
PRED, PLATES, PREFIX, GQ, OUT = sys.argv[1:6]; VMAX = K.clock_vmax
A = pd.concat([pd.read_parquet(f) for f in glob.glob(f"{PRED}/*.parquet")]); A = A[A.V < VMAX]
PL = pd.read_csv(PLATES).set_index("plate")
fs = s3fs.S3FileSystem(client_kwargs={"region_name": "eu-north-1"})
GQF = sorted(glob.glob(f"{GQ}/*_gq_*.npz")); GQC = SkyCoord([float(f.rsplit("_", 2)[1]) for f in GQF], [float(f.rsplit("_", 1)[1][:-4]) for f in GQF], unit="deg")


def gaia_saved(c):
    """Gaia brighter than clock_gaia_gmax at the plate epoch from the saved query centred nearest c (deg arrays ra, dec)."""
    z = np.load(GQF[int(np.argmin(c.separation(GQC).deg))]); k = z["G"] < K.clock_gaia_gmax
    return z["ra"][k], z["dec"][k]


k_ = Gaussian2DKernel(D_.smooth_sigma_px, x_size=D_.smooth_size_px, y_size=D_.smooth_size_px); k_.normalize()
rows, no_scan = [], []
for plate, a in A.groupby("plate"):
    if plate not in PL.index or not fs.exists(f"{PREFIX}/{plate}.fits.fz"):
        no_scan.append(plate); continue
    with fs.open(f"{PREFIX}/{plate}.fits.fz", "rb", block_size=2 ** 20) as fh, fits.open(fh, lazy_load_hdus=True) as hl:
        h = hl[1].header; nx, ny = h["ZNAXIS1"] if "ZNAXIS1" in h else h["NAXIS1"], h["ZNAXIS2"] if "ZNAXIS2" in h else h["NAXIS2"]
        if PL.loc[plate, "wcs"] == "meridian":
            w = WCS(naxis=2); w.wcs.ctype = ["RA---TAN", "DEC--TAN"]; w.wcs.crval = [PL.loc[plate, "ra_c"], PL.loc[plate, "dec_c"]]
            w.wcs.crpix = [nx / 2 + 1, ny / 2 + 1]; w.wcs.cd = config.A.plates.pixel_arcsec / 3600 * np.array([[-1.0, 0], [0, -1.0]])   # the scan's scale, east left and north down
        else:
            h2 = h.copy(); h2["CTYPE1"], h2["CTYPE2"] = "RA---TAN-SIP", "DEC--TAN-SIP"; w = WCS(h2)
        epoch = Time(f"{h['DATE-OBS']}T{h['TM_START']}").jyear; EXP = h.get("EXPTIME", 0) / 60
        for d, g in a.groupby("desig"):
            g = g.set_index("when")
            if not {"start", "mid", "end"} <= set(g.index):
                continue
            s, m, e = (g.loc[t, ["ra", "dec"]].to_numpy(float) for t in ("start", "mid", "end"))
            px, py = w.world_to_pixel(SkyCoord(m[0], m[1], unit="deg"))
            H = K.clock_cutout_half_px
            if not (H < px < nx - H and H < py < ny - H):
                continue
            x0, y0 = int(px - H), int(py - H)
            try:
                im = np.asarray(hl[1].section[y0:y0 + 2 * H, x0:x0 + 2 * H], np.float32)
                bx, fx = K.clock_bg_box_px, K.clock_bg_filter_px
                bk = Background2D(im, (bx, bx), filter_size=(fx, fx), bkg_estimator=MedianBackground()); r = im - bk.background
                cv = convolve_fft(r, k_); sd = sigma_clipped_stats(cv[::2, ::2], sigma=D_.clip_sigma)[2]
                seg = detect_sources(cv, D_.detect_sigma * sd, n_pixels=D_.detect_npix)
                T = SourceCatalog(r, seg, convolved_data=cv).to_table(columns=["label", "x_centroid", "y_centroid", "segment_flux", "semimajor_axis", "semiminor_axis"])
                X, Y = np.asarray(T["x_centroid"], float) + x0, np.asarray(T["y_centroid"], float) + y0
                c = w.pixel_to_world(px, py)
                gra, gde = gaia_saved(c); keep = c.separation(SkyCoord(gra, gde, unit="deg")).arcsec < H * K.clock_gaia_margin; gra, gde = gra[keep], gde[keep]
                fr = c.skyoffset_frame(); tg = SkyCoord(gra, gde, unit="deg").transform_to(fr); GX = np.c_[tg.lon.arcsec, tg.lat.arcsec]
                td = w.pixel_to_world(X, Y).transform_to(fr); DX = np.c_[td.lon.arcsec, td.lat.arcsec]
                br = np.argsort(-np.asarray(T["segment_flux"], float))[:K.clock_n_bright]
                ro = K.clock_offset_arcsec; off = (GX[None] - DX[br][:, None]).reshape(-1, 2); off = off[(np.abs(off) < ro).all(1)]
                Hh, xe, ye = np.histogram2d(off[:, 0], off[:, 1], bins=K.clock_offset_bins, range=[[-ro, ro], [-ro, ro]]); kk = np.unravel_index(Hh.argmax(), Hh.shape)
                Dp = DX + [(xe[kk[0]] + xe[kk[0] + 1]) / 2, (ye[kk[1]] + ye[kk[1] + 1]) / 2]; tree = cKDTree(GX)
                for rad in K.clock_match_arcsec:
                    dd, jj = tree.query(Dp, distance_upper_bound=rad); mm = np.isfinite(dd)
                    Am = np.c_[np.ones(mm.sum()), X[mm], Y[mm]]; cx = np.linalg.lstsq(Am, GX[jj[mm], 0], rcond=None)[0]; cy = np.linalg.lstsq(Am, GX[jj[mm], 1], rcond=None)[0]
                    Dp = np.c_[np.c_[np.ones(len(X)), X, Y] @ cx, np.c_[np.ones(len(X)), X, Y] @ cy]
                dd, jj = tree.query(Dp, distance_upper_bound=K.clock_match_arcsec[-1]); nfit = int(np.isfinite(dd).sum()); res = float(np.median(dd[np.isfinite(dd)])) if nfit else np.nan
                tp = SkyCoord(np.r_[s[0], m[0], e[0]], np.r_[s[1], m[1], e[1]], unit="deg").transform_to(fr); Pp = np.c_[tp.lon.arcsec, tp.lat.arcsec]
                L = np.linalg.norm(Pp[2] - Pp[0]); uu = (Pp[2] - Pp[0]) / max(L, 1e-9); nn = np.array([-uu[1], uu[0]])
                al, ac = (Dp - Pp[1]) @ uu, (Dp - Pp[1]) @ nn
                near = np.flatnonzero((np.abs(ac) < K.clock_search_across_arcsec) & (np.abs(al) < L + K.clock_search_along_arcsec) & ~np.isfinite(dd))
                lab = np.asarray(T["label"]); fl = np.asarray(T["segment_flux"], float)
                best = near[np.argmax(fl[near])] if len(near) else -1
                for i in near:
                    pix = np.argwhere(seg.data == lab[i]); sp = np.c_[np.ones(len(pix)), pix[:, 1] + x0, pix[:, 0] + y0]
                    sal = np.c_[sp @ cx, sp @ cy] @ uu - Pp[1] @ uu; lo, hi = np.percentile(sal, K.clock_end_pct)
                    rows.append(dict(plate=plate, epoch=round(epoch, 3), exp_min=round(EXP, 1), desig=d, V=round(float(g.V.iloc[0]), 2), rate=L / EXP,
                                     along=round(float(al[i]), 2), across=round(float(ac[i]), 2), dt_min=round(float(al[i] / (L / EXP)), 2),
                                     end_lo=round(float(lo), 1), end_hi=round(float(hi), 1), length_pred=round(float(L), 1), flux=float(fl[i]),
                                     ab=round(float(T["semimajor_axis"][i].value / max(T["semiminor_axis"][i].value, 1e-3)), 2),
                                     fit_stars=nfit, fit_res=round(res, 2), trail=bool(i == best)))
            except Exception as ex:
                rows.append(dict(plate=plate, desig=d, error=str(ex)[:80]))
    pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"{len(rows)} rows -> {OUT}; {len(no_scan)} plates without a kept scan" + (f": {' '.join(no_scan)}" if no_scan else ""))

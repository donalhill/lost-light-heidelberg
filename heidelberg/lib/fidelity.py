"""Injection fidelity from a pair's saved products: star-like glints (no sharpening, no duration) against Gaia stars of
the Sun's colour (BP-RP 0.6-1.1, as the glints' donors are; on these blue-sensitive plates a red star of the same G
records much fainter) and of the same G (within 0.25 mag; aperture significance saturates for bright sources, so it
cannot be the matching variable), both paired and measured as candidates are. Three comparisons: the twin agreement
test, the catalogue minor axis on each plate (overall and by G), and the background noise (9-12 px) around them.
Used by heidelberg/stages/measure_pair.py; can be rerun on saved products. The criteria (config/analysis.json, section
fidelity) are the spec's injection fidelity criteria (notes/pipeline_spec.md)."""
import glob

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from heidelberg.config import A
from heidelberg.lib import selection as SEL

F_ = A.fidelity
PASS_PTS, PAIR_PTS, AXIS_DEX, AXIS_BIN, NOISE = F_.fid_shortfall, F_.fid_pair_shortfall, F_.fid_axis_dex, F_.fid_axis_bin_dex, F_.fid_noise
DG, SUN, ISO, NSIG = F_.fid_match_dg, F_.fid_sun_bprp, F_.fid_iso_widths, F_.fid_pair_sigmas
STAMP = A.detection.stamp_half_px


def unit(r, d):
    r, d = np.radians(r), np.radians(d); return np.c_[np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]


def report(base, gq=None):
    """The fidelity report for the pair whose products are at base (base_a.npz, base_b.npz, base_inj.parquet,
    base_inj_stamps.npz, base_stamps_a.npz, and the Gaia query base_gq_*.npz or gq)."""
    za, zb = np.load(base + "_a.npz"), np.load(base + "_b.npz"); INJ = pd.read_parquet(base + "_inj.parquet").reset_index(drop=True)
    if not len(INJ) or "found_a" not in INJ:
        return dict(passed=None, note="no injections")
    INJST = np.load(base + "_inj_stamps.npz"); SA = np.load(base + "_stamps_a.npz")
    q = np.load(gq or sorted(glob.glob(base + "_gq_*.npz"))[-1]); tq = cKDTree(unit(q["ra"], q["dec"]))
    def g_of(ra_, de_):
        d_, j_ = tq.query(unit(np.asarray(ra_, float), np.asarray(de_, float))); ok = np.degrees(d_) * 3600 < F_.fid_gaia_match_arcsec
        return np.where(ok, q["G"][j_], np.nan), np.where(ok, q["bprp"][j_], np.nan)
    def nearest_G(Gt, Gs, k=F_.fid_match_k):
        out_ = []
        for v in Gt:
            o = np.argsort(np.abs(Gs - v))[:k]; out_ += list(o[np.abs(Gs[o] - v) < DG])
        return np.unique(np.array(out_, int))
    sun = lambda c: (c > SUN[0]) & (c < SUN[1])
    feats, score, CLIM, WMAX = SEL.consistency(za, zb)
    def agree_rate(T):
        sep = np.degrees(np.linalg.norm(unit(T.ra, T.dec) - unit(T.ra_b, T.dec_b), axis=1)) * 3600
        X, Mn = feats(T.flux_a.values, T.flux_b.values, T.fwhm_a.values, T.fwhm_b.values, T.a_a.values / np.maximum(T.b_a.values, 1e-3), T.a_b.values / np.maximum(T.b_b.values, 1e-3), sep,
                      T.x_a.values, T.y_a.values, T.x_b.values, T.y_b.values, T.theta_a.values, T.theta_b.values, T.b_a.values, T.b_b.values, T.a_a.values, T.a_b.values)
        okx = np.isfinite(X).all(1); sc = np.full(len(T), np.inf); sc[okx] = score(X[okx])
        return (sc <= CLIM) & (np.nan_to_num(Mn[:, 0], nan=9) <= WMAX[0]) & (np.nan_to_num(Mn[:, 1], nan=9) <= WMAX[1]) & (np.nan_to_num(Mn[:, 2], nan=9) <= 1) & (np.nan_to_num(Mn[:, 3], nan=9) <= 1)
    star_like = (INJ.get("dur", 0) == 0) & (INJ.get("sharpen", 1) == 1)
    J, Ti = SEL.injection_table(INJ[star_like]); keep = SEL.injection_pairing(J.ra_f_a, J.dec_f_a, J.ra_f_b, J.dec_f_b, za, zb); J, Ti = J[keep], Ti[keep]   # paired as candidates are
    zr = lambda T: (T.zo_a >= SEL.Z_MIN) & (T.zo_b >= SEL.Z_MIN) & (T.zf_a >= SEL.Z_MIN) & (T.zf_b >= SEL.Z_MIN)
    # isolated only, glints and stars alike: no other detection within ISO image widths on either plate (a glint landing
    # next to a source blends with it, as a real flash would, and that loss is in the efficiency; the Gaia comparison
    # stars under-represent blends, whose shifted centroid can break the Gaia match)
    ta_, tb_ = cKDTree(np.c_[za["x"], za["y"]]), cKDTree(np.c_[zb["x"], zb["y"]])
    def isolated(xa_, ya_, wa_, xb_, yb_, wb_, self_=False):
        k_ = 2 if self_ else 1
        da_ = ta_.query(np.c_[xa_, ya_], k=k_)[0]; db_ = tb_.query(np.c_[xb_, yb_], k=k_)[0]
        if self_:
            da_, db_ = da_[:, 1], db_[:, 1]
        return (da_ > ISO * np.asarray(wa_)) & (db_ > ISO * np.asarray(wb_))
    Ti = Ti[zr(Ti)]; Ti = Ti[isolated(Ti.x_a, Ti.y_a, Ti.fwhm_a, Ti.x_b, Ti.y_b, Ti.fwhm_b)]; Gi = INJ.loc[Ti.index, "G"].values
    ia, jb = SEL.candidates(za, zb); St = SEL.det_table(za, zb, ia, jb); St = St[za["gaia"][ia] & (za["gsep"][ia] < A.selection.ref_gaia_arcsec) & zb["gaia"][jb] & (zb["gsep"][jb] < A.selection.ref_gaia_arcsec) & zr(St).values]   # Gaia on both plates
    St = St[isolated(St.x_a, St.y_a, St.fwhm_a, St.x_b, St.y_b, St.fwhm_b, self_=True)]
    Gst, Cst = g_of(St.ra, St.dec); k_ = sun(Cst) & np.isfinite(Gst); St, Gst = St[k_], Gst[k_]
    pick = nearest_G(Gi, Gst)
    ai = agree_rate(Ti) if len(Ti) else np.zeros(0, bool); as_ = agree_rate(St.iloc[pick]) if len(pick) else np.zeros(0, bool)
    r0 = dict(n_inj=int(len(ai)), n_star=int(len(as_)), agree_inj=float(ai.mean()) if len(ai) else np.nan, agree_star=float(as_.mean()) if len(as_) else np.nan)
    ax = []
    for tag, z in (("a", za), ("b", zb)):
        s_ = z["inside"] & z["gaia"] & (z["gsep"] < F_.fid_gaia_match_arcsec); Gz, Cz = g_of(z["ra"][s_], z["dec"][s_]); bm = z["b"][s_]; ks = sun(Cz)
        for (_, c), Gt in zip(Ti.iterrows(), Gi):
            k = ks & (np.abs(Gz - Gt) < DG)
            if k.sum() >= F_.fid_axis_min_stars:
                ax.append((Gt, np.log10(c[f"b_{tag}"] / np.median(bm[k]))))
    ax = np.array(ax).reshape(-1, 2); bins = []
    for lo, hi in F_.fid_axis_g_bins:
        qq = (ax[:, 0] >= lo) & (ax[:, 0] < hi)
        if qq.sum() >= F_.fid_bin_min:
            bins.append(dict(G=f"{lo}-{hi}", n=int(qq.sum()), minor_axis_offset_dex=float(np.median(ax[qq, 1])), se=float(1.2533 * np.std(ax[qq, 1]) / np.sqrt(qq.sum()))))
    yy_, xx_ = np.mgrid[-STAMP:STAMP + 1, -STAMP:STAMP + 1]; rr = np.hypot(xx_, yy_); ann = rr >= F_.fid_noise_inner_px; ring = np.floor(rr).astype(int)
    def bg_noise(st):
        """The spread in the ring from fid_noise_inner_px to the stamp's edge (9-12 px) after removing each stamp's own azimuthal profile (1 px rings), so a bright
        image's falling wing is not read as noise (its shape is judged by the minor-axis and agreement tests)."""
        st = st.astype(float); st = st[np.isfinite(st).all((1, 2))]; out_ = []
        for s_ in st:
            prof_ = np.array([np.median(s_[ring == r]) for r in range(ring.max() + 1)]); res = (s_ - prof_[ring])[ann]
            out_.append(1.4826 * np.median(np.abs(res - np.median(res))))
        return float(np.median(out_)) if out_ else np.nan
    # paired: each glint's ring noise against the median of its own matched stars (5 nearest in G within DG), then the
    # median over glints. A median over the unique matched stars would weight G differently from the glints (the stars
    # picked are 0.3-0.7 mag fainter on the tuning pairs) while ring noise depends on G.
    # Both plates alike: each plate's glints against that plate's stars, pooled.
    per = []
    for tag, z_ in (("a", za), ("b", zb)):
        S_ = SA if tag == "a" else np.load(base + "_stamps_b.npz")
        sl = (star_like & INJ[f"found_{tag}"].fillna(False).astype(bool)).values
        gi_ = S_["idx"][S_["gaia"]]; Gsa, Csa = g_of(z_["ra"][gi_], z_["dec"][gi_]); Gsa = np.where(sun(Csa), Gsa, 99.0); Gsa = np.nan_to_num(Gsa, nan=99.0)
        gst = INJST[tag][np.flatnonzero(sl)]; Gg = INJ.G.values[sl]; ok_ = np.isfinite(gst.astype(float)).all((1, 2)); gst, Gg = gst[ok_], Gg[ok_]
        stS = S_["stamp"][S_["gaia"]]; ns_cache = {}
        def ns_of(k):
            if k not in ns_cache:
                ns_cache[k] = bg_noise(stS[k:k + 1])
            return ns_cache[k]
        for s_, g_ in zip(gst, Gg):
            o = np.argsort(np.abs(Gsa - g_))[:F_.fid_match_k]; o = o[np.abs(Gsa[o] - g_) < DG]
            if len(o):
                per.append(bg_noise(s_[None]) / np.median([ns_of(int(k)) for k in o]))
    per = np.array(per)
    noise_ratio = float(np.median(per)) if len(per) else np.nan
    rng_ = np.random.default_rng(F_.fid_boot_seed)                       # its counting noise (bootstrap over glints)
    noise_se = float(np.std([np.median(per[rng_.integers(0, len(per), len(per))]) for _ in range(F_.fid_boot_n)])) if len(per) > 1 else np.nan
    short = r0["agree_star"] - r0["agree_inj"]; n_ = max(r0["n_inj"], 1)
    sig = float(np.sqrt(max(r0["agree_inj"] * (1 - r0["agree_inj"]), 1 / n_) / n_)) if np.isfinite(r0["agree_inj"]) else np.nan
    axis = float(np.median(ax[:, 1])) if len(ax) else np.nan; axis_se = float(1.2533 * np.std(ax[:, 1]) / np.sqrt(len(ax))) if len(ax) else np.nan
    bin_ok = all(abs(b_["minor_axis_offset_dex"]) <= AXIS_BIN + NSIG * b_["se"] for b_ in bins) and bool(bins)   # each bin, within its counting noise
    return dict(pooled=r0, agree_shortfall=short, shortfall_sigma=sig, minor_axis_offset_dex=axis, minor_axis_by_G=bins, n_axis=int(len(ax)),
                background_noise_ratio=noise_ratio, background_noise_statistic=f"paired: median over glints of glint / median of its {F_.fid_match_k} G-matched stars", matched_in=f"Gaia G within {DG} mag, stars of the Sun's colour (BP-RP {SUN[0]}-{SUN[1]}), isolated (no detection within {ISO:g} widths) for agreement",
                minor_axis_se=axis_se, background_noise_se=noise_se,
                # per pair, gross failures only, each limit with its own counting noise; the acceptance test is pooled over
                # all pairs (heidelberg/stages/summarise.py) with the strict limits
                passed=bool(np.isfinite(short) and short <= PAIR_PTS + NSIG * sig and np.isfinite(axis) and abs(axis) <= AXIS_DEX + NSIG * axis_se and bin_ok
                            and n_ >= F_.fid_min_glints and np.isfinite(noise_ratio) and noise_ratio <= NOISE + NSIG * np.nan_to_num(noise_se)),
                criteria=dict(per_pair=f"shortfall <= {PAIR_PTS} + {NSIG} sigma; |minor axis| <= {AXIS_DEX} dex + {NSIG} s.e. overall and <= {AXIS_BIN} dex + {NSIG} s.e. in each G bin "
                                       f"(n >= {F_.fid_bin_min}); background noise ratio <= {NOISE} + {NSIG} s.e.; >= {F_.fid_min_glints} isolated star-like glints",
                              pooled=f"shortfall <= {PASS_PTS}; |minor axis| <= {AXIS_DEX} dex; background noise ratio <= {NOISE}"))

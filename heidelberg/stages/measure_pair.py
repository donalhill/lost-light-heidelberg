"""Stage 1-2, one Heidelberg twin pair (EC2, both scans): detect on each plate, solve each against Gaia DR3 (header
pre-correction; the blind meridian search when the header is wrong or the solve is not accepted), forced photometry of
each plate at its own detections and at the twin's (and at the twin's moved 20" in eight directions, for the chance
pairs), the dark ring, and injected glints measured exactly as candidates are (own and forced on each plate). Writes the
pair's products to OUT_DIR/PAIR*, with a manifest of input and output hashes (the injection fidelity report is
heidelberg/stages/fidelity.py's, from these products). Its settings are in config/analysis.json (sections detection, solve, and injections).
python -I -m heidelberg.stages.measure_pair A.fits B.fits OUT_DIR"""
import hashlib
import json
import os
import sys

import astropy.units as u
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from heidelberg import config
from heidelberg.lib import measure as M

FA, FB, OD = sys.argv[1:4]
os.makedirs(OD, exist_ok=True)
PAIR = os.path.basename(FA).replace(".fits", ""); OUT = f"{OD}/{PAIR}"
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.environ.get("HD_GAIA", os.path.join(os.path.dirname(HERE), "hd_gaia"))
DET, SOL = config.A.detection, config.A.solve
N_INJ = config.A.injections.n_inj                    # about 40 isolated star-like glints per pair for the fidelity check


def sha(f, chunk=1 << 24):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


PA, PB = M.detect(FA), M.detect(FB)
for P in (PA, PB):
    ok, alt = M.meridian.header_ok(P["hd"], P["w"]); P["wcs"] = dict(source="header", alt=round(float(alt), 1))
    br = np.argsort(-P["flux"])[:DET.trail_n_bright]; P["trail"] = round(float(np.median(P["a"][br] / P["b"][br])), 2)
    if not ok:
        M.blind(P, CACHE)
for attempt in range(SOL.solve_attempts):
    c0 = PA["w"].pixel_to_world(PA["shape"][1] / 2, PA["shape"][0] / 2)
    gq = f"{OUT}_gq_{c0.ra.deg:.1f}_{c0.dec.deg:.1f}.npz"                       # the query, kept for reruns
    if os.path.exists(gq) and "tiled" in np.load(gq).files:
        G4 = tuple(np.load(gq)[k] for k in ("ra", "dec", "G", "bprp"))
    else:
        cB = PB["w"].pixel_to_world(PB["shape"][1] / 2, PB["shape"][0] / 2)                # both boxes covered, whatever the plates' offset
        G4 = M.gaia(c0, (PA["W"] * M.SCALE / 3600 * SOL.gaia_width_factor + SOL.gaia_width_pad_deg + 2 * c0.separation(cB).deg) * u.deg, PA["epoch"]); np.savez(gq, ra=G4[0], dec=G4[1], G=G4[2], bprp=G4[3], tiled=True)
    G, BPRP = G4[:3], G4[3]
    M.precorrect(PA, G); M.precorrect(PB, G)
    ra_a, de_a, ga, _, sa = M.solve(PA, G); ra_b, de_b, gb, _, sb = M.solve(PB, G)
    for P_, s_, r_, d_ in ((PA, sa, ra_a, de_a), (PB, sb, ra_b, de_b)):
        s_["g50_astrometric"], s_["match_grid"] = M.match_grid(P_, r_, d_, G)
        s_["accepted"] = s_["gaia_matched"] >= SOL.solve_min_stars and M.grid_ok(s_["match_grid"]) and s_["resid_median_arcsec"] < SOL.solve_resid_arcsec   # the plan's limits
    bad = [P for P, s_ in ((PA, sa), (PB, sb)) if not s_["accepted"] and P["wcs"]["source"] == "header"]
    if not bad:
        break
    for P in bad:
        M.blind(P, CACHE)
for P_, s_ in ((PA, sa), (PB, sb)):                   # the solved centre of the plate and of the test box (for the roles)
    rc, dc = P_["to_sky"](np.array([P_["shape"][1] / 2]), np.array([P_["shape"][0] / 2])); s_["centre"] = [round(float(rc[0]), 5), round(float(dc[0]), 5)]
    rb, db_ = P_["to_sky"](np.array([P_["x0"] + P_["W"] / 2]), np.array([P_["y0"] + P_["W"] / 2])); s_["box_centre"] = [round(float(rb[0]), 5), round(float(db_[0]), 5)]
    s_["precorrect"], s_["wcs"], s_["trail_axis_ratio"] = P_.get("precorrect"), P_["wcs"], P_["trail"]
ua, ub = M.unit(ra_a, de_a), M.unit(ra_b, de_b); ta, tb = cKDTree(ua), cKDTree(ub)
inside_a = M.in_box(PB, ra_a, de_a) & M.in_box(PA, ra_a, de_a); inside_b = M.in_box(PA, ra_b, de_b) & M.in_box(PB, ra_b, de_b)
R = np.radians(config.A.pairing.pair_radius_arcsec / 3600)
da, _ = tb.query(ua, distance_upper_bound=R); db, _ = ta.query(ub, distance_upper_bound=R)
PA["apnoise"], PB["apnoise"] = M.aperture_noise(PA), M.aperture_noise(PB)
zself_a = M.forced(PA, PA["x"], PA["y"])[1]; zself_b = M.forced(PB, PB["x"], PB["y"])[1]
ring_a, ring_b = M.ring(PA, PA["x"], PA["y"], PA["fwhm"]), M.ring(PB, PB["x"], PB["y"], PB["fwhm"])
xa_in_b = PB["to_pix"](ra_a, de_a); xb_in_a = PA["to_pix"](ra_b, de_b)
ztw_a = M.forced(PB, *xa_in_b)[1]; ztw_b = M.forced(PA, *xb_in_a)[1]
zsh_a = M.shifted_forced(PB, ra_a, de_a, zself_a, -1); zsh_b = M.shifted_forced(PA, ra_b, de_b, zself_b, +1)
zshx_a = np.stack([M.shifted_forced(PB, ra_a, de_a, zself_a, -1, r_) for r_ in M.CH_XRADII], axis=2)   # (n, 8, radii)
zshx_b = np.stack([M.shifted_forced(PA, ra_b, de_b, zself_b, +1, r_) for r_ in M.CH_XRADII], axis=2)
INJ, INJST = M.inject(PA, PB, G, BPRP, FA, FB, n=N_INJ) if N_INJ > 0 else (pd.DataFrame(), {t: np.zeros((0, 2 * M.STAMP + 1, 2 * M.STAMP + 1), np.float32) for t in "ab"})
INJ.to_parquet(f"{OUT}_inj.parquet"); np.savez_compressed(f"{OUT}_inj_stamps.npz", a=INJST["a"], b=INJST["b"])
# stamps for offline shape tests: every non-Gaia detection with z >= stamp_zmin inside both boxes, and some Gaia stars
for name, P, inside, dg, zs in (("a", PA, inside_a, ga, zself_a), ("b", PB, inside_b, gb, zself_b)):
    (z0, z1) = DET.stamp_gaia_z; k = np.flatnonzero(inside & ~np.isfinite(dg) & (np.nan_to_num(zs) >= DET.stamp_zmin))
    ks = np.flatnonzero(inside & np.isfinite(dg) & (np.nan_to_num(zs) > z0) & (np.nan_to_num(zs) < z1)); ks = np.sort(np.random.default_rng(DET.stamp_seed).permutation(ks)[:DET.stamp_gaia_n])
    idx = np.r_[k, ks]; st_, off_ = M.stamps(P["img"], P["x"][idx], P["y"][idx], P["x0"], P["y0"])
    np.savez_compressed(f"{OUT}_stamps_{name}.npz", idx=idx, gaia=np.r_[np.zeros(len(k), bool), np.ones(len(ks), bool)], stamp=st_, offset=off_)
out = {"pair": PAIR, "aperture_noise": [PA["apnoise"], PB["apnoise"]], "plates": [os.path.basename(FA), os.path.basename(FB)], "box_px": PA["W"],
       "epoch": PA["epoch"], "gaia_stars": len(G[0]), "detections": [len(ra_a), len(ra_b)], "solve": [sa, sb],
       "common": [int(inside_a.sum()), int(inside_b.sum())], "gq": os.path.basename(gq)}
for name, P, inside, dtw, dg in (("a", PA, inside_a, da, ga), ("b", PB, inside_b, db, gb)):
    both = np.isfinite(dtw) & inside; one = ~np.isfinite(dtw) & inside; star = np.isfinite(dg)
    out[name] = {"both": int(both.sum()), "one_only": int(one.sum()), "one_only_gaia": int((one & star).sum()), "one_only_no_gaia": int((one & ~star).sum()),
                 "both_no_gaia": int((both & ~star).sum())}
    zs, zt, zsh, zshx = (zself_a, ztw_a, zsh_a, zshx_a) if name == "a" else (zself_b, ztw_b, zsh_b, zshx_b)
    xt, yt = xa_in_b if name == "a" else xb_in_a
    np.savez(f"{OUT}_{name}.npz", zself=zs, ztwin=zt, zshift=zsh, zshift_x=zshx, xshift_radii=np.array(M.CH_XRADII), ra=[ra_a, ra_b][name == "b"], dec=[de_a, de_b][name == "b"], flux=P["flux"], inside=inside,
             twin=np.isfinite(dtw), gaia=star, x=P["x"], y=P["y"], fwhm=P["fwhm"], peak=P["peak"], gsep=dg, a=P["a"], b=P["b"], npix=P["npix"], xt=xt, yt=yt,
             theta=P["theta"], ring=[ring_a, ring_b][name == "b"])
gra, gde, gm = G; inbox = np.ones(len(gra), bool); det = {}; bm = DET.box_margin_px
for name, P, ra_, de_ in (("a", PA, ra_a, de_a), ("b", PB, ra_b, de_b)):
    px, py = P["to_pix"](gra, gde)
    inbox &= (px > P["x0"] + bm) & (py > P["y0"] + bm) & (px < P["x0"] + P["W"] - bm) & (py < P["y0"] + P["W"] - bm)
    dd, jj = cKDTree(M.unit(ra_, de_)).query(M.unit(gra, gde), distance_upper_bound=np.radians(DET.gaia_det_arcsec / 3600))
    det[name] = np.isfinite(dd); det[name + "_flux"] = np.where(np.isfinite(dd), P["flux"][np.minimum(jj, len(P["flux"]) - 1)], np.nan)
np.savez(f"{OUT}_gaia.npz", G=gm, ra=gra, dec=gde, inbox=inbox, det_a=det["a"], det_b=det["b"], flux_a=det["a_flux"], flux_b=det["b_flux"])
out["gaia_in_both_boxes"] = int(inbox.sum())

out["injections_requested"], out["injections_placed"] = N_INJ, int(len(INJ))
import resource                                              # peak memory of this process (Linux: KiB; macOS: bytes)
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss; out["max_rss_gb"] = round(rss / (1e9 if sys.platform == "darwin" else 1e6), 2)
products = sorted(f for f in os.listdir(OD) if f.startswith(PAIR) and not f.endswith((".json", "_manifest.json")))
out["manifest"] = dict(inputs={os.path.basename(FA): sha(FA), os.path.basename(FB): sha(FB)}, gaia_cache=os.path.basename(gq),
                       outputs={f: sha(f"{OD}/{f}") for f in products})
json.dump(out, open(f"{OUT}.json", "w"), indent=1, default=float)
print(json.dumps({k: out[k] for k in ("pair", "detections", "common")}, indent=1, default=float))

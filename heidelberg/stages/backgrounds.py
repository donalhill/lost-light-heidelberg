"""Backgrounds that can appear on both twins, which the shifted-twin chance estimate cannot see.
1. Fixed-offset copies (a double exposure or a guiding jump moves the whole mount, so every star gets a copy at one
   vector on both twins): two-plate sources with no Gaia star, stacked by their offset (arcsec, east and north) from
   every Gaia star within 120", against the same stack with the offsets rotated at random (the radial profile kept);
   bins whose excess is above 5 standard deviations are reported with their vector.
2. Lens ghosts (reflections inside the twin lenses put a faint image of a bright star near its position reflected
   through the optical axis): two-plate sources with no Gaia star, stacked around the reflected pixel positions of the
   stars brighter than G 7 (scale 0.9-1.1 about the plate centre), against random positions.
3. The scanner-frame map: sharp single-plate sources (>= 10 sigma, no twin, no Gaia star, FWHM below the stars') that
   recur at the same scan pixel (within 2 px) on three or more different scans: marks on the scanner, vetoed by pixel.
Its settings are in config/analysis.json (section backgrounds; the scanner map's scans in section selection).
python -E -s -m heidelberg.stages.backgrounds RESULTS_DIR [OUT_DIR]  (RESULTS_DIR/hd_res/*_a.npz, ... or a measure_pair.py products dir)"""
import glob, json, os, sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from heidelberg import config

K = config.A.backgrounds
S = sys.argv[1]; O = sys.argv[2] if len(sys.argv) > 2 else S
def unit(r, d):
    r, d = np.radians(r), np.radians(d); return np.c_[np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]
rng = np.random.default_rng(K.bg_seed); fixed, ghosts, marks = [], [], []
FILES = sorted(glob.glob(f"{S}/hd_res/*.json")) or sorted(f for f in glob.glob(f"{S}/*.json") if "_" not in os.path.basename(f))   # old layout, or a products dir
for f in FILES:
    base = f[:-5]; pair = os.path.basename(base)
    try:
        o = json.load(open(f)); za, zb, g = np.load(f"{base}_a.npz"), np.load(f"{base}_b.npz"), np.load(f"{base}_gaia.npz")
    except Exception as e:
        print(f"{pair}: skipped ({e})", flush=True); continue
    if min(x["gaia_matched"] for x in o["solve"]) < config.A.solve.solve_min_stars:
        continue
    src = za["inside"] & za["twin"] & ~za["gaia"] & (za["zself"] >= K.bg_zmin)
    # 1. fixed-offset stack
    st = g["inbox"] & (g["G"] < K.bg_gaia_g); c0 = np.median(za["dec"][za["inside"]])
    S_ = np.c_[za["ra"][src] * np.cos(np.radians(c0)), za["dec"][src]] * 3600; G_ = np.c_[g["ra"][st] * np.cos(np.radians(c0)), g["dec"][st]] * 3600
    RO, NB = K.bg_offset_arcsec, K.bg_offset_bins; pr = cKDTree(G_).query_ball_point(S_, RO)
    off = np.array([S_[i] - G_[j] for i, js in enumerate(pr) for j in js]) if any(len(x) for x in pr) else np.zeros((0, 2))
    if len(off) > K.bg_min_offsets:
        r = np.hypot(*off.T); t = rng.uniform(0, 2 * np.pi, len(r)); rot = np.c_[r * np.cos(t), r * np.sin(t)]
        H, xe, ye = np.histogram2d(off[:, 0], off[:, 1], bins=NB, range=[[-RO, RO], [-RO, RO]])
        Hr = np.mean([np.histogram2d(*(np.c_[r * np.cos(tt), r * np.sin(tt)]).T, bins=NB, range=[[-RO, RO], [-RO, RO]])[0]
                      for tt in [rng.uniform(0, 2 * np.pi, len(r)) for _ in range(K.bg_rotations)]], 0)
        z = (H - Hr) / np.sqrt(Hr + 1); ks = np.argwhere(z > K.bg_peak_z)
        for kx, ky in ks:
            fixed.append(dict(pair=pair, east=round((xe[kx] + xe[kx + 1]) / 2, 1), north=round((ye[ky] + ye[ky + 1]) / 2, 1), n=int(H[kx, ky]), expected=round(float(Hr[kx, ky]), 1), z=round(float(z[kx, ky]), 1)))
    # 2. lens ghosts (pixel space of plate a)
    br = g["G"] < K.ghost_star_g
    if br.any() and "xt" in za.files:
        # the bright stars' pixels: nearest detection, or skip
        tz = cKDTree(unit(za["ra"], za["dec"])); d, j = tz.query(unit(g["ra"][br], g["dec"][br])); d = np.degrees(d) * 3600
        bx, by = za["x"][j[d < K.ghost_star_arcsec]], za["y"][j[d < K.ghost_star_arcsec]]
        if len(bx):
            cx, cy = np.median(za["x"]), np.median(za["y"]); P = np.c_[za["x"][src], za["y"][src]]; tp = cKDTree(P)
            near = lambda q: sum(len(v) for v in tp.query_ball_point(q, K.ghost_r_px))
            refl = np.concatenate([np.c_[cx - k * (bx - cx), cy - k * (by - cy)] for k in np.arange(*K.ghost_scale)])
            nr = K.ghost_random_n; rnd = np.c_[rng.uniform(P[:, 0].min(), P[:, 0].max(), nr), rng.uniform(P[:, 1].min(), P[:, 1].max(), nr)]
            ghosts.append(dict(pair=pair, bright=len(bx), at_reflection=near(refl) / len(refl), at_random=near(rnd) / len(rnd)))
    # 3. scanner marks: single-plate sharp sources in scan pixels, per plate
    for z_, tag in ((za, o["plates"][0]), (zb, o["plates"][1])):
        sing = z_["inside"] & ~z_["twin"] & ~z_["gaia"] & (z_["zself"] >= K.mark_zmin)
        wst = np.median(z_["fwhm"][z_["inside"] & z_["gaia"]]); sh = sing & (z_["fwhm"] < wst)
        marks.append((os.path.basename(tag).replace(".fits", ""), np.c_[z_["x"][sh], z_["y"][sh]]))
pd.DataFrame(fixed).to_csv(f"{O}/bg_fixed_offsets.csv", index=False)
pd.DataFrame(ghosts).to_csv(f"{O}/bg_lens_ghosts.csv", index=False)
# scanner map: positions recurring on >= 3 different scans
allp = np.concatenate([m[1] for m in marks]); lab = np.concatenate([[i] * len(m[1]) for i, m in enumerate(marks)])
t = cKDTree(allp); grp = t.query_ball_point(allp, K.mark_r_px); NS = config.A.selection.scanner_scans
rec = np.array([len(set(lab[g_])) for g_ in grp]); keep = rec >= NS
M = pd.DataFrame({"x": allp[keep, 0], "y": allp[keep, 1], "scans": rec[keep]}).round(1).drop_duplicates()
M.to_csv(f"{O}/scanner_marks.csv", index=False)
print(f"fixed-offset peaks (> {K.bg_peak_z} sd): {len(fixed)} in {len(set(x['pair'] for x in fixed))} pairs")
print("lens ghosts: ", pd.DataFrame(ghosts).describe().round(3).to_string() if ghosts else "no bright stars")
print(f"scanner marks: {len(M)} positions recurring on >= {NS} of {len(marks)} scans")

"""The sample and coverage map registered before any analysis (the challenge: teams register "their sample and coverage
map, or a checksum of them, with a timestamp before any analysis"). From the measurement alone (no selection, no
catalogue): per pair, the two plates, the twins' common exposure (dateObs +- exposure / 2 on both plates: its
midpoint and length; times from the retimed plate table, stages/plate_times.py), the footprint (the convex hull, on the
sky, of the injections introduced on both plates, which cover the common region uniformly), its area, whether both
solves were accepted, each plate's solved centre and its offset from GAVO's listed centre (requirement S8: GAVO's
centres are wrong on many plates, and the record set was drawn from them), and the role from the record set; and per
measured file its sha256, which fixes the sample. Plates a and b are treated alike throughout.
Writes OUT_DIR/coverage_map.csv, OUT_DIR/coverage_footprints.json, OUT_DIR/sample_manifest.csv, and OUT_DIR/registration.txt
(the sha256 of each of the three; the record run's analysis part checks the products against sample_manifest.csv).
python -E -s -m heidelberg.stages.coverage_map PRODUCTS_DIR LSW_PLATES.csv RECORD_SET.csv OUT_DIR [ROLE]
(ROLE: every pair's role, e.g. different-night for the controls, whose plates share no exposure: their common exposure
is zero and their time plate a's mid-exposure)"""
import glob
import hashlib
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from astropy.time import Time
from scipy.spatial import ConvexHull

from heidelberg import config

warnings.filterwarnings("ignore", message=".*dubious year")       # UTC before 1960 has no leap-second table
OFFSET_REPORT = config.A.roles.offset_report_deg                        # deg: pairs whose solved centre is further from GAVO's are counted by role (S8)


def sha(f):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def footprint(ra, dec):
    """Convex hull of (ra, dec) in the tangent plane at their centre: vertices (deg, anticlockwise) and area (deg^2)."""
    c = np.c_[np.cos(np.radians(dec)) * np.cos(np.radians(ra)), np.cos(np.radians(dec)) * np.sin(np.radians(ra)), np.sin(np.radians(dec))].mean(0)
    ra0, de0 = np.degrees(np.arctan2(c[1], c[0])) % 360, np.degrees(np.arcsin(c[2] / np.linalg.norm(c)))
    r, d, r0, d0 = np.radians(ra), np.radians(dec), np.radians(ra0), np.radians(de0)
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0)
    xi = np.degrees(np.cos(d) * np.sin(r - r0) / cosc); eta = np.degrees((np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(r - r0)) / cosc)
    H = ConvexHull(np.c_[xi, eta])
    return [[round(float(ra[k]), 5), round(float(dec[k]), 5)] for k in H.vertices], float(H.volume), (float(ra0), float(de0))


def offset_deg(ra1, dec1, ra2, dec2):
    """Angle between two positions (deg)."""
    u = lambda r, d: np.array([np.cos(np.radians(d)) * np.cos(np.radians(r)), np.cos(np.radians(d)) * np.sin(np.radians(r)), np.sin(np.radians(d))])
    return float(np.degrees(np.arccos(np.clip(u(ra1, dec1) @ u(ra2, dec2), -1, 1))))


def main():
    PD, LSW, RS, OD = sys.argv[1:5]; ROLE = sys.argv[5] if len(sys.argv) > 5 else None
    os.makedirs(OD, exist_ok=True)
    L = pd.read_csv(LSW).drop_duplicates("plateID").set_index("plateID"); R = pd.read_csv(RS).set_index("a")
    rows, polys, files = [], {}, []
    for jf in sorted(f for f in glob.glob(f"{PD}/*.json") if "_" not in os.path.basename(f)):
        p = os.path.basename(jf)[:-5]; o = json.load(open(jf))
        pa, pb = (os.path.splitext(x)[0] for x in o["plates"]); A, B = L.loc[pa], L.loc[pb]
        lo = max(A.dateObs - A.exposure / 172800, B.dateObs - B.exposure / 172800); hi = min(A.dateObs + A.exposure / 172800, B.dateObs + B.exposure / 172800)
        if hi < lo:                                           # different nights: no common exposure
            lo = hi = A.dateObs
        J = pd.read_parquet(f"{PD}/{p}_inj.parquet"); J = J[J.placed_a.astype(bool) & J.placed_b.astype(bool)]
        verts, area, (ra0, de0) = footprint(J.ra.to_numpy(float), J.dec.to_numpy(float))
        rows.append(dict(pair=p, plate_a=pa, plate_b=pb, role=ROLE or R.role.get(p, "control"),
                         mid_utc=Time(0.5 * (lo + hi), format="mjd", scale="utc").isot[:19], mid_mjd=round(0.5 * (lo + hi), 6),
                         common_exposure_s=round((hi - lo) * 86400, 1), centre_ra=round(ra0, 5), centre_dec=round(de0, 5),
                         area_deg2=round(area, 4), accepted=all(s.get("accepted") for s in o["solve"]), injections_both=len(J),
                         **{f"{k}_{x}": v for x, P_, s in zip("ab", (A, B), o["solve"]) for k, v in (
                             ("solved_ra", round(s["centre"][0], 5)), ("solved_dec", round(s["centre"][1], 5)),
                             ("gavo_offset_deg", round(offset_deg(s["centre"][0], s["centre"][1], P_.centerAlpha, P_.centerDelta), 4)))}))
        polys[p] = verts
    for f in sorted(glob.glob(f"{PD}/*")):                  # every measured file (logs aside), as the analysis part will find them
        if os.path.isfile(f) and not f.endswith(".log"):
            b = os.path.basename(f); files.append(dict(pair=b.split("_")[0].split(".")[0], file=b, bytes=os.path.getsize(f), sha256=sha(f)))
    C = pd.DataFrame(rows); C.to_csv(f"{OD}/coverage_map.csv", index=False)
    json.dump(polys, open(f"{OD}/coverage_footprints.json", "w"), separators=(",", ":"), sort_keys=True)
    pd.DataFrame(files).to_csv(f"{OD}/sample_manifest.csv", index=False)
    reg = "".join(f"{sha(f'{OD}/{n}')}  {n}\n" for n in ("coverage_map.csv", "coverage_footprints.json", "sample_manifest.csv"))
    open(f"{OD}/registration.txt", "w").write(reg)
    print(f"{len(C)} pairs ({', '.join(f'{v} {k}' for k, v in C.role.value_counts().items())}); {int(C.accepted.sum())} with both solves accepted;"
          f" {C.area_deg2.sum():.0f} deg2 in all; {len(files)} measured files")
    far = C[np.maximum(C.gavo_offset_deg_a, C.gavo_offset_deg_b) > OFFSET_REPORT]
    print(f"solved centre more than {OFFSET_REPORT:g} deg from GAVO's on either plate: {len(far)} pairs "
          f"({', '.join(f'{v} {k}' for k, v in far.role.value_counts().items()) or 'none'}); roles and membership unchanged")
    print(reg, end="")


if __name__ == "__main__":
    main()

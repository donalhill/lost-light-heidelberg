"""Check 3: the products measured here (NEW_ROOT/PAIR/) against a reference run's (REF_DIR, flat), pair by pair: the same
files, arrays, shapes, and columns; integers, booleans, and text exactly equal; float64 within RTOL relative plus ATOL absolute,
float32 within SCALE32 of the array's largest value (NaN equal to NaN), since machines of different sizes sum in a different order (numpy's linear algebra splits
its sums across the machine's threads); the largest relative difference is printed. Every PAIR.json likewise, except
max_rss_gb and the manifest's input and code hashes, which differ by design (the inputs are scans restored from their
compressed copies, whose integrity cards differ; the code hash covers every library file). The fidelity report
(PAIR_fidelity.json) is the fidelity step's, not the measurement's, and is not compared. Then every RESTORE file must
say the copy restores.
python tools/same_products.py NEW_ROOT REF_DIR RESTORE..."""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

RTOL, ATOL, SCALE32 = 1e-5, 1e-8, 1e-6     # float64 relative (rounding passes through arctangents and the like); float32
# (the flash stamps) within SCALE32 of the array's largest value, since their rounding is a fixed quantum (0.004 counts)
# from sums at the image's scale, which a relative test fails on pixels near zero
new, ref, restore = sys.argv[1], sys.argv[2], sys.argv[3:]
bad, notes, n = [], [], 0
worst = [0.0, ""]


def close(a, b, what):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or (a.dtype.kind in "fc") != (b.dtype.kind in "fc"):
        return False
    if a.dtype.kind not in "fc":
        return bool(pd.Series(a.ravel()).equals(pd.Series(b.ravel())))
    f32 = a.dtype == np.float32
    a, b = a.astype(float), b.astype(float)
    if not (np.array_equal(np.isnan(a), np.isnan(b)) and np.array_equal(a[np.isinf(a)], b[np.isinf(a)])):
        return False
    m = np.isfinite(a) & np.isfinite(b)
    if m.any():
        r = float(np.max(np.abs(a[m] - b[m]) / (np.abs(b[m]) + ATOL)))
        if r > worst[0]:
            worst[:] = [r, what]
    if f32:
        return bool(np.all(np.abs(a[m] - b[m]) <= SCALE32 * np.max(np.abs(b[m]))))
    return bool(np.allclose(a[m], b[m], rtol=RTOL, atol=ATOL))


def same_json(a, b, what):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same_json(a[k], b[k], what) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same_json(x, y, what) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        return isinstance(a, (int, float)) and isinstance(b, (int, float)) and close([a], [b], what)
    return a == b
for d in sorted(glob.glob(f"{new}/*/")):
    pair = os.path.basename(d.rstrip("/"))
    want = sorted(f for f in os.listdir(ref) if f.startswith(pair) and not f.endswith((".log", "_fidelity.json")))
    if not want:
        bad.append(f"{pair}: no reference products"); continue
    for f in want:
        a, b = os.path.join(d, f), os.path.join(ref, f)
        if not os.path.exists(a):
            bad.append(f"{pair}: {f} missing"); continue
        if f.endswith(".npz"):
            A, B = np.load(a, allow_pickle=False), np.load(b, allow_pickle=False)
            ok = sorted(A.files) == sorted(B.files) and all(close(A[k], B[k], f"{f}:{k}") for k in A.files)
        elif f.endswith(".parquet"):
            A, B = pd.read_parquet(a), pd.read_parquet(b)
            ok = list(A.columns) == list(B.columns) and len(A) == len(B) and all(close(A[c].to_numpy(), B[c].to_numpy(), f"{f}:{c}") for c in A.columns)
        elif f.endswith(".json"):
            A, B = json.load(open(a)), json.load(open(b))
            for J in (A, B):
                J.pop("max_rss_gb", None); m = J.get("manifest", {}); m.pop("inputs", None); m.pop("code", None)
            if A.get("manifest", {}).get("outputs") != B.get("manifest", {}).get("outputs"):
                notes.append(f"{pair}: product file hashes differ (contents compared above)")
            for J in (A, B):
                J.get("manifest", {}).pop("outputs", None)
            ok = same_json(A, B, f)
        else:
            ok = open(a, "rb").read() == open(b, "rb").read()
        n += 1
        if not ok:
            bad.append(f"{pair}: {f} differs")
for r in restore:
    if "restores" not in open(r).read():
        bad.append(f"{os.path.basename(r)}: the copy does not restore")
for x in bad + notes:
    print(x)
print(f"largest relative difference {worst[0]:.2g} ({worst[1]})")
print(f"{n} product files compared, {len(restore)} copies checked: " + (f"equal within rounding" if not bad else f"{len(bad)} differences"))
sys.exit(1 if bad else 0)

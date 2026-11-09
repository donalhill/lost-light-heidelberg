"""Whether every file of a rebuilt output folder equals the reference's: the same file names, parquet tables equal by
value (floats within tools/same.py's 1e-9, since x86 and arm64 differ in the last digits), JSON equal by value except run_key (a cache key hashing every library file, so it changes
with the set of modules, not with the result), CSV and .npy by tools/same.py, anything else byte for byte. Log files are
skipped. Prints each difference and a summary line.
python tools/same_dir.py NEW REFERENCE"""
import filecmp
import json
import os
import sys

import pandas as pd

from same import close, same

new, ref = sys.argv[1:3]
files = lambda d: {os.path.relpath(os.path.join(r, f), d) for r, _, fs in os.walk(d) for f in fs if not f.endswith(".log") and not f.startswith(".snakemake")}
A, B = files(new), files(ref)
bad = [f"only in {new}: {f}" for f in sorted(A - B)] + [f"only in {ref}: {f}" for f in sorted(B - A)]
for f in sorted(A & B):
    a, b = os.path.join(new, f), os.path.join(ref, f)
    if f.endswith(".parquet"):
        Ta, Tb = pd.read_parquet(a), pd.read_parquet(b)
        ok = list(Ta.columns) == list(Tb.columns) and len(Ta) == len(Tb) and all(close(Ta[c].to_numpy(), Tb[c].to_numpy())[0] for c in Ta.columns)
    elif f.endswith(".json"):
        ja, jb = json.load(open(a)), json.load(open(b))
        for j in (ja, jb):
            if isinstance(j, dict):
                j.pop("run_key", None)
        ok = ja == jb
    elif f.endswith((".csv", ".npy")):
        ok = same(a, b)[0]
    else:
        ok = filecmp.cmp(a, b, shallow=False)
    if not ok:
        bad.append(f"differs: {f}")
for b in bad[:40]:
    print(b)
print(f"{len(A & B)} files compared: " + ("all equal" if not bad else f"{len(bad)} differences"))
sys.exit(1 if bad else 0)

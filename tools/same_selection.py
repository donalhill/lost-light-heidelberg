"""The offline repeat of a selection against the online run: every pair's flag table and sealed candidates identical
(every row, flag, and value), and every PAIR_select.json equal apart from the fields that say how it was run (offline,
cache_used, run_key). Differences are listed; the last line says whether the repeat is reproducible.
python tools/same_selection.py ONLINE_DIR OFFLINE_DIR"""
import glob
import json
import os
import sys

import pandas as pd

RUN = {"offline", "cache_used", "run_key"}
a, b = sys.argv[1:3]; bad = []
tables = sorted(glob.glob(f"{a}/*_flags.parquet") + glob.glob(f"{a}/sealed/*_candidates.parquet"))
for f in tables:
    g = f.replace(a, b, 1)
    if not os.path.exists(g) or not pd.read_parquet(f).equals(pd.read_parquet(g)):
        bad.append(f"differs: {os.path.relpath(f, a)}")
for f in sorted(glob.glob(f"{a}/*_select.json") + glob.glob(f"{a}/sealed/*_select.json")):
    g = f.replace(a, b, 1)
    if not os.path.exists(g):
        bad.append(f"missing: {os.path.relpath(g, b)}"); continue
    x, y = json.load(open(f)), json.load(open(g))
    keys = sorted(k for k in set(x) | set(y) if k not in RUN and x.get(k) != y.get(k))
    if keys:
        bad.append(f"{os.path.relpath(f, a)}: {keys}")
for x in bad[:40]:
    print(x)
print(f"{len(tables)} tables compared: " + ("reproducible" if not bad else f"{len(bad)} differences"))
sys.exit(1 if bad else 0)

"""The injection fidelity report (heidelberg/lib/fidelity.py) recomputed from a pair's saved products with the current
code, written beside them as PAIR_fidelity.json, which heidelberg/stages/select.py requires for every pair whose solve
was accepted. The measured products are not touched.
python -E -s -m heidelberg.stages.fidelity PRODUCTS_DIR [PAIRS...]"""
import glob
import json
import os
import sys

from heidelberg.lib import fidelity as FID

PD = sys.argv[1]; ONLY = set(sys.argv[2:])
for f in sorted(glob.glob(f"{PD}/*.json")):
    pair = os.path.basename(f)[:-5]
    if "_" in pair or (ONLY and pair not in ONLY):
        continue
    o = json.load(open(f))
    if not all(s.get("accepted") for s in o["solve"]):
        continue
    r = FID.report(f"{PD}/{pair}", f"{PD}/{o['gq']}")
    json.dump(r, open(f"{PD}/{pair}_fidelity.json", "w"), indent=1, default=float)
    print(pair, "passed" if r["passed"] else "FAILED", f"noise {r['background_noise_ratio']:.3f}", flush=True)

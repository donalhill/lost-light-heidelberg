"""The asteroid clock's plate list (requirement N2) from the measured products: per plate its solved centre and whether
the measurement used the scan's header WCS or the meridian search (heidelberg/lib/meridian.py), as clock_batch.py reads it.
python -E -s -m heidelberg.stages.clock_plates PRODUCTS_DIR OUT.csv"""
import glob
import json
import os
import sys

import pandas as pd

PD, OUT = sys.argv[1:3]; rows = []
for f in sorted(glob.glob(f"{PD}/*.json")):
    if "_" in os.path.basename(f):
        continue
    o = json.load(open(f))
    for pl, s in zip(o["plates"], o["solve"]):
        if s.get("accepted"):
            rows.append(dict(plate=os.path.splitext(pl)[0], ra_c=s["centre"][0], dec_c=s["centre"][1],
                             wcs="meridian" if s.get("wcs", {}).get("source") == "meridian" else "header"))
pd.DataFrame(rows).to_csv(OUT, index=False); print(f"{len(rows)} plates with accepted solves -> {OUT}")

"""Part two's first check: every registered product file present with its registered SHA-256, and none extra.
python tools/check_manifest.py SAMPLE_MANIFEST.csv PRODUCTS_DIR"""
import glob
import hashlib
import os
import sys

import pandas as pd

M, d = pd.read_csv(sys.argv[1]), sys.argv[2]
have = {os.path.basename(f) for f in glob.glob(f"{d}/*") if not f.endswith(".log")}
bad = [r.file for r in M.itertuples() if not os.path.exists(f"{d}/{r.file}") or hashlib.sha256(open(f"{d}/{r.file}", "rb").read()).hexdigest() != r.sha256]
bad += [f"not registered: {x}" for x in sorted(have - set(M.file))]
for b in bad[:20]:
    print("DIFFERS or missing:", b)
print(f"{len(M)} files checked: " + ("as registered" if not bad else f"{len(bad)} differ")); sys.exit(1 if bad else 0)
